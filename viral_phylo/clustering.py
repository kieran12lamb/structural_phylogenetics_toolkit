#!/usr/bin/env python3
"""Similarity-driven clustering of structures prior to alignment and tree inference.

The previous strategy aligned every structure together and then clustered greedily
on that master alignment, defining "coverage" as the count of mutually non-gapped
columns. That is circular - it judges similarity through an alignment which is
itself poor when the inputs are heterogeneous - and being greedy and
longest-representative-first it produced one dominant cluster plus a long tail of
singletons (226 clusters from 1193 taxa: 104 singletons, only 22 with >= 10 taxa,
from a master alignment that was 91.6% gaps).

This module instead measures similarity directly:

  * ``mmseqs cluster``   - amino-acid sequence identity
  * ``foldseek cluster`` - structural/3Di identity

A sweep over identity and coverage thresholds yields a family of candidate
partitions. Each is scored by silhouette in two independent dense metric spaces
(ESM-2 embeddings, and an all-vs-all identity matrix) and the selection is
*constrained* rather than maximised blind: silhouette alone rewards splitting and
would happily reproduce the singleton pathology it is meant to solve.
"""

import glob
import json
import os
import shutil
import subprocess
from typing import Any, Dict, List, Optional, Sequence, Tuple

import numpy as np

from viral_phylo.binaries import find_foldseek_bin, find_mmseqs_bin, tool_is_available

# Amino-acid clustering sweeps sequence identity; structural clustering sweeps
# TM-score, where 0.5 is the conventional "same fold" boundary. These are different
# quantities and must not share a grid.
DEFAULT_MIN_SEQ_IDS: Tuple[float, ...] = (0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9)
DEFAULT_TMSCORE_THRESHOLDS: Tuple[float, ...] = (0.4, 0.5, 0.6, 0.7)
DEFAULT_COVERAGES: Tuple[float, ...] = (0.5, 0.7, 0.8)

# Default structural operating point: foldseek's native 3Di+AA alignment at
# e-value 1.0 and 50% coverage. On the 1,135-taxa Nipah binder set this yields 85
# clusters, 21 of them with 10+ taxa, retaining 87% of taxa.
#
# Measured against TM-align at its best setting (TM>=0.5, c=0.7), 3Di+AA is better
# on every axis - 87% vs 72% retention, +0.131 vs +0.102 ESM-2 silhouette, +0.284
# vs +0.233 silhouette in TM-score space itself, 85 vs 190 clusters, and a 54s
# sweep against 152s. It is also the mode foldseek is designed around, and is
# consistent with the rest of this toolkit, which works in 3Di throughout.
#
# Caveat worth remembering: an e-value depends on database size, so it is NOT
# comparable across datasets of different sizes the way a TM-score is. The taxa
# count is recorded alongside it in the manifest for that reason. Use
# alignment_type=1 with a tmscore threshold when cross-dataset comparability
# matters more than quality on this one.
DEFAULT_EVALUE = 1.0
DEFAULT_CLUSTER_COVERAGE = 0.50
DEFAULT_EVALUES: Tuple[float, ...] = (0.001, 0.01, 1.0, 10.0)

# Retained for the optional TM-align mode (alignment_type=1).
DEFAULT_TMSCORE = 0.50

# How much extra retention a more inclusive operating point must buy before it is
# worth putting in front of the user.
SUGGESTION_MIN_RETENTION_GAIN = 0.10

# Thresholds for splitting an already-formed cluster further. These are deliberately
# looser than the whole-set gate: that gate decides "can this entire heterogeneous
# input be one alignment", whereas a cluster alignment at 40-55% gaps and ~2x
# expansion is ordinary for a diverse protein family. Measured on the Nipah binder
# set, the whole-set gate's 50%/1.5x would have flagged 19 of 21 usable clusters and
# fragmented the partition; 65%/3.0x flags the two genuinely heterogeneous ones
# (244 and 88 taxa, at 65%/3.3x and 62%/3.0x).
RECURSE_MAX_GAP_FRACTION = 0.65
RECURSE_MAX_EXPANSION = 3.0

# A cluster needs at least this many taxa for a tree (IQ-TREE refuses bootstrap below 4).
TREE_MIN_TAXA = 4

# Finer e-values tried when splitting a cluster, gentlest first. The original sweep
# grid jumped straight from 1.0 to 0.01; re-clustering an 88-taxon Nipah cluster at
# 0.01 gave 55 clusters, 41 of them singletons, keeping only 28% of its taxa in
# tree-eligible clusters. At 0.3 it kept 88%, at 0.1 69%.
RECURSE_EVALUE_STEPS: Tuple[float, ...] = (10.0, 1.0, 0.3, 0.1, 0.03, 0.01, 0.001)
RECURSE_MAX_ATTEMPTS = 3

# A split is only accepted if it keeps at least this fraction of the parent's taxa in
# tree-eligible clusters. Otherwise it is shattering the cluster rather than refining
# it, and the cluster is kept whole and reported as not divisible.
RECURSE_MIN_KEPT = 0.70

# Retention differences within this band are treated as equivalent when ranking
# suggestions, so the tie is broken on how many usable clusters the setting yields.
RETENTION_TIE_BAND = 0.03


class ClusteringToolMissing(RuntimeError):
    """Raised when a required external clustering tool is not installed."""


# --------------------------------------------------------------------------- QC


def parse_fasta(path: str) -> Dict[str, str]:
    """Read a FASTA file into an ordered ``{id: sequence}`` mapping."""
    records: Dict[str, str] = {}
    name: Optional[str] = None
    chunks: List[str] = []
    with open(path, "r", encoding="utf-8", errors="ignore") as handle:
        for line in handle:
            line = line.strip()
            if line.startswith(">"):
                if name is not None:
                    records[name] = "".join(chunks)
                name = line[1:].split()[0] if len(line) > 1 else ""
                chunks = []
            elif name is not None:
                chunks.append(line)
    if name is not None:
        records[name] = "".join(chunks)
    return records


def qc_filter_sequences(
    sequences: Dict[str, str],
    min_length: int = 50,
    min_length_frac: float = 0.4,
) -> Tuple[Dict[str, str], List[Dict[str, Any]]]:
    """Drop sequences too short to align or cluster meaningfully.

    Length heterogeneity is what inflates a global alignment: a 9-residue fragment
    and a 259-residue protein cannot share meaningful columns, and forcing them
    into one alignment produces a staircase of gaps. Both an absolute floor and a
    floor relative to the median are applied; every exclusion is returned with its
    reason so the caller can record it rather than dropping it silently.
    """
    ungapped = {name: seq.replace("-", "") for name, seq in sequences.items()}
    lengths = sorted(len(s) for s in ungapped.values())
    if not lengths:
        return {}, []

    median = lengths[len(lengths) // 2]
    relative_floor = int(median * min_length_frac)
    floor = max(int(min_length), relative_floor)

    kept: Dict[str, str] = {}
    excluded: List[Dict[str, Any]] = []
    for name, seq in sequences.items():
        length = len(ungapped[name])
        if length < floor:
            excluded.append({
                "taxon_id": name,
                "length": length,
                "reason": f"length {length} < floor {floor} (max of min_length={min_length}, "
                          f"{min_length_frac:g} x median {median})",
            })
        else:
            kept[name] = seq
    return kept, excluded


def write_fasta(sequences: Dict[str, str], path: str, ungap: bool = True) -> str:
    """Write sequences to FASTA, stripping alignment gaps by default."""
    os.makedirs(os.path.dirname(os.path.abspath(path)) or ".", exist_ok=True)
    with open(path, "w", encoding="utf-8") as handle:
        for name, seq in sequences.items():
            payload = seq.replace("-", "") if ungap else seq
            handle.write(f">{name}\n")
            for i in range(0, len(payload), 60):
                handle.write(payload[i:i + 60] + "\n")
    return path


def alignment_gap_fraction(alignment: Dict[str, str]) -> Dict[str, float]:
    """Summarise how gappy an alignment is, and how far it is stretched.

    ``expansion`` is alignment width divided by median ungapped length: the factor
    by which the aligner had to stretch a typical sequence to fit.
    """
    if not alignment:
        return {"taxa": 0, "columns": 0, "gap_fraction": 0.0, "expansion": 0.0}
    columns = len(next(iter(alignment.values())))
    taxa = len(alignment)
    gaps = sum(seq.count("-") for seq in alignment.values())
    ungapped = sorted(len(seq.replace("-", "")) for seq in alignment.values())
    median = ungapped[len(ungapped) // 2] or 1
    return {
        "taxa": taxa,
        "columns": columns,
        "gap_fraction": round(gaps / (taxa * columns), 4) if taxa and columns else 0.0,
        "median_ungapped_length": median,
        "expansion": round(columns / median, 2),
    }


def recursion_schedule(threshold_kind: str, current: float,
                       max_attempts: int = RECURSE_MAX_ATTEMPTS) -> List[float]:
    """Thresholds to try, in order, when splitting one cluster.

    For e-values the current value comes first: an e-value depends on database size,
    so the same nominal value on a subset is already stricter - on the 88-taxon
    cluster it split it into 6 while keeping 94% of taxa. Then progressively finer
    values, in small steps rather than the 100x jump of the sweep grid. TM-score and
    sequence identity do not depend on database size, so those start one step stricter.
    """
    if threshold_kind == "evalue":
        finer = [e for e in RECURSE_EVALUE_STEPS if e < current]
        return ([current] + finer)[:max_attempts]
    schedule = []
    value = current
    for _ in range(max_attempts):
        value = stricter_threshold(threshold_kind, value)
        if value is None:
            break
        schedule.append(value)
    return schedule


def fraction_kept(sub_clusters: Dict[str, List[str]], parent_size: int,
                  min_taxa: int = TREE_MIN_TAXA) -> float:
    """Fraction of a parent cluster's taxa that land in tree-eligible sub-clusters."""
    if not parent_size:
        return 0.0
    kept = sum(len(m) for m in sub_clusters.values() if len(m) >= min_taxa)
    return kept / parent_size


def stricter_threshold(threshold_kind: str, current: float) -> Optional[float]:
    """Next value along the sweep grid that clusters more strictly.

    For e-values stricter means smaller; for TM-score and sequence identity it means
    larger. Returns None when the grid is exhausted.
    """
    if threshold_kind == "evalue":
        candidates = [e for e in DEFAULT_EVALUES if e < current]
        return max(candidates) if candidates else None
    grid = DEFAULT_TMSCORE_THRESHOLDS if threshold_kind == "tmscore" else DEFAULT_MIN_SEQ_IDS
    candidates = [t for t in grid if t > current]
    return min(candidates) if candidates else None


def alignment_needs_splitting(
    stats: Dict[str, float],
    max_gap_fraction: float = RECURSE_MAX_GAP_FRACTION,
    max_expansion: float = RECURSE_MAX_EXPANSION,
) -> Tuple[bool, str]:
    """Whether a cluster's own alignment is poor enough to warrant re-clustering."""
    gap = stats.get("gap_fraction", 0.0)
    expansion = stats.get("expansion", 0.0)
    reasons = []
    if gap >= max_gap_fraction:
        reasons.append(f"{gap:.0%} gaps (limit {max_gap_fraction:.0%})")
    if expansion > max_expansion:
        reasons.append(f"{expansion:.2f}x expansion (limit {max_expansion:g}x)")
    return (bool(reasons), "; ".join(reasons) if reasons else "alignment quality acceptable")


def cluster_needs_split(
    stats: Optional[Dict[str, Any]],
    recurse_on: str = "3di",
    max_gap_fraction: float = RECURSE_MAX_GAP_FRACTION,
    max_expansion: float = RECURSE_MAX_EXPANSION,
) -> Tuple[bool, str]:
    """Decide whether a cluster should be split, from its alignment statistics.

    ``stats`` holds the 3Di alignment's numbers at the top level and the AA
    alignment's under ``"aa"``. ``recurse_on="3di"`` (the default) judges only the
    3Di alignment, since structural alignment is the purpose of this tool.
    ``recurse_on="both"`` also splits when the AA alignment is over threshold - useful
    with ``--tree-type both`` under MAFFT, where the AA alignment is built separately
    and is typically gappier.
    """
    if not stats:
        return False, "no alignment statistics"
    split_3di, why_3di = alignment_needs_splitting(stats, max_gap_fraction, max_expansion)
    if recurse_on != "both":
        return split_3di, (f"3Di: {why_3di}" if split_3di else why_3di)
    aa = stats.get("aa") or {}
    split_aa, why_aa = alignment_needs_splitting(aa, max_gap_fraction, max_expansion) if aa else (False, "")
    reasons = [f"3Di: {why_3di}"] if split_3di else []
    if split_aa:
        reasons.append(f"AA: {why_aa}")
    return (split_3di or split_aa), ("; ".join(reasons) if reasons else "alignment quality acceptable")


def full_set_alignment_is_usable(
    stats: Dict[str, float],
    max_gap_fraction: float = 0.50,
    max_expansion: float = 1.5,
) -> Tuple[bool, str]:
    """Decide whether the whole input can reasonably be aligned as one set.

    Automates the judgement that previously had to be made by eye. On the 25 Sep
    Nipah run this returns False: 91.6% gaps and a 12.15x expansion.
    """
    gap = stats.get("gap_fraction", 1.0)
    expansion = stats.get("expansion", 99.0)
    reasons = []
    if gap >= max_gap_fraction:
        reasons.append(f"gap fraction {gap:.1%} >= {max_gap_fraction:.0%}")
    if expansion > max_expansion:
        reasons.append(f"alignment is {expansion:.2f}x the median sequence length (limit {max_expansion:g}x)")
    if reasons:
        return False, "; ".join(reasons)
    return True, f"gap fraction {gap:.1%}, expansion {expansion:.2f}x"


# ------------------------------------------------------------------- clustering


def _run(cmd: Sequence[str], **kwargs) -> subprocess.CompletedProcess:
    return subprocess.run(list(cmd), capture_output=True, text=True, check=True, **kwargs)


def _parse_cluster_tsv(path: str) -> Dict[str, List[str]]:
    """Parse an mmseqs/foldseek ``*_cluster.tsv`` (representative<TAB>member)."""
    clusters: Dict[str, List[str]] = {}
    with open(path, "r", encoding="utf-8", errors="ignore") as handle:
        for line in handle:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 2:
                continue
            rep, member = parts[0].strip(), parts[1].strip()
            # Foldseek reports structure filenames; normalise to the taxon id.
            rep = _strip_structure_suffix(rep)
            member = _strip_structure_suffix(member)
            clusters.setdefault(rep, []).append(member)
    return clusters


def _strip_structure_suffix(name: str) -> str:
    for suffix in (".pdb", ".cif", ".mmcif", ".pdb.gz", ".cif.gz"):
        if name.endswith(suffix):
            return name[: -len(suffix)]
    return name


def mmseqs_cluster(
    fasta_path: str,
    workdir: str,
    min_seq_id: float,
    coverage: float,
    cov_mode: int = 0,
    threads: int = 0,
    mmseqs_bin: Optional[str] = None,
) -> Dict[str, List[str]]:
    """Cluster amino-acid sequences with ``mmseqs easy-cluster``."""
    if not tool_is_available("mmseqs", override=mmseqs_bin):
        raise ClusteringToolMissing(
            "mmseqs is not installed. Install with: conda install -c bioconda mmseqs2")
    binary = find_mmseqs_bin(mmseqs_bin)
    os.makedirs(workdir, exist_ok=True)
    prefix = os.path.join(workdir, "mm")
    tmp = os.path.join(workdir, "tmp")
    cmd = [
        binary, "easy-cluster", fasta_path, prefix, tmp,
        "--min-seq-id", str(min_seq_id),
        "-c", str(coverage),
        "--cov-mode", str(cov_mode),
        "-v", "1",
    ]
    if threads:
        cmd += ["--threads", str(threads)]
    _run(cmd)
    tsv = f"{prefix}_cluster.tsv"
    clusters = _parse_cluster_tsv(tsv)
    shutil.rmtree(tmp, ignore_errors=True)
    return clusters


def foldseek_cluster(
    structure_dir: str,
    workdir: str,
    evalue: float = DEFAULT_EVALUE,
    coverage: float = DEFAULT_CLUSTER_COVERAGE,
    cov_mode: int = 0,
    alignment_type: int = 2,
    tmscore_threshold: Optional[float] = None,
    min_seq_id: float = 0.0,
    threads: int = 0,
    foldseek_bin: Optional[str] = None,
) -> Dict[str, List[str]]:
    """Cluster structures with ``foldseek easy-cluster`` on structural similarity.

    Deliberately does NOT impose a sequence-identity threshold. Structural
    clustering exists to group things whose folds match even when their sequences
    have diverged, so requiring sequence identity as well defeats the purpose - on
    the Nipah binder set, adding ``--min-seq-id 0.3`` inflated the result from 345
    clusters to 480 and cut retention from 55% to 42%.

    Defaults to ``alignment_type=2``, foldseek's native 3Di+amino-acid alignment,
    where ``evalue`` is the control. Passing ``alignment_type=1`` selects TM-align
    instead, in which case ``tmscore_threshold`` applies; that is slower and scored
    worse here, but a TM-score is comparable across datasets whereas an e-value is
    not. See the DEFAULT_EVALUE comment for the measured comparison.
    """
    if not tool_is_available("foldseek", override=foldseek_bin):
        raise ClusteringToolMissing(
            "foldseek is not installed. Install with: conda install -c bioconda foldseek")
    binary = find_foldseek_bin(foldseek_bin)
    os.makedirs(workdir, exist_ok=True)
    prefix = os.path.join(workdir, "fs")
    tmp = os.path.join(workdir, "tmp")
    cmd = [
        binary, "easy-cluster", structure_dir, prefix, tmp,
        "--min-seq-id", str(min_seq_id),
        "-c", str(coverage),
        "--cov-mode", str(cov_mode),
        "--alignment-type", str(alignment_type),
        "-v", "1",
    ]
    if alignment_type == 1:
        cmd += ["--tmscore-threshold", str(tmscore_threshold
                                           if tmscore_threshold is not None else DEFAULT_TMSCORE)]
    else:
        cmd += ["-e", str(evalue)]
    if threads:
        cmd += ["--threads", str(threads)]
    _run(cmd)
    tsv = f"{prefix}_cluster.tsv"
    clusters = _parse_cluster_tsv(tsv)
    shutil.rmtree(tmp, ignore_errors=True)
    return clusters


# ------------------------------------------------------------- distance spaces


def embedding_distance_matrix(npz_path: str) -> Tuple[List[str], np.ndarray]:
    """Dense cosine distance matrix from a saved ESM-2 embedding bundle.

    Reuses embeddings the pipeline has already paid to compute. Because this space
    is independent of the mmseqs/foldseek identity used to form the clusters, it
    validates a partition rather than grading its own homework.
    """
    data = np.load(npz_path, allow_pickle=True)
    taxa = [str(t) for t in data["taxa"]]
    embeddings = np.asarray(data["embeddings"], dtype=np.float64)
    norms = np.linalg.norm(embeddings, axis=1, keepdims=True)
    norms = np.where(norms == 0, 1e-12, norms)
    unit = embeddings / norms
    similarity = np.clip(unit @ unit.T, -1.0, 1.0)
    distance = 1.0 - similarity
    np.fill_diagonal(distance, 0.0)
    return taxa, distance


def identity_distance_matrix(
    query: str,
    taxa: Sequence[str],
    workdir: str,
    tool: str = "mmseqs",
    binary: Optional[str] = None,
    threads: int = 0,
    max_seqs: int = 5000,
    score: str = "fident",
    evalue: str = "10000",
    alignment_type: Optional[int] = None,
) -> Tuple[List[str], np.ndarray]:
    """Dense distance matrix, ``1 - score``, from an all-vs-all search.

    ``score="fident"`` gives sequence-identity distance (mmseqs). For structures,
    foldseek with ``score="alntmscore"`` gives structural distance, 1 - TM-score.
    Foldseek's native 3Di+AA search (``alignment_type=2``) reports TM-scores too, and
    on the 1,135-structure Nipah set it took 9 s against 131 s for a TM-align search,
    ranked all 12 sweep partitions the same way (Spearman rho 0.993) and selected
    the same one - so it is the default for the structural space.

    mmseqs and foldseek only emit hits above threshold, so the matrix is sparse by
    nature; unreported pairs are set to the maximum distance of 1.0. That
    completion biases silhouette values, which is why this space is reported
    alongside the embedding space rather than used to select on its own.
    """
    if not tool_is_available(tool, override=binary):
        raise ClusteringToolMissing(f"{tool} is not installed.")
    resolved = find_mmseqs_bin(binary) if tool == "mmseqs" else find_foldseek_bin(binary)
    os.makedirs(workdir, exist_ok=True)
    out_tsv = os.path.join(workdir, f"{tool}_allvall.tsv")
    tmp = os.path.join(workdir, "tmp_search")
    cmd = [
        resolved, "easy-search", query, query, out_tsv, tmp,
        "--format-output", f"query,target,{score}",
        "-e", str(evalue),
        "--max-seqs", str(max_seqs),
        "-v", "1",
    ]
    if alignment_type is not None:
        cmd += ["--alignment-type", str(alignment_type)]
    if threads:
        cmd += ["--threads", str(threads)]
    _run(cmd)

    order = {name: i for i, name in enumerate(taxa)}
    distance = np.ones((len(taxa), len(taxa)), dtype=np.float64)
    np.fill_diagonal(distance, 0.0)
    with open(out_tsv, "r", encoding="utf-8", errors="ignore") as handle:
        for line in handle:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 3:
                continue
            q = _strip_structure_suffix(parts[0].strip())
            t = _strip_structure_suffix(parts[1].strip())
            if q not in order or t not in order:
                continue
            try:
                fident = float(parts[2])
            except ValueError:
                continue
            d = max(0.0, min(1.0, 1.0 - fident))
            i, j = order[q], order[t]
            # Keep the best (smallest) distance seen for a pair, and symmetrise.
            if d < distance[i, j]:
                distance[i, j] = d
                distance[j, i] = d
    np.fill_diagonal(distance, 0.0)
    shutil.rmtree(tmp, ignore_errors=True)
    return list(taxa), distance


# ----------------------------------------------------------------------- scoring


def clusters_to_labels(clusters: Dict[str, List[str]], taxa: Sequence[str]) -> np.ndarray:
    """Map ``{representative: members}`` onto an integer label per taxon.

    Taxa absent from the clustering are labelled -1.
    """
    assignment = {}
    for index, (_rep, members) in enumerate(sorted(clusters.items())):
        for member in members:
            assignment[member] = index
    return np.array([assignment.get(name, -1) for name in taxa], dtype=int)


def partition_shape(labels: np.ndarray, min_cluster_size: int) -> Dict[str, Any]:
    """Describe a partition's size distribution and how much of the data it retains."""
    valid = labels[labels >= 0]
    sizes: Dict[int, int] = {}
    for label in valid:
        sizes[int(label)] = sizes.get(int(label), 0) + 1
    size_list = sorted(sizes.values(), reverse=True)
    retained = sum(s for s in size_list if s >= min_cluster_size)
    total = int(labels.size)
    return {
        "n_clusters": len(size_list),
        "sizes": size_list[:25],
        "largest": size_list[0] if size_list else 0,
        "singletons": sum(1 for s in size_list if s == 1),
        "clusters_at_or_above_min": sum(1 for s in size_list if s >= min_cluster_size),
        "taxa_total": total,
        "taxa_retained": retained,
        "retained_fraction": round(retained / total, 4) if total else 0.0,
        "unassigned": int((labels < 0).sum()),
    }


def silhouette_for_partition(
    labels: np.ndarray,
    distance: np.ndarray,
    min_cluster_size: int,
) -> Optional[float]:
    """Silhouette over the taxa that fall in clusters of at least ``min_cluster_size``.

    Small and singleton clusters are excluded from the score rather than allowed to
    inflate it; a partition of mostly singletons would otherwise look excellent.
    """
    from sklearn.metrics import silhouette_score

    sizes: Dict[int, int] = {}
    for label in labels:
        if label >= 0:
            sizes[int(label)] = sizes.get(int(label), 0) + 1
    keep = np.array([lab >= 0 and sizes.get(int(lab), 0) >= min_cluster_size for lab in labels])
    if keep.sum() < 3:
        return None
    kept_labels = labels[keep]
    if len(set(kept_labels.tolist())) < 2:
        return None
    sub = distance[np.ix_(keep, keep)]
    try:
        return float(silhouette_score(sub, kept_labels, metric="precomputed"))
    except Exception:
        return None


def score_partition(
    clusters: Dict[str, List[str]],
    taxa: Sequence[str],
    spaces: Dict[str, np.ndarray],
    min_cluster_size: int,
) -> Dict[str, Any]:
    """Score one candidate partition in every available distance space."""
    labels = clusters_to_labels(clusters, taxa)
    result = partition_shape(labels, min_cluster_size)
    result["silhouette"] = {
        name: silhouette_for_partition(labels, matrix, min_cluster_size)
        for name, matrix in spaces.items()
    }
    return result


# --------------------------------------------------------------- sweep + select


def sweep_clusterings(
    taxa: Sequence[str],
    spaces: Dict[str, np.ndarray],
    cluster_fn,
    workdir: str,
    thresholds: Sequence[float] = DEFAULT_MIN_SEQ_IDS,
    coverages: Sequence[float] = DEFAULT_COVERAGES,
    min_cluster_size: int = 10,
    label: str = "mmseqs",
    threshold_kind: str = "min_seq_id",
    verbose: bool = True,
) -> List[Dict[str, Any]]:
    """Cluster across a grid of thresholds and coverages, scoring each result.

    ``threshold_kind`` records what the first axis means - sequence identity for
    mmseqs, TM-score for foldseek - so downstream readers cannot confuse the two.
    """
    candidates: List[Dict[str, Any]] = []
    for threshold in thresholds:
        for coverage in coverages:
            run_dir = os.path.join(workdir, f"{label}_t{int(threshold * 100)}_c{int(coverage * 100)}")
            try:
                clusters = cluster_fn(run_dir, threshold, coverage)
            except ClusteringToolMissing:
                raise
            except subprocess.CalledProcessError as exc:
                if verbose:
                    tail = (exc.stderr or "").strip().splitlines()[-1:] or [str(exc)]
                    print(f"[Cluster] [!] {label} failed at {threshold_kind}={threshold} "
                          f"c={coverage}: {tail[0]}")
                continue

            scored = score_partition(clusters, taxa, spaces, min_cluster_size)
            scored.update({"tool": label, "threshold": threshold, "threshold_kind": threshold_kind,
                           "coverage": coverage, "clusters": clusters})
            candidates.append(scored)
            if verbose:
                sil = scored["silhouette"]
                pretty = ", ".join(
                    f"{k}={'None' if v is None else f'{v:.3f}'}" for k, v in sil.items())
                print(f"[Cluster] {label} {threshold_kind}={threshold:.2f} c={coverage:.2f} -> "
                      f"{scored['n_clusters']:>4} clusters, "
                      f"{scored['clusters_at_or_above_min']:>3} >= {min_cluster_size}, "
                      f"retained {scored['retained_fraction']:.0%}, {pretty}")
    return candidates


def select_partition(
    candidates: List[Dict[str, Any]],
    min_cluster_size: int = 10,
    min_retained: float = 0.70,
    primary_space: str = "esm2",
) -> Dict[str, Any]:
    """Pick the best partition subject to size and coverage constraints.

    Silhouette is maximised only *within* the feasible set. Maximising it directly
    would favour many tight clusters - exactly the fragmentation this replaces.
    """
    feasible = [
        c for c in candidates
        if c["retained_fraction"] >= min_retained and c["clusters_at_or_above_min"] >= 2
    ]
    pool = feasible or candidates
    relaxed = not feasible

    def silhouette_of(candidate: Dict[str, Any]) -> float:
        # The primary space first; otherwise the first other space that has a value.
        sil = candidate["silhouette"]
        value = sil.get(primary_space)
        if value is None:
            value = next((v for k, v in sil.items() if k != primary_space and v is not None), None)
        return value if value is not None else -1.0

    if not pool:
        return {"selected": None, "reason": "no candidate partitions were produced",
                "relaxed_constraints": relaxed, "candidates": []}

    if feasible:
        # Inside the feasible set every candidate already keeps enough data together,
        # so silhouette is a safe thing to maximise. Ties - including the case where
        # no silhouette could be computed at all, e.g. a sweep on a fresh run before
        # any embeddings exist - fall back to retention, then usable cluster count,
        # rather than to whichever candidate happened to be first in the grid.
        best = max(pool, key=lambda c: (silhouette_of(c), c["retained_fraction"],
                                        c["clusters_at_or_above_min"]))
    else:
        # Nothing met the constraints. Maximising silhouette here would be actively
        # harmful: the most fragmented partition scores best precisely because
        # near-singleton clusters are trivially well separated. Prefer the partition
        # that keeps the most data in usable clusters, and break ties on silhouette.
        best = max(pool, key=lambda c: (c["retained_fraction"], silhouette_of(c)))

    # Report where the two spaces disagree rather than hiding it behind one number.
    disagreement = None
    other_spaces = [s for s in best["silhouette"] if s != primary_space]
    for space in other_spaces:
        best_other = max(pool, key=lambda c: (c["silhouette"].get(space) if c["silhouette"].get(space) is not None else -1.0))
        if (best_other.get("threshold"), best_other.get("coverage"), best_other.get("tool")) != \
           (best.get("threshold"), best.get("coverage"), best.get("tool")):
            disagreement = {
                "space": space,
                "primary_choice": {"tool": best["tool"], "threshold": best["threshold"],
                                   "threshold_kind": best.get("threshold_kind"),
                                   "coverage": best["coverage"]},
                "alternative_choice": {"tool": best_other["tool"], "threshold": best_other["threshold"],
                                       "threshold_kind": best_other.get("threshold_kind"),
                                       "coverage": best_other["coverage"]},
            }

    recommendation = None
    if relaxed:
        best_retention = max(c["retained_fraction"] for c in candidates)
        recommendation = (
            f"No threshold kept {min_retained:.0%} of taxa in clusters of {min_cluster_size}+ "
            f"(best achievable was {best_retention:.0%}). This means the input does not "
            f"contain groups of that size at this similarity level - the sequences are too "
            f"divergent for the requested granularity. Consider lowering --min-cluster-size, "
            f"clustering on structure (foldseek) instead of sequence, or accepting that the "
            f"set has no meaningful subfamily structure."
        )

    return {
        "selected": best,
        "relaxed_constraints": relaxed,
        "reason": (
            f"constraints unmet; selected the partition retaining the most taxa "
            f"({best['retained_fraction']:.0%}) rather than the best silhouette, which would "
            f"have chosen the most fragmented partition"
        ) if relaxed else (
            f"best {primary_space} silhouette among partitions retaining "
            f">= {min_retained:.0%} of taxa in clusters of {min_cluster_size}+"
        ),
        "recommendation": recommendation,
        "disagreement": disagreement,
        "alternative": suggest_alternative(candidates, best, min_cluster_size=min_cluster_size),
        "constraints": {"min_cluster_size": min_cluster_size, "min_retained": min_retained},
    }


def _flags_for(candidate: Dict[str, Any]) -> str:
    """Render the CLI flags that reproduce a candidate's operating point."""
    kind = candidate.get("threshold_kind")
    coverage = candidate["coverage"]
    threshold = candidate["threshold"]
    if kind == "tmscore":
        return f"--cluster-tmscore {threshold} --cluster-coverage {coverage}"
    if kind == "min_seq_id":
        return f"--cluster-min-seq-id {threshold} --cluster-coverage {coverage}"
    return f"--cluster-evalue {threshold} --cluster-coverage {coverage}"


def suggest_alternative(
    candidates: List[Dict[str, Any]],
    selected: Optional[Dict[str, Any]],
    min_retention_gain: float = SUGGESTION_MIN_RETENTION_GAIN,
    min_cluster_size: int = 10,
) -> Optional[Dict[str, Any]]:
    """Offer a more inclusive operating point than the selected one, if a good one exists.

    Selection maximises silhouette, which favours cleanly separated groups. That is
    the right default, but it is not the only reasonable goal: if the aim is to get
    as much of the data as possible into alignable groups, a looser threshold can be
    a better trade. Rather than making that choice silently, surface the option.

    Returns the candidate retaining the most taxa, provided it beats the selection by
    ``min_retention_gain`` and still yields a comparable number of usable clusters.
    """
    if not selected or not candidates:
        return None

    baseline = selected["retained_fraction"]
    usable_floor = max(2, int(selected["clusters_at_or_above_min"] * 0.75))
    viable = [
        c for c in candidates
        if c["retained_fraction"] >= baseline + min_retention_gain
        and c["clusters_at_or_above_min"] >= usable_floor
        and c is not selected
    ]
    if not viable:
        return None

    # Among the most inclusive options, a percentage point or two of retention is
    # noise, whereas the number of usable clusters is not: a setting that keeps
    # essentially the same taxa but resolves them into more alignable groups is the
    # better suggestion. So take everything within a small band of the best
    # retention, then prefer the one yielding the most clusters of usable size.
    best_retention = max(c["retained_fraction"] for c in viable)
    band = [c for c in viable if c["retained_fraction"] >= best_retention - RETENTION_TIE_BAND]
    best = max(band, key=lambda c: (c["clusters_at_or_above_min"], c["retained_fraction"]))
    return {
        "tool": best["tool"],
        "threshold": best["threshold"],
        "threshold_kind": best.get("threshold_kind"),
        "coverage": best["coverage"],
        "n_clusters": best["n_clusters"],
        "clusters_at_or_above_min": best["clusters_at_or_above_min"],
        "retained_fraction": best["retained_fraction"],
        "silhouette": best["silhouette"],
        "rationale": (
            f"retains {best['retained_fraction']:.0%} of taxa vs {baseline:.0%} "
            f"(+{(best['retained_fraction'] - baseline) * 100:.0f} points) in "
            f"{best['clusters_at_or_above_min']} clusters of {min_cluster_size}+ taxa. "
            f"Lower silhouette, so groups are less cleanly separated - prefer this when the "
            f"goal is to align as much of the set as possible rather than to find the "
            f"crispest subfamilies."
        ),
        "how_to_use": _flags_for(best),
    }


def summarise_sweep(candidates: List[Dict[str, Any]]) -> List[Dict[str, Any]]:
    """Strip cluster membership from a sweep so it can be stored in the manifest."""
    summary = []
    for candidate in candidates:
        row = {k: v for k, v in candidate.items() if k != "clusters"}
        summary.append(row)
    return summary


def save_partition(clusters: Dict[str, List[str]], path: str) -> str:
    os.makedirs(os.path.dirname(os.path.abspath(path)) or ".", exist_ok=True)
    with open(path, "w", encoding="utf-8") as handle:
        json.dump({"clusters": {k: sorted(v) for k, v in clusters.items()}}, handle, indent=2)
    return path


# ------------------------------------------------------- end-to-end partitioning


def _cluster_tag(threshold_kind: str, threshold: float) -> str:
    """Short directory suffix encoding how a partition was produced."""
    return {"evalue": f"e{threshold:g}", "tmscore": f"tm{threshold:g}",
            "min_seq_id": f"id{threshold:g}"}.get(threshold_kind, f"t{threshold:g}")


def partition_structures_by_similarity(
    pdb_dir: str,
    output_dir: str,
    mode: str = "structural",
    evalue: float = DEFAULT_EVALUE,
    coverage: float = DEFAULT_CLUSTER_COVERAGE,
    min_seq_id: float = 0.5,
    tmscore: Optional[float] = None,
    alignment_type: int = 2,
    sweep: bool = False,
    min_cluster_size: int = 10,
    min_retained: float = 0.70,
    embeddings_npz: Optional[str] = None,
    min_length: int = 50,
    min_length_frac: float = 0.4,
    recurse: bool = False,
    recurse_max_depth: int = 2,
    recurse_max_gap: float = RECURSE_MAX_GAP_FRACTION,
    recurse_max_expansion: float = RECURSE_MAX_EXPANSION,
    recurse_min_kept: float = RECURSE_MIN_KEPT,
    recurse_on: str = "3di",
    mafft_leavegappyregion: bool = False,
    aligner: str = "foldmason",
    foldmason_bin: Optional[str] = None,
    mafft_bin: Optional[str] = None,
    mafft_matrix: str = "matrices/mat3di.out",
    threads: int = 0,
    manifest: Any = None,
) -> Dict[str, Any]:
    """Partition structures by measured similarity, then align each cluster.

    Replaces the greedy coverage partitioner, which clustered on mutual non-gap
    columns of a master alignment - a circular measure, since that alignment is
    itself poor when the inputs are heterogeneous.

    Emits the same ``multi_alignment_summary.json`` contract as before so that the
    tree stage and the web viewer keep working unchanged.
    """
    from viral_phylo.alignment import align_structures, extract_unaligned_sequences

    os.makedirs(output_dir, exist_ok=True)
    work = os.path.join(output_dir, "_clustering")
    os.makedirs(work, exist_ok=True)

    structure_files: List[str] = []
    for pattern in ("*.pdb", "*.cif", "*.mmcif"):
        structure_files.extend(glob.glob(os.path.join(pdb_dir, pattern)))
    structure_files = sorted(set(structure_files))
    if len(structure_files) < 2:
        raise ValueError(f"Need at least 2 structures to partition, found {len(structure_files)} in '{pdb_dir}'.")
    taxon_to_file = {os.path.splitext(os.path.basename(p))[0]: p for p in structure_files}

    # 1. Raw, genuinely unaligned sequences (FoldMason is used only as an extractor).
    print(f"[Cluster] Extracting sequences from {len(structure_files)} structures...")
    unaligned_aa, _unaligned_3di = extract_unaligned_sequences(
        structure_files, work, foldmason_bin=foldmason_bin, threads=threads)
    sequences = parse_fasta(unaligned_aa)

    # 2. QC. Short fragments cannot be meaningfully aligned or clustered and are
    #    what stretch a global alignment into a staircase of gaps.
    kept, excluded = qc_filter_sequences(sequences, min_length=min_length,
                                         min_length_frac=min_length_frac)
    print(f"[Cluster] QC: {len(sequences)} sequences -> {len(kept)} kept, {len(excluded)} excluded "
          f"(min length {min_length}, or {min_length_frac:g} x median)")
    if manifest is not None:
        manifest.record_qc(excluded, retained=len(kept))

    qc_fasta = write_fasta(kept, os.path.join(work, "qc_sequences.fasta"))
    qc_structures = os.path.join(work, "qc_structures")
    os.makedirs(qc_structures, exist_ok=True)
    for taxon in kept:
        source = taxon_to_file.get(taxon)
        if not source:
            continue
        link = os.path.join(qc_structures, os.path.basename(source))
        if not os.path.exists(link):
            try:
                os.symlink(os.path.abspath(source), link)
            except OSError:
                shutil.copy2(source, link)

    taxa = [t for t in kept if t in taxon_to_file]

    # 3. Cluster, either at the configured operating point or across a sweep.
    use_structural = mode in ("structural", "both")
    threshold_kind = "evalue" if (use_structural and alignment_type != 1) else (
        "tmscore" if use_structural else "min_seq_id")

    def _structural(run_dir, threshold, cov):
        return foldseek_cluster(qc_structures, run_dir, evalue=threshold, coverage=cov,
                                alignment_type=alignment_type, tmscore_threshold=tmscore,
                                threads=threads)

    def _sequence(run_dir, threshold, cov):
        return mmseqs_cluster(qc_fasta, run_dir, threshold, cov, threads=threads)

    cluster_fn = _structural if use_structural else _sequence
    selection: Dict[str, Any] = {}

    if sweep:
        spaces: Dict[str, np.ndarray] = {}
        if embeddings_npz and os.path.isfile(embeddings_npz):
            emb_taxa, emb_dist = embedding_distance_matrix(embeddings_npz)
            index = {t: i for i, t in enumerate(emb_taxa)}
            usable = [t for t in taxa if t in index]
            if len(usable) >= 3:
                order = np.array([index[t] for t in usable])
                taxa = usable
                spaces["esm2"] = emb_dist[np.ix_(order, order)]
        # Second, independent space from the data itself: structural distance
        # (1 - TM-score) for structural clustering, sequence-identity distance for
        # sequence clustering. Computed once per sweep, not per candidate.
        second = "structural" if use_structural else "identity"
        try:
            if use_structural:
                sp_taxa, sp_dist = identity_distance_matrix(
                    qc_structures, taxa, os.path.join(work, "space_search"), tool="foldseek",
                    threads=threads, score="alntmscore", evalue="10", alignment_type=2,
                    max_seqs=max(2000, len(taxa)))
            else:
                sp_taxa, sp_dist = identity_distance_matrix(
                    qc_fasta, taxa, os.path.join(work, "space_search"), tool="mmseqs",
                    threads=threads, score="fident", max_seqs=max(2000, len(taxa)))
            spaces[second] = sp_dist
        except Exception as exc:
            print(f"[Cluster] [!] Notice: could not build the {second} distance space ({exc}); "
                  "scoring without it.")
        if "esm2" not in spaces:
            print("[Cluster] [!] Notice: no ESM-2 embeddings for this run (add --embed to score "
                  f"in the ESM-2 space too); scoring by {second} silhouette only." if spaces else
                  "[Cluster] [!] Notice: no distance space available; the sweep will be scored "
                  "on retention and cluster size only.")

        thresholds = (DEFAULT_TMSCORE_THRESHOLDS if threshold_kind == "tmscore"
                      else DEFAULT_EVALUES if threshold_kind == "evalue" else DEFAULT_MIN_SEQ_IDS)
        candidates = sweep_clusterings(
            taxa, spaces, cluster_fn, work, thresholds=thresholds,
            coverages=DEFAULT_COVERAGES, min_cluster_size=min_cluster_size,
            label="foldseek" if use_structural else "mmseqs", threshold_kind=threshold_kind,
        )
        primary = "esm2" if "esm2" in spaces else second
        selection = select_partition(candidates, min_cluster_size=min_cluster_size,
                                     min_retained=min_retained, primary_space=primary)
        selection["spaces"] = {"used": sorted(spaces), "primary": primary if spaces else None}
        chosen = selection.get("selected")
        if not chosen:
            raise RuntimeError("clustering sweep produced no usable partition")
        clusters = chosen["clusters"]
        threshold_used, coverage_used = chosen["threshold"], chosen["coverage"]
        selection["sweep"] = summarise_sweep(candidates)
        sil_text = ", ".join(f"{k} silhouette {v:+.3f}" for k, v in chosen["silhouette"].items() if v is not None)
        print(f"[Cluster] Selected {threshold_kind}={threshold_used} coverage={coverage_used}: "
              f"{chosen['n_clusters']} clusters, {chosen['clusters_at_or_above_min']} with "
              f"{min_cluster_size}+ taxa, {chosen['retained_fraction']:.0%} of taxa retained"
              + (f" ({sil_text})." if sil_text else "."))
        if selection.get("disagreement"):
            dis = selection["disagreement"]
            alt = dis["alternative_choice"]
            print(f"[Cluster] [!] Notice: the {dis['space']} space would have chosen "
                  f"{alt['threshold_kind']}={alt['threshold']} coverage={alt['coverage']} instead; "
                  "both are recorded in the manifest.")
        if selection.get("alternative"):
            alt = selection["alternative"]
            print(f"[Cluster] A more inclusive option is available: {alt['how_to_use']}")
            print(f"[Cluster]   {alt['rationale']}")
    else:
        threshold_used = (tmscore if threshold_kind == "tmscore" else
                          evalue if threshold_kind == "evalue" else min_seq_id)
        coverage_used = coverage
        print(f"[Cluster] Clustering at {threshold_kind}={threshold_used} coverage={coverage_used} "
              f"({'foldseek structural' if use_structural else 'mmseqs sequence'})...")
        clusters = cluster_fn(os.path.join(work, "operating_point"), threshold_used, coverage_used)
        labels = clusters_to_labels(clusters, taxa)
        shape = partition_shape(labels, min_cluster_size)
        selection = {"selected": {**shape, "tool": "foldseek" if use_structural else "mmseqs",
                                  "threshold": threshold_used, "threshold_kind": threshold_kind,
                                  "coverage": coverage_used}}
        print(f"[Cluster] {shape['n_clusters']} clusters, {shape['clusters_at_or_above_min']} with "
              f"{min_cluster_size}+ taxa, {shape['retained_fraction']:.0%} of taxa retained.")

    # An e-value shifts with database size, so record what it was measured against.
    selection["operating_point"] = {
        "mode": mode, "threshold_kind": threshold_kind, "threshold": threshold_used,
        "coverage": coverage_used, "alignment_type": alignment_type,
        "taxa_clustered": len(taxa),
        "note": ("e-values depend on database size and are not comparable across datasets "
                 "of different sizes" if threshold_kind == "evalue" else None),
    }

    # 4. Materialise each cluster and align it. Optionally split any cluster whose
    #    own alignment is still poor, at a stricter threshold - the symptom being
    #    corrected is a gappy alignment, so that is what drives the decision rather
    #    than an arbitrary size cap.
    def _materialise(name, members):
        c_dir = os.path.join(output_dir, name)
        c_structs = os.path.join(c_dir, "structures")
        os.makedirs(c_structs, exist_ok=True)
        for taxon in members:
            source = taxon_to_file.get(taxon)
            if not source:
                continue
            link = os.path.join(c_structs, os.path.basename(source))
            if not os.path.exists(link):
                try:
                    os.symlink(os.path.abspath(source), link)
                except OSError:
                    shutil.copy2(source, link)
        return c_dir, c_structs

    def _align_and_measure(c_dir, c_structs, members):
        if len(members) < 2:
            return None, None, None
        try:
            align_structures(c_structs, c_dir, foldmason_bin=foldmason_bin, aligner=aligner,
                             mafft_bin=mafft_bin, mafft_matrix=mafft_matrix,
                             mafft_leavegappyregion=mafft_leavegappyregion, threads=threads)
        except Exception as exc:
            print(f"  [!] Notice: alignment failed for {os.path.basename(c_dir)}: {exc}")
            return None, None, None
        aln_3di = os.path.join(c_dir, "foldmason.fasta_3di.fa")
        aln_aa = os.path.join(c_dir, "foldmason.fasta_aa.fa")
        # Record which aligner really wrote these files: under MAFFT the
        # foldmason.fasta_*.fa names are compatibility copies of mafft.fasta_*.fa.
        try:
            from viral_phylo.manifest import write_alignment_info
            write_alignment_info(
                c_dir, aligner=aligner,
                produced={os.path.basename(aln_3di): aligner, os.path.basename(aln_aa): aligner},
                mafft_matrix=mafft_matrix if aligner == "mafft" else None,
                tool_paths={"mafft": mafft_bin},
                extra={"mafft_leavegappyregion": mafft_leavegappyregion if aligner == "mafft" else None,
                       "taxa": len(members)})
        except Exception as exc:
            print(f"  [!] Notice: could not write alignment_info.json for {os.path.basename(c_dir)}: {exc}")
        # Report quality with the 3Di alignment at the top level: structural
        # alignment is the purpose of this tool, so by default splits are judged on
        # it (see cluster_needs_split). The AA alignment's numbers are kept under
        # "aa". Under --aligner mafft the AA alignment is built separately from
        # sequence alone, is typically ~5 points gappier and is not column-paired
        # with the 3Di one, so it is not used for the decision unless asked for.
        stats = alignment_gap_fraction(parse_fasta(aln_3di)) if os.path.isfile(aln_3di) else None
        stats_aa = alignment_gap_fraction(parse_fasta(aln_aa)) if os.path.isfile(aln_aa) else None
        if stats is not None:
            stats = {**stats, "measured_on": "3di", "aa": stats_aa}
        return aln_3di, aln_aa, stats

    # (representative, members, threshold, depth, parent). Names are never reused: a
    # split parent's index is retired rather than handed to its first child, so every
    # name in the recursion log refers to exactly one cluster.
    pending = [(rep, members, threshold_used, 0, None) for rep, members in
               sorted(clusters.items(), key=lambda kv: (-len(kv[1]), kv[0]))]
    cluster_results: List[Dict[str, Any]] = []
    recursion_log: List[Dict[str, Any]] = []
    aligned = skipped_small = 0
    index = 0

    while pending:
        representative, members, threshold_at, depth, parent = pending.pop(0)
        index += 1
        name = f"cluster_{index}_{_cluster_tag(threshold_kind, threshold_at)}"
        c_dir, c_structs = _materialise(name, members)
        aln_3di, aln_aa, stats = _align_and_measure(c_dir, c_structs, members)

        if aln_3di:
            aligned += 1
        elif len(members) < 2:
            skipped_small += 1

        # Should this cluster be split further?
        not_divisible = False
        if (recurse and stats and depth < recurse_max_depth
                and len(members) >= 2 * max(2, min_cluster_size)):
            needs_split, why = cluster_needs_split(
                stats, recurse_on=recurse_on,
                max_gap_fraction=recurse_max_gap, max_expansion=recurse_max_expansion)
            if needs_split:
                attempts = []
                accepted = None
                for candidate in recursion_schedule(threshold_kind, threshold_at):
                    sub_dir = os.path.join(work, f"recurse_{name}_{_cluster_tag(threshold_kind, candidate)}")
                    try:
                        if use_structural:
                            sub = foldseek_cluster(
                                c_structs, sub_dir, evalue=candidate, coverage=coverage_used,
                                alignment_type=alignment_type,
                                tmscore_threshold=candidate if threshold_kind == "tmscore" else tmscore,
                                threads=threads)
                        else:
                            sub_fasta = write_fasta({t: kept[t] for t in members if t in kept},
                                                    os.path.join(sub_dir, "sub.fasta"))
                            sub = mmseqs_cluster(sub_fasta, sub_dir, candidate, coverage_used, threads=threads)
                    except Exception as exc:
                        print(f"  [!] Notice: re-clustering {name} at {threshold_kind}={candidate:g} failed: {exc}")
                        continue
                    kept_frac = fraction_kept(sub, len(members))
                    attempts.append({"threshold": candidate, "sub_clusters": len(sub),
                                     "kept_fraction": round(kept_frac, 3)})
                    if len(sub) > 1 and kept_frac >= recurse_min_kept:
                        accepted = (candidate, sub, kept_frac)
                        break

                if accepted is None:
                    # Either nothing divided it, or every division shattered it. Keep the
                    # cluster whole, with a known-poor alignment, and say so.
                    tried = ", ".join(f"{a['threshold']:g}->{a['sub_clusters']} ({a['kept_fraction']:.0%} kept)"
                                      for a in attempts) or "none"
                    print(f"[Cluster] [!] Notice: {name} ({len(members)} taxa) has a poor alignment ({why}) "
                          f"but no split kept >= {recurse_min_kept:.0%} of its taxa in clusters of "
                          f"{TREE_MIN_TAXA}+; keeping it whole. Tried: {tried}")
                    recursion_log.append({
                        "cluster": name, "taxa": len(members), "depth": depth, "reason": why,
                        "from_threshold": threshold_at, "attempts": attempts,
                        "outcome": "not_divisible",
                    })
                    not_divisible = True
                else:
                    tighter, sub, kept_frac = accepted
                    print(f"[Cluster] Splitting {name} ({len(members)} taxa): {why} "
                          f"-> {len(sub)} sub-clusters at {threshold_kind}={tighter:g}, "
                          f"{kept_frac:.0%} of taxa kept in clusters of {TREE_MIN_TAXA}+")
                    recursion_log.append({
                        "cluster": name, "taxa": len(members), "depth": depth,
                        "reason": why, "from_threshold": threshold_at,
                        "to_threshold": tighter, "sub_clusters": len(sub),
                        "kept_fraction": round(kept_frac, 3), "attempts": attempts,
                        "outcome": "split",
                    })
                    shutil.rmtree(c_dir, ignore_errors=True)
                    if aln_3di:
                        aligned -= 1
                    pending = [(r, m, tighter, depth + 1, name) for r, m in
                               sorted(sub.items(), key=lambda kv: (-len(kv[1]), kv[0]))] + pending
                    continue

        over_3di = bool(stats and alignment_needs_splitting(stats, recurse_max_gap, recurse_max_expansion)[0])
        meta = {
            "cluster_id": index,
            "name": name,
            "taxa": sorted(members),
            "taxa_count": len(members),
            "directory": c_dir,
            "aln_3di": aln_3di,
            "aln_aa": aln_aa,
            "mean_coverage": round(1.0 - stats["gap_fraction"], 4) if stats else 0.0,
            "alignment_stats": stats,
            "representative": representative,
            "aligner": aligner,
            "threshold": threshold_at,
            "recursion_depth": depth,
            "parent": parent,
            "tree_eligible": len(members) >= 4,
            "aa_over_threshold": bool(stats and stats.get("aa") and alignment_needs_splitting(
                stats["aa"], recurse_max_gap, recurse_max_expansion)[0]),
            # Clusters kept although their 3Di alignment is over the recursion thresholds,
            # and why - so a tree built on a gappy alignment is never unexplained.
            "over_threshold_3di": over_3di,
            "kept_over_threshold_because": (None if not over_3di else
                                            "recursion disabled" if not recurse else
                                            "no split kept enough taxa in clusters of 4+" if not_divisible else
                                            "recursion depth limit reached" if depth >= recurse_max_depth else
                                            "too small to split"),
        }
        cluster_results.append(meta)
        with open(os.path.join(c_dir, "cluster_info.json"), "w", encoding="utf-8") as handle:
            json.dump(meta, handle, indent=2, default=str)

    print(f"[Cluster] Aligned {aligned} cluster(s); {skipped_small} singleton(s) left unaligned.")

    summary = {
        "dataset": pdb_dir,
        "method": "similarity",
        "mode": mode,
        "operating_point": selection["operating_point"],
        "qc": {"input": len(sequences), "kept": len(kept), "excluded": len(excluded)},
        "min_coverage": coverage_used,
        "total_taxa": len(taxa),
        "num_clusters": len(cluster_results),
        "multi_taxa_clusters": sum(1 for c in cluster_results if c["taxa_count"] >= min_cluster_size),
        "tree_eligible_clusters": sum(1 for c in cluster_results if c["tree_eligible"]),
        # Tree-eligible clusters whose AA alignment is over the recursion thresholds.
        # Only 3Di drives splitting by default, so these are where an AA tree
        # (--tree-type both) would rest on a gappier alignment than the 3Di tree.
        "aa_over_threshold_clusters": [c["name"] for c in cluster_results
                                       if c["tree_eligible"] and c.get("aa_over_threshold")],
        "over_threshold_3di_clusters": {c["name"]: c["kept_over_threshold_because"] for c in cluster_results
                                        if c["tree_eligible"] and c.get("over_threshold_3di")},
        "clusters": cluster_results,
        "recursion": {"enabled": recurse, "max_depth": recurse_max_depth,
                      "max_gap": recurse_max_gap, "max_expansion": recurse_max_expansion,
                      "min_kept": recurse_min_kept, "judged_on": recurse_on,
                      "splits": recursion_log},
        "selection": {k: v for k, v in selection.items() if k != "selected"},
        "master_alignment": {"aln_3di": None, "aln_aa": None},
    }

    summary_path = os.path.join(output_dir, "multi_alignment_summary.json")
    with open(summary_path, "w", encoding="utf-8") as handle:
        json.dump(summary, handle, indent=2, default=str)
    print(f"[Cluster] Saved partition summary to '{summary_path}'.\n")

    if manifest is not None:
        manifest.record_clustering({k: v for k, v in summary.items() if k != "clusters"})

    return summary
