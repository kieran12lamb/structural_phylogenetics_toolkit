"""Multiple sequence and structure alignment (FoldMason / MAFFT), coverage calculation, and partitioning."""

import os
import glob
import json
import shutil
import subprocess
from pathlib import Path
from typing import Dict, List, Optional, Union, Any

from viral_phylo.binaries import find_foldmason_bin, find_mafft_bin
from viral_phylo.matrices import ensure_3di_matrix

def parse_alignment_fasta(fasta_path: str) -> dict:
    """Parse a FASTA file into a dictionary {taxon_id: sequence_string}."""
    seqs = {}
    cur = None
    if not os.path.isfile(fasta_path):
        return seqs
    with open(fasta_path, "r", encoding="utf-8", errors="ignore") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                cur = line[1:].split()[0]
                seqs[cur] = []
            elif cur:
                seqs[cur].append(line)
    return {k: "".join(v) for k, v in seqs.items()}


def write_alignment_fasta(aln_dict: dict, out_path: str):
    """Write an alignment dictionary to a FASTA file."""
    os.makedirs(os.path.dirname(os.path.abspath(out_path)), exist_ok=True)
    with open(out_path, "w", encoding="utf-8") as f:
        for tid, seq in aln_dict.items():
            f.write(f">{tid}\n{seq}\n")


def strip_all_gap_columns_dict(alignment_dict: dict) -> dict:
    """Strip columns where all sequences have gaps."""
    if not alignment_dict:
        return {}
    first_seq = next(iter(alignment_dict.values()))
    seq_len = len(first_seq)
    keep_indices = [c for c in range(seq_len) if any(seq[c] != "-" for seq in alignment_dict.values())]
    return {tid: "".join(seq[i] for i in keep_indices) for tid, seq in alignment_dict.items()}


def compute_alignment_coverage(aln_input, mode: str = "both") -> dict:
    """Compute per-sequence coverage, column occupancy, and pairwise mutual coverage.
    
    Coverage Definitions:
    - sequence_coverage: (non-gap residues in sequence) / (total alignment length)
    - domain_coverage: (non-gap residues in sequence) / (maximum sequence length in alignment)
    - column_occupancy: fraction of sequences with non-gap residues at each column
    - pairwise_coverage: mutual non-gap columns / min(lenA, lenB) and / max(lenA, lenB)
    """
    if isinstance(aln_input, (str, Path)):
        aln = parse_alignment_fasta(str(aln_input))
    else:
        aln = dict(aln_input)

    if not aln:
        return {"taxa": {}, "aln_length": 0, "max_seq_len": 0, "column_occupancy": [], "pairwise": {}}

    taxa = list(aln.keys())
    first_seq = next(iter(aln.values()))
    aln_len = len(first_seq)
    taxa_lens = {t: sum(1 for c in aln[t] if c != "-") for t in taxa}
    max_seq_len = max(taxa_lens.values()) if taxa_lens else 1

    taxa_stats = {}
    for t in taxa:
        ng = taxa_lens[t]
        seq_cov = round(ng / aln_len, 4) if aln_len > 0 else 0.0
        dom_cov = round(ng / max_seq_len, 4) if max_seq_len > 0 else 0.0
        taxa_stats[t] = {
            "non_gaps": ng,
            "aln_length": aln_len,
            "sequence_coverage": seq_cov,
            "domain_coverage": dom_cov,
            "coverage": seq_cov,
        }

    # Column occupancy
    col_occupancy = []
    num_taxa = len(taxa)
    for c in range(aln_len):
        non_gaps_col = sum(1 for t in taxa if aln[t][c] != "-")
        col_occupancy.append(round(non_gaps_col / num_taxa, 4))

    # Pairwise mutual coverage
    pairwise = {}
    for i in range(len(taxa)):
        for j in range(i + 1, len(taxa)):
            t1, t2 = taxa[i], taxa[j]
            s1, s2 = aln[t1], aln[t2]
            l1, l2 = taxa_lens[t1], taxa_lens[t2]
            mutual = sum(1 for c1, c2 in zip(s1, s2) if c1 != "-" and c2 != "-")
            cov_min = round(mutual / min(l1, l2), 4) if min(l1, l2) > 0 else 0.0
            cov_max = round(mutual / max(l1, l2), 4) if max(l1, l2) > 0 else 0.0
            pairwise[f"{t1}_vs_{t2}"] = {
                "mutual": mutual,
                "cov_min": cov_min,
                "cov_max": cov_max,
            }

    return {
        "taxa": taxa_stats,
        "aln_length": aln_len,
        "max_seq_len": max_seq_len,
        "column_occupancy": col_occupancy,
        "pairwise": pairwise,
    }


def filter_alignment_by_coverage(
    aln_3di_path: str,
    aln_aa_path: str = None,
    min_coverage: float = 0.70,
    output_dir: str = None,
    use_domain_cov: bool = False,
) -> dict:
    """Filter an alignment to retain only sequences with coverage >= min_coverage."""
    aln_3di = parse_alignment_fasta(aln_3di_path)
    aln_aa = parse_alignment_fasta(aln_aa_path) if aln_aa_path and os.path.isfile(aln_aa_path) else None

    cov_info = compute_alignment_coverage(aln_3di)
    metric_key = "domain_coverage" if use_domain_cov else "sequence_coverage"

    passed_taxa = [t for t, st in cov_info["taxa"].items() if st[metric_key] >= min_coverage]
    removed_taxa = [t for t, st in cov_info["taxa"].items() if st[metric_key] < min_coverage]

    filtered_3di = {t: aln_3di[t] for t in passed_taxa}
    filtered_3di = strip_all_gap_columns_dict(filtered_3di)

    filtered_aa = None
    if aln_aa:
        filtered_aa = {t: aln_aa[t] for t in passed_taxa if t in aln_aa}
        filtered_aa = strip_all_gap_columns_dict(filtered_aa)

    out_3di = None
    out_aa = None
    if output_dir:
        os.makedirs(output_dir, exist_ok=True)
        cov_pct = int(min_coverage * 100)
        out_3di = os.path.join(output_dir, f"foldmason_cov{cov_pct}_filtered.fasta_3di.fa")
        write_alignment_fasta(filtered_3di, out_3di)
        if filtered_aa:
            out_aa = os.path.join(output_dir, f"foldmason_cov{cov_pct}_filtered.fasta_aa.fa")
            write_alignment_fasta(filtered_aa, out_aa)

        if "mafft" in os.path.basename(aln_3di_path).lower():
            m_3di = os.path.join(output_dir, f"mafft_cov{cov_pct}_filtered.fasta_3di.fa")
            shutil.copyfile(out_3di, m_3di)
            if out_aa and os.path.isfile(out_aa):
                m_aa = os.path.join(output_dir, f"mafft_cov{cov_pct}_filtered.fasta_aa.fa")
                shutil.copyfile(out_aa, m_aa)

    return {
        "passed_taxa": passed_taxa,
        "removed_taxa": removed_taxa,
        "filtered_3di": filtered_3di,
        "filtered_aa": filtered_aa,
        "out_3di": out_3di,
        "out_aa": out_aa,
        "min_coverage": min_coverage,
    }


def partition_structures_by_coverage(
    pdb_dir: str,
    output_dir: str,
    min_coverage: float = 0.70,
    aligner: str = "foldmason",
    foldmason_bin: str = None,
    mafft_bin: str = None,
    mafft_matrix: str = "matrices/mat3di.out",
    min_cluster_size: int = 2,
) -> dict:
    """Partition a single structure dataset into multiple sub-alignments where members share >= min_coverage.

    1. Runs initial alignment on all structures to compute pairwise mutual structural coverage.
    2. Clusters structures greedily (longest representative first, adding members with >= min_coverage mutual overlap).
    3. For each cluster with >= min_cluster_size structures, aligns the cluster into its own MSTA.
    4. Emits cluster metadata and summary JSON.
    """
    os.makedirs(output_dir, exist_ok=True)

    extensions = ("*.pdb", "*.cif", "*.mmcif")
    all_pdb_files = []
    for ext in extensions:
        all_pdb_files.extend(glob.glob(os.path.join(pdb_dir, ext)))
    all_pdb_files = sorted(set(all_pdb_files))

    if len(all_pdb_files) < 2:
        raise ValueError(f"Need at least 2 structure files to partition, found {len(all_pdb_files)} in '{pdb_dir}'.")

    # Map taxon_id -> file path
    taxon_to_file = {}
    for p in all_pdb_files:
        tid = os.path.splitext(os.path.basename(p))[0]
        taxon_to_file[tid] = p

    # Step 1: Master alignment on all structures
    master_aln_dir = os.path.join(output_dir, "master_alignment")
    os.makedirs(master_aln_dir, exist_ok=True)
    master_3di = os.path.join(master_aln_dir, "foldmason.fasta_3di.fa")

    if not (os.path.isfile(master_3di) and os.path.getsize(master_3di) > 0):
        print(f"[Multi-Alignment] Running initial master alignment across all {len(all_pdb_files)} structures...")
        align_structures(
            pdb_dir,
            master_aln_dir,
            foldmason_bin=foldmason_bin,
            aligner=aligner,
            mafft_bin=mafft_bin,
            mafft_matrix=mafft_matrix,
        )

    # Step 2: Compute pairwise coverage from master alignment
    aln = parse_alignment_fasta(master_3di)
    taxa = list(aln.keys())
    taxa_lens = {t: sum(1 for c in aln[t] if c != "-") for t in taxa}
    sorted_taxa = sorted(taxa, key=lambda t: taxa_lens[t], reverse=True)

    clusters_taxa = []
    assigned = set()
    for rep in sorted_taxa:
        if rep in assigned:
            continue
        cluster = [rep]
        assigned.add(rep)
        rep_seq = aln[rep]
        rep_len = taxa_lens[rep]
        for other in sorted_taxa:
            if other in assigned:
                continue
            other_seq = aln[other]
            other_len = taxa_lens[other]
            mutual = sum(1 for c1, c2 in zip(rep_seq, other_seq) if c1 != "-" and c2 != "-")
            # Coverage relative to representative (longer) sequence
            cov = mutual / max(rep_len, other_len) if max(rep_len, other_len) > 0 else 0.0
            if cov >= min_coverage:
                cluster.append(other)
                assigned.add(other)
        clusters_taxa.append(cluster)

    cov_pct = int(min_coverage * 100)
    print(f"[Multi-Alignment] Partitioned {len(taxa)} taxa into {len(clusters_taxa)} cluster(s) at {cov_pct}% coverage threshold.")

    # Step 3: Align each cluster with >= min_cluster_size
    cluster_results = []
    for c_idx, c_taxa in enumerate(clusters_taxa, 1):
        c_name = f"cluster_{c_idx}_cov{cov_pct}"
        c_dir = os.path.join(output_dir, c_name)
        c_struct_dir = os.path.join(c_dir, "structures")
        os.makedirs(c_struct_dir, exist_ok=True)

        for t in c_taxa:
            src_p = taxon_to_file.get(t)
            if src_p and os.path.isfile(src_p):
                dst_p = os.path.join(c_struct_dir, os.path.basename(src_p))
                if not os.path.exists(dst_p):
                    try:
                        os.symlink(os.path.abspath(src_p), dst_p)
                    except (OSError, AttributeError):
                        shutil.copy2(src_p, dst_p)

        c_aln_3di = None
        c_aln_aa = None
        if len(c_taxa) >= min_cluster_size:
            print(f"[Multi-Alignment] Aligning Cluster {c_idx} ({len(c_taxa)} taxa: {', '.join(c_taxa[:4])}{'...' if len(c_taxa)>4 else ''}) -> '{c_dir}'...")
            align_structures(
                c_struct_dir,
                c_dir,
                foldmason_bin=foldmason_bin,
                aligner=aligner,
                mafft_bin=mafft_bin,
                mafft_matrix=mafft_matrix,
            )
            c_aln_3di = os.path.join(c_dir, "foldmason.fasta_3di.fa")
            c_aln_aa = os.path.join(c_dir, "foldmason.fasta_aa.fa")

        # Calculate mean mutual coverage for cluster members
        mutual_covs = []
        for i in range(len(c_taxa)):
            for j in range(i + 1, len(c_taxa)):
                t1, t2 = c_taxa[i], c_taxa[j]
                s1, s2 = aln[t1], aln[t2]
                l1, l2 = taxa_lens[t1], taxa_lens[t2]
                mut = sum(1 for c1, c2 in zip(s1, s2) if c1 != "-" and c2 != "-")
                if max(l1, l2) > 0:
                    mutual_covs.append(mut / max(l1, l2))
        mean_cov = round(sum(mutual_covs) / len(mutual_covs), 4) if mutual_covs else 1.0

        c_meta = {
            "cluster_id": c_idx,
            "name": c_name,
            "taxa": c_taxa,
            "taxa_count": len(c_taxa),
            "directory": c_dir,
            "aln_3di": c_aln_3di,
            "aln_aa": c_aln_aa,
            "mean_coverage": mean_cov,
            "representative": c_taxa[0],
        }
        cluster_results.append(c_meta)

        # Save individual cluster metadata
        with open(os.path.join(c_dir, "cluster_info.json"), "w", encoding="utf-8") as f:
            json.dump(c_meta, f, indent=2)

    summary = {
        "dataset": pdb_dir,
        "min_coverage": min_coverage,
        "total_taxa": len(taxa),
        "num_clusters": len(clusters_taxa),
        "multi_taxa_clusters": sum(1 for c in cluster_results if c["taxa_count"] >= min_cluster_size),
        "clusters": cluster_results,
        "master_alignment": {
            "aln_3di": master_3di,
            "aln_aa": os.path.join(master_aln_dir, "foldmason.fasta_aa.fa"),
        }
    }

    summary_path = os.path.join(output_dir, "multi_alignment_summary.json")
    with open(summary_path, "w", encoding="utf-8") as f:
        json.dump(summary, f, indent=2)
    print(f"[Multi-Alignment] Saved multi-alignment summary index to '{summary_path}'.\n")

    return summary


def count_fasta_seqs(aln_path: str) -> int:

    """Count number of sequences in a FASTA file."""
    if not os.path.isfile(aln_path):
        return 0
    with open(aln_path, "r", errors="ignore") as f:
        return sum(1 for line in f if line.startswith(">"))


def align_structures(
    pdb_dir: str,
    output_dir: str,
    foldmason_bin: str = None,
    report_mode: int = None,
    aligner: str = "foldmason",
    mafft_bin: str = None,
    mafft_matrix: str = "matrices/mat3di.out",
    min_coverage: float = None,
    multi_alignment: bool = False,
    filter_coverage: bool = False,
):
    """Run FoldMason or MAFFT multiple sequence/structure alignment on PDB files."""
    if not foldmason_bin:
        foldmason_bin = find_foldmason_bin()

    os.makedirs(output_dir, exist_ok=True)
    tmp_folder = os.path.join(output_dir, "tmpFolder")
    os.makedirs(tmp_folder, exist_ok=True)

    extensions = ("*.pdb", "*.cif", "*.mmcif")
    pdb_files = []
    for ext in extensions:
        pdb_files.extend(glob.glob(os.path.join(pdb_dir, ext)))
    pdb_files = sorted(set(pdb_files))

    if len(pdb_files) < 2:
        raise ValueError(f"Need at least 2 structure files in '{pdb_dir}' to align, found {len(pdb_files)}.")

    aligner_norm = (aligner or "foldmason").lower().strip()

    if aligner_norm == "mafft":
        actual_mafft = find_mafft_bin(mafft_bin)
        actual_mat3di = ensure_3di_matrix(mafft_matrix)

        print(f"[MAFFT] Extracting unaligned 3Di and AA sequences for {len(pdb_files)} structures via FoldMason createdb...")
        sdb_prefix = os.path.join(tmp_folder, "sDB")
        cmd_createdb = [foldmason_bin, "createdb", *pdb_files, sdb_prefix]
        subprocess.run(cmd_createdb, check=True)

        # Prepare header database for sDB_ss conversion
        for ext in ("", ".index", ".dbtype"):
            src = f"{sdb_prefix}_h{ext}"
            dst = f"{sdb_prefix}_ss_h{ext}"
            if os.path.isfile(src) and not os.path.isfile(dst):
                try:
                    shutil.copyfile(src, dst)
                except Exception:
                    pass

        unaligned_aa = os.path.join(tmp_folder, "unaligned_aa.fa")
        unaligned_3di = os.path.join(tmp_folder, "unaligned_3di.fa")

        # Convert to FASTA
        subprocess.run([foldmason_bin, "convert2fasta", sdb_prefix, unaligned_aa], check=True)
        subprocess.run([foldmason_bin, "convert2fasta", f"{sdb_prefix}_ss", unaligned_3di], check=True)

        # Output alignment paths: explicitly save MAFFT alignment files and also maintain foldmason names for downstream compatibility
        mafft_3di = os.path.join(output_dir, "mafft.fasta_3di.fa")
        mafft_aa = os.path.join(output_dir, "mafft.fasta_aa.fa")
        aln_3di = os.path.join(output_dir, "foldmason.fasta_3di.fa")
        aln_aa = os.path.join(output_dir, "foldmason.fasta_aa.fa")

        print(f"[MAFFT] Aligning 3Di structural sequences with substitution matrix '{actual_mat3di}'...")
        with open(mafft_3di, "w") as f_out:
            subprocess.run([actual_mafft, "--aamatrix", actual_mat3di, "--auto", unaligned_3di], stdout=f_out, check=True)

        print(f"[MAFFT] Aligning amino acid sequences...")
        with open(mafft_aa, "w") as f_out:
            subprocess.run([actual_mafft, "--auto", unaligned_aa], stdout=f_out, check=True)

        # Mirror alignments to foldmason.fasta_*.fa for full backward compatibility with downstream pipeline tools
        shutil.copyfile(mafft_3di, aln_3di)
        shutil.copyfile(mafft_aa, aln_aa)

        if not os.path.isfile(mafft_3di) or os.path.getsize(mafft_3di) == 0:
            raise FileNotFoundError(f"Expected MAFFT 3Di alignment '{mafft_3di}' was not generated.")
        print(f"[MAFFT] Alignment complete: generated '{mafft_3di}' and '{mafft_aa}' (also linked as '{aln_3di}' and '{aln_aa}').\n")

    else:
        # Default FoldMason easy-msa
        if report_mode is None:
            report_mode = 0 if len(pdb_files) > 40 else 1

        out_prefix = os.path.join(output_dir, "foldmason.fasta")
        cmd = [foldmason_bin, "easy-msa", *pdb_files, out_prefix, tmp_folder, "--report-mode", str(report_mode)]

        print(f"[FoldMason] Aligning {len(pdb_files)} structures from '{pdb_dir}'...")
        print("  " + " ".join(cmd[:10]) + (f" ... [{len(pdb_files)-10} more files]" if len(pdb_files) > 10 else ""))
        subprocess.run(cmd, check=True)

        aln_3di = os.path.join(output_dir, "foldmason.fasta_3di.fa")
        aln_aa = os.path.join(output_dir, "foldmason.fasta_aa.fa")
        if not os.path.isfile(aln_3di):
            raise FileNotFoundError(f"Expected 3Di alignment '{aln_3di}' was not generated.")
        print(f"[FoldMason] Alignment complete: generated '{aln_3di}'.\n")

    # Multi-Alignment Partitioning or Single Coverage Filtering (supported for both FoldMason and MAFFT)
    if multi_alignment and min_coverage is not None and min_coverage > 0:
        multi_summary = partition_structures_by_coverage(
            pdb_dir,
            output_dir,
            min_coverage=min_coverage,
            aligner=aligner,
            foldmason_bin=foldmason_bin,
            mafft_bin=mafft_bin,
            mafft_matrix=mafft_matrix,
        )
        return aln_3di
    elif filter_coverage and min_coverage is not None and min_coverage > 0:
        filter_res = filter_alignment_by_coverage(aln_3di, aln_aa, min_coverage=min_coverage, output_dir=output_dir)
        print(f"[Coverage Filter] Filtered to {len(filter_res['passed_taxa'])} taxa (>= {int(min_coverage*100)}% coverage). Removed {len(filter_res['removed_taxa'])} taxa.")
        return filter_res['out_3di'] or aln_3di

    return aln_3di
