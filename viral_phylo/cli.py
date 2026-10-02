"""Command Line Interface for the Viral Structural Phylogenetics Toolkit."""

import argparse
import glob
import json
import os
import shutil
import subprocess
import sys
import time
from pathlib import Path
from typing import Optional

from viral_phylo.binaries import (
    find_iqtree_bin,
    find_foldmason_bin,
    find_mafft_bin,
    find_torch_python,
)
from viral_phylo.matrices import ensure_matrix_file, ensure_3di_matrix
from viral_phylo.fetch import fetch_structures, fetch_alphafold_structures
from viral_phylo.alignment import align_structures
from viral_phylo.clustering import (
    alignment_gap_fraction,
    full_set_alignment_is_usable,
    parse_fasta as cluster_parse_fasta,
    partition_structures_by_similarity,
)
from viral_phylo.manifest import RunManifest, write_alignment_info
from viral_phylo.structures import extract_ca_dict, write_ca_bundle
from viral_phylo.tree import build_tree
from viral_phylo.parallel_trees import plan_tree_jobs, run_tree_jobs
from viral_phylo.metadata import parse_metadata
from viral_phylo.web.report import build_tree_report


def pipeline_preflight(args, manifest) -> None:
    """Check everything this run will need before doing any work.

    Several failures used to surface only hours in - a missing tool at the tree
    stage, or an embedding interpreter without ``transformers`` after every tree was
    built. This checks the tools implied by the chosen flags, the input and
    reference files, and the embedding interpreter, and stops immediately with a
    clear list if anything is missing. Results are recorded in the manifest.
    """
    from viral_phylo.binaries import (find_fasttree_bin, find_torch_python,
                                      probe_torch_python, tool_is_available)

    blockers, warnings, checks = [], [], {}

    needed = {"foldmason": "structure database / 3Di extraction"}
    if args.aligner == "mafft":
        needed["mafft"] = "--aligner mafft"
    if getattr(args, "multi_alignment", False):
        mode = getattr(args, "cluster_mode", "structural")
        if mode == "structural":
            needed["foldseek"] = "--cluster-mode structural"
        elif mode == "sequence":
            needed["mmseqs"] = "--cluster-mode sequence"
    methods = {getattr(args, "master_method", "fasttree"), args.method}
    if "iqtree" in methods:
        needed["iqtree"] = "tree method iqtree"
    for tool, why in needed.items():
        ok = tool_is_available(tool)
        checks[tool] = ok
        if not ok:
            blockers.append(f"'{tool}' not found (needed for {why})")
    if getattr(args, "guide_trees", False) and "fasttree" not in methods and "iqtree" in methods:
        # Starting trees are an optimisation: without FastTree IQ-TREE simply runs unguided,
        # so this must not block the run.
        checks["guide_trees_fasttree"] = bool(find_fasttree_bin())
        if not checks["guide_trees_fasttree"]:
            warnings.append("no FastTree/VeryFastTree found: IQ-TREE starting trees are disabled and large "
                            "clusters will run unguided")
    if "fasttree" in methods:
        ft = find_fasttree_bin()
        checks["fasttree"] = bool(ft)
        if not ft:
            blockers.append("no FastTree/VeryFastTree binary found (needed for tree method fasttree)")

    for label, path in (("input folder", getattr(args, "input_folder", None)),
                        ("metadata", getattr(args, "metadata", None))):
        if path and not os.path.exists(path):
            blockers.append(f"{label} '{path}' does not exist")
    if args.aligner == "mafft" and not os.path.isfile(getattr(args, "mafft_matrix", "")):
        blockers.append(f"MAFFT 3Di matrix '{args.mafft_matrix}' not found - launch from the repository root")

    if getattr(args, "embed", False):
        py = find_torch_python()
        probe = probe_torch_python(py)
        checks["embedding_interpreter"] = {"python": py, **{k: probe[k] for k in
                                           ("usable", "cuda", "torch", "transformers", "missing")}}
        if not probe["usable"]:
            blockers.append(f"no Python interpreter can run the embeddings; best candidate {py} is "
                            f"missing {', '.join(probe['missing'])}")
        elif not probe["cuda"]:
            warnings.append(f"embeddings will run on CPU ({py}: torch {probe['torch']} cannot use the GPU "
                            "here); expect this stage to take much longer than on a GPU")

    manifest.data["preflight"] = {"checks": checks, "blockers": blockers, "warnings": warnings}
    for w in warnings:
        print(f"[Preflight] [!] {w}")
        manifest.add_warning(w)
    if blockers:
        for b in blockers:
            print(f"[Preflight] BLOCKER: {b}")
        manifest.record_stage("preflight", "failed", reason="; ".join(blockers))
        manifest.finish("failed")
        raise SystemExit(f"[Preflight] {len(blockers)} problem(s) must be fixed before this run can start.")
    manifest.record_stage("preflight", "ok", reason=(f"{len(warnings)} warning(s)" if warnings else None))
    print("[Preflight] All required tools, inputs and interpreters are present.")

def build_cli_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Viral Structural Phylogenetics & PLM Embedding Suite.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    subparsers = parser.add_subparsers(dest="subcommand", help="Subcommand to run")

    # Subcommand: fetch
    p_fetch = subparsers.add_parser("fetch", help="Download viral PDB structures from Viro3D or AlphaFold Database.")
    p_fetch.add_argument("--source", choices=["viro3d", "alphafold", "afdb"], default="viro3d", help="Structure repository: 'viro3d' or 'alphafold' (default: viro3d)")
    p_fetch.add_argument("-q", "--qualifier", default="glycoprotein", help="Protein search term (default: glycoprotein)")
    p_fetch.add_argument("-u", "--uniprot", default=None, help="UniProt ID(s) for AlphaFold DB, comma-separated (e.g. 'P00520,P04637')")
    p_fetch.add_argument("--uniprot-file", default=None, help="Path to text file with one UniProt ID per line")
    p_fetch.add_argument("-m", "--max-sequences", type=int, default=6, help="Max structures to download (default: 6)")
    p_fetch.add_argument("--format", choices=["pdb", "cif"], default="pdb", help="File format for AlphaFold DB (default: pdb)")
    p_fetch.add_argument("--download-pae", action="store_true", help="Also download Predicted Aligned Error (PAE) JSON from AlphaFold DB")
    p_fetch.add_argument("-o", "--output-dir", default="viro_3d_structures", help="Output directory for PDBs")

    # Subcommand: align
    p_align = subparsers.add_parser("align", help="Align PDB structures with FoldMason or MAFFT to generate 3Di alignments.")
    p_align.add_argument("-i", "--folder", default="viro_3d_structures", help="Folder containing PDB files (local or fetched)")
    p_align.add_argument("-o", "--output-dir", default="foldmason_alignments", help="Output directory for MSA")
    p_align.add_argument("--aligner", choices=["foldmason", "mafft"], default="foldmason", help="Multiple sequence alignment engine (default: foldmason)")
    p_align.add_argument("--foldmason-bin", default=None, help="Path to foldmason binary")
    p_align.add_argument("--mafft-bin", default=None, help="Path to mafft binary (if using --aligner mafft)")
    p_align.add_argument("--mafft-matrix", default="matrices/mat3di.out", help="3Di substitution matrix for MAFFT (default: matrices/mat3di.out)")
    p_align.add_argument("--min-coverage", type=float, default=0.70, help="Minimum alignment/sequence coverage threshold (0.0 to 1.0, default: 0.70)")
    p_align.add_argument("--multi-alignment", action="store_true", help="Partition dataset into multiple high-coverage sub-alignments (>= min-coverage)")
    p_align.add_argument("--filter-coverage", action="store_true", help="Filter alignment to only retain sequences meeting >= min-coverage")
    p_align.add_argument("-t", "--threads", default="8",
                         help="Maximum CPU threads for FoldMason/MAFFT (default: 8; 'AUTO' = every core)")

    # Subcommand: tree
    p_tree = subparsers.add_parser("tree", help="Build structural and/or amino acid phylogeny using IQ-TREE or FoldMason.")
    p_tree.add_argument("-a", "--alignment", default="foldmason_alignments/foldmason.fasta_3di.fa", help="3Di alignment file or directory")
    p_tree.add_argument("--alignment-aa", default=None, help="Amino acid alignment file (default: auto-detected)")
    p_tree.add_argument("--tree-type", choices=["3di", "aa", "both", "tanglegram"], default="both", help="Phylogeny type: '3di' (structural), 'aa' (sequence), or 'both'/'tanglegram' (default: both)")
    p_tree.add_argument("-m", "--method", choices=["iqtree", "fasttree", "foldmason"], default="iqtree",
                        help="Tree method (default: iqtree). 'fasttree' uses VeryFastTree/FastTree, with the --matrix "
                             "3Di model applied via -trans for 3Di trees.")
    p_tree.add_argument("--matrix", choices=["alphafold", "af", "esmfold", "llm", "both", "auto"], default="both", help="3Di substitution matrix")
    p_tree.add_argument("--rate-heterogeneity", default="auto", help="Rate heterogeneity: '+G4', '+I+G4', '+R', or 'auto' (default: auto)")
    p_tree.add_argument("--criterion", choices=["BIC", "AIC", "AICc"], default="BIC", help="Model selection criterion (default: BIC)")
    p_tree.add_argument("-b", "--bootstrap", type=int, default=1000, help="Ultrafast bootstrap replicates (default: 1000)")
    p_tree.add_argument("--alrt", type=int, default=1000, help="SH-aLRT replicates (default: 1000)")
    p_tree.add_argument("-t", "--threads", default="8", help="CPU threads for IQ-TREE (default: 8). 'AUTO' lets IQ-TREE benchmark thread counts, which is very slow on many-core machines.")
    p_tree.add_argument("--fast", action="store_true", help="Enable fast heuristic search mode (auto-enabled for large datasets)")
    p_tree.add_argument("-o", "--output-dir", default="phylogeny_results", help="Output directory for trees")
    p_tree.add_argument("-p", "--prefix", default="viral_tree", help="Tree output filename prefix")
    p_tree.add_argument("--iqtree-bin", default=None, help="Path to iqtree binary")
    p_tree.add_argument("--guide-trees", action=argparse.BooleanOptionalAction, default=True,
                        help="Start IQ-TREE from a validated VeryFastTree/FastTree tree on alignments of at least "
                             "--guide-tree-min-taxa taxa (default: on). A guided run that fails is rerun unguided.")
    p_tree.add_argument("--guide-tree-min-taxa", type=int, default=80,
                        help="Smallest alignment given a starting tree (default: 80)")

    # Subcommand: embed (PLM Embeddings & Hierarchical Clustering)
    p_embed = subparsers.add_parser("embed", help="Extract PLM embeddings (ESM-2 / ESM-C) and construct hierarchical clustering trees.")
    p_embed.add_argument("-i", "--input", "--fasta", dest="input", required=True, help="Path to input FASTA file (e.g. foldmason.fasta_aa.fa) or directory of PDB structures.")
    p_embed.add_argument("-o", "--output-dir", default="plm_results", help="Directory to store resulting treefile and embeddings (default: plm_results).")
    p_embed.add_argument("-m", "--model", default="esm2", choices=["esm2", "esmc", "facebook/esm2_t33_650M_UR50D", "biohub/ESMC-600M"], help="Protein language model (default: esm2 [650M]).")
    p_embed.add_argument("-c", "--clustering", default="upgma", choices=["upgma", "nj", "average", "complete", "single", "ward"], help="Hierarchical clustering method (default: upgma).")
    p_embed.add_argument("--metric", default="cosine", choices=["cosine", "euclidean", "l1", "cityblock", "manhattan"], help="Pairwise distance metric (default: cosine).")
    p_embed.add_argument("-p", "--prefix", default="viral_plm", help="Output file prefix (default: viral_plm).")
    p_embed.add_argument("--batch-size", type=int, default=1, help="Inference batch size (default: 1).")
    p_embed.add_argument("--max-length", type=int, default=1024, help="Max sequence length for truncation (default: 1024).")
    p_embed.add_argument("-t", "--threads", default="8",
                         help="Maximum CPU threads for embedding and clustering (default: 8; 'AUTO' = every core)")

    # Subcommand: report (visualise an existing tree file)
    p_report = subparsers.add_parser("report", help="Build a standalone interactive HTML report for an existing phylogenetic tree file.")
    p_report.add_argument("-t", "--tree", required=True, help="Tree file in Newick or NEXUS format (e.g. .treefile, .nwk, .tre, .nex)")
    p_report.add_argument("-meta", "--metadata", default=None, help="Optional metadata table (.xlsx, .csv, .tsv, .json) with one row per leaf")
    p_report.add_argument("-o", "--output", default=None, help="Output HTML path (default: <tree>_report.html next to the tree file)")
    p_report.add_argument("--title", default=None, help="Report title (default: tree file name)")

    # Subcommand: pipeline (all-in-one)
    p_pipe = subparsers.add_parser("pipeline", help="Run full pipeline: fetch/local -> align -> tree (+ optional PLM embed).")
    p_pipe.add_argument("-i", "--input-folder", default=None, help="Local directory of PDB structures (skips Viro3D fetch if specified)")
    p_pipe.add_argument("--source", choices=["viro3d", "alphafold", "afdb"], default="viro3d", help="Structure repository: 'viro3d' or 'alphafold' (default: viro3d)")
    p_pipe.add_argument("-q", "--qualifier", default="glycoprotein", help="Protein search term if fetching (default: glycoprotein)")
    p_pipe.add_argument("-u", "--uniprot", default=None, help="UniProt ID(s) for AlphaFold DB (comma-separated)")
    p_pipe.add_argument("--uniprot-file", default=None, help="Path to text file with one UniProt ID per line")
    p_pipe.add_argument("-c", "--count", type=int, default=6, help="Number of structures to fetch (default: 6)")
    p_pipe.add_argument("--format", choices=["pdb", "cif"], default="pdb", help="Format for AlphaFold DB: pdb or cif (default: pdb)")
    p_pipe.add_argument("--download-pae", action="store_true", help="Download PAE error matrix from AlphaFold DB")
    p_pipe.add_argument("--aligner", choices=["foldmason", "mafft"], default="foldmason", help="Multiple sequence alignment engine (default: foldmason)")
    p_pipe.add_argument("--mafft-bin", default=None, help="Path to mafft binary (if using --aligner mafft)")
    p_pipe.add_argument("--mafft-matrix", default="matrices/mat3di.out", help="3Di substitution matrix for MAFFT (default: matrices/mat3di.out)")

    # Similarity clustering. Replaces the greedy coverage partitioner, which measured
    # similarity through the master alignment it was supposed to be improving.
    p_pipe.add_argument("--cluster-mode", choices=["structural", "sequence", "coverage"], default="structural",
                        help="How to partition before aligning: 'structural' (foldseek, 3Di+AA), "
                             "'sequence' (mmseqs, amino acid), or 'coverage' (legacy greedy partitioner). "
                             "Default: structural.")
    p_pipe.add_argument("--cluster-evalue", type=float, default=1.0,
                        help="foldseek e-value for structural clustering (default: 1.0). Note an e-value "
                             "depends on database size and is not comparable across datasets of different sizes.")
    p_pipe.add_argument("--cluster-coverage", type=float, default=0.5,
                        help="Alignment coverage required to join a cluster (default: 0.5)")
    p_pipe.add_argument("--cluster-min-seq-id", type=float, default=0.5,
                        help="mmseqs sequence identity for --cluster-mode sequence (default: 0.5)")
    p_pipe.add_argument("--cluster-tmscore", type=float, default=None,
                        help="Use TM-align structural clustering at this TM-score instead of 3Di+AA. "
                             "Slower and scored worse on the Nipah binder set, but a TM-score is "
                             "comparable across datasets whereas an e-value is not.")
    p_pipe.add_argument("--cluster-sweep", action="store_true",
                        help="Sweep thresholds and select the best partition by silhouette, subject to "
                             "--cluster-min-size and --cluster-min-retained, instead of using one operating point.")
    p_pipe.add_argument("--cluster-min-size", type=int, default=10,
                        help="Cluster size counted as 'usable' when scoring a sweep (default: 10)")
    p_pipe.add_argument("--cluster-min-retained", type=float, default=0.70,
                        help="Fraction of taxa a sweep partition must keep in usable clusters (default: 0.70)")
    p_pipe.add_argument("--min-seq-length", type=int, default=50,
                        help="Drop sequences shorter than this before clustering (default: 50)")
    p_pipe.add_argument("--min-length-frac", type=float, default=0.4,
                        help="Also drop sequences shorter than this fraction of the median (default: 0.4)")
    p_pipe.add_argument("--full-set-max-gap", type=float, default=0.50,
                        help="Reject aligning the whole set as one when its gap fraction reaches this (default: 0.50)")
    p_pipe.add_argument("--full-set-max-expansion", type=float, default=1.5,
                        help="Reject the whole-set alignment when it exceeds this multiple of the median "
                             "ungapped sequence length (default: 1.5)")
    p_pipe.add_argument("--cluster-recurse", action="store_true",
                        help="Split any cluster whose own alignment is still poor, re-clustering it at a "
                             "stricter threshold. Driven by alignment quality rather than cluster size.")
    p_pipe.add_argument("--recurse-max-depth", type=int, default=2,
                        help="How many times a cluster may be split (default: 2)")
    p_pipe.add_argument("--recurse-max-gap", type=float, default=0.65,
                        help="Split a cluster whose alignment reaches this gap fraction (default: 0.65). "
                             "Deliberately looser than --full-set-max-gap: 40-55%% gaps is ordinary for a "
                             "diverse family, and the whole-set threshold would split almost every cluster.")
    p_pipe.add_argument("--recurse-on", choices=["3di", "both"], default="3di",
                        help="Which alignment decides whether a cluster is split (default: 3di, since "
                             "structural alignment is the purpose of this tool). 'both' also splits when the "
                             "AA alignment is over threshold, which matters with --aligner mafft and "
                             "--tree-type both.")
    p_pipe.add_argument("--recurse-min-kept", type=float, default=0.70,
                        help="Only accept a split that keeps at least this fraction of the cluster's taxa in "
                             "clusters of 4+ (default: 0.70). Guards against a split shattering a cluster "
                             "into singletons.")
    p_pipe.add_argument("--mafft-leavegappyregion", action="store_true",
                        help="Pass --leavegappyregion to MAFFT's amino-acid alignment, which inserts fewer "
                             "gaps into gap-rich regions. Only affects --aligner mafft.")
    p_pipe.add_argument("--recurse-max-expansion", type=float, default=3.0,
                        help="Split a cluster whose alignment exceeds this multiple of its median ungapped "
                             "sequence length (default: 3.0)")
    p_pipe.add_argument("--min-coverage", type=float, default=0.70, help="Minimum alignment/sequence coverage threshold (0.0 to 1.0, default: 0.70)")
    p_pipe.add_argument("--multi-alignment", action="store_true", help="Partition dataset into multiple high-coverage sub-alignments (>= min-coverage)")
    p_pipe.add_argument("--filter-coverage", action="store_true", help="Filter alignment to only retain sequences meeting >= min-coverage")
    p_pipe.add_argument("--tree-type", choices=["3di", "aa", "both", "tanglegram"], default="both", help="Phylogeny type (default: both)")
    p_pipe.add_argument("-m", "--method", choices=["iqtree", "fasttree", "foldmason"], default="iqtree",
                        help="Tree method for the per-cluster subtrees (default: iqtree). Set to "
                             "fasttree to use approximate-ML everywhere.")
    p_pipe.add_argument("--guide-trees", action=argparse.BooleanOptionalAction, default=True,
                        help="Start IQ-TREE from a VeryFastTree/FastTree tree on large alignments (see "
                             "--guide-tree-min-taxa). On by default: on the 210-taxon Nipah cluster all three guided "
                             "runs reached higher-likelihood trees than all three unguided ones (~30 log-lik units), "
                             "at similar runtime. The starting tree is repaired and validated for IQ-TREE, and a "
                             "guided run that fails is rerun unguided, so it can never lose a tree. "
                             "Use --no-guide-trees to disable.")
    p_pipe.add_argument("--guide-tree-min-taxa", type=int, default=80,
                        help="Smallest alignment given a FastTree starting tree with --guide-trees (default: 80, "
                             "the size at which IQ-TREE switches to a fixed model, so no model selection is affected).")
    p_pipe.add_argument("--master-method", choices=["fasttree", "iqtree", "foldmason"], default="fasttree",
                        help="Tree method for the whole-set tree (default: fasttree). IQ-TREE on a large, "
                             "gap-rich whole-set alignment can take days - a 1,193-taxon run was ~78 hours - "
                             "for a tree that is an overview rather than the analysis. FastTree gives the "
                             "overview in about a minute and uses the --matrix 3Di model via -trans.")
    p_pipe.add_argument("--matrix", choices=["alphafold", "af", "esmfold", "llm", "both", "auto"], default="both", help="3Di matrix (default: both)")
    p_pipe.add_argument("--rate-heterogeneity", default="auto", help="Rate heterogeneity model (default: auto)")
    p_pipe.add_argument("-b", "--bootstrap", type=int, default=1000, help="Ultrafast bootstrap replicates (default: 1000)")
    p_pipe.add_argument("--alrt", type=int, default=1000, help="SH-aLRT replicates (default: 1000)")
    p_pipe.add_argument("-t", "--threads", default="8",
                        help="Maximum CPU threads for the whole pipeline (default: 8). Every stage runs one at a "
                             "time and each tool is capped at this count: FoldMason, MAFFT, Foldseek/MMseqs2, "
                             "IQ-TREE, VeryFastTree/FastTree, ESM embeddings and NumPy/BLAS. 'AUTO' removes the cap "
                             "(tools use every core, and IQ-TREE benchmarks thread counts, which is very slow on "
                             "many-core machines).")
    p_pipe.add_argument("--tree-threads", default="8",
                        help="Threads for each per-cluster tree job (default: 8). IQ-TREE gains little from more on "
                             "short cluster alignments; spare threads run more clusters at once instead.")
    p_pipe.add_argument("--tree-jobs", type=int, default=None,
                        help="Cluster trees to build at once (default: --threads / --tree-threads, e.g. 32/8 = 4). "
                             "Never more than fit within --threads.")
    p_pipe.add_argument("-meta", "--metadata", default=None, help="Path to metadata file (.xlsx, .csv, .tsv, .json). If omitted, inferred from structures or Viro3D API.")
    p_pipe.add_argument("--fast", action="store_true", help="Enable fast search mode")
    p_pipe.add_argument("--embed", action="store_true", help="Extract PLM embeddings and build hierarchical clustering tree (ESM-2 / ESM-C)")
    p_pipe.add_argument("--embed-model", default="esm2", choices=["esm2", "esmc"], help="PLM model for embeddings (default: esm2)")
    p_pipe.add_argument("--embed-clustering", default="upgma", choices=["upgma", "nj", "average", "complete", "single", "ward"], help="Clustering method for PLM tree (default: upgma). Note: 'nj' produces no linkage matrix, so no silhouette profile is computed.")
    p_pipe.add_argument("--embed-metric", default="cosine", choices=["cosine", "euclidean", "l1", "cityblock", "manhattan"], help="Pairwise distance metric (default: cosine)")
    p_pipe.add_argument("-o", "--output-dir", default="glycoprotein_workflow", help="Parent output directory")

    return parser


# Thread pools that OpenMP and the BLAS libraries (and so NumPy, SciPy, scikit-learn,
# PyTorch and OpenMP builds of FastTree) size from the environment. Left unset they
# use every core.
THREAD_ENV_VARS = ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS",
                   "BLIS_NUM_THREADS", "VECLIB_MAXIMUM_THREADS", "NUMEXPR_NUM_THREADS")


def apply_thread_cap(threads):
    """Cap every thread pool this process and its children can open at ``threads``.

    Tools that take an explicit count are passed it at their call sites; this covers
    the rest - child processes through the environment (a lower value the user
    already set is kept), and this process's already-loaded BLAS through
    threadpoolctl. Returns the cap, or None when ``threads`` is AUTO or unset.
    """
    try:
        cap = int(threads)
    except (TypeError, ValueError):
        return None
    if cap <= 0:
        return None
    for var in THREAD_ENV_VARS:
        try:
            current = int(os.environ.get(var, ""))
        except ValueError:
            current = None
        os.environ[var] = str(min(cap, current) if current and current > 0 else cap)
    try:
        from threadpoolctl import threadpool_limits
        threadpool_limits(limits=cap)
    except Exception:
        pass
    return cap


def main():
    parser = build_cli_parser()
    args = parser.parse_args()
    thread_cap = apply_thread_cap(getattr(args, "threads", None))

    if args.subcommand == "fetch":
        source = getattr(args, "source", "viro3d")
        if source in ("alphafold", "afdb") or getattr(args, "uniprot", None) or getattr(args, "uniprot_file", None):
            u_ids = []
            if getattr(args, "uniprot", None):
                u_ids.extend([x.strip() for x in args.uniprot.split(",") if x.strip()])
            if getattr(args, "uniprot_file", None) and os.path.isfile(args.uniprot_file):
                with open(args.uniprot_file, "r") as uf:
                    u_ids.extend([line.strip() for line in uf if line.strip() and not line.startswith("#")])
            fetch_alphafold_structures(
                uniprot_ids=u_ids if u_ids else None,
                query=args.qualifier if not u_ids else None,
                max_sequences=args.max_sequences,
                output_dir=args.output_dir,
                file_format=getattr(args, "format", "pdb"),
                download_pae=getattr(args, "download_pae", False),
            )
        else:
            fetch_structures(args.qualifier, args.max_sequences, args.output_dir)
    elif args.subcommand == "align":
        align_structures(
            args.folder,
            args.output_dir,
            foldmason_bin=args.foldmason_bin,
            aligner=args.aligner,
            mafft_bin=args.mafft_bin,
            mafft_matrix=args.mafft_matrix,
            min_coverage=getattr(args, "min_coverage", 0.70),
            multi_alignment=getattr(args, "multi_alignment", False),
            filter_coverage=getattr(args, "filter_coverage", False),
            threads=args.threads,
        )
    elif args.subcommand == "tree":
        aln_file = args.alignment
        if os.path.isdir(aln_file):
            for cand in ["mafft.fasta_3di.fa", "foldmason.fasta_3di.fa", "foldmason_3di.fa"]:
                if os.path.isfile(os.path.join(aln_file, cand)):
                    aln_file = os.path.join(aln_file, cand)
                    break
            else:
                aln_file = os.path.join(aln_file, "foldmason.fasta_3di.fa")
        build_tree(
            alignment_file=aln_file,
            method=args.method,
            matrix=args.matrix,
            rate_het=args.rate_heterogeneity,
            criterion=args.criterion,
            bootstrap=args.bootstrap,
            alrt=args.alrt,
            threads=args.threads,
            output_dir=args.output_dir,
            prefix=args.prefix,
            iqtree_bin=args.iqtree_bin,
            tree_type=args.tree_type,
            alignment_aa=args.alignment_aa,
            fast=args.fast,
            guide_trees=args.guide_trees,
            guide_min_taxa=args.guide_tree_min_taxa,
        )
    elif args.subcommand == "embed":
        embed_script = os.path.join(os.path.dirname(__file__), "embed_and_cluster.py")
        if not os.path.isfile(embed_script):
            embed_script = "scripts/embed_and_cluster.py"
        target_py = find_torch_python()
        cmd = [
            target_py, embed_script,
            "-i", args.input,
            "-o", args.output_dir,
            "-m", args.model,
            "-c", args.clustering,
            "--metric", args.metric,
            "-p", args.prefix,
            "--batch-size", str(args.batch_size),
            "--max-length", str(args.max_length),
        ]
        print(f"[PLM] Running embedding extraction using interpreter: {target_py}")
        sub_env = dict(os.environ)
        sub_env["KMP_DUPLICATE_LIB_OK"] = "TRUE"
        subprocess.run(cmd, check=True, env=sub_env)
    elif args.subcommand == "report":
        out = build_tree_report(args.tree, metadata_path=args.metadata, output_path=args.output, title=args.title)
        print(f"[Report] Interactive tree report written to: file://{out.resolve()}")
    elif args.subcommand == "pipeline":
        out_base = args.output_dir
        aln_dir = os.path.join(out_base, "alignment")
        phy_dir = os.path.join(out_base, "phylogeny")

        # Start the run manifest before anything else, so that even a run which dies
        # in the first stage leaves a record of how it was invoked. It is rewritten
        # after every stage rather than only at the end.
        manifest = RunManifest(out_base)
        manifest.install_exit_handlers()
        manifest.data["resources"] = {
            "thread_cap": thread_cap,
            "thread_env": {v: os.environ.get(v) for v in THREAD_ENV_VARS},
            "cpu_count": os.cpu_count(),
        }
        if thread_cap is None:
            manifest.add_warning("--threads AUTO: no thread cap; tools may use all "
                                 f"{os.cpu_count()} cores")
        # Line-buffer stdout so a nohup'd log shows progress as it happens; when
        # redirected to a file Python otherwise block-buffers and the log lags badly.
        try:
            sys.stdout.reconfigure(line_buffering=True)
        except (AttributeError, ValueError):
            pass
        manifest.record_invocation(argv=sys.argv, args=args)
        manifest.record_tool_versions(overrides={"mafft": getattr(args, "mafft_bin", None)})
        print(f"[Pipeline] Run manifest: '{manifest.json_path}'")
        pipeline_preflight(args, manifest)

        # 1. Fetch or Local Input
        _t_stage = time.time()
        source = getattr(args, "source", "viro3d")
        if args.input_folder:
            if not os.path.isdir(args.input_folder):
                raise FileNotFoundError(f"Specified local structures directory '{args.input_folder}' does not exist.")
            pdb_dir = args.input_folder
            print(f"[Pipeline] Using local structure folder: '{pdb_dir}'")
        elif source in ("alphafold", "afdb") or getattr(args, "uniprot", None) or getattr(args, "uniprot_file", None):
            pdb_dir = os.path.join(out_base, "structures")
            u_ids = []
            if getattr(args, "uniprot", None):
                u_ids.extend([x.strip() for x in args.uniprot.split(",") if x.strip()])
            if getattr(args, "uniprot_file", None) and os.path.isfile(args.uniprot_file):
                with open(args.uniprot_file, "r") as uf:
                    u_ids.extend([line.strip() for line in uf if line.strip() and not line.startswith("#")])
            fetch_alphafold_structures(
                uniprot_ids=u_ids if u_ids else None,
                query=args.qualifier if not u_ids else None,
                max_sequences=args.count,
                output_dir=pdb_dir,
                file_format=getattr(args, "format", "pdb"),
                download_pae=getattr(args, "download_pae", False),
            )
        else:
            pdb_dir = os.path.join(out_base, "structures")
            fetch_structures(args.qualifier, args.count, pdb_dir)

        manifest.record_input(pdb_dir)
        manifest.record_stage(
            "input" if args.input_folder else "fetch",
            "ok" if os.path.isdir(pdb_dir) else "failed",
            reason=None if os.path.isdir(pdb_dir) else f"structure directory '{pdb_dir}' missing",
            duration_s=time.time() - _t_stage, outputs=[pdb_dir],
        )
        _t_stage = time.time()

        # 2. Align (with optional Coverage Threshold & Multi-Alignment Partitioning)
        min_cov = getattr(args, "min_coverage", 0.70)
        multi_aln = getattr(args, "multi_alignment", False)
        filter_cov = getattr(args, "filter_coverage", False)

        cluster_mode = getattr(args, "cluster_mode", "structural")
        use_similarity = multi_aln and cluster_mode != "coverage"

        # Always build the whole-set alignment first: the master tree needs it, and it
        # is what the full-set gate below is judged on. Under similarity clustering it
        # is a diagnostic only - it no longer drives the partitioning.
        aln_3di = align_structures(
            pdb_dir,
            aln_dir,
            aligner=args.aligner,
            mafft_bin=getattr(args, "mafft_bin", None),
            mafft_matrix=getattr(args, "mafft_matrix", "matrices/mat3di.out"),
            min_coverage=min_cov,
            multi_alignment=multi_aln and not use_similarity,
            filter_coverage=filter_cov,
            mafft_leavegappyregion=getattr(args, "mafft_leavegappyregion", False),
            threads=args.threads,
        )

        # Full-set gate: can the whole input reasonably be aligned as one set?
        # Judged on the 3Di alignment, consistent with cluster recursion: structural
        # alignment is the purpose of this tool. The AA alignment is recorded alongside.
        full_set_stats = full_set_usable = full_set_reason = None
        _full_3di = os.path.join(aln_dir, "foldmason.fasta_3di.fa")
        _full_aa = os.path.join(aln_dir, "foldmason.fasta_aa.fa")
        if os.path.isfile(_full_3di):
            full_set_stats = alignment_gap_fraction(cluster_parse_fasta(_full_3di))
            full_set_usable, full_set_reason = full_set_alignment_is_usable(
                full_set_stats,
                max_gap_fraction=getattr(args, "full_set_max_gap", 0.50),
                max_expansion=getattr(args, "full_set_max_expansion", 1.5),
            )
            verdict = "USABLE" if full_set_usable else "NOT usable"
            print(f"\n[Pipeline] Whole-set 3Di alignment: {full_set_stats['taxa']} taxa x "
                  f"{full_set_stats['columns']} columns, {full_set_stats['gap_fraction']:.1%} gaps, "
                  f"{full_set_stats['expansion']:.2f}x median length -> {verdict} ({full_set_reason})")
            if not full_set_usable and not multi_aln:
                print("[Pipeline] [!] Notice: the whole set does not align well as one group. "
                      "Re-run with --multi-alignment to partition it first.")
            manifest.record_alignment({"full_set": {
                **full_set_stats, "measured_on": "3di", "usable": full_set_usable, "reason": full_set_reason,
                "aa": alignment_gap_fraction(cluster_parse_fasta(_full_aa)) if os.path.isfile(_full_aa) else None,
            }})

        if args.input_folder:
            prefix = "local_tree"
        elif source in ("alphafold", "afdb"):
            prefix = "alphafold_tree"
        else:
            prefix = f"{args.qualifier}_tree"

        # 2b. Optional PLM embeddings. They run before clustering, so a --cluster-sweep
        # in the same run scores its partitions with this run's embeddings, never
        # stale ones left in the output directory by an earlier run.
        _t_align_start, _t_stage = _t_stage, time.time()
        if getattr(args, "embed", False):
            print(f"\n[Pipeline] Extracting PLM embeddings ({args.embed_model}) & constructing hierarchical tree ({args.embed_clustering})...")
            embed_script = os.path.join(os.path.dirname(__file__), "embed_and_cluster.py")
            if not os.path.isfile(embed_script):
                embed_script = "scripts/embed_and_cluster.py"
            aa_aln = os.path.join(aln_dir, "foldmason.fasta_aa.fa")
            for cand in ["mafft.fasta_aa.fa", "foldmason.fasta_aa.fa", "foldmason_aa.fa"]:
                if os.path.isfile(os.path.join(aln_dir, cand)):
                    aa_aln = os.path.join(aln_dir, cand)
                    break
            embed_src = aa_aln if os.path.isfile(aa_aln) else pdb_dir
            target_py = find_torch_python()
            cmd = [
                target_py, embed_script,
                "-i", embed_src,
                "-o", phy_dir,
                "-m", args.embed_model,
                "-c", args.embed_clustering,
                "--metric", args.embed_metric,
                "-p", prefix,
                "--batch-size", "1",
            ]
            sub_env = dict(os.environ)
            sub_env["KMP_DUPLICATE_LIB_OK"] = "TRUE"
            _embed_error = None
            try:
                subprocess.run(cmd, check=True, env=sub_env)
            except (subprocess.CalledProcessError, OSError) as exc:
                _embed_error = f"{type(exc).__name__}: {exc}"
                print(f"[Pipeline] [!] Notice: embedding stage failed: {_embed_error}. "
                      "Continuing to build the dashboard without it.")

            _embed_outputs = sorted(
                glob.glob(os.path.join(phy_dir, f"{prefix}_*_{args.embed_model}.*"))
                + glob.glob(os.path.join(phy_dir, f"{prefix}_{args.embed_model}.*"))
            )
            _sil = [f for f in _embed_outputs if "_silhouette_" in f]
            manifest.record_stage(
                "embed", "failed" if _embed_error else ("ok" if _embed_outputs else "skipped"),
                reason=_embed_error or (None if _sil else
                        f"no silhouette profile: clustering '{args.embed_clustering}' yields no linkage matrix"),
                duration_s=time.time() - _t_stage, command=cmd, outputs=_embed_outputs,
            )
            if not _sil and not _embed_error:
                manifest.add_warning(
                    f"no silhouette profile produced (--embed-clustering {args.embed_clustering}); "
                    "use upgma|average|complete|ward to get one")
        else:
            manifest.record_stage("embed", "skipped", reason="--embed not requested")
        # The align stage is recorded after clustering; keep embedding time out of it.
        _t_stage = _t_align_start + (time.time() - _t_stage)

        multi_summary = None
        if use_similarity:
            _emb_npz = None
            if getattr(args, "embed", False) and args.embed_model == "esm2":
                _cand = os.path.join(phy_dir, f"{prefix}_embeddings_esm2.npz")
                _emb_npz = _cand if os.path.isfile(_cand) else None
            multi_summary = partition_structures_by_similarity(
                pdb_dir, aln_dir,
                mode=cluster_mode,
                evalue=getattr(args, "cluster_evalue", 1.0),
                coverage=getattr(args, "cluster_coverage", 0.5),
                min_seq_id=getattr(args, "cluster_min_seq_id", 0.5),
                tmscore=getattr(args, "cluster_tmscore", None),
                alignment_type=1 if getattr(args, "cluster_tmscore", None) else 2,
                sweep=getattr(args, "cluster_sweep", False),
                min_cluster_size=getattr(args, "cluster_min_size", 10),
                min_retained=getattr(args, "cluster_min_retained", 0.70),
                embeddings_npz=_emb_npz,
                min_length=getattr(args, "min_seq_length", 50),
                min_length_frac=getattr(args, "min_length_frac", 0.4),
                recurse=getattr(args, "cluster_recurse", False),
                recurse_max_depth=getattr(args, "recurse_max_depth", 2),
                recurse_max_gap=getattr(args, "recurse_max_gap", 0.65),
                recurse_max_expansion=getattr(args, "recurse_max_expansion", 3.0),
                recurse_min_kept=getattr(args, "recurse_min_kept", 0.70),
                recurse_on=getattr(args, "recurse_on", "3di"),
                mafft_leavegappyregion=getattr(args, "mafft_leavegappyregion", False),
                aligner=args.aligner,
                mafft_bin=getattr(args, "mafft_bin", None),
                mafft_matrix=getattr(args, "mafft_matrix", "matrices/mat3di.out"),
                threads=int(args.threads) if str(args.threads).isdigit() else 0,
                manifest=manifest,
            )
            shutil.copy2(os.path.join(aln_dir, "multi_alignment_summary.json"),
                         os.path.join(out_base, "multi_alignment_summary.json"))
            _over3 = multi_summary.get("over_threshold_3di_clusters") or {}
            if _over3:
                from collections import Counter
                why = ", ".join(f"{n} {r}" for r, n in Counter(_over3.values()).items())
                msg = (f"{len(_over3)} tree-eligible cluster(s) kept with a 3Di alignment over the recursion "
                       f"thresholds ({why})")
                print(f"[Pipeline] [!] Notice: {msg}.")
                manifest.add_warning(msg + "; see clustering.over_threshold_3di_clusters.")
            _aa_flagged = multi_summary.get("aa_over_threshold_clusters") or []
            if _aa_flagged and args.tree_type in ("aa", "both", "tanglegram"):
                if getattr(args, "recurse_on", "3di") == "both":
                    why = ("AA was judged too, but these could not be split further (no split kept enough "
                           "taxa, or the depth limit was reached)")
                else:
                    why = "splits were judged on 3Di only (use --recurse-on both to also split on AA)"
                msg = (f"{len(_aa_flagged)} tree-eligible cluster(s) have an AA alignment over the "
                       f"recursion thresholds, so their AA trees rest on a gappier alignment than their "
                       f"3Di trees; {why}.")
                print(f"[Pipeline] [!] Notice: {msg}")
                manifest.add_warning(msg + " See clustering.aa_over_threshold_clusters.")
        elif multi_aln:
            summary_f = os.path.join(aln_dir, "multi_alignment_summary.json")
            if os.path.isfile(summary_f):
                with open(summary_f, "r", encoding="utf-8") as f:
                    multi_summary = json.load(f)
                # Copy summary to parent output directory
                shutil.copy2(summary_f, os.path.join(out_base, "multi_alignment_summary.json"))

        aln_aa = os.path.join(aln_dir, "foldmason.fasta_aa.fa")
        write_alignment_info(
            aln_dir,
            aligner=args.aligner,
            produced={
                os.path.basename(aln_3di): args.aligner,
                os.path.basename(aln_aa): args.aligner,
            },
            mafft_matrix=getattr(args, "mafft_matrix", None) if args.aligner == "mafft" else None,
            tool_paths={"mafft": getattr(args, "mafft_bin", None)},
            extra={"multi_alignment": multi_aln, "filter_coverage": filter_cov, "min_coverage": min_cov,
                   "mafft_leavegappyregion": getattr(args, "mafft_leavegappyregion", False)},
        )
        manifest.record_alignment({
            "aligner": args.aligner,
            "mafft_matrix": getattr(args, "mafft_matrix", None) if args.aligner == "mafft" else None,
            "alignment_dir": aln_dir,
            "note": "foldmason.fasta_*.fa is a legacy compatibility name; see 'aligner' for the real producer.",
        })
        manifest.record_stage(
            "align", "ok" if os.path.isfile(aln_3di) else "skipped",
            reason=None if os.path.isfile(aln_3di) else "no 3Di alignment was produced",
            duration_s=time.time() - _t_stage, outputs=[aln_3di, aln_aa],
        )
        _t_stage = time.time()

        # 3. Tree

        _master_method = getattr(args, "master_method", "fasttree")
        _master_error = None
        try:
            build_tree(
                alignment_file=aln_3di,
                method=_master_method,
                matrix=args.matrix,
                rate_het=args.rate_heterogeneity,
                criterion="BIC",
                bootstrap=args.bootstrap,
                alrt=args.alrt,
                threads=args.threads,
                output_dir=phy_dir,
                prefix=prefix,
                tree_type=args.tree_type,
                fast=args.fast,
                guide_trees=getattr(args, "guide_trees", False),
                guide_min_taxa=getattr(args, "guide_tree_min_taxa", 80),
            )
        except Exception as exc:
            # The whole-set tree is an overview; the per-cluster trees are the analysis.
            # A failure here must not cancel them.
            _master_error = f"{type(exc).__name__}: {exc}"
            print(f"[Pipeline] [!] Notice: whole-set tree ({_master_method}) failed: {_master_error}. "
                  "Continuing with the per-cluster trees.")

        _master_trees = sorted(glob.glob(os.path.join(phy_dir, f"{prefix}*.treefile")))
        _master_info = [json.load(open(t + ".info.json")) for t in _master_trees if os.path.isfile(t + ".info.json")]
        manifest.record_stage(
            "tree-master", "failed" if _master_error else ("ok" if _master_trees else "skipped"),
            reason=_master_error or (None if _master_trees else f"{_master_method} produced no .treefile"),
            duration_s=time.time() - _t_stage, outputs=_master_trees,
            command="; ".join(i["command"] for i in _master_info) or _master_method,
        )
        if _master_info:
            manifest.data.setdefault("trees", {})["master"] = _master_info
        _t_stage = time.time()

        # 3b. Build trees for coverage-partitioned clusters
        if multi_summary and "clusters" in multi_summary:
            print(f"\n[Pipeline] Building phylogenetic trees for {len(multi_summary['clusters'])} coverage-partitioned cluster(s)...")
            _cluster_ok, _cluster_failed, _cluster_skipped = [], [], []
            _n_jobs, _job_threads, _plan_note = plan_tree_jobs(
                args.threads, getattr(args, "tree_threads", "8"), getattr(args, "tree_jobs", None))
            if _plan_note:
                print(f"[Pipeline] [!] Notice: {_plan_note}.")
                manifest.add_warning(_plan_note)
            _jobs = []
            for c in multi_summary["clusters"]:
                if c.get("taxa_count", 0) >= 4 and c.get("aln_3di") and os.path.isfile(c["aln_3di"]):
                    c_phy_dir = os.path.join(phy_dir, c["name"])
                    cmd = [
                        sys.executable, "-c", "import sys; from viral_phylo.cli import main; sys.exit(main())",
                        "tree",
                        "--alignment", c["aln_3di"],
                        "--method", args.method,
                        "--matrix", args.matrix,
                        "--rate-heterogeneity", args.rate_heterogeneity,
                        "--criterion", "BIC",
                        "--bootstrap", str(args.bootstrap),
                        "--alrt", str(args.alrt),
                        "--threads", str(_job_threads),
                        "--output-dir", c_phy_dir,
                        "--prefix", f"{prefix}_{c['name']}",
                        "--tree-type", args.tree_type,
                        "--guide-trees" if getattr(args, "guide_trees", False) else "--no-guide-trees",
                        "--guide-tree-min-taxa", str(getattr(args, "guide_tree_min_taxa", 80)),
                    ]
                    if args.fast:
                        cmd.append("--fast")
                    _jobs.append({"name": c["name"], "taxa": c.get("taxa_count", 0), "cmd": cmd,
                                  "log": os.path.join(c_phy_dir, "tree_job.log")})
                else:
                    _cluster_skipped.append({"cluster": c["name"], "taxa": c.get("taxa_count"),
                                             "reason": ("fewer than 4 taxa - too few for a meaningful tree"
                                                        if c.get("taxa_count", 0) < 4
                                                        else "no 3Di alignment was produced")})
            if _jobs:
                print(f"[Pipeline] {len(_jobs)} cluster tree(s) to build: {_n_jobs} at a time, "
                      f"{_job_threads} threads each (cap --threads {args.threads}), largest first. "
                      "Each job's full output goes to <cluster>/tree_job.log.")
                # Children size OpenMP/BLAS pools from the environment: give them the
                # per-job count, so concurrent jobs together stay within the cap.
                _job_env = dict(os.environ)
                if str(_job_threads).isdigit():
                    for _var in THREAD_ENV_VARS:
                        _job_env[_var] = str(_job_threads)
                _results = run_tree_jobs(_jobs, _n_jobs, env=_job_env)
                for r in _results:
                    if r.get("returncode") == 0:
                        _cluster_ok.append(r["name"])
                    else:
                        _cluster_failed.append({"cluster": r["name"], "taxa": r["taxa"], "error": r.get("error"),
                                                "log": r["log"], "log_tail": r.get("log_tail")})
                manifest.data.setdefault("trees", {})["cluster_jobs"] = {
                    "parallel_jobs": _n_jobs, "threads_per_job": _job_threads, "thread_cap": args.threads,
                    "jobs": [{"cluster": r["name"], "taxa": r["taxa"], "returncode": r.get("returncode"),
                              "seconds": r.get("seconds"), "log": r["log"],
                              "command": "viral-phylo " + " ".join(r["cmd"][3:])}
                             for r in _results],
                }

            manifest.record_stage(
                "tree-clusters",
                "ok" if _cluster_ok else "skipped",
                reason=(f"{len(_cluster_ok)} ok, {len(_cluster_failed)} failed, "
                        f"{sum(1 for x in _cluster_skipped if (x.get('taxa') or 0) < 4)} skipped as too small, "
                        f"{sum(1 for x in _cluster_skipped if (x.get('taxa') or 0) >= 4)} skipped for lack of an alignment"),
                duration_s=time.time() - _t_stage,
            )
            manifest.data["clustering"]["per_cluster_trees"] = {
                "ok": _cluster_ok, "failed": _cluster_failed, "skipped": _cluster_skipped,
            }
            if _cluster_failed:
                manifest.add_warning(f"{len(_cluster_failed)} cluster tree(s) failed; see clustering.per_cluster_trees")
        _t_stage = time.time()

        if getattr(args, "guide_trees", False):
            _st = []
            for side in sorted(glob.glob(os.path.join(phy_dir, "**", "*.start_tree.json"), recursive=True)):
                try:
                    rec = json.load(open(side))
                except Exception:
                    continue
                rec["tree"] = os.path.relpath(side[:-len(".start_tree.json")], phy_dir)
                _st.append(rec)
            manifest.data.setdefault("trees", {})["start_trees"] = _st
            _count = {k: sum(1 for r in _st if r.get("outcome") == k) for k in ("guided", "fell_back", "not_used")}
            print(f"[Pipeline] FastTree starting trees: {_count['guided']} guided, {_count['fell_back']} fell back "
                  f"to unguided, {_count['not_used']} not used.")
            if _count["fell_back"] or _count["not_used"]:
                manifest.add_warning(f"starting trees: {_count['fell_back']} guided IQ-TREE run(s) failed and were "
                                     f"rerun unguided, {_count['not_used']} starting tree(s) were rejected; "
                                     "see trees.start_trees")

        # 4. Universal Metadata Processing
        print("\n[Pipeline] Analyzing and indexing metadata...")
        meta_src = getattr(args, "metadata", None)
        if not meta_src and not args.input_folder:
            cand_meta = os.path.join(pdb_dir, "taxa_metadata.json")
            if os.path.isfile(cand_meta):
                meta_src = cand_meta
        meta_res = parse_metadata(meta_src, structure_dir=pdb_dir)
        meta_out_path = os.path.join(out_base, "taxa_metadata.json")
        with open(meta_out_path, "w", encoding="utf-8") as f:
            json.dump(meta_res, f, indent=2)
        print(f"[Pipeline] Saved standardized metadata schema to: '{meta_out_path}' ({len(meta_res['taxa'])} taxa, {len(meta_res['columns'])} columns).")

        manifest.record_stage(
            "metadata", "ok", duration_s=time.time() - _t_stage, outputs=[meta_out_path],
        )
        _t_stage = time.time()

        # 5. PLM embeddings ran before clustering (see above).

        # 6. Generate Standalone Interactive HTML Visualizations & Datasets
        print("\n[Pipeline] Compiling interactive HTML suite and supporting assets...")

        try:
            # Structure parsing lives in viral_phylo.structures, which handles both PDB
            # fixed-column records and mmCIF _atom_site loops. The previous inline
            # parser applied PDB column slicing to every input, so mmCIF structures
            # silently yielded no backbone at all.
            ca_dict = extract_ca_dict(pdb_dir)
            if ca_dict:
                write_ca_bundle(ca_dict, os.path.join(out_base, "ca_structures.js"))
            else:
                print(f"[Pipeline] [!] Notice: no C-alpha backbones extracted from '{pdb_dir}'; "
                      "the 3D structure viewer will have no data for this run.")

            # Also ensure all cohort structure packs exist in out_base so all datasets load cleanly
            for sf in ["ca_500_structures.js", "ca_1193_structures.js", "ca_100_structures.js", "alignments_data.js"]:
                results_sf = os.path.join("results", sf)
                dst_sf = os.path.join(out_base, sf)
                if os.path.isfile(results_sf) and not os.path.isfile(dst_sf):
                    try:
                        shutil.copy2(results_sf, dst_sf)
                    except Exception:
                        pass
            manifest.record_stage(
                "ca-extraction", "ok" if ca_dict else "skipped",
                reason=None if ca_dict else f"no C-alpha atoms parsed from any structure in '{pdb_dir}'",
                duration_s=time.time() - _t_stage,
                outputs=[os.path.join(out_base, "ca_structures.js")] if ca_dict else [],
            )
        except Exception as e:
            manifest.record_stage("ca-extraction", "failed", reason=f"{type(e).__name__}: {e}",
                                  duration_s=time.time() - _t_stage)
            print(f"[Pipeline] Notice: C-alpha extraction encountered: {e}")
        _t_stage = time.time()

        # Build the alignment bundle and the interactive HTML in-process. These were
        # previously shelled out to scripts looked up under viral_phylo/, where they
        # have never lived, so os.path.isfile() was always False and both steps were
        # skipped without any warning while the pipeline still reported success.
        wf_html = os.path.join(out_base, "interactive_tree.html")
        root_html = os.path.abspath("interactive_tree.html")
        _html_error = None
        try:
            from viral_phylo.web.alignments import build_alignments_data
            from viral_phylo.web.builder import build_interactive_tree

            build_alignments_data()
            build_interactive_tree()
        except Exception as e:
            _html_error = f"{type(e).__name__}: {e}"
            print(f"[Pipeline] [!] Notice: Automated HTML compilation failed: {_html_error}")

        # Count only files written by this run: the repo-root dashboard and older
        # per-run copies already exist from earlier runs, so mere existence proves nothing.
        def _fresh(path):
            return os.path.isfile(path) and os.path.getmtime(path) >= _t_stage - 1
        _html_outputs = [p for p in (wf_html, root_html, os.path.join(out_base, "alignments_data.js"))
                         if _fresh(p)]
        manifest.record_stage(
            "html-compile",
            "failed" if _html_error else ("ok" if _html_outputs else "skipped"),
            reason=_html_error or (None if _html_outputs else "no interactive HTML or alignment bundle was written"),
            duration_s=time.time() - _t_stage, outputs=_html_outputs,
        )

        # The banner must agree with the manifest: a dashboard existing does not make a
        # run successful if stages failed or were skipped along the way.
        _failed = [st["name"] for st in manifest.data["stages"] if st["status"] == "failed"]
        _skipped = [st["name"] for st in manifest.data["stages"] if st["status"] == "skipped"]
        print("\n[Pipeline] ========================================================")
        if _failed:
            print(f"[Pipeline] ❌ Pipeline finished with {len(_failed)} FAILED stage(s): {', '.join(_failed)}")
            print("[Pipeline]    See the reasons in the run manifest below.")
        elif not _fresh(wf_html):
            print("[Pipeline] ⚠️  Pipeline finished, but NO interactive HTML was produced.")
            print("[Pipeline]     Re-run `build-alignments` then `build-tree-view`, or check the notices above.")
        elif _skipped:
            print(f"[Pipeline] ✅ Pipeline complete, with skipped stage(s): {', '.join(_skipped)}")
        else:
            print("[Pipeline] ✅ Full Pipeline Execution & Visualization Complete!")
        if _fresh(wf_html):
            print(f"[Pipeline] 📊 Standalone Visualizer: file://{os.path.abspath(wf_html)}")
        if _fresh(root_html):
            print(f"[Pipeline] 🌐 Dashboard Visualizer:  file://{root_html}")

        # Finalise the manifest: the inventory of what this run actually produced,
        # and an overall status that reflects any failed or skipped stage.
        manifest.record_outputs([
            wf_html, root_html, meta_out_path,
            os.path.join(out_base, "ca_structures.js"),
            os.path.join(out_base, "alignments_data.js"),
            os.path.join(out_base, "multi_alignment_summary.json"),
        ] + sorted(glob.glob(os.path.join(phy_dir, "**", "*.treefile"), recursive=True)))
        _statuses = [st["status"] for st in manifest.data["stages"]]
        manifest.finish("failed" if "failed" in _statuses
                        else ("completed_with_skips" if "skipped" in _statuses else "completed"))
        print(f"[Pipeline] 🧾 Run manifest:          file://{os.path.abspath(manifest.json_path)}")
        print("[Pipeline] ========================================================\n")




if __name__ == "__main__":
    main()
