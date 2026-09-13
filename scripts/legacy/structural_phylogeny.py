#!/usr/bin/env python3
"""Build structural phylogenetic trees using FoldMason guide trees or IQ-TREE with 3Di substitution matrices."""

import argparse
import os
import shutil
import subprocess
import sys
import urllib.request

# Persistent Edmond Dataverse URLs for 3Di substitution matrices (doi:10.17617/3.1MJJBH)
EDMOND_MATRICES = {
    "alphafold": {
        "filename": "Q.3Di.AF",
        "file_id": 311466,
        "url": "https://edmond.mpg.de/api/access/datafile/311466",
        "desc": "AlphaFold-derived 3Di substitution matrix (Q.3Di.AF)",
    },
    "esmfold": {
        "filename": "Q.3Di.LLM",
        "file_id": 311467,
        "url": "https://edmond.mpg.de/api/access/datafile/311467",
        "desc": "ESMFold / LLM-derived 3Di substitution matrix (Q.3Di.LLM)",
    },
}


def ensure_matrix_file(matrix_key_or_path: str, matrices_dir: str = "matrices") -> str:
    """Ensure the selected substitution matrix is available locally, downloading if necessary."""
    norm_key = matrix_key_or_path.lower().strip()
    if norm_key in ("both", "auto"):
        p_af = ensure_matrix_file("alphafold", matrices_dir)
        p_llm = ensure_matrix_file("esmfold", matrices_dir)
        return f"{p_af},{p_llm}"
    elif norm_key in ("af", "alphafold", "q.3di.af"):
        meta = EDMOND_MATRICES["alphafold"]
    elif norm_key in ("llm", "esm", "esmfold", "q.3di.llm"):
        meta = EDMOND_MATRICES["esmfold"]
    else:
        if os.path.isfile(matrix_key_or_path):
            return matrix_key_or_path
        raise FileNotFoundError(f"Custom matrix path '{matrix_key_or_path}' not found.")

    os.makedirs(matrices_dir, exist_ok=True)
    target_path = os.path.join(matrices_dir, meta["filename"])

    if os.path.isfile(target_path) and os.path.getsize(target_path) > 0:
        return target_path

    print(f"Downloading {meta['desc']} from Edmond...")
    try:
        req = urllib.request.Request(meta["url"], headers={"User-Agent": "Mozilla/5.0"})
        with urllib.request.urlopen(req) as resp, open(target_path, "wb") as f:
            f.write(resp.read())
        print(f"Successfully downloaded to {target_path} ({os.path.getsize(target_path)} bytes).")
        return target_path
    except Exception as e:
        raise RuntimeError(
            f"Failed to download matrix from {meta['url']}: {e}\n"
            f"Please verify internet connectivity or manually place '{meta['filename']}' in '{matrices_dir}'."
        )


def resolve_alignment_file(input_path: str) -> str:
    """Resolve alignment file from path or directory, prioritizing 3Di fasta."""
    if os.path.isdir(input_path):
        # Look for 3di fasta in directory
        candidates = [
            os.path.join(input_path, "foldmason.fasta_3di.fa"),
            os.path.join(input_path, "foldmason_3di.fa"),
            os.path.join(input_path, "alignment_3di.fa"),
        ]
        for c in candidates:
            if os.path.isfile(c):
                return c
        # Fallback to any .fa / .fasta
        for f in os.listdir(input_path):
            if f.endswith(("_3di.fa", "_3di.fasta", ".fasta", ".fa")):
                return os.path.join(input_path, f)
        raise FileNotFoundError(f"No suitable FASTA alignment found in directory '{input_path}'.")

    if os.path.isfile(input_path):
        # If user passed base fasta prefix (e.g. foldmason.fasta), check if _3di.fa exists
        cand_3di = input_path + "_3di.fa"
        if os.path.isfile(cand_3di):
            return cand_3di
        return input_path

    raise FileNotFoundError(f"Alignment path '{input_path}' not found.")


def parse_args():
    parser = argparse.ArgumentParser(
        description="Construct structural phylogenetic trees from FoldMason alignments using IQ-TREE or FoldMason."
    )
    parser.add_argument(
        "-a",
        "--alignment",
        default="foldmason_alignments/foldmason.fasta_3di.fa",
        help="Path to 3Di alignment file or FoldMason output directory (default: foldmason_alignments/foldmason.fasta_3di.fa)",
    )
    parser.add_argument(
        "-m",
        "--method",
        choices=["iqtree", "foldmason"],
        default="iqtree",
        help="Tree inference method: 'iqtree' (Maximum Likelihood with 3Di matrix) or 'foldmason' (guide tree) (default: iqtree)",
    )
    parser.add_argument(
        "--matrix",
        choices=["alphafold", "af", "esmfold", "llm", "both", "auto"],
        default="alphafold",
        help="3Di substitution matrix for IQ-TREE: 'alphafold' (Q.3Di.AF), 'esmfold' (Q.3Di.LLM), or 'both'/'auto' to let ModelFinder compare both (default: alphafold)",
    )
    parser.add_argument(
        "--rate-heterogeneity",
        default="auto",
        help="Rate heterogeneity model: '+G4', '+I+G4', '+R4', '' (none), or 'auto'/'MFP' to have ModelFinder automatically select the best model (default: auto)",
    )
    parser.add_argument(
        "--criterion",
        choices=["BIC", "AIC", "AICc"],
        default="BIC",
        help="Selection criterion for ModelFinder (default: BIC)",
    )
    parser.add_argument(
        "-b",
        "--bootstrap",
        type=int,
        default=1000,
        help="Ultrafast bootstrap replicates for IQ-TREE (e.g. 1000, or 0 to disable) (default: 1000)",
    )
    parser.add_argument(
        "--alrt",
        type=int,
        default=1000,
        help="SH-like approximate likelihood ratio test (SH-aLRT) replicates (e.g. 1000, or 0 to disable) (default: 1000)",
    )
    parser.add_argument(
        "-t",
        "--threads",
        default="AUTO",
        help="Number of CPU threads for IQ-TREE (e.g. 2, 4, or AUTO) (default: AUTO)",
    )
    parser.add_argument(
        "-o",
        "--output-dir",
        default="phylogeny_results",
        help="Output directory for phylogenetic trees and reports (default: phylogeny_results)",
    )
    parser.add_argument(
        "-p",
        "--prefix",
        default="structural_tree",
        help="Prefix for output tree files (default: structural_tree)",
    )
    parser.add_argument(
        "--iqtree-bin",
        default="iqtree",
        help="Name or path of the iqtree binary (default: iqtree)",
    )
    return parser.parse_args()


def run_foldmason_tree(alignment_file: str, output_dir: str, prefix: str):
    """Retrieve or copy the FoldMason guide tree."""
    os.makedirs(output_dir, exist_ok=True)
    out_tree = os.path.join(output_dir, f"{prefix}_foldmason.nwk")

    # FoldMason typically generates a .nw file alongside the alignment
    dir_path = os.path.dirname(alignment_file)
    cand_nw = os.path.join(dir_path, "foldmason.fasta.nw")
    if not os.path.isfile(cand_nw):
        base = os.path.splitext(alignment_file)[0]
        if base.endswith("_3di"):
            base = base[:-4]
        cand_nw = base + ".nw"

    if os.path.isfile(cand_nw):
        shutil.copyfile(cand_nw, out_tree)
        print(f"Retrieved FoldMason guide tree from '{cand_nw}'.")
        print(f"Saved tree to '{out_tree}'.")
        with open(out_tree) as f:
            print("\nNewick Tree:")
            print(f.read().strip())
    else:
        print(
            f"Error: Could not locate FoldMason guide tree (.nw file) near '{alignment_file}'.",
            file=sys.stderr,
        )
        sys.exit(1)


def run_iqtree(args, alignment_file: str):
    """Execute IQ-TREE with the designated 3Di substitution matrix."""
    matrix_path = ensure_matrix_file(args.matrix)
    os.makedirs(args.output_dir, exist_ok=True)
    output_prefix = os.path.join(args.output_dir, f"{args.prefix}_{args.matrix}")

    # Model definition
    is_auto = args.rate_heterogeneity.lower() in ("auto", "mfp", "mf", "test")
    if is_auto:
        model_str = f"ModelFinder (evaluating {matrix_path})"
        cmd = [
            args.iqtree_bin,
            "-s",
            alignment_file,
            "-mset",
            matrix_path,
            "-m",
            "MFP",
            "-merit",
            args.criterion,
            "-pre",
            output_prefix,
            "-nt",
            str(args.threads),
        ]
    else:
        rate_het = args.rate_heterogeneity if args.rate_heterogeneity else ""
        if rate_het and not rate_het.startswith("+"):
            rate_het = "+" + rate_het
        model_str = f"{matrix_path}{rate_het}"
        cmd = [
            args.iqtree_bin,
            "-s",
            alignment_file,
            "-m",
            model_str,
            "-pre",
            output_prefix,
            "-nt",
            str(args.threads),
        ]

    if args.bootstrap and args.bootstrap > 0:
        cmd.extend(["-B", str(args.bootstrap)])
    if args.alrt and args.alrt > 0:
        cmd.extend(["-alrt", str(args.alrt)])

    # Overwrite if exists
    cmd.append("-redo")

    print("=" * 70)
    print("RUNNING STRUCTURAL PHYLOGENY WITH IQ-TREE (3Di ALPHABET)")
    print(f"Alignment:          {alignment_file}")
    print(f"3Di Matrix:         {matrix_path} ({args.matrix})")
    print(f"Model Selection:    {model_str}")
    if is_auto:
        print(f"Selection Criterion:{args.criterion}")
    if args.bootstrap > 0:
        print(f"UFBoot replicates:  {args.bootstrap}")
    if args.alrt > 0:
        print(f"SH-aLRT replicates: {args.alrt}")
    print(f"Output Prefix:      {output_prefix}")
    print(f"Command:            {' '.join(cmd)}")
    print("=" * 70)

    # Check iqtree executable
    if not shutil.which(args.iqtree_bin):
        # Try conda env if not found in current PATH
        print(f"Note: '{args.iqtree_bin}' not found directly in current shell PATH.")
        print("Looking for conda environment 'skills_hackathon'...")
        conda_iqtree = os.path.expanduser("~/miniconda3/envs/skills_hackathon/bin/iqtree")
        if os.path.isfile(conda_iqtree):
            cmd[0] = conda_iqtree
        else:
            print(
                f"Error: IQ-TREE binary '{args.iqtree_bin}' not found. Please install via: conda install -c bioconda iqtree",
                file=sys.stderr,
            )
            sys.exit(127)

    try:
        subprocess.run(cmd, check=True)
    except subprocess.CalledProcessError as e:
        print(f"\nError: IQ-TREE execution failed with return code {e.returncode}.", file=sys.stderr)
        sys.exit(e.returncode)

    treefile = f"{output_prefix}.treefile"
    contree = f"{output_prefix}.contree"
    reportfile = f"{output_prefix}.iqtree"

    print("\n" + "=" * 70)
    print("PHYLOGENY INFERENCE COMPLETE")
    print("=" * 70)
    if os.path.isfile(reportfile):
        with open(reportfile) as f:
            for line in f:
                if "Best-fit model" in line:
                    print(f"Selected Model:              {line.strip()}")
                    break
        print(f"Detailed IQ-TREE report:     {reportfile}")
    if os.path.isfile(treefile):
        print(f"Maximum-likelihood tree:     {treefile}")
        with open(treefile) as f:
            nwk = f.read().strip()
            print(f"\nML Newick string:\n{nwk}\n")
    if os.path.isfile(contree):
        print(f"Consensus tree with support: {contree}")


def main():
    args = parse_args()
    alignment_file = resolve_alignment_file(args.alignment)
    print(f"Using alignment: {alignment_file}")

    if args.method == "foldmason":
        run_foldmason_tree(alignment_file, args.output_dir, args.prefix)
    else:
        run_iqtree(args, alignment_file)


if __name__ == "__main__":
    main()
