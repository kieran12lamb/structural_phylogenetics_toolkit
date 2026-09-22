"""Phylogenetic tree inference using IQ-TREE with empirical 3Di substitution matrices and FoldMason."""

import os
import shutil
import subprocess
from typing import Optional, Dict, Any

from viral_phylo.binaries import find_iqtree_bin
from viral_phylo.matrices import ensure_matrix_file
from viral_phylo.alignment import count_fasta_seqs

def build_tree(
    alignment_file: str,
    method: str = "iqtree",
    matrix: str = "both",
    rate_het: str = "auto",
    criterion: str = "BIC",
    bootstrap: int = 1000,
    alrt: int = 1000,
    threads: str = "AUTO",
    output_dir: str = "phylogeny_results",
    prefix: str = "viral_tree",
    iqtree_bin: str = None,
    tree_type: str = "both",
    alignment_aa: str = None,
    fast: bool = False,
):
    """Build structural and/or amino acid phylogeny using FoldMason or IQ-TREE."""
    os.makedirs(output_dir, exist_ok=True)

    # 0. Path Resolution & Sanitization
    if alignment_file.startswith("Users/") or alignment_file.startswith("home/"):
        abs_cand = "/" + alignment_file
        if os.path.exists(abs_cand):
            print(f"[Path] Auto-corrected missing leading slash: '{alignment_file}' -> '{abs_cand}'")
            alignment_file = abs_cand

    if os.path.isdir(alignment_file):
        candidates = [
            os.path.join(alignment_file, "mafft.fasta_3di.fa"),
            os.path.join(alignment_file, "foldmason.fasta_3di.fa"),
            os.path.join(alignment_file, "foldmason_3di.fa"),
        ]
        found = None
        for c in candidates:
            if os.path.isfile(c):
                found = c
                break
        if not found:
            import glob
            cands = glob.glob(os.path.join(alignment_file, "*3di*.fa*"))
            if cands:
                found = cands[0]
        if found:
            print(f"[Path] Auto-resolved directory to 3Di alignment: '{found}'")
            alignment_file = found
        else:
            raise FileNotFoundError(f"Alignment directory '{alignment_file}' does not contain 'mafft.fasta_3di.fa', 'foldmason.fasta_3di.fa', or any 3Di alignment file.")

    if not os.path.isfile(alignment_file):
        raise FileNotFoundError(f"Alignment file '{alignment_file}' not found. Please provide a path to a FASTA alignment file (e.g., 300_rdrp/alignment/foldmason.fasta_3di.fa).")

    if alignment_aa and (alignment_aa.startswith("Users/") or alignment_aa.startswith("home/")):
        abs_cand = "/" + alignment_aa
        if os.path.exists(abs_cand):
            alignment_aa = abs_cand

    if method == "foldmason":
        # Extract FoldMason guide tree
        dir_path = os.path.dirname(alignment_file)
        cand_nw = os.path.join(dir_path, "foldmason.fasta.nw")
        out_tree = os.path.join(output_dir, f"{prefix}_foldmason.nwk")
        if os.path.isfile(cand_nw):
            shutil.copyfile(cand_nw, out_tree)
            print(f"[Phylogeny] Saved FoldMason guide tree to: {out_tree}")
            return out_tree
        raise FileNotFoundError(f"FoldMason guide tree not found at '{cand_nw}'.")

    # IQ-TREE
    if not iqtree_bin:
        iqtree_bin = find_iqtree_bin()

    seq_count = count_fasta_seqs(alignment_file)
    is_large = seq_count >= 80

    # IQ-TREE strictly forbids combining '--fast' with '-B' (Ultrafast bootstrap) or '-alrt'
    # If user explicitly passed --fast, prioritize --fast and omit bootstrap/alrt.
    # Otherwise, if bootstrap is requested, omit --fast to preserve bootstrap support.
    if fast:
        use_fast = True
        if bootstrap > 0 or alrt > 0:
            print("[IQ-TREE] Notice: '--fast' mode was specified. Disabling ultrafast bootstrap (-B) and SH-aLRT (-alrt) as IQ-TREE does not support bootstrapping in fast mode.")
            bootstrap = 0
            alrt = 0
    else:
        # If user did not specify --fast, do not auto-enable --fast if bootstrapping is active
        use_fast = is_large and (bootstrap <= 0 and alrt <= 0)
        if is_large and (bootstrap > 0 or alrt > 0):
            print(f"[IQ-TREE] Large dataset ({seq_count} taxa) with bootstrap/aLRT enabled: using standard tree search with model optimization (omitting '--fast' for bootstrap compatibility).")

    results = {}

    # 1. Structural 3Di Tree
    if tree_type in ("3di", "both", "tanglegram"):
        matrix_path = ensure_matrix_file(matrix)
        out_prefix = os.path.join(output_dir, f"{prefix}_{matrix}")

        cmd = [iqtree_bin, "-s", alignment_file, "-pre", out_prefix, "-nt", str(threads), "-redo"]
        if use_fast:
            cmd.append("--fast")

        is_auto = rate_het.lower() in ("auto", "mfp", "mf", "test")
        if is_auto:
            if is_large or use_fast:
                first_mat = matrix_path.split(",")[0]
                cmd.extend(["-m", f"{first_mat}+G4"])
                print(f"[IQ-TREE 3Di] Large dataset ({seq_count} taxa): using verified model '{first_mat}+G4'.")
            else:
                cmd.extend(["-mset", matrix_path, "-m", "MFP", "-merit", criterion])
        else:
            if rate_het and not rate_het.startswith("+"):
                rate_het = "+" + rate_het
            cmd.extend(["-m", f"{matrix_path}{rate_het}"])


        if bootstrap > 0:
            cmd.extend(["-B", str(bootstrap)])
        if alrt > 0:
            cmd.extend(["-alrt", str(alrt)])

        print(f"\n[IQ-TREE 3Di] Running structural phylogeny inference:")
        print("  " + " ".join(cmd))
        subprocess.run(cmd, check=True)

        treefile = f"{out_prefix}.treefile"
        results["3di"] = treefile
        print(f"[IQ-TREE 3Di] Tree saved to: {treefile}")
        if os.path.isfile(treefile):
            with open(treefile) as f:
                print(f"3Di Newick: {f.read().strip()}")

    # 2. Amino Acid Sequence Tree
    if tree_type in ("aa", "both", "tanglegram"):
        # Auto-resolve AA alignment if not given
        if not alignment_aa:
            dir_p = os.path.dirname(alignment_file)
            cand_aa = [
                os.path.join(dir_p, "mafft.fasta_aa.fa"),
                os.path.join(dir_p, "foldmason.fasta_aa.fa"),
                os.path.join(dir_p, "foldmason_aa.fa"),
                alignment_file.replace("_3di.fa", "_aa.fa"),
            ]
            for c in cand_aa:
                if os.path.isfile(c):
                    alignment_aa = c
                    break

        if alignment_aa and os.path.isfile(alignment_aa):
            aa_prefix = os.path.join(output_dir, f"{prefix}_aa")
            cmd_aa = [
                iqtree_bin,
                "-s",
                alignment_aa,
                "-pre",
                aa_prefix,
                "-nt",
                str(threads),
                "-redo",
            ]
            if use_fast:
                cmd_aa.append("--fast")

            if is_large or use_fast:
                cmd_aa.extend(["-m", "LG+G4"])
                print(f"[IQ-TREE AA] Large dataset ({seq_count} taxa): using verified model 'LG+G4'.")
            else:
                cmd_aa.extend(["-m", "MFP"])

            if bootstrap > 0:
                cmd_aa.extend(["-B", str(bootstrap)])
            if alrt > 0:
                cmd_aa.extend(["-alrt", str(alrt)])

            print(f"\n[IQ-TREE AA] Running amino acid sequence phylogeny inference:")
            print("  " + " ".join(cmd_aa))
            subprocess.run(cmd_aa, check=True)

            aa_treefile = f"{aa_prefix}.treefile"
            results["aa"] = aa_treefile
            print(f"[IQ-TREE AA] Tree saved to: {aa_treefile}")
            if os.path.isfile(aa_treefile):
                with open(aa_treefile) as f:
                    print(f"AA Newick: {f.read().strip()}")
        else:
            print(f"[Warning] AA alignment not found. Skipping amino acid tree inference.")

    return results.get("3di") or results.get("aa")
