"""Phylogenetic tree inference using IQ-TREE with empirical 3Di substitution matrices and FoldMason."""

import json
import os
import shutil
import subprocess
from typing import Optional, Dict, Any

from viral_phylo.binaries import find_iqtree_bin, find_fasttree_bin
from viral_phylo.matrices import ensure_matrix_file, paml_to_fasttree_trans, uppercase_matrix_alias
from viral_phylo.alignment import count_fasta_seqs

def resolve_aa_alignment(alignment_file: str, alignment_aa: str = None):
    """Find the amino-acid alignment that accompanies a 3Di alignment.

    Used by every tree method, so that ``--tree-type both`` builds an AA tree
    whichever method is chosen. Previously only the IQ-TREE path did this lookup;
    the FastTree path relied on an explicitly passed AA file, which the pipeline
    never passes, so the AA tree was silently skipped.
    """
    if alignment_aa and os.path.isfile(alignment_aa):
        return alignment_aa
    dir_p = os.path.dirname(alignment_file)
    for candidate in (
        os.path.join(dir_p, "mafft.fasta_aa.fa"),
        os.path.join(dir_p, "foldmason.fasta_aa.fa"),
        os.path.join(dir_p, "foldmason_aa.fa"),
        alignment_file.replace("_3di.fa", "_aa.fa"),
    ):
        if candidate != alignment_file and os.path.isfile(candidate):
            return candidate
    return None


def _run_iqtree_maybe_guided(cmd, alignment, out_prefix, alphabet, guide, matrix_path=None,
                             threads="8", seed=1):
    """Run an IQ-TREE command, optionally from a FastTree starting tree.

    Records what happened in ``<out_prefix>.start_tree.json`` whenever guidance was
    requested: whether a starting tree was used, how many polytomies were resolved,
    and - if the guided run failed - that it fell back to an unguided run and why.
    """
    if not guide:
        subprocess.run(cmd, check=True)
        return None

    from viral_phylo.start_trees import make_start_tree

    record = make_start_tree(alignment, out_prefix + ".start_tree.nwk", alphabet,
                             matrix_path=matrix_path, threads=threads, seed=seed)
    if record["used"]:
        # -keep-ident: IQ-TREE otherwise drops duplicate sequences and then rejects a
        # starting tree that still contains them.
        guided = cmd + ["-t", record["path"], "-keep-ident"]
        print(f"[IQ-TREE {alphabet}] Starting from a FastTree tree "
              f"({record['polytomies_resolved']} polytomies resolved):")
        print("  " + " ".join(guided))
        try:
            subprocess.run(guided, check=True)
            record["outcome"] = "guided"
            record["iqtree_command"] = " ".join(guided)
        except subprocess.CalledProcessError as exc:
            record["outcome"] = "fell_back"
            record["fallback_reason"] = f"guided IQ-TREE exited {exc.returncode}"
            print(f"[IQ-TREE {alphabet}] [!] Notice: guided run failed (exit {exc.returncode}); "
                  "rerunning without the starting tree.")
            subprocess.run(cmd, check=True)
            record["iqtree_command"] = " ".join(cmd)
    else:
        record["outcome"] = "not_used"
        print(f"[IQ-TREE {alphabet}] [!] Notice: no starting tree used ({record['reason']}); running unguided.")
        subprocess.run(cmd, check=True)
        record["iqtree_command"] = " ".join(cmd)

    import json as _json
    with open(out_prefix + ".start_tree.json", "w", encoding="utf-8") as handle:
        _json.dump(record, handle, indent=2)
    return record


def build_tree(
    alignment_file: str,
    method: str = "iqtree",
    matrix: str = "both",
    rate_het: str = "auto",
    criterion: str = "BIC",
    bootstrap: int = 1000,
    alrt: int = 1000,
    threads: str = "8",
    output_dir: str = "phylogeny_results",
    prefix: str = "viral_tree",
    iqtree_bin: str = None,
    tree_type: str = "both",
    alignment_aa: str = None,
    fast: bool = False,
    guide_trees: bool = False,
    guide_min_taxa: int = 80,
    guide_seed: int = 1,
):
    """Build structural and/or amino acid phylogeny using FoldMason or IQ-TREE.

    With ``guide_trees``, IQ-TREE runs on alignments of at least ``guide_min_taxa``
    taxa start from a repaired and validated FastTree tree (see start_trees.py). A
    guided run that fails is rerun unguided, so this can never lose a tree.
    """
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

    if method == "fasttree":
        # Approximate-maximum-likelihood tree. Orders of magnitude faster than IQ-TREE
        # on large alignments: the 1,193-taxon whole-set tree that IQ-TREE was ~78
        # hours into runs here in about a minute.
        #
        # The 3Di tree uses the requested 3Di matrix via FastTree's -trans option, after
        # converting it from PAML format (see matrices.paml_to_fasttree_trans). On the
        # 1,193-taxon whole-set alignment that raised the log-likelihood from -201,283
        # under LG to -151,272, for ~10 s more. FastTree cannot choose between several
        # matrices, so with --matrix both the first is used and this is recorded.
        fasttree_bin = find_fasttree_bin(getattr(build_tree, "_fasttree_bin", None))
        if not fasttree_bin:
            raise FileNotFoundError(
                "No FastTree-compatible binary found (looked for veryfasttree, VeryFastTree, "
                "FastTreeMP, FastTree, fasttree). Install with: "
                "conda install -c bioconda veryfasttree")
        os.makedirs(output_dir, exist_ok=True)
        results = {}
        targets = []
        # Named for the method so a FastTree tree is never mistaken for an IQ-TREE one;
        # the "3di"/"aa" tokens keep the dashboard's tree-type detection working.
        if tree_type in ("3di", "both", "tanglegram"):
            targets.append(("3di", alignment_file, os.path.join(output_dir, f"{prefix}_3di_fasttree.treefile")))
        if tree_type in ("aa", "both", "tanglegram"):
            alignment_aa = resolve_aa_alignment(alignment_file, alignment_aa)
            if alignment_aa:
                targets.append(("aa", alignment_aa, os.path.join(output_dir, f"{prefix}_aa_fasttree.treefile")))
            else:
                print(f"[FastTree aa] [!] Notice: tree type '{tree_type}' asks for an AA tree, but no "
                      f"AA alignment was found next to '{alignment_file}'; skipping it.")

        # Resolve the 3Di model once: convert the requested matrix for -trans.
        trans_file, trans_source, trans_note = None, None, None
        if any(t[0] == "3di" for t in targets):
            try:
                resolved = ensure_matrix_file(matrix).split(",")
                trans_source = resolved[0].strip()
                if len(resolved) > 1:
                    trans_note = (f"--matrix {matrix} names {len(resolved)} matrices; FastTree cannot select "
                                  f"between them, so {os.path.basename(trans_source)} was used")
                    print(f"[FastTree 3di] [!] Notice: {trans_note}.")
                trans_file = paml_to_fasttree_trans(
                    trans_source, os.path.join(output_dir, os.path.basename(trans_source) + ".fasttree"))
            except Exception as exc:
                trans_note = f"could not use matrix '{matrix}' ({type(exc).__name__}: {exc}); fell back to LG"
                print(f"[FastTree 3di] [!] Notice: {trans_note}.")
                trans_file = None
        is_very = "veryfasttree" in os.path.basename(fasttree_bin).lower()

        for label, aln, out_tree in targets:
            if label == "3di" and trans_file:
                cmd = [fasttree_bin, "-trans", trans_file, "-gamma"]
                model = f"{os.path.basename(trans_source)} (converted for -trans) + Gamma"
            else:
                cmd = [fasttree_bin, "-lg", "-gamma"]
                model = "LG + Gamma (-lg -gamma)"
            if is_very:
                # VeryFastTree defaults to single precision, which cost ~0.75 log-lik units
                # on a 32-taxon tree and ~280 on the 1,193-taxon one; double precision
                # matches IQ-TREE exactly for ~10% more time.
                cmd.append("-double-precision")
                if str(threads).isdigit():
                    cmd += ["-threads", str(threads)]
            cmd += [aln]
            print(f"[FastTree {label}] Inferring approximate-ML tree ({count_fasta_seqs(aln)} taxa):")
            print("  " + " ".join(cmd))
            with open(out_tree, "w", encoding="utf-8") as handle:
                subprocess.run(cmd, stdout=handle, check=True)
            if not os.path.isfile(out_tree) or os.path.getsize(out_tree) == 0:
                raise FileNotFoundError(f"FastTree produced no tree at '{out_tree}'.")
            print(f"[FastTree {label}] Tree saved to: {out_tree}")
            # Sidecar recording exactly what produced the tree, since the Newick file
            # itself carries no model or program information.
            try:
                from viral_phylo.binaries import tool_version
                with open(out_tree + ".info.json", "w", encoding="utf-8") as handle:
                    json.dump({
                        "method": "fasttree", "binary": fasttree_bin,
                        "version": tool_version("fasttree"),
                        "command": " ".join(cmd), "alignment": aln, "alphabet": label,
                        "model": model,
                        "matrix_source": trans_source if (label == "3di" and trans_file) else None,
                        "support": "SH-like local support (single value per branch, 0-1)",
                        "note": trans_note if label == "3di" else None,
                    }, handle, indent=2)
            except Exception as exc:
                print(f"[FastTree {label}] [!] Notice: could not write tree provenance: {exc}")
            results[label] = out_tree
        return results.get("3di") or results.get("aa")

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
    use_guide = bool(guide_trees) and seq_count >= guide_min_taxa and not fast

    # IQ-TREE refuses ultrafast bootstrap and SH-aLRT below 4 sequences ("It makes no
    # sense to perform bootstrap with less than 4 sequences") and exits non-zero, which
    # previously surfaced as an unexplained failure for every small cluster. Drop the
    # support values and still build the tree.
    if seq_count < 4 and (bootstrap > 0 or alrt > 0):
        print(f"[IQ-TREE] Notice: {seq_count} taxa - disabling bootstrap (-B) and SH-aLRT "
              "(-alrt), which IQ-TREE requires at least 4 sequences for.")
        bootstrap = 0
        alrt = 0

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
                # IQ-TREE upper-cases -mset, so hand it an upper-case path that
                # survives that transformation (see uppercase_matrix_alias).
                cmd.extend(["-mset", uppercase_matrix_alias(matrix_path),
                            "-m", "MFP", "-merit", criterion])
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
        _run_iqtree_maybe_guided(cmd, alignment_file, out_prefix, "3di", use_guide,
                                 matrix_path=matrix_path.split(",")[0], threads=threads, seed=guide_seed)

        treefile = f"{out_prefix}.treefile"
        results["3di"] = treefile
        print(f"[IQ-TREE 3Di] Tree saved to: {treefile}")
        if os.path.isfile(treefile):
            with open(treefile) as f:
                print(f"3Di Newick: {f.read().strip()}")

    # 2. Amino Acid Sequence Tree
    if tree_type in ("aa", "both", "tanglegram"):
        alignment_aa = resolve_aa_alignment(alignment_file, alignment_aa)

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
            _run_iqtree_maybe_guided(cmd_aa, alignment_aa, aa_prefix, "aa", use_guide,
                                     threads=threads, seed=guide_seed)

            aa_treefile = f"{aa_prefix}.treefile"
            results["aa"] = aa_treefile
            print(f"[IQ-TREE AA] Tree saved to: {aa_treefile}")
            if os.path.isfile(aa_treefile):
                with open(aa_treefile) as f:
                    print(f"AA Newick: {f.read().strip()}")
        else:
            print(f"[Warning] AA alignment not found. Skipping amino acid tree inference.")

    return results.get("3di") or results.get("aa")
