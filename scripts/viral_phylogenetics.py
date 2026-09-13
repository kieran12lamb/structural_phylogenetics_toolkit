#!/usr/bin/env python3
"""Unified CLI tool for the Viral Structural Phylogenetics Skill.

Subcommands:
  fetch     - Download viral protein structures (PDB) from Viro3D.
  align     - Run FoldMason structural multiple alignment (MSTA) to generate 3Di alignments.
  tree      - Infer structural phylogenetic tree with IQ-TREE using 3Di substitution matrices.
  pipeline  - Run full end-to-end workflow: fetch -> align -> tree.
"""

import argparse
import glob
import json
import os
import shutil
import subprocess
import sys
import urllib.parse
import urllib.request
from concurrent.futures import ThreadPoolExecutor, as_completed

try:
    from scripts.metadata_handler import parse_metadata
except ImportError:
    try:
        from metadata_handler import parse_metadata
    except ImportError:
        parse_metadata = None

BASE_URL = "https://viro3d.cvr.gla.ac.uk/api"


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


def find_iqtree_bin():
    """Find iqtree binary in PATH or conda environments."""
    if shutil.which("iqtree"):
        return "iqtree"
    conda_prefix = os.environ.get("CONDA_PREFIX", "")
    candidates = [
        os.path.join(conda_prefix, "bin/iqtree") if conda_prefix else "",
        os.path.expanduser("~/miniconda3/envs/spt/bin/iqtree"),
        os.path.expanduser("~/miniconda3/bin/iqtree"),
        "/opt/homebrew/bin/iqtree",
    ]
    for c in candidates:
        if c and os.path.isfile(c) and os.access(c, os.X_OK):
            return c
    return "iqtree"


def find_foldmason_bin():
    """Find foldmason binary in PATH or conda environments."""
    if shutil.which("foldmason"):
        return "foldmason"
    conda_prefix = os.environ.get("CONDA_PREFIX", "")
    candidates = [
        os.path.join(conda_prefix, "bin/foldmason") if conda_prefix else "",
        os.path.expanduser("~/miniconda3/envs/spt/bin/foldmason"),
        os.path.expanduser("~/miniconda3/bin/foldmason"),
        "/opt/homebrew/bin/foldmason",
    ]
    for c in candidates:
        if c and os.path.isfile(c) and os.access(c, os.X_OK):
            return c
    return "foldmason"


def find_mafft_bin(mafft_bin: str = None) -> str:
    """Find mafft binary in specified path, PATH, or conda environments."""
    if mafft_bin and os.path.isfile(mafft_bin) and os.access(mafft_bin, os.X_OK):
        return mafft_bin
    if shutil.which("mafft"):
        return "mafft"
    conda_prefix = os.environ.get("CONDA_PREFIX", "")
    candidates = [
        os.path.join(conda_prefix, "bin/mafft") if conda_prefix else "",
        os.path.expanduser("~/miniconda3/envs/spt/bin/mafft"),
        os.path.expanduser("~/miniconda3/bin/mafft"),
        "/opt/homebrew/bin/mafft",
        "/usr/local/bin/mafft",
    ]
    for c in candidates:
        if c and os.path.isfile(c) and os.access(c, os.X_OK):
            return c
    return "mafft"


def ensure_3di_matrix(matrix_path: str = "matrices/mat3di.out") -> str:
    """Ensure the Foldseek 3Di substitution matrix (mat3di.out) is present locally."""
    if os.path.isfile(matrix_path) and os.path.getsize(matrix_path) > 100:
        return matrix_path

    os.makedirs(os.path.dirname(os.path.abspath(matrix_path)), exist_ok=True)
    url = "https://raw.githubusercontent.com/steineggerlab/foldseek/master/data/mat3di.out"
    print(f"Downloading Foldseek 3Di substitution matrix (mat3di.out) from {url}...")
    req = urllib.request.Request(url, headers={"User-Agent": "Mozilla/5.0"})
    with urllib.request.urlopen(req) as resp, open(matrix_path, "wb") as f:
        f.write(resp.read())
    return matrix_path


def ensure_matrix_file(matrix_key: str, matrices_dir: str = "matrices") -> str:
    """Ensure substitution matrices are locally present, downloading if necessary."""
    norm = matrix_key.lower().strip()
    if norm in ("both", "auto"):
        af = ensure_matrix_file("alphafold", matrices_dir)
        llm = ensure_matrix_file("esmfold", matrices_dir)
        return f"{af},{llm}"
    elif norm in ("af", "alphafold", "q.3di.af"):
        meta = EDMOND_MATRICES["alphafold"]
    elif norm in ("llm", "esm", "esmfold", "q.3di.llm"):
        meta = EDMOND_MATRICES["esmfold"]
    else:
        if os.path.isfile(matrix_key):
            return matrix_key
        raise FileNotFoundError(f"Matrix file '{matrix_key}' not found.")

    os.makedirs(matrices_dir, exist_ok=True)
    target = os.path.join(matrices_dir, meta["filename"])
    if os.path.isfile(target) and os.path.getsize(target) > 0:
        return target

    print(f"Downloading {meta['desc']} from Edmond...")
    req = urllib.request.Request(meta["url"], headers={"User-Agent": "Mozilla/5.0"})
    with urllib.request.urlopen(req) as resp, open(target, "wb") as f:
        f.write(resp.read())
    return target


def fetch_structures(qualifier: str, max_sequences: int, output_dir: str, workers: int = 16):
    """Download viral structures from Viro3D with concurrent worker threads."""
    os.makedirs(output_dir, exist_ok=True)
    print(f"\n[Viro3D] Querying API for '{qualifier}' (target: {max_sequences} structures)...")
    
    url = f"{BASE_URL}/proteins/protein_name/"
    query = urllib.parse.urlencode({"qualifier": qualifier, "page_size": max_sequences, "page_num": 1})
    req = urllib.request.Request(f"{url}?{query}", headers={"User-Agent": "Mozilla/5.0"})
    with urllib.request.urlopen(req) as resp:
        data = json.loads(resp.read().decode())

    items = data.get("protein_structures", [])[:max_sequences]
    if not items:
        print(f"No structures found for '{qualifier}'.")
        return []

    # Save metadata JSON for downstream visualizers and tree annotations
    meta_path = os.path.join(output_dir, "taxa_metadata.json")
    with open(meta_path, "w", encoding="utf-8") as f:
        json.dump(items, f, indent=2)
    print(f"[Viro3D] Saved metadata index for {len(items)} structures to '{meta_path}'.")

    downloaded = []
    print(f"[Viro3D] Downloading {len(items)} PDB structures concurrently ({workers} workers)...")

    def _dl_one(item):
        rid = item["record_id"]
        pdb_url = f"{BASE_URL}/pdb/CF-{rid}_relaxed.pdb"
        out_file = os.path.join(output_dir, f"{rid}.pdb")
        if os.path.isfile(out_file) and os.path.getsize(out_file) > 500:
            return rid, out_file, True
        try:
            p_req = urllib.request.Request(pdb_url, headers={"User-Agent": "Mozilla/5.0"})
            with urllib.request.urlopen(p_req, timeout=30) as r, open(out_file, "wb") as f:
                f.write(r.read())
            return rid, out_file, True
        except Exception as e:
            return rid, str(e), False

    with ThreadPoolExecutor(max_workers=workers) as pool:
        futures = [pool.submit(_dl_one, it) for it in items]
        completed = 0
        for fut in as_completed(futures):
            rid, path_or_err, success = fut.result()
            completed += 1
            if success:
                downloaded.append(path_or_err)
            if completed % 50 == 0 or completed == len(items):
                print(f"  Progress: {completed}/{len(items)} downloaded ({len(downloaded)} successful)")

    print(f"[Viro3D] Download complete: {len(downloaded)}/{len(items)} structure(s) saved to '{output_dir}'.\n")
    return downloaded



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

        # Output alignment paths (matching foldmason conventions for downstream compatibility)
        aln_3di = os.path.join(output_dir, "foldmason.fasta_3di.fa")
        aln_aa = os.path.join(output_dir, "foldmason.fasta_aa.fa")

        print(f"[MAFFT] Aligning 3Di structural sequences with substitution matrix '{actual_mat3di}'...")
        with open(aln_3di, "w") as f_out:
            subprocess.run([actual_mafft, "--aamatrix", actual_mat3di, "--auto", unaligned_3di], stdout=f_out, check=True)

        print(f"[MAFFT] Aligning amino acid sequences...")
        with open(aln_aa, "w") as f_out:
            subprocess.run([actual_mafft, "--auto", unaligned_aa], stdout=f_out, check=True)

        if not os.path.isfile(aln_3di) or os.path.getsize(aln_3di) == 0:
            raise FileNotFoundError(f"Expected MAFFT 3Di alignment '{aln_3di}' was not generated.")
        print(f"[MAFFT] Alignment complete: generated '{aln_3di}' and '{aln_aa}'.\n")
        return aln_3di

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
        if not os.path.isfile(aln_3di):
            raise FileNotFoundError(f"Expected 3Di alignment '{aln_3di}' was not generated.")
        print(f"[FoldMason] Alignment complete: generated '{aln_3di}'.\n")
        return aln_3di


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
            raise FileNotFoundError(f"Alignment directory '{alignment_file}' does not contain 'foldmason.fasta_3di.fa' or any 3Di alignment file.")

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



def find_torch_python() -> str:
    """Find a Python interpreter to use for ML/embeddings.

    Priority:
    1. The install script environment 'spt' (e.g. envs/spt/bin/python)
    2. The user's active Conda environment ($CONDA_PREFIX/bin/python)
    3. The currently executing Python interpreter (sys.executable)
    4. Fallback candidate environments
    """
    def is_valid_python(path: str) -> bool:
        return bool(path and os.path.isfile(path) and os.access(path, os.X_OK))

    # 1. Primary default: The install script environment 'spt'
    spt_candidates = []
    if os.environ.get("CONDA_DEFAULT_ENV") == "spt" and os.environ.get("CONDA_PREFIX"):
        spt_candidates.append(os.path.join(os.environ["CONDA_PREFIX"], "bin", "python"))

    conda_base_hints = [
        os.environ.get("CONDA_PREFIX", ""),
        os.path.expanduser("~/miniconda3"),
        os.path.expanduser("~/anaconda3"),
        os.path.expanduser("~/miniforge3"),
        os.path.expanduser("~/micromamba"),
        "/opt/conda",
        "/opt/homebrew/Caskroom/miniconda/base"
    ]
    for hint in conda_base_hints:
        if not hint:
            continue
        if "/envs/" in hint:
            base_dir = hint.split("/envs/")[0]
            spt_candidates.append(os.path.join(base_dir, "envs", "spt", "bin", "python"))
        else:
            spt_candidates.append(os.path.join(hint, "envs", "spt", "bin", "python"))

    for cand in spt_candidates:
        if is_valid_python(cand):
            return cand

    # 2. Secondary default: Assume the user's active Conda environment
    if os.environ.get("CONDA_PREFIX"):
        active_py = os.path.join(os.environ["CONDA_PREFIX"], "bin", "python")
        if is_valid_python(active_py):
            return active_py

    # 3. Third priority: Currently executing Python interpreter
    try:
        import torch  # noqa: F401
        import transformers  # noqa: F401
        import scipy  # noqa: F401
        return sys.executable
    except ImportError:
        pass

    # 4. Fallback candidate environments
    fallbacks = [
        os.path.expanduser("~/miniconda3/envs/nipah/bin/python"),
        os.path.expanduser("~/miniconda3/bin/python"),
        "/usr/local/bin/python3",
        sys.executable,
    ]
    for cand in fallbacks:
        if is_valid_python(cand):
            return cand

    return sys.executable


def build_cli_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Viral Structural Phylogenetics & PLM Embedding Suite.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    subparsers = parser.add_subparsers(dest="subcommand", help="Subcommand to run")

    # Subcommand: fetch
    p_fetch = subparsers.add_parser("fetch", help="Download viral PDB structures from Viro3D.")
    p_fetch.add_argument("-q", "--qualifier", default="glycoprotein", help="Protein search term (default: glycoprotein)")
    p_fetch.add_argument("-m", "--max-sequences", type=int, default=6, help="Max structures to download (default: 6)")
    p_fetch.add_argument("-o", "--output-dir", default="viro_3d_structures", help="Output directory for PDBs")

    # Subcommand: align
    p_align = subparsers.add_parser("align", help="Align PDB structures with FoldMason or MAFFT to generate 3Di alignments.")
    p_align.add_argument("-i", "--folder", default="viro_3d_structures", help="Folder containing PDB files (local or fetched)")
    p_align.add_argument("-o", "--output-dir", default="foldmason_alignments", help="Output directory for MSA")
    p_align.add_argument("--aligner", choices=["foldmason", "mafft"], default="foldmason", help="Multiple sequence alignment engine (default: foldmason)")
    p_align.add_argument("--foldmason-bin", default=None, help="Path to foldmason binary")
    p_align.add_argument("--mafft-bin", default=None, help="Path to mafft binary (if using --aligner mafft)")
    p_align.add_argument("--mafft-matrix", default="matrices/mat3di.out", help="3Di substitution matrix for MAFFT (default: matrices/mat3di.out)")

    # Subcommand: tree
    p_tree = subparsers.add_parser("tree", help="Build structural and/or amino acid phylogeny using IQ-TREE or FoldMason.")
    p_tree.add_argument("-a", "--alignment", default="foldmason_alignments/foldmason.fasta_3di.fa", help="3Di alignment file or directory")
    p_tree.add_argument("--alignment-aa", default=None, help="Amino acid alignment file (default: auto-detected)")
    p_tree.add_argument("--tree-type", choices=["3di", "aa", "both", "tanglegram"], default="both", help="Phylogeny type: '3di' (structural), 'aa' (sequence), or 'both'/'tanglegram' (default: both)")
    p_tree.add_argument("-m", "--method", choices=["iqtree", "foldmason"], default="iqtree", help="Tree method (default: iqtree)")
    p_tree.add_argument("--matrix", choices=["alphafold", "af", "esmfold", "llm", "both", "auto"], default="both", help="3Di substitution matrix")
    p_tree.add_argument("--rate-heterogeneity", default="auto", help="Rate heterogeneity: '+G4', '+I+G4', '+R', or 'auto' (default: auto)")
    p_tree.add_argument("--criterion", choices=["BIC", "AIC", "AICc"], default="BIC", help="Model selection criterion (default: BIC)")
    p_tree.add_argument("-b", "--bootstrap", type=int, default=1000, help="Ultrafast bootstrap replicates (default: 1000)")
    p_tree.add_argument("--alrt", type=int, default=1000, help="SH-aLRT replicates (default: 1000)")
    p_tree.add_argument("-t", "--threads", default="AUTO", help="CPU threads (default: AUTO)")
    p_tree.add_argument("--fast", action="store_true", help="Enable fast heuristic search mode (auto-enabled for large datasets)")
    p_tree.add_argument("-o", "--output-dir", default="phylogeny_results", help="Output directory for trees")
    p_tree.add_argument("-p", "--prefix", default="viral_tree", help="Tree output filename prefix")
    p_tree.add_argument("--iqtree-bin", default=None, help="Path to iqtree binary")

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

    # Subcommand: pipeline (all-in-one)
    p_pipe = subparsers.add_parser("pipeline", help="Run full pipeline: fetch/local -> align -> tree (+ optional PLM embed).")
    p_pipe.add_argument("-i", "--input-folder", default=None, help="Local directory of PDB structures (skips Viro3D fetch if specified)")
    p_pipe.add_argument("-q", "--qualifier", default="glycoprotein", help="Protein search term if fetching (default: glycoprotein)")
    p_pipe.add_argument("-c", "--count", type=int, default=6, help="Number of structures to fetch (default: 6)")
    p_pipe.add_argument("--aligner", choices=["foldmason", "mafft"], default="foldmason", help="Multiple sequence alignment engine (default: foldmason)")
    p_pipe.add_argument("--mafft-bin", default=None, help="Path to mafft binary (if using --aligner mafft)")
    p_pipe.add_argument("--mafft-matrix", default="matrices/mat3di.out", help="3Di substitution matrix for MAFFT (default: matrices/mat3di.out)")
    p_pipe.add_argument("--tree-type", choices=["3di", "aa", "both", "tanglegram"], default="both", help="Phylogeny type (default: both)")
    p_pipe.add_argument("-m", "--method", choices=["iqtree", "foldmason"], default="iqtree", help="Tree method (default: iqtree)")
    p_pipe.add_argument("--matrix", choices=["alphafold", "af", "esmfold", "llm", "both", "auto"], default="both", help="3Di matrix (default: both)")
    p_pipe.add_argument("--rate-heterogeneity", default="auto", help="Rate heterogeneity model (default: auto)")
    p_pipe.add_argument("-b", "--bootstrap", type=int, default=1000, help="Ultrafast bootstrap replicates (default: 1000)")
    p_pipe.add_argument("--alrt", type=int, default=1000, help="SH-aLRT replicates (default: 1000)")
    p_pipe.add_argument("-t", "--threads", default="AUTO", help="CPU threads (default: AUTO)")
    p_pipe.add_argument("-meta", "--metadata", default=None, help="Path to metadata file (.xlsx, .csv, .tsv, .json). If omitted, inferred from structures or Viro3D API.")
    p_pipe.add_argument("--fast", action="store_true", help="Enable fast search mode")
    p_pipe.add_argument("--embed", action="store_true", help="Extract PLM embeddings and build hierarchical clustering tree (ESM-2 / ESM-C)")
    p_pipe.add_argument("--embed-model", default="esm2", choices=["esm2", "esmc"], help="PLM model for embeddings (default: esm2)")
    p_pipe.add_argument("--embed-clustering", default="upgma", choices=["upgma", "nj", "average", "complete"], help="Clustering method for PLM tree (default: upgma)")
    p_pipe.add_argument("--embed-metric", default="cosine", choices=["cosine", "euclidean", "l1", "cityblock", "manhattan"], help="Pairwise distance metric (default: cosine)")
    p_pipe.add_argument("-o", "--output-dir", default="glycoprotein_workflow", help="Parent output directory")

    return parser


def main():
    parser = build_cli_parser()
    args = parser.parse_args()

    if args.subcommand == "fetch":
        fetch_structures(args.qualifier, args.max_sequences, args.output_dir)
    elif args.subcommand == "align":
        align_structures(
            args.folder,
            args.output_dir,
            foldmason_bin=args.foldmason_bin,
            aligner=args.aligner,
            mafft_bin=args.mafft_bin,
            mafft_matrix=args.mafft_matrix,
        )
    elif args.subcommand == "tree":
        aln_file = args.alignment
        if os.path.isdir(aln_file):
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
        subprocess.run(cmd, check=True)
    elif args.subcommand == "pipeline":
        out_base = args.output_dir
        aln_dir = os.path.join(out_base, "alignment")
        phy_dir = os.path.join(out_base, "phylogeny")

        # 1. Fetch or Local Input
        if args.input_folder:
            if not os.path.isdir(args.input_folder):
                raise FileNotFoundError(f"Specified local structures directory '{args.input_folder}' does not exist.")
            pdb_dir = args.input_folder
            print(f"[Pipeline] Using local structure folder: '{pdb_dir}'")
        else:
            pdb_dir = os.path.join(out_base, "structures")
            fetch_structures(args.qualifier, args.count, pdb_dir)

        # 2. Align
        aln_3di = align_structures(
            pdb_dir,
            aln_dir,
            aligner=args.aligner,
            mafft_bin=getattr(args, "mafft_bin", None),
            mafft_matrix=getattr(args, "mafft_matrix", "matrices/mat3di.out"),
        )
        # 3. Tree
        prefix = "local_tree" if args.input_folder else f"{args.qualifier}_tree"
        build_tree(
            alignment_file=aln_3di,
            method=args.method,
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
        )

        # 4. Universal Metadata Processing
        if parse_metadata:
            print(f"\n[Pipeline] Analyzing and indexing metadata...")
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

        # 5. Optional PLM Embedding & Hierarchical Clustering Tree
        if getattr(args, "embed", False):
            print(f"\n[Pipeline] Extracting PLM embeddings ({args.embed_model}) & constructing hierarchical tree ({args.embed_clustering})...")
            embed_script = os.path.join(os.path.dirname(__file__), "embed_and_cluster.py")
            if not os.path.isfile(embed_script):
                embed_script = "scripts/embed_and_cluster.py"
            aa_aln = os.path.join(aln_dir, "foldmason.fasta_aa.fa")
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
            subprocess.run(cmd, check=True)




if __name__ == "__main__":
    main()
