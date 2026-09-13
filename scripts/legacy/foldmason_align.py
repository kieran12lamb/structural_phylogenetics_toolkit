#!/usr/bin/env python3
"""Run FoldMason multiple structural alignment on protein structure files."""

import argparse
import glob
import os
import shutil
import subprocess
import sys


def parse_args():
    parser = argparse.ArgumentParser(
        description="Run FoldMason structural alignment on PDB/mmCIF files."
    )
    parser.add_argument(
        "-t",
        "--tool",
        default="easy-msa",
        help="FoldMason subtool to run (default: easy-msa)",
    )
    parser.add_argument(
        "pdb_dir_pos",
        nargs="?",
        default=None,
        help="Folder containing PDB/mmCIF files (positional shorthand, default: viro_3d_structures)",
    )
    parser.add_argument(
        "-i",
        "--folder",
        "--pdb-dir",
        dest="pdb_dir_flag",
        default=None,
        help="Folder containing PDB/mmCIF files (default: viro_3d_structures)",
    )
    parser.add_argument(
        "-o",
        "--output-dir",
        default="foldmason_alignments",
        help="Output directory for alignment results (default: foldmason_alignments)",
    )
    parser.add_argument(
        "-s",
        "--output-sequences",
        "--output-fasta",
        dest="output_sequences",
        default="foldmason.fasta",
        help="Output sequences filename / prefix (default: foldmason.fasta)",
    )
    parser.add_argument(
        "-r",
        "--report-mode",
        type=int,
        default=1,
        help="Report mode for FoldMason (default: 1)",
    )
    parser.add_argument(
        "--tmp-dir",
        default=None,
        help="Temporary directory (default: <output_dir>/tmpFolder)",
    )
    parser.add_argument(
        "--foldmason-bin",
        default="foldmason",
        help="Path or name of foldmason executable (default: foldmason)",
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Print the generated command without executing it",
    )
    return parser.parse_args()


def find_structure_files(pdb_dir: str):
    """Find all PDB and mmCIF structure files in the given directory."""
    extensions = ("*.pdb", "*.cif", "*.mmcif", "*.ent", "*.pdb.gz", "*.cif.gz")
    files = []
    for ext in extensions:
        files.extend(glob.glob(os.path.join(pdb_dir, ext)))
    return sorted(set(files))


def main():
    args = parse_args()

    tool = args.tool.strip().replace(" ", "-")
    pdb_dir = args.pdb_dir_flag or args.pdb_dir_pos or "viro_3d_structures"
    output_dir = args.output_dir
    output_sequences = args.output_sequences
    report_mode = args.report_mode
    foldmason_bin = args.foldmason_bin

    if not os.path.exists(pdb_dir):
        print(f"Error: Input directory '{pdb_dir}' does not exist.", file=sys.stderr)
        sys.exit(1)

    structure_files = find_structure_files(pdb_dir)
    if not structure_files:
        print(
            f"Error: No PDB/mmCIF files found in '{pdb_dir}'. "
            "Download structures first (e.g. using viro3d_structures.py).",
            file=sys.stderr,
        )
        sys.exit(1)

    print(f"Found {len(structure_files)} structure file(s) in '{pdb_dir}'.")
    if len(structure_files) < 2:
        print(
            "Warning: FoldMason typically requires at least 2 structures for alignment.",
            file=sys.stderr,
        )

    # Prepare directories
    os.makedirs(output_dir, exist_ok=True)
    tmp_folder = args.tmp_dir or os.path.join(output_dir, "tmpFolder")
    os.makedirs(tmp_folder, exist_ok=True)

    result_path = os.path.join(output_dir, output_sequences)

    # Construct FoldMason command:
    # foldmason <tool> <PDB/mmCIF files> <result> <tmpFolder> --report-mode <report_mode>
    cmd = [
        foldmason_bin,
        tool,
        *structure_files,
        result_path,
        tmp_folder,
        "--report-mode",
        str(report_mode),
    ]

    print("\nExecuting FoldMason command:")
    print(" ".join(cmd))
    print("-" * 60)

    if args.dry_run:
        print("Dry-run requested. Skipping execution.")
        return

    # Check executable existence
    if not shutil.which(foldmason_bin):
        print(
            f"\nError: '{foldmason_bin}' executable not found in PATH.",
            file=sys.stderr,
        )
        print(
            "Please ensure FoldMason is installed and on your PATH (e.g. conda install -c bioconda foldmason).",
            file=sys.stderr,
        )
        sys.exit(127)

    try:
        completed = subprocess.run(cmd, check=True)
        print(f"\nFoldMason completed successfully (exit code {completed.returncode}).")
        print(f"Results saved to '{output_dir}'.")
    except subprocess.CalledProcessError as e:
        print(f"\nError: FoldMason execution failed with code {e.returncode}.", file=sys.stderr)
        sys.exit(e.returncode)


if __name__ == "__main__":
    main()
