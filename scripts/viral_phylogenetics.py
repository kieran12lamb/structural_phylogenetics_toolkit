#!/usr/bin/env python3
import sys
from pathlib import Path
_repo_dir = str(Path(__file__).resolve().parent.parent)
if _repo_dir not in sys.path:
    sys.path.insert(0, _repo_dir)
"""Viral Structural Phylogenetics, PLM Embeddings, and Interactive Tanglegrams.

Backward-compatible facade delegating to the modular viral_phylo package.
"""

from viral_phylo.binaries import (
    find_iqtree_bin,
    find_foldmason_bin,
    find_mafft_bin,
    find_torch_python,
)
from viral_phylo.matrices import (
    ensure_3di_matrix,
    ensure_matrix_file,
    EDMOND_MATRICES,
)
from viral_phylo.fetch import (
    BASE_URL,
    fetch_structures,
    fetch_alphafold_structures,
    query_uniprot_for_accessions,
)
from viral_phylo.alignment import (
    parse_alignment_fasta,
    write_alignment_fasta,
    strip_all_gap_columns_dict,
    compute_alignment_coverage,
    filter_alignment_by_coverage,
    partition_structures_by_coverage,
    count_fasta_seqs,
    align_structures,
)
from viral_phylo.tree import build_tree
from viral_phylo.metadata import parse_metadata
from viral_phylo.cli import build_cli_parser, main

if __name__ == "__main__":
    main()
