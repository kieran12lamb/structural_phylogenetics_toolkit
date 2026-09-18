"""Viral Structural Phylogenetics, PLM Embeddings, and Interactive Cophylogenetic Tanglegrams.

A unified framework for structural phylogenetics (FoldMason, IQ-TREE with 3Di substitution matrices),
Protein Language Model embeddings (ESM-2 / ESM-C), 2D UMAP projections, coverage-based sequence filtering,
and zero-dependency interactive visualization suites.
"""

__version__ = "1.0.0"

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
from viral_phylo.tree import (
    build_tree,
)
from viral_phylo.embeddings import (
    generate_esm_embeddings,
    build_hierarchical_tree,
    compute_silhouette_profile,
    compute_umap_projection,
    run_embedding_pipeline,
    extract_sequences_from_fasta,
    extract_sequences_from_pdb_dir,
)
from viral_phylo.metadata import (
    parse_metadata,
    classify_columns,
    choose_default_color_column,
    EXPANDED_PALETTE,
    KNOWN_VALUE_COLORS,
)

__all__ = [
    "find_iqtree_bin",
    "find_foldmason_bin",
    "find_mafft_bin",
    "find_torch_python",
    "ensure_3di_matrix",
    "ensure_matrix_file",
    "EDMOND_MATRICES",
    "fetch_structures",
    "fetch_alphafold_structures",
    "query_uniprot_for_accessions",
    "parse_alignment_fasta",
    "write_alignment_fasta",
    "strip_all_gap_columns_dict",
    "compute_alignment_coverage",
    "filter_alignment_by_coverage",
    "partition_structures_by_coverage",
    "count_fasta_seqs",
    "align_structures",
    "build_tree",
    "generate_esm_embeddings",
    "build_hierarchical_tree",
    "compute_silhouette_profile",
    "compute_umap_projection",
    "run_embedding_pipeline",
    "extract_sequences_from_fasta",
    "extract_sequences_from_pdb_dir",
    "parse_metadata",
    "classify_columns",
    "choose_default_color_column",
    "EXPANDED_PALETTE",
    "KNOWN_VALUE_COLORS",
]
