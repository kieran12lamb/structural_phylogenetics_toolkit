#!/usr/bin/env python3
import sys
from pathlib import Path
_repo_dir = str(Path(__file__).resolve().parent.parent)
if _repo_dir not in sys.path:
    sys.path.insert(0, _repo_dir)
"""Backward-compatible facade for viral_phylo.web.alignments."""
from viral_phylo.web.alignments import (
    parse_fasta,
    compute_cov_dict,
    find_cluster_partitions,
    build_alignments_data,
    main,
)

if __name__ == "__main__":
    main()
