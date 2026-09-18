#!/usr/bin/env python3
"""Dynamic metadata interactive tree builder.

Backward-compatible facade delegating to the modular viral_phylo.web.builder.
"""

import sys
from pathlib import Path

_repo_dir = str(Path(__file__).resolve().parent.parent)
if _repo_dir not in sys.path:
    sys.path.insert(0, _repo_dir)

from viral_phylo.web.builder import (
    build_interactive_tree,
    assemble_html,
    load_or_compute_umap,
    DATASETS,
    DATASET_CONGRUENCE,
    main,
)

if __name__ == "__main__":
    main()
