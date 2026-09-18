#!/usr/bin/env python3
"""Backward-compatible facade for viral_phylo.web.server."""
import sys
from pathlib import Path

_repo_dir = str(Path(__file__).resolve().parent.parent)
if _repo_dir not in sys.path:
    sys.path.insert(0, _repo_dir)

from viral_phylo.web.server import main

if __name__ == "__main__":
    main()
