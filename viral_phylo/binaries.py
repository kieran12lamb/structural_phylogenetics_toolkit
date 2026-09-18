"""Binary and interpreter discovery utilities for structural phylogenetics and ML tools."""

import os
import shutil
import sys
from typing import Optional


def find_iqtree_bin() -> str:
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


def find_foldmason_bin() -> str:
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


def find_mafft_bin(mafft_bin: Optional[str] = None) -> str:
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
