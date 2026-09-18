"""3Di substitution matrices and model parameter management."""

import os
import urllib.request
from typing import Dict, Any

EDMOND_MATRICES: Dict[str, Dict[str, Any]] = {
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

    print(f"Downloading {meta[desc]} from Edmond...")
    req = urllib.request.Request(meta["url"], headers={"User-Agent": "Mozilla/5.0"})
    with urllib.request.urlopen(req) as resp, open(target, "wb") as f:
        f.write(resp.read())
    return target
