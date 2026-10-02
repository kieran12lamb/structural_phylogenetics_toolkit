"""3Di substitution matrices and model parameter management."""

import os
import shutil
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

    print(f"Downloading {meta['desc']} from Edmond...")
    req = urllib.request.Request(meta["url"], headers={"User-Agent": "Mozilla/5.0"})
    with urllib.request.urlopen(req) as resp, open(target, "wb") as f:
        f.write(resp.read())
    return target


def uppercase_matrix_alias(matrix_path: str) -> str:
    """Return a path to the same matrix that survives IQ-TREE's uppercasing of -mset.

    IQ-TREE upper-cases the whole ``-mset`` argument, so on a case-sensitive
    filesystem ``-mset matrices/Q.3Di.AF`` is looked up as ``MATRICES/Q.3DI.AF``
    and fails with "File not found". That is why every small-dataset tree in the
    25 Sep run exited 2 while the large master tree - which takes the ``-m`` branch,
    and ``-m`` is not upper-cased - succeeded.

    Upper-casing is idempotent, so an already-upper-case path passes through
    unchanged. The copy is always made at ``MATRICES/<BASENAME>`` relative to the
    working directory IQ-TREE runs in, which works for absolute inputs too (an
    absolute ``/raid/...`` path would be upper-cased to ``/RAID/...``, which cannot
    exist). Accepts a comma-separated list, as ``ensure_matrix_file('both')`` returns.
    """
    import hashlib

    def _digest(path):
        with open(path, "rb") as handle:
            return hashlib.sha256(handle.read()).hexdigest()

    aliases = []
    for part in str(matrix_path).split(","):
        part = part.strip()
        if not part:
            continue
        if part == part.upper() or not os.path.isfile(part):
            aliases.append(part)
            continue
        target = os.path.join("MATRICES", os.path.basename(part).upper())
        try:
            os.makedirs("MATRICES", exist_ok=True)
            if os.path.exists(target) and os.path.samefile(part, target):
                aliases.append(target)          # case-insensitive filesystem
                continue
            if os.path.isfile(target) and _digest(target) != _digest(part):
                # A different matrix already holds this name; do not overwrite it
                # (another run may be using it) - use a content-addressed name.
                target = f"{target}.{_digest(part)[:8].upper()}"
            if not os.path.isfile(target):
                # Write-then-rename: concurrent tree jobs may create the same alias,
                # and none may ever see a half-written matrix.
                tmp = f"{target}.tmp.{os.getpid()}"
                shutil.copyfile(part, tmp)
                os.replace(tmp, target)
            aliases.append(target)
        except OSError as exc:
            print(f"[IQ-TREE] [!] Notice: could not create an upper-case alias for '{part}' ({exc}); "
                  "IQ-TREE -mset is likely to fail to find it.")
            aliases.append(part)
    return ",".join(aliases)


# Amino-acid order used by PAML-format matrices (and so by IQ-TREE) and by FastTree.
# The 20-letter 3Di alphabet is written with these same letters, so a 3Di matrix
# converts exactly like an amino-acid one.
PAML_AA_ORDER = "ARNDCQEGHILKMFPSTWYV"


def read_paml_matrix(paml_path: str):
    """Read a PAML-format exchangeability matrix: 190 lower-triangle values, then 20
    stationary frequencies. Returns ``(S, pi)`` with ``S`` a full symmetric 20x20."""
    import re

    with open(paml_path, "r", encoding="utf-8") as handle:
        numbers = [float(x) for x in re.findall(
            r"[-+]?\d*\.\d+(?:[eE][-+]?\d+)?|[-+]?\d+(?:[eE][-+]?\d+)?", handle.read())]
    if len(numbers) < 210:
        raise ValueError(f"'{paml_path}' has {len(numbers)} numbers; a PAML amino-acid matrix needs 210 "
                         "(190 exchangeabilities followed by 20 frequencies).")
    tri, pi = numbers[:190], numbers[190:210]
    if abs(sum(pi) - 1.0) > 1e-3:
        raise ValueError(f"stationary frequencies in '{paml_path}' sum to {sum(pi):.6f}, not 1.")
    S = [[0.0] * 20 for _ in range(20)]
    k = 0
    for i in range(1, 20):
        for j in range(i):
            S[i][j] = S[j][i] = tri[k]
            k += 1
    return S, pi


def paml_to_fasttree_trans(paml_path: str, out_path: str) -> str:
    """Convert a PAML exchangeability matrix into a FastTree ``-trans`` model file.

    FastTree (and VeryFastTree) can use a custom amino-acid model, but expect a
    different representation from PAML. Worked out against FastTree's own
    validation messages and confirmed by likelihood on a fixed tree:

      * the table holds the rate matrix Q, not the exchangeabilities S:
        Q_ij = S_ij * pi_j off the diagonal, with the diagonal making each row of Q
        sum to zero; normalised to one expected substitution per unit time, as
        IQ-TREE does;
      * it is written transposed - row i, column j is the rate from j to i - so it
        is the *columns* that sum to zero;
      * tab-delimited, header ``A R N ... V *`` with no leading cell, one row per
        letter, and a final ``*`` column holding the stationary frequencies.

    Validation: on a 32-taxon 3Di cluster with branch lengths optimised on a fixed
    topology, IQ-TREE ``-m Q.3Di.AF`` gave log-likelihood -3830.308 and FastTree
    ``-trans`` with this file gave -3830.305 (VeryFastTree with -double-precision
    likewise); the LG control agreed to 0.002.
    """
    S, pi = read_paml_matrix(paml_path)
    Q = [[S[i][j] * pi[j] if i != j else 0.0 for j in range(20)] for i in range(20)]
    for i in range(20):
        Q[i][i] = -sum(Q[i][j] for j in range(20) if j != i)
    rate = -sum(pi[i] * Q[i][i] for i in range(20))
    Q = [[q / rate for q in row] for row in Q]

    os.makedirs(os.path.dirname(os.path.abspath(out_path)) or ".", exist_ok=True)
    with open(out_path, "w", encoding="utf-8") as handle:
        handle.write("\t".join(PAML_AA_ORDER) + "\t*\n")
        for i, letter in enumerate(PAML_AA_ORDER):
            handle.write(letter + "\t" + "\t".join(repr(Q[j][i]) for j in range(20)) + f"\t{pi[i]!r}\n")
    return out_path
