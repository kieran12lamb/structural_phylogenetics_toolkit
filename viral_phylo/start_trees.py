#!/usr/bin/env python3
"""FastTree starting trees for IQ-TREE on large clusters.

Large clusters dominate the pipeline's runtime, and giving IQ-TREE a good
starting topology can shorten its search: on a 175-taxon 3Di cluster the guided
run finished in 793 s against 1052 s, at a marginally better likelihood.

A FastTree tree cannot be handed to IQ-TREE as-is, though. Two incompatibilities
were hit on real clusters:

  * IQ-TREE drops duplicate sequences before searching, then rejects a starting
    tree that still contains them ("Tree taxon ... does not appear in the
    alignment"). Your binder set has 27 duplicates. Guided runs therefore pass
    ``-keep-ident`` so IQ-TREE keeps every taxon the tree contains.
  * FastTree joins identical sequences in a multifurcation. IQ-TREE's NNI search
    asserts every node is bifurcating and aborts ("Assertion node1->degree() == 3").
    Polytomies are therefore resolved into random binary splits of near-zero
    length, seeded so the result is reproducible.

Every starting tree is then validated - leaf set identical to the alignment,
strictly bifurcating, no labels - and anything that fails validation is simply not
used. The caller also falls back to an unguided IQ-TREE run if a guided one fails,
so a starting tree can cost time but never lose a tree.
"""

import io
import os
import random
import subprocess
from typing import Any, Dict, List, Optional

MIN_BRANCH = 1e-6


def _alignment_names(path: str) -> List[str]:
    with open(path, "r", encoding="utf-8", errors="ignore") as handle:
        return [line[1:].split()[0] for line in handle if line.startswith(">") and len(line) > 1]


def resolve_to_binary(tree, seed: int = 1) -> int:
    """Make a Bio.Phylo tree strictly bifurcating, unrooted-style, in place.

    Non-root clades end up with exactly two children and the root with three, which
    is the form IQ-TREE's NNI search requires. Multifurcations are split at random
    (seeded) with near-zero branch lengths; unary nodes are collapsed. Internal
    labels and support values are removed, and missing or non-positive branch
    lengths are set to a small positive value. Returns the number of polytomies
    resolved.
    """
    from Bio.Phylo.BaseTree import Clade

    rng = random.Random(seed)
    resolved = 0

    def _fix(clade, is_root):
        nonlocal resolved
        # Collapse unary internal nodes first.
        changed = True
        while changed:
            changed = False
            for i, child in enumerate(list(clade.clades)):
                if len(child.clades) == 1:
                    grandchild = child.clades[0]
                    grandchild.branch_length = (grandchild.branch_length or 0.0) + (child.branch_length or 0.0)
                    clade.clades[i] = grandchild
                    changed = True
        for child in clade.clades:
            if child.clades:
                _fix(child, False)
        target = 3 if is_root else 2
        if len(clade.clades) > target:
            resolved += 1
            children = list(clade.clades)
            rng.shuffle(children)
            while len(children) > target:
                a, b = children.pop(), children.pop()
                children.insert(rng.randrange(len(children) + 1),
                                Clade(branch_length=MIN_BRANCH, clades=[a, b]))
            clade.clades = children

    root = tree.root
    _fix(root, True)
    # An explicitly rooted binary tree (root with two children) is unrooted by
    # merging one internal child into the root, giving the trifurcating root IQ-TREE
    # expects for an unrooted tree.
    if len(root.clades) == 2:
        internal = [c for c in root.clades if c.clades]
        if internal:
            merge = internal[0]
            other = [c for c in root.clades if c is not merge][0]
            other.branch_length = (other.branch_length or 0.0) + (merge.branch_length or 0.0)
            root.clades = [other] + list(merge.clades)

    for clade in tree.find_clades():
        if clade is not root and (clade.branch_length is None or clade.branch_length <= 0):
            clade.branch_length = MIN_BRANCH
        clade.confidence = None
        if clade.clades:
            clade.name = None
    root.branch_length = None
    return resolved


def validate_start_tree(tree, alignment_names: List[str]) -> Optional[str]:
    """Return None if the tree is safe to give IQ-TREE, else the reason it is not."""
    leaves = [t.name for t in tree.get_terminals()]
    if len(leaves) != len(set(leaves)):
        return "tree has duplicate leaf names"
    if set(leaves) != set(alignment_names) or len(leaves) != len(alignment_names):
        missing = set(alignment_names) - set(leaves)
        extra = set(leaves) - set(alignment_names)
        return (f"leaf set does not match the alignment ({len(missing)} missing, {len(extra)} extra; "
                f"e.g. {sorted(missing or extra)[:3]})")
    root = tree.root
    if len(root.clades) != 3:
        return f"root has {len(root.clades)} children, expected 3"
    for clade in tree.find_clades():
        if clade is not root and clade.clades and len(clade.clades) != 2:
            return f"a node has {len(clade.clades)} children; tree is not bifurcating"
    for bad in ("(", ")", ":", ",", ";", " ", "'", "[", "]"):
        if any(bad in name for name in leaves):
            return f"a leaf name contains the Newick metacharacter {bad!r}"
    return None


def make_start_tree(
    alignment: str,
    out_path: str,
    alphabet: str,
    matrix_path: Optional[str] = None,
    threads: str = "8",
    seed: int = 1,
) -> Dict[str, Any]:
    """Build, repair and validate a FastTree starting tree for one IQ-TREE run.

    ``alphabet`` is "3di" (uses ``matrix_path`` via -trans, matching IQ-TREE's
    Q.3Di model) or "aa" (LG, matching IQ-TREE's LG). Returns a record with
    ``used`` True only when the tree passed validation.
    """
    from viral_phylo.binaries import find_fasttree_bin
    from viral_phylo.matrices import paml_to_fasttree_trans

    record: Dict[str, Any] = {"used": False, "path": None, "reason": None, "command": None,
                              "polytomies_resolved": 0, "alphabet": alphabet, "seed": seed}
    try:
        from Bio import Phylo
    except ImportError:
        # Without Biopython the tree cannot be repaired, so it must not be used; the
        # caller then runs IQ-TREE unguided rather than failing.
        record["reason"] = "Biopython is not installed, so a starting tree cannot be repaired"
        return record
    binary = find_fasttree_bin()
    if not binary:
        record["reason"] = "no FastTree/VeryFastTree binary found"
        return record

    os.makedirs(os.path.dirname(os.path.abspath(out_path)) or ".", exist_ok=True)
    if alphabet == "3di" and matrix_path:
        trans = paml_to_fasttree_trans(matrix_path, out_path + ".trans")
        model = ["-trans", trans]
    else:
        model = ["-lg"]
    cmd = [binary, *model, "-gamma"]
    if "veryfasttree" in os.path.basename(binary).lower():
        cmd.append("-double-precision")
        if str(threads).isdigit():
            cmd += ["-threads", str(threads)]
    cmd.append(alignment)
    record["command"] = " ".join(cmd)

    try:
        proc = subprocess.run(cmd, capture_output=True, text=True, check=True, timeout=6 * 3600)
    except (subprocess.CalledProcessError, subprocess.TimeoutExpired, OSError) as exc:
        record["reason"] = f"FastTree failed: {type(exc).__name__}"
        return record
    if not proc.stdout.strip():
        record["reason"] = "FastTree produced no tree"
        return record

    try:
        tree = Phylo.read(io.StringIO(proc.stdout), "newick")
        record["polytomies_resolved"] = resolve_to_binary(tree, seed=seed)
    except Exception as exc:
        record["reason"] = f"could not parse or repair the FastTree tree: {type(exc).__name__}: {exc}"
        return record

    problem = validate_start_tree(tree, _alignment_names(alignment))
    if problem:
        record["reason"] = f"starting tree rejected: {problem}"
        return record

    buffer = io.StringIO()
    Phylo.write(tree, buffer, "newick")
    with open(out_path, "w", encoding="utf-8") as handle:
        handle.write(buffer.getvalue().strip() + "\n")
    record.update({"used": True, "path": out_path})
    return record
