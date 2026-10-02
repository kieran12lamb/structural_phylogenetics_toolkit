#!/usr/bin/env python3
"""Backbone (C-alpha) extraction from PDB and mmCIF structure files.

The interactive viewer renders a C-alpha trace per taxon and colours it by the
per-residue confidence score (pLDDT for predicted structures). This module is
the single place that turns a structure file into that trace.

mmCIF is parsed by reading the ``_atom_site`` loop header and mapping field
names to column indices. Fixed-column slicing - which is correct for PDB - is
wrong for mmCIF, whose rows are whitespace-delimited, and silently matches
nothing: the pipeline previously produced an empty C-alpha bundle for every
mmCIF input without reporting an error.
"""

import glob
import os
from typing import Dict, List, Optional

# One residue: [x, y, z, confidence, residue_number, residue_name]
Residue = List[object]

_CIF_PREFIX = "_atom_site."


def _select_confidence(occupancies: List[float], b_factors: List[float]) -> List[float]:
    """Choose which column carries the per-residue confidence score.

    AlphaFold and ESMFold conventionally write pLDDT into the B-factor column and
    leave occupancy at 1.00, but not every writer does: the ESMFold mmCIF files in
    this project declare ``occupancy`` first and put the varying pLDDT there, with
    ``B_iso_or_equiv`` pinned at 1.00. Reading the spec-correct field alone would
    yield a flat trace and silently lose the confidence colouring.

    So prefer whichever column actually varies, falling back to the B-factor.
    """
    if len(set(b_factors)) > 1:
        return b_factors
    if len(set(occupancies)) > 1:
        return occupancies
    return b_factors


def _parse_cif(path: str) -> List[Residue]:
    """Parse C-alpha atoms from an mmCIF file via its ``_atom_site`` loop header."""
    fields: List[str] = []
    in_loop_header = False
    rows: List[List[str]] = []

    with open(path, "r", encoding="utf-8", errors="ignore") as handle:
        for line in handle:
            stripped = line.strip()
            if stripped.startswith(_CIF_PREFIX):
                fields.append(stripped[len(_CIF_PREFIX):].split()[0])
                in_loop_header = True
                continue
            if in_loop_header and (stripped.startswith("ATOM") or stripped.startswith("HETATM")):
                rows.append(stripped.split())
            elif in_loop_header and stripped in ("#", "loop_", ""):
                # End of this loop's data block; keep the header we collected.
                continue

    if not fields or not rows:
        return []

    index = {name: i for i, name in enumerate(fields)}
    atom_key = "label_atom_id" if "label_atom_id" in index else "auth_atom_id"
    comp_key = "label_comp_id" if "label_comp_id" in index else "auth_comp_id"
    seq_key = "label_seq_id" if "label_seq_id" in index else "auth_seq_id"

    required = [atom_key, "Cartn_x", "Cartn_y", "Cartn_z"]
    if any(key not in index for key in required):
        return []

    coords: List[List[float]] = []
    occupancies: List[float] = []
    b_factors: List[float] = []
    meta: List[List[object]] = []

    for row in rows:
        if len(row) < len(fields):
            continue
        if row[index[atom_key]] != "CA":
            continue
        try:
            x = round(float(row[index["Cartn_x"]]), 1)
            y = round(float(row[index["Cartn_y"]]), 1)
            z = round(float(row[index["Cartn_z"]]), 1)
        except (ValueError, IndexError):
            continue

        def _number(key: str, default: float = 0.0) -> float:
            if key not in index:
                return default
            try:
                return round(float(row[index[key]]), 1)
            except (ValueError, IndexError):
                return default

        resnum_raw = row[index[seq_key]] if seq_key in index else "0"
        try:
            resnum = int(resnum_raw)
        except ValueError:
            resnum = len(meta) + 1

        coords.append([x, y, z])
        occupancies.append(_number("occupancy"))
        b_factors.append(_number("B_iso_or_equiv"))
        meta.append([resnum, row[index[comp_key]] if comp_key in index else "UNK"])

    if not coords:
        return []

    confidence = _select_confidence(occupancies, b_factors)
    return [
        [c[0], c[1], c[2], conf, m[0], m[1]]
        for c, conf, m in zip(coords, confidence, meta)
    ]


def _parse_pdb(path: str) -> List[Residue]:
    """Parse C-alpha atoms from a PDB file using the fixed-column record format."""
    coords: List[List[float]] = []
    occupancies: List[float] = []
    b_factors: List[float] = []
    meta: List[List[object]] = []

    with open(path, "r", encoding="utf-8", errors="ignore") as handle:
        for line in handle:
            if not line.startswith("ATOM"):
                continue
            if line[12:16].strip() != "CA":
                continue
            try:
                resname = line[17:20].strip()
                resnum = int(line[22:26].strip())
                x = round(float(line[30:38].strip()), 1)
                y = round(float(line[38:46].strip()), 1)
                z = round(float(line[46:54].strip()), 1)
            except (ValueError, IndexError):
                continue

            def _column(start: int, end: int) -> float:
                try:
                    return round(float(line[start:end].strip()), 1)
                except (ValueError, IndexError):
                    return 0.0

            coords.append([x, y, z])
            occupancies.append(_column(54, 60))
            b_factors.append(_column(60, 66))
            meta.append([resnum, resname])

    if not coords:
        return []

    confidence = _select_confidence(occupancies, b_factors)
    return [
        [c[0], c[1], c[2], conf, m[0], m[1]]
        for c, conf, m in zip(coords, confidence, meta)
    ]


def parse_structure_ca(path: str) -> List[Residue]:
    """Extract the C-alpha trace from a PDB or mmCIF structure file.

    Returns a list of ``[x, y, z, confidence, residue_number, residue_name]``,
    or an empty list if the file holds no parseable C-alpha atoms.
    """
    suffix = os.path.splitext(path)[1].lower()
    if suffix in (".cif", ".mmcif"):
        return _parse_cif(path)
    return _parse_pdb(path)


def extract_ca_dict(structure_dir: str, verbose: bool = True) -> Dict[str, List[Residue]]:
    """Extract C-alpha traces for every structure in a directory, keyed by taxon id.

    Reports how many structures yielded no backbone rather than returning a
    silently empty bundle.
    """
    paths: List[str] = []
    for pattern in ("*.pdb", "*.cif", "*.mmcif"):
        paths.extend(glob.glob(os.path.join(structure_dir, pattern)))
    paths = sorted(set(paths))

    ca_dict: Dict[str, List[Residue]] = {}
    failed: List[str] = []
    for path in paths:
        taxon_id = os.path.splitext(os.path.basename(path))[0]
        trace = parse_structure_ca(path)
        if trace:
            ca_dict[taxon_id] = trace
        else:
            failed.append(os.path.basename(path))

    if verbose and paths:
        print(f"[Structures] Extracted C-alpha backbones for {len(ca_dict)}/{len(paths)} structures.")
        if failed:
            preview = ", ".join(failed[:5]) + (f" ... (+{len(failed) - 5} more)" if len(failed) > 5 else "")
            print(f"[Structures] [!] Notice: no C-alpha atoms parsed from {len(failed)} file(s): {preview}")

    return ca_dict


def write_ca_bundle(ca_dict: Dict[str, List[Residue]], output_path: str) -> Optional[str]:
    """Write a C-alpha bundle as the JS asset the interactive viewer loads."""
    import json

    if not ca_dict:
        return None
    payload = "window.CA_STRUCTURES = Object.assign(window.CA_STRUCTURES || {}, " + json.dumps(ca_dict) + ");\n"
    with open(output_path, "w", encoding="utf-8") as handle:
        handle.write(payload)
    return output_path
