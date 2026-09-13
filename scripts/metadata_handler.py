#!/usr/bin/env python3
"""Universal metadata ingestion and schema analysis utility.

Handles:
- Formats: .xlsx, .xls, .csv, .tsv, .json, or direct structure folder extraction
- Automatic ID detection matching against structure file stems or tree leaves
- Categorical vs continuous column classification
- Color palette generation for discrete and gradient variables
- Export to structured JSON for interactive visualizations
"""

import os
import sys
import glob
import json
import re
import csv
import zipfile
import xml.etree.ElementTree as ET
from pathlib import Path

# High-contrast 40-color categorical palette (accessible, distinct, visually pleasant)
EXPANDED_PALETTE = [
    "#38bdf8", "#f97316", "#a855f7", "#10b981", "#ec4899",
    "#eab308", "#06b6d4", "#8b5cf6", "#f43f5e", "#14b8a6",
    "#6366f1", "#84cc16", "#e11d48", "#0284c7", "#ca8a04",
    "#d946ef", "#059669", "#b45309", "#4f46e5", "#f59e0b",
    "#22c55e", "#0ea5e9", "#d97706", "#9333ea", "#2dd4bf",
    "#fb7185", "#3b82f6", "#16a34a", "#c026d3", "#64748b",
    "#e2e8f0", "#78716c", "#a3e635", "#34d399", "#818cf8",
    "#c084fc", "#f472b6", "#fb923c", "#facc15", "#4ade80"
]

# Standard biological domain maps
KNOWN_VALUE_COLORS = {
    # Structural Classes
    "mainly alpha": "#38bdf8",
    "alpha beta": "#a855f7",
    "mainly beta": "#f97316",
    "few secondary structures": "#10b981",
    "unclassified": "#94a3b8",
    
    # Binding Strengths
    "strong": "#10b981",
    "medium": "#38bdf8",
    "weak": "#f59e0b",
    "none": "#64748b",
    
    # Viral Families
    "rhabdoviridae": "#38bdf8",
    "orthoherpesviridae": "#ec4899",
    "phenuiviridae": "#10b981",
    "hantaviridae": "#f97316",
    "peribunyaviridae": "#8b5cf6",
    "arenaviridae": "#06b6d4",
    "nairoviridae": "#eab308",
    "togaviridae": "#a855f7",
    "arteriviridae": "#14b8a6",
    "phasmaviridae": "#f43f5e",
    "chuviridae": "#6366f1",
    "matonaviridae": "#84cc16",
    "orthomyxoviridae": "#e11d48",
    "nyamiviridae": "#10b981",
    "lispiviridae": "#ca8a04",
    "poxviridae": "#d946ef",
    "aliusviridae": "#0284c7",
    "tobaniviridae": "#f59e0b",
    "flaviviridae": "#e11d48",
    "xinmoviridae": "#4f46e5",
    "baculoviridae": "#059669",
    "bornaviridae": "#b45309",
    "unknown": "#94a3b8"
}


def read_xlsx_native(xlsx_path):
    """Read an Excel .xlsx file using Python's standard zipfile and xml libraries."""
    with zipfile.ZipFile(xlsx_path, "r") as z:
        shared_strings = []
        if "xl/sharedStrings.xml" in z.namelist():
            tree = ET.fromstring(z.read("xl/sharedStrings.xml"))
            for si in tree.findall("{http://schemas.openxmlformats.org/spreadsheetml/2006/main}si"):
                text = "".join(t.text for t in si.findall(".//{http://schemas.openxmlformats.org/spreadsheetml/2006/main}t") if t.text)
                shared_strings.append(text)

        # Pick first worksheet
        ws_names = [n for n in z.namelist() if n.startswith("xl/worksheets/sheet") and n.endswith(".xml")]
        ws_name = sorted(ws_names)[0] if ws_names else "xl/worksheets/sheet1.xml"
        
        ws_xml = z.read(ws_name)
        ws_tree = ET.fromstring(ws_xml)
        ns = {"s": "http://schemas.openxmlformats.org/spreadsheetml/2006/main"}

        raw_rows = []
        for row in ws_tree.findall(".//s:row", ns):
            row_dict = {}
            for c in row.findall("s:c", ns):
                r_ref = c.get("r", "")
                col_letters = "".join(ch for ch in r_ref if ch.isalpha())
                t_attr = c.get("t")
                v_el = c.find("s:v", ns)
                val = v_el.text if v_el is not None else None
                if val is not None:
                    if t_attr == "s":
                        try:
                            val = shared_strings[int(val)]
                        except Exception:
                            pass
                    row_dict[col_letters] = val
            if row_dict:
                raw_rows.append(row_dict)

    if not raw_rows:
        return []

    # Map column letters to header labels
    header_row = raw_rows[0]
    col_map = {col: str(val).strip() for col, val in header_row.items() if val is not None}
    
    records = []
    for r in raw_rows[1:]:
        rec = {}
        for col_letter, val in r.items():
            if col_letter in col_map:
                rec[col_map[col_letter]] = val
        if rec:
            records.append(rec)
    return records


def read_delimited(file_path):
    """Read CSV, TSV, or delimited file with dialect sniffing."""
    with open(file_path, "r", encoding="utf-8", errors="ignore") as f:
        sample = f.read(8192)
        f.seek(0)
        try:
            dialect = csv.Sniffer().sniff(sample, delimiters=",\t;|")
        except Exception:
            dialect = csv.excel
        reader = csv.DictReader(f, dialect=dialect)
        return [row for row in reader]


def read_json_file(file_path):
    """Read JSON file (dict of dicts, or list of dicts)."""
    with open(file_path, "r", encoding="utf-8") as f:
        data = json.load(f)
    if isinstance(data, list):
        return data
    elif isinstance(data, dict):
        if "taxa" in data and isinstance(data["taxa"], dict):
            # Already normalized
            out = []
            for tid, fields in data["taxa"].items():
                item = dict(fields)
                item["id"] = tid
                out.append(item)
            return out
        else:
            out = []
            for k, v in data.items():
                if isinstance(v, dict):
                    item = dict(v)
                    item["id"] = k
                    out.append(item)
                else:
                    out.append({"id": k, "value": v})
            return out
    return []


def extract_metadata_from_structures(structure_dir):
    """Fallback: extract basic metadata (length, avg pLDDT/B-factor) directly from structure files."""
    extensions = ("*.pdb", "*.cif", "*.mmcif", "*.ent")
    files = []
    for ext in extensions:
        files.extend(glob.glob(os.path.join(structure_dir, ext)))
    files = sorted(set(files))

    records = []
    for fpath in files:
        stem = Path(fpath).stem
        # Clean common prefixes/suffixes
        clean_id = stem.replace("CF-", "").replace("_relaxed", "").replace(".cif", "").replace(".pdb", "")
        
        ca_count = 0
        b_factors = []
        is_cif = fpath.endswith(".cif") or fpath.endswith(".mmcif")
        
        try:
            with open(fpath, "r", errors="ignore") as f:
                for line in f:
                    if is_cif:
                        if line.startswith("ATOM") or line.startswith("HETATM"):
                            parts = line.split()
                            # Check for CA atom (usually col 3 in CIF)
                            if len(parts) > 10 and ("CA" in parts[3] or parts[3] == "CA"):
                                ca_count += 1
                                try:
                                    # CIF B-factor is often near column 14 or 15
                                    b_val = float(parts[-4]) if len(parts) >= 15 else float(parts[14])
                                    b_factors.append(b_val)
                                except Exception:
                                    pass
                    else:
                        if line.startswith("ATOM") and len(line) >= 54:
                            atom_name = line[12:16].strip()
                            if atom_name == "CA":
                                ca_count += 1
                                try:
                                    bf = float(line[60:66].strip())
                                    b_factors.append(bf)
                                except Exception:
                                    pass
        except Exception:
            pass

        avg_plddt = round(sum(b_factors) / len(b_factors), 1) if b_factors else 75.0
        records.append({
            "id": clean_id,
            "filename": os.path.basename(fpath),
            "length": ca_count if ca_count > 0 else 300,
            "plddt": avg_plddt,
            "source": "Local Structure"
        })
    return records


def find_best_id_column(records, candidate_ids=None):
    """Determine which column contains the unique structure/taxon identifiers."""
    if not records:
        return "id"

    cols = list(records[0].keys())
    if not cols:
        return "id"

    # If candidate IDs are provided (e.g. from tree leaves or folder files)
    if candidate_ids:
        cand_set = set(str(c).strip().lower() for c in candidate_ids)
        best_col = None
        best_overlap = 0
        for c in cols:
            matches = sum(1 for r in records if str(r.get(c, "")).strip().lower() in cand_set)
            if matches > best_overlap:
                best_overlap = matches
                best_col = c
        if best_col and best_overlap > len(records) * 0.3:
            return best_col

    # Standard naming heuristics
    priority_names = [
        "record_id", "recordid", "id", "taxon_id", "taxon", "name",
        "structure_id", "sequence_id", "seq_id", "genbank", "accession",
        "filename", "file", "design_id"
    ]
    col_lower = {c.lower().replace("_", "").replace("-", "").strip(): c for c in cols}
    for p in priority_names:
        p_clean = p.replace("_", "").replace("-", "")
        if p_clean in col_lower:
            return col_lower[p_clean]

    # Return the first column with high uniqueness (>70%)
    for c in cols:
        vals = [str(r.get(c, "")) for r in records if r.get(c) is not None]
        if len(set(vals)) >= len(vals) * 0.7:
            return c

    return cols[0]


def classify_columns(records, id_col):
    """Analyze and classify each column as categorical or continuous, with generated palettes."""
    if not records:
        return []

    cols = [c for c in records[0].keys() if c != id_col]
    column_specs = []

    for col in cols:
        vals = [r.get(col) for r in records if r.get(col) is not None and str(r.get(col)).strip() != ""]
        if not vals:
            continue

        # Check if numerical
        is_num = True
        num_vals = []
        for v in vals:
            try:
                num_vals.append(float(v))
            except (ValueError, TypeError):
                is_num = False
                break

        clean_label = col.replace("_", " ").strip().title()

        if is_num and len(set(num_vals)) > 8:
            # Continuous numerical variable
            min_val = min(num_vals)
            max_val = max(num_vals)
            column_specs.append({
                "key": col,
                "label": clean_label,
                "type": "continuous",
                "min": round(min_val, 2),
                "max": round(max_val, 2),
                "cardinality": len(set(num_vals))
            })
        else:
            # Categorical variable
            str_vals = [str(v).strip() for v in vals]
            unique_vals = sorted(list(set(str_vals)))
            cardinality = len(unique_vals)

            if cardinality <= 65:
                # Assign distinct palette colors
                color_map = {}
                pal_idx = 0
                for u in unique_vals:
                    u_lower = u.lower().strip()
                    if u_lower in KNOWN_VALUE_COLORS:
                        color_map[u] = KNOWN_VALUE_COLORS[u_lower]
                    else:
                        color_map[u] = EXPANDED_PALETTE[pal_idx % len(EXPANDED_PALETTE)]
                        pal_idx += 1

                column_specs.append({
                    "key": col,
                    "label": clean_label,
                    "type": "categorical",
                    "values": unique_vals,
                    "colors": color_map,
                    "cardinality": cardinality
                })
            else:
                # High cardinality text field (e.g. description, sequence, note)
                column_specs.append({
                    "key": col,
                    "label": clean_label,
                    "type": "text",
                    "cardinality": cardinality
                })

    return column_specs


def choose_default_color_column(column_specs):
    """Intelligently select the most biologically relevant initial color column."""
    cat_specs = [c for c in column_specs if c["type"] == "categorical"]
    if not cat_specs:
        cont_specs = [c for c in column_specs if c["type"] == "continuous"]
        return cont_specs[0]["key"] if cont_specs else None

    # Priority 1: Family / Taxonomy
    for c in cat_specs:
        k = c["key"].lower()
        if "family" in k or "taxon" in k or "lineage" in k or "clade" in k:
            return c["key"]

    # Priority 2: Structural classification
    for c in cat_specs:
        k = c["key"].lower()
        if "class" in k or "fold" in k or "struct" in k or "cluster" in k:
            return c["key"]

    # Priority 3: Binding or phenotype
    for c in cat_specs:
        k = c["key"].lower()
        if "bind" in k or "affinity" in k or "host" in k or "phenotype" in k:
            return c["key"]

    # Priority 4: First categorical with 2 to 30 categories
    for c in cat_specs:
        if 2 <= c["cardinality"] <= 30:
            return c["key"]

    return cat_specs[0]["key"]


def parse_metadata(metadata_source=None, structure_dir=None, candidate_ids=None):
    """Main entry point: loads metadata from any format or extracts it from structures.
    
    Returns:
      {
        "taxa": { taxon_id: { ...fields... } },
        "columns": [ ...column_specs... ],
        "default_color_col": "..."
      }
    """
    raw_records = []

    if metadata_source and os.path.isfile(metadata_source):
        ext = os.path.splitext(metadata_source)[1].lower()
        if ext in (".xlsx", ".xls"):
            raw_records = read_xlsx_native(metadata_source)
        elif ext in (".csv", ".tsv", ".tab", ".txt"):
            raw_records = read_delimited(metadata_source)
        elif ext == ".json":
            raw_records = read_json_file(metadata_source)
        else:
            raise ValueError(f"Unsupported metadata file extension: '{ext}'")

    if not raw_records and structure_dir and os.path.isdir(structure_dir):
        print(f"[Metadata] No metadata file provided. Extracting structural metrics from '{structure_dir}'...")
        raw_records = extract_metadata_from_structures(structure_dir)

    if not raw_records:
        return {"taxa": {}, "columns": [], "default_color_col": None}

    # Identify ID column
    id_col = find_best_id_column(raw_records, candidate_ids)
    print(f"[Metadata] Identified ID column: '{id_col}' ({len(raw_records)} records).")

    # Classify schema
    column_specs = classify_columns(raw_records, id_col)
    default_col = choose_default_color_column(column_specs)

    # Format taxa index
    taxa_map = {}
    for r in raw_records:
        tid = str(r.get(id_col, "")).strip()
        if not tid:
            continue
        # Format values cleanly
        clean_record = {}
        for k, v in r.items():
            if k == id_col:
                continue
            if v is None:
                continue
            # Try parsing float / int
            try:
                fv = float(v)
                clean_record[k] = int(fv) if fv.is_integer() else round(fv, 3)
            except (ValueError, TypeError):
                clean_record[k] = str(v).strip()
        taxa_map[tid] = clean_record

    return {
        "taxa": taxa_map,
        "columns": column_specs,
        "default_color_col": default_col
    }


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description="Universal metadata parser and schema classifier.")
    parser.add_argument("-i", "--input", default=None, help="Path to .xlsx, .csv, .tsv, or .json metadata file")
    parser.add_argument("-s", "--structures", default=None, help="Path to structure files folder")
    parser.add_argument("-o", "--output", default="taxa_metadata_index.json", help="Output JSON path")
    args = parser.parse_args()

    result = parse_metadata(args.input, args.structures)
    with open(args.output, "w", encoding="utf-8") as f:
        json.dump(result, f, indent=2)
    print(f"[Metadata] Saved index to '{args.output}' with {len(result['taxa'])} taxa and {len(result['columns'])} columns.")
