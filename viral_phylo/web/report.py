#!/usr/bin/env python3
"""Build a standalone interactive report for a single user-supplied phylogenetic tree.

The tree (Newick or NEXUS) and optional metadata table are embedded as a
``custom_tree`` dataset; the viewer normalises the tree text and matches
metadata rows to leaf names in the browser, so the same code path serves both
this CLI report and the in-browser "Load tree file" control.
"""

import html
import json
import re
from pathlib import Path
from typing import Optional

from viral_phylo import __version__
from viral_phylo.metadata import parse_metadata

TEMPLATE_DIR = Path(__file__).resolve().parent / "template"


def assemble_html(datasets_json: str, tool_version_str: str, extra_scale_options: str) -> str:
    """Assemble interactive tree HTML from modular templates."""
    index_html = (TEMPLATE_DIR / "index.html").read_text(encoding="utf-8")
    viewer_css = (TEMPLATE_DIR / "viewer.css").read_text(encoding="utf-8")
    viewer_js = (TEMPLATE_DIR / "viewer.js").read_text(encoding="utf-8")

    js_filled = viewer_js.replace("/*__DATASETS_JSON__*/", datasets_json)
    html_filled = index_html.replace("/*__INLINE_CSS__*/", viewer_css).replace("/*__INLINE_JS__*/", js_filled)
    html_filled = html_filled.replace("<!-- TOOL_VERSION -->", tool_version_str)
    html_filled = html_filled.replace("<!-- EXTRA_SCALE_OPTIONS -->", extra_scale_options)
    return html_filled


def newick_leaf_names(tree_text: str) -> list:
    """Approximate leaf names from Newick/NEXUS text (used only as an ID-column hint for metadata)."""
    text = re.sub(r"\[[^\]]*\]", "", tree_text)
    translate = {}
    m = re.search(r"\btranslate\b(.*?);", text, flags=re.IGNORECASE | re.DOTALL)
    if m:
        for key, val in re.findall(r"(\S+)\s+('(?:[^']|'')*'|\"[^\"]*\"|[^,\s]+)\s*(?:,|$)", m.group(1).strip()):
            translate[key] = val[1:-1].replace("''", "'") if val.startswith("'") else val.strip('"')
    tree_m = re.search(r"\btree\b[^=]*=\s*(.*?;)", text, flags=re.IGNORECASE | re.DOTALL)
    newick = tree_m.group(1) if tree_m else text
    names = re.findall(r"[(,]\s*('(?:[^']|'')*'|[^\s:,();]+)", newick)
    return [translate.get(n, n[1:-1].replace("''", "'") if n.startswith("'") else n) for n in names]


def build_tree_report(tree_path: str, metadata_path: Optional[str] = None,
                      output_path: Optional[str] = None, title: Optional[str] = None) -> Path:
    """Render an interactive HTML report for one tree file with optional metadata."""
    tree_file = Path(tree_path)
    if not tree_file.is_file():
        raise FileNotFoundError(f"Tree file '{tree_path}' does not exist.")
    tree_text = tree_file.read_text(encoding="utf-8", errors="ignore")
    if "(" not in tree_text:
        raise ValueError(f"'{tree_path}' does not look like a Newick or NEXUS tree file.")

    raw_metadata = None
    if metadata_path:
        if not Path(metadata_path).is_file():
            raise FileNotFoundError(f"Metadata file '{metadata_path}' does not exist.")
        raw_metadata = parse_metadata(metadata_path, candidate_ids=newick_leaf_names(tree_text))

    title = title or tree_file.stem
    datasets = {
        "custom": {
            "custom_tree": True,
            "title": f"📂 {title}",
            "tree_filename": tree_file.name,
            "tree_text": tree_text,
            "raw_metadata": raw_metadata,
            "metadata_filename": Path(metadata_path).name if metadata_path else None,
        }
    }
    # Escape "</" so embedded labels can never terminate the inline <script> block
    datasets_json = json.dumps(datasets).replace("</", "<\\/")
    option = f'              <option value="custom" selected>📂 {html.escape(title)}</option>'
    out_html = assemble_html(datasets_json, f"v{__version__}", option)
    out_html = out_html.replace("<!-- DATA_SCRIPTS_PLACEHOLDER -->", "")
    out_html = out_html.replace('let currentScale = "1193";', 'let currentScale = "custom";')

    out = Path(output_path) if output_path else tree_file.with_name(f"{tree_file.stem}_report.html")
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(out_html, encoding="utf-8")
    return out
