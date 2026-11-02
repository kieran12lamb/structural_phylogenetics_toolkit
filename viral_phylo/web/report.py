#!/usr/bin/env python3
"""Build a standalone interactive report for a single user-supplied phylogenetic tree.

The tree (Newick or NEXUS) and optional metadata table are embedded as a
``custom_tree`` dataset; the viewer normalises the tree text and matches
metadata rows to leaf names in the browser, so the same code path serves both
this CLI report and the in-browser "Load tree file" control.
"""

import html
import json
import os
import re
import zipfile
from pathlib import Path
from typing import Optional, Union

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


_DEFAULT_SCALE_OPTIONS_BLOCK = (
    '              <option value="1193" selected>🧬 1,193 Nipah ESMFold Structures</option>\n'
    '              <option value="500">🌐 500 Viral Glycoproteins</option>\n'
    '              <option value="100">🧬 100 RNA-dependent RNA Polymerases</option>\n'
    '              <option value="6">🔬 6 Benchmark Glycoproteins</option>\n'
)


def build_run_dashboard(run_dir: Union[str, Path]) -> Optional[Path]:
    """Generate a self-contained ``interactive_tree.html`` and ``alignments_data.js`` inside ``run_dir``."""
    from viral_phylo.web.alignments import write_workflow_alignments_bundle
    from viral_phylo.web.builder import build_workflow_dataset_entry, tool_version_str

    d = Path(run_dir).resolve()
    if not d.is_dir():
        return None

    write_workflow_alignments_bundle(d)
    entry = build_workflow_dataset_entry(d)
    if not entry:
        return None

    key = d.name
    datasets = {key: entry}
    datasets_json = json.dumps(datasets).replace("</", "<\\/")
    option = f'              <option value="{html.escape(key)}" selected>{html.escape(entry["title"])}</option>'

    out_html = assemble_html(datasets_json, tool_version_str, option)
    out_html = out_html.replace(_DEFAULT_SCALE_OPTIONS_BLOCK, "")

    script_lines = []
    if (d / "ca_structures.js").is_file():
        script_lines.append('  <script src="ca_structures.js"></script>')
    if (d / "alignments_data.js").is_file():
        script_lines.append('  <script src="alignments_data.js"></script>')
    out_html = out_html.replace("<!-- DATA_SCRIPTS_PLACEHOLDER -->", "\n".join(script_lines))
    out_html = out_html.replace('let currentScale = "1193";', f'let currentScale = {json.dumps(key)};')

    out_file = d / "interactive_tree.html"
    out_file.write_text(out_html, encoding="utf-8")
    return out_file


_ALIGNMENT_SCRATCH_DIRS = {"tmpFolder", "_clustering_work", "_qc_structures", "structures", "_search"}
_ALIGNMENT_REPORT_EXTS = {".fa", ".fasta", ".faa", ".nw", ".nwk", ".html", ".json", ".tsv", ".csv", ".txt"}


def package_run_report(run_dir: Union[str, Path], zip_path: Optional[Union[str, Path]] = None) -> Path:
    """Package all report data for a completed pipeline run into a portable ``.zip`` archive.

    Produces a self-contained directory inside the archive (``<run_name>_report/``)
    containing:
      - ``interactive_tree.html`` (and ``index.html``) with embedded trees, metadata,
        clustering silhouette & UMAP projections, and run diagnostics
      - ``ca_structures.js`` (3D C-alpha backbones & confidence tracks)
      - ``alignments_data.js`` (whole-set and per-cluster 3Di & AA alignments)
      - ``run_manifest.json``, ``run_manifest.txt``, ``taxa_metadata.json``,
        ``multi_alignment_summary.json``
      - ``phylogeny/`` (all whole-set and per-cluster Newick trees, IQ-TREE reports,
        guide/start trees, PLM embeddings, UMAP coordinates, and silhouette profiles)
      - ``alignment/`` (all 3Di/AA FASTA alignments, FoldMason guide trees & HTML
        reports, and cluster summaries, excluding intermediate scratch databases)
    """
    d = Path(run_dir).resolve()
    if not d.is_dir():
        raise FileNotFoundError(f"Run directory '{run_dir}' does not exist.")

    # Ensure standalone per-run interactive HTML and alignments_data.js are up to date
    try:
        build_run_dashboard(d)
    except Exception as exc:
        print(f"[Report] [!] Notice: standalone run dashboard generation encountered: {exc}")

    out_zip = Path(zip_path).resolve() if zip_path else (d / f"{d.name}_report.zip")
    out_zip.parent.mkdir(parents=True, exist_ok=True)
    arc_root = Path(f"{d.name}_report")

    top_level_files = [
        "interactive_tree.html",
        "ca_structures.js",
        "alignments_data.js",
        "run_manifest.json",
        "run_manifest.txt",
        "taxa_metadata.json",
        "multi_alignment_summary.json",
        "search_summary.json",
    ]

    with zipfile.ZipFile(out_zip, "w", compression=zipfile.ZIP_DEFLATED, compresslevel=6) as zf:
        for fname in top_level_files:
            fpath = d / fname
            if fpath.is_file():
                zf.write(fpath, arcname=str(arc_root / fname))
                if fname == "interactive_tree.html":
                    zf.write(fpath, arcname=str(arc_root / "index.html"))

        # Include all phylogeny outputs (master & cluster trees, IQ-TREE logs/reports, embeddings)
        phy_dir = d / "phylogeny"
        if phy_dir.is_dir():
            for root, dirs, files in os.walk(phy_dir):
                dirs[:] = sorted(dn for dn in dirs if not dn.startswith("."))
                for fn in sorted(files):
                    if fn.startswith(".") or fn.endswith(".zip"):
                        continue
                    fp = Path(root) / fn
                    if fp.is_file():
                        rel = fp.relative_to(d)
                        zf.write(fp, arcname=str(arc_root / rel))

        # Include alignment FASTA/NW/HTML/JSON reports, skipping FoldMason/Foldseek scratch dirs
        aln_dir = d / "alignment"
        if aln_dir.is_dir():
            for root, dirs, files in os.walk(aln_dir):
                dirs[:] = sorted(
                    dn for dn in dirs
                    if not dn.startswith(".") and dn not in _ALIGNMENT_SCRATCH_DIRS
                )
                for fn in sorted(files):
                    if fn.startswith(".") or fn.startswith("sDB"):
                        continue
                    fp = Path(root) / fn
                    if fp.is_file() and fp.suffix.lower() in _ALIGNMENT_REPORT_EXTS:
                        rel = fp.relative_to(d)
                        zf.write(fp, arcname=str(arc_root / rel))

        readme_text = (
            f"Viral Structural Phylogenetics Toolkit — Run Report Archive ({d.name})\n"
            f"=====================================================================\n\n"
            f"Quick Start:\n"
            f"  1. Extract this .zip archive completely to a local folder.\n"
            f"  2. Open 'interactive_tree.html' (or 'index.html') in any modern web browser,\n"
            f"     or serve locally via:\n"
            f"       python3 -m http.server 8000\n"
            f"     and visit http://localhost:8000\n\n"
            f"Contents:\n"
            f"  - interactive_tree.html        Standalone interactive tree, 3D fold & MSA viewer\n"
            f"  - ca_structures.js             3D C-alpha backbone coordinates & confidence tracks\n"
            f"  - alignments_data.js           Whole-set & per-cluster 3Di / AA alignment bundle\n"
            f"  - taxa_metadata.json           Standardized taxon metadata & color schemas\n"
            f"  - multi_alignment_summary.json Similarity clustering & partition metrics\n"
            f"  - run_manifest.json / .txt     Stage-by-stage provenance, warnings & parameters\n"
            f"  - phylogeny/                   Whole-set & cluster Newick trees, IQ-TREE reports, UMAP\n"
            f"  - alignment/                   Whole-set & cluster 3Di/AA FASTA alignments & HTML reports\n"
        )
        zf.writestr(str(arc_root / "README.txt"), readme_text)

    return out_zip

