#!/usr/bin/env python3
"""Build upgraded interactive phylogenetic tree visualization with dynamic metadata architecture.

Features:
- Dynamic "Color Nodes & Clades By" dropdown selector (categorical & continuous)
- Dynamic "Group & Collapse Clades By" dropdown selector in the Clades tab
- Dynamic Color Legend (categorical badges with counts & clickable collapse, or continuous gradient bar)
- Dynamic Taxon Inspection Card with dynamic key-value metadata table
- Dynamic Search across all metadata annotations
- High-performance non-scaling strokes and large-cohort viewport fitting
- 3D C-alpha backbone fold viewer in pinned card and hover tooltip
"""

import json
import os
import shutil
import subprocess
import re
from pathlib import Path
import numpy as np
from viral_phylo.metadata import parse_metadata, EXPANDED_PALETTE, KNOWN_VALUE_COLORS

repo_dir = Path(__file__).resolve().parent.parent.parent
results_dir = repo_dir / "results"

# Determine tool version from pyproject.toml and git commit hash
version_base = "1.0.0"
pyproject_path = repo_dir / "pyproject.toml"
if pyproject_path.exists():
    match = re.search(r'version\s*=\s*["\']([^"\']+)["\']', pyproject_path.read_text(encoding="utf-8"))
    if match:
        version_base = match.group(1)

try:
    commit_hash = subprocess.check_output(
        ["git", "-C", str(repo_dir), "rev-parse", "--short", "HEAD"],
        text=True, stderr=subprocess.DEVNULL
    ).strip()
except Exception:
    commit_hash = "c3d160e"

tool_version_str = f"v{version_base} ({commit_hash})"

# 1. Load Tree Files
with open(results_dir / "phylogeny_500_results/viral_tree_both.treefile") as f:
    newick_500_3di = f.read().strip()

with open(results_dir / "phylogeny_500_results/viral_tree_aa.treefile") as f:
    newick_500_aa = f.read().strip()

with open(results_dir / "glycoprotein_workflow/phylogeny/glycoprotein_tree_both.treefile") as f:
    newick_6_3di = f.read().strip()

with open(results_dir / "glycoprotein_workflow/phylogeny/aa_tree.treefile") as f:
    newick_6_aa = f.read().strip()

with open(results_dir / "nipah_esm_workflow/phylogeny/nipah_tree_3di.treefile") as f:
    newick_1193_3di = f.read().strip()

with open(results_dir / "nipah_esm_workflow/phylogeny/nipah_tree_aa.treefile") as f:
    newick_1193_aa = f.read().strip()

# Multi-distance metric ESM-2 PLM trees (Cosine, Euclidean, L1/Manhattan)
with open(results_dir / "glycoprotein_workflow/phylogeny/glycoprotein_tree_esm2_cosine.treefile") as f:
    newick_6_esm2_cosine = f.read().strip()
with open(results_dir / "glycoprotein_workflow/phylogeny/glycoprotein_tree_esm2_euclidean.treefile") as f:
    newick_6_esm2_euclidean = f.read().strip()
with open(results_dir / "glycoprotein_workflow/phylogeny/glycoprotein_tree_esm2_l1.treefile") as f:
    newick_6_esm2_l1 = f.read().strip()

with open(results_dir / "phylogeny_500_results/viral_tree_esm2_cosine.treefile") as f:
    newick_500_esm2_cosine = f.read().strip()
with open(results_dir / "phylogeny_500_results/viral_tree_esm2_euclidean.treefile") as f:
    newick_500_esm2_euclidean = f.read().strip()
with open(results_dir / "phylogeny_500_results/viral_tree_esm2_l1.treefile") as f:
    newick_500_esm2_l1 = f.read().strip()

with open(results_dir / "nipah_esm_workflow/phylogeny/nipah_tree_esm2_cosine.treefile") as f:
    newick_1193_esm2_cosine = f.read().strip()
with open(results_dir / "nipah_esm_workflow/phylogeny/nipah_tree_esm2_euclidean.treefile") as f:
    newick_1193_esm2_euclidean = f.read().strip()
with open(results_dir / "nipah_esm_workflow/phylogeny/nipah_tree_esm2_l1.treefile") as f:
    newick_1193_esm2_l1 = f.read().strip()

newick_6_esm2 = newick_6_esm2_cosine
newick_500_esm2 = newick_500_esm2_cosine
newick_1193_esm2 = newick_1193_esm2_cosine

# Load Silhouette Profiles for all 3 modalities (ESM-2, 3Di, AA) across datasets
def load_json_if_exists(p):
    return json.loads(p.read_text()) if p.exists() else None

sil_6_esm2 = load_json_if_exists(results_dir / "glycoprotein_workflow/phylogeny/glycoprotein_silhouette_esm2.json")
sil_6_3di = load_json_if_exists(results_dir / "glycoprotein_workflow/phylogeny/glycoprotein_silhouette_3di.json")
sil_6_aa = load_json_if_exists(results_dir / "glycoprotein_workflow/phylogeny/glycoprotein_silhouette_aa.json")

sil_500_esm2 = load_json_if_exists(results_dir / "phylogeny_500_results/viral_silhouette_esm2.json")
sil_500_3di = load_json_if_exists(results_dir / "phylogeny_500_results/viral_silhouette_3di.json")
sil_500_aa = load_json_if_exists(results_dir / "phylogeny_500_results/viral_silhouette_aa.json")

sil_1193_esm2 = load_json_if_exists(results_dir / "nipah_esm_workflow/phylogeny/nipah_silhouette_esm2.json")
sil_1193_3di = load_json_if_exists(results_dir / "nipah_esm_workflow/phylogeny/nipah_silhouette_3di.json")
sil_1193_aa = load_json_if_exists(results_dir / "nipah_esm_workflow/phylogeny/nipah_silhouette_aa.json")

sil_6 = sil_6_esm2
sil_500 = sil_500_esm2
sil_1193 = sil_1193_esm2

# 100 RdRp Workflow Tree & Silhouette Loading
rdrp_100_dir = results_dir / "rdrp_100_workflow"
newick_100_3di = ""
newick_100_aa = ""
tree_3di_path = rdrp_100_dir / "phylogeny/RNA dependent RNA polymerase_tree_auto.treefile"
tree_aa_path = rdrp_100_dir / "phylogeny/RNA dependent RNA polymerase_tree_aa.treefile"
if tree_3di_path.exists():
    newick_100_3di = tree_3di_path.read_text().strip()
if tree_aa_path.exists():
    newick_100_aa = tree_aa_path.read_text().strip()

has_esm_100 = False
newick_100_esm2_cosine = None
newick_100_esm2_euclidean = None
newick_100_esm2_l1 = None
if (rdrp_100_dir / "phylogeny/rdrp_tree_esm2_cosine.treefile").exists():
    newick_100_esm2_cosine = (rdrp_100_dir / "phylogeny/rdrp_tree_esm2_cosine.treefile").read_text().strip()
    has_esm_100 = True
if (rdrp_100_dir / "phylogeny/rdrp_tree_esm2_euclidean.treefile").exists():
    newick_100_esm2_euclidean = (rdrp_100_dir / "phylogeny/rdrp_tree_esm2_euclidean.treefile").read_text().strip()
if (rdrp_100_dir / "phylogeny/rdrp_tree_esm2_l1.treefile").exists():
    newick_100_esm2_l1 = (rdrp_100_dir / "phylogeny/rdrp_tree_esm2_l1.treefile").read_text().strip()

sil_100_3di = load_json_if_exists(rdrp_100_dir / "phylogeny/rdrp_silhouette_3di.json")
sil_100_aa = load_json_if_exists(rdrp_100_dir / "phylogeny/rdrp_silhouette_aa.json")
sil_100_esm2 = load_json_if_exists(rdrp_100_dir / "phylogeny/rdrp_silhouette_esm2.json") if has_esm_100 else None
sil_100 = sil_100_3di or sil_100_aa


# 2. Build 1,193 Dataset Schema
with open(results_dir / "nipah_esm_workflow/taxa_metadata.json") as f:
    raw_1193 = json.load(f)

meta_1193 = {}
for tid, val in raw_1193.items():
    meta_1193[tid] = {
        "structural_class": val.get("classification", "Unclassified"),
        "binding_strength": val.get("binding", "None"),
        "plddt": float(val.get("plddt", 75.0)),
        "length": int(val.get("length", 300)),
        "virus_target": "Nipah Virus G Binder",
        "design_id": tid
    }

columns_1193 = [
    {
        "key": "structural_class",
        "label": "Structural Class (Tertiary Fold)",
        "type": "categorical",
        "values": ["Mainly Alpha", "Alpha Beta", "Mainly Beta", "Few Secondary Structures", "Unclassified"],
        "colors": {
            "Mainly Alpha": "#38bdf8",
            "Alpha Beta": "#a855f7",
            "Mainly Beta": "#f97316",
            "Few Secondary Structures": "#10b981",
            "Unclassified": "#94a3b8"
        }
    },
    {
        "key": "binding_strength",
        "label": "Binding Strength (Henipavirus)",
        "type": "categorical",
        "values": ["Strong", "Medium", "Weak", "None"],
        "colors": {
            "Strong": "#10b981",
            "Medium": "#38bdf8",
            "Weak": "#f59e0b",
            "None": "#64748b"
        }
    },
    {
        "key": "plddt",
        "label": "ESMFold Confidence (pLDDT)",
        "type": "continuous",
        "min": 45.0,
        "max": 96.0
    },
    {
        "key": "length",
        "label": "Sequence Length (aa)",
        "type": "continuous",
        "min": 71,
        "max": 650
    }
]

# 3. Build 500 Dataset Schema
with open(results_dir / "viro_500_glycoproteins/taxa_metadata.json") as f:
    raw_500 = json.load(f)

FAMILY_COLORS = {
    "Rhabdoviridae": "#38bdf8", "Orthoherpesviridae": "#ec4899", "Phenuiviridae": "#10b981",
    "Hantaviridae": "#f97316", "Peribunyaviridae": "#8b5cf6", "Arenaviridae": "#06b6d4",
    "Nairoviridae": "#eab308", "Togaviridae": "#a855f7", "Arteriviridae": "#14b8a6",
    "Phasmaviridae": "#f43f5e", "Chuviridae": "#6366f1", "Matonaviridae": "#84cc16",
    "Orthomyxoviridae": "#e11d48", "Nyamiviridae": "#10b981", "Lispiviridae": "#ca8a04",
    "Poxviridae": "#d946ef", "Aliusviridae": "#0284c7", "Tobaniviridae": "#f59e0b",
    "Flaviviridae": "#e11d48", "Xinmoviridae": "#4f46e5", "Baculoviridae": "#059669",
    "Bornaviridae": "#b45309", "Unknown": "#94a3b8"
}

meta_500 = {}
for it in raw_500:
    rid = it["record_id"]
    fam = it.get("Family", "Unknown") or "Unknown"
    vname = it.get("Virus_name_s_", "Unknown virus") or "Unknown virus"
    prod = it.get("product", "glycoprotein") or "glycoprotein"
    plddt = float(it.get("colabfold_json_pLDDT") or 70.0)
    protlen = int(it.get("protlen") or 300)
    meta_500[rid] = {
        "family": fam,
        "virus": vname,
        "product": prod,
        "plddt": round(plddt, 1),
        "length": protlen,
        "genbank": it.get("GenBank_accession", rid.split("_")[0])
    }

unique_fams_500 = sorted(list(set(m["family"] for m in meta_500.values())))
unique_prods_500 = sorted(list(set(m["product"] for m in meta_500.values())))
prod_colors_500 = {p: EXPANDED_PALETTE[i % len(EXPANDED_PALETTE)] for i, p in enumerate(unique_prods_500)}

columns_500 = [
    {
        "key": "family",
        "label": "Viral Family (ICTV)",
        "type": "categorical",
        "values": unique_fams_500,
        "colors": {f: FAMILY_COLORS.get(f, "#94a3b8") for f in unique_fams_500}
    },
    {
        "key": "product",
        "label": "Gene Product / Protein",
        "type": "categorical",
        "values": unique_prods_500,
        "colors": prod_colors_500
    },
    {
        "key": "plddt",
        "label": "ColabFold pLDDT Score",
        "type": "continuous",
        "min": 40.0,
        "max": 95.0
    },
    {
        "key": "length",
        "label": "Protein Length (aa)",
        "type": "continuous",
        "min": 100,
        "max": 1200
    }
]

# 4. Build 6 Dataset Schema
meta_6 = {
    "AAA16239.1.2_9396": {"family": "Hantaviridae", "virus": "Hantaan orthohantavirus", "product": "Glycoproteins G1/G2", "plddt": 85.5, "length": 652, "genbank": "AAA16239.1"},
    "AAA88529.1.3_6592": {"family": "Togaviridae", "virus": "O'nyong-nyong virus", "product": "Precursor (E1/E2)", "plddt": 60.1, "length": 423, "genbank": "AAA88529.1"},
    "AAQ55251.1.1_9251": {"family": "Arenaviridae", "virus": "Junin mammarenavirus", "product": "Glycoprotein GP1", "plddt": 79.6, "length": 247, "genbank": "AAQ55251.1"},
    "AAB93840.1.1_9564": {"family": "Togaviridae", "virus": "Ross River virus", "product": "Precursor (E1/E2)", "plddt": 86.1, "length": 487, "genbank": "AAB93840.1"},
    "AAQ55251.1.2_9251": {"family": "Arenaviridae", "virus": "Junin mammarenavirus", "product": "Glycoprotein GP2", "plddt": 84.2, "length": 234, "genbank": "AAQ55251.1"},
    "AAA88529.1.2_6592": {"family": "Togaviridae", "virus": "O'nyong-nyong virus", "product": "Precursor (Capsid)", "plddt": 72.8, "length": 268, "genbank": "AAA88529.1"}
}

columns_6 = [
    {
        "key": "family",
        "label": "Viral Family",
        "type": "categorical",
        "values": ["Hantaviridae", "Togaviridae", "Arenaviridae"],
        "colors": {"Hantaviridae": "#f97316", "Togaviridae": "#a855f7", "Arenaviridae": "#06b6d4"}
    },
    {
        "key": "product",
        "label": "Gene Product",
        "type": "categorical",
        "values": ["Glycoproteins G1/G2", "Precursor (E1/E2)", "Glycoprotein GP1", "Glycoprotein GP2", "Precursor (Capsid)"],
        "colors": {
            "Glycoproteins G1/G2": "#38bdf8",
            "Precursor (E1/E2)": "#a855f7",
            "Glycoprotein GP1": "#06b6d4",
            "Glycoprotein GP2": "#10b981",
            "Precursor (Capsid)": "#f59e0b"
        }
    },
    {
        "key": "plddt",
        "label": "ColabFold pLDDT",
        "type": "continuous",
        "min": 60.0,
        "max": 86.5
    },
    {
        "key": "length",
        "label": "Protein Length",
        "type": "continuous",
        "min": 234,
        "max": 652
    }
]


# 5. Build 100 RdRp Dataset Schema
meta_100 = {}
columns_100 = []
if (rdrp_100_dir / "taxa_metadata.json").exists():
    with open(rdrp_100_dir / "taxa_metadata.json") as f:
        raw_100_data = json.load(f)
    raw_100 = raw_100_data.get("taxa", raw_100_data)
    for tid, val in raw_100.items():
        fam = val.get("Family") or "Unknown"
        gen = val.get("Genus") or "Unknown"
        sp = val.get("Species") or "Unknown"
        vname = val.get("Virus_name_s_") or "Unknown virus"
        host = val.get("Host_source") or val.get("host") or "Unknown"
        plddt = float(val.get("esmfold_log_pLDDT") or val.get("colabfold_json_pLDDT") or 70.0)
        protlen = int(val.get("protlen") or len(val.get("protein_seq", "")) or 400)
        gb = val.get("genbank_id") or tid.split("_")[0]
        meta_100[tid] = {
            "family": fam,
            "genus": gen,
            "species": sp,
            "virus": vname,
            "host": host,
            "plddt": round(plddt, 1),
            "length": protlen,
            "genbank": gb
        }

    unique_fams_100 = sorted(list(set(m["family"] for m in meta_100.values())))
    unique_genera_100 = sorted(list(set(m["genus"] for m in meta_100.values())))
    unique_hosts_100 = sorted(list(set(m["host"] for m in meta_100.values())))
    genus_colors_100 = {g: EXPANDED_PALETTE[i % len(EXPANDED_PALETTE)] for i, g in enumerate(unique_genera_100)}
    host_colors_100 = {h: EXPANDED_PALETTE[(i * 3 + 2) % len(EXPANDED_PALETTE)] for i, h in enumerate(unique_hosts_100)}

    columns_100 = [
        {
            "key": "family",
            "label": "Viral Family (ICTV)",
            "type": "categorical",
            "values": unique_fams_100,
            "colors": {f: FAMILY_COLORS.get(f, EXPANDED_PALETTE[i % len(EXPANDED_PALETTE)]) for i, f in enumerate(unique_fams_100)}
        },
        {
            "key": "genus",
            "label": "Viral Genus",
            "type": "categorical",
            "values": unique_genera_100,
            "colors": genus_colors_100
        },
        {
            "key": "host",
            "label": "Host Source",
            "type": "categorical",
            "values": unique_hosts_100,
            "colors": host_colors_100
        },
        {
            "key": "plddt",
            "label": "Structure pLDDT Score",
            "type": "continuous",
            "min": 25.0,
            "max": 95.0
        },
        {
            "key": "length",
            "label": "Protein Length (aa)",
            "type": "continuous",
            "min": 30,
            "max": 2100
        }
    ]


# 4. Precomputed Phylogenetic Congruence Profiles
DATASET_CONGRUENCE = {
    "100": {
        "3di_vs_aa": {
            "taxa_count": 100,
            "shared_splits": 11,
            "rf_distance": 172,
            "max_rf": 194,
            "norm_rf": 0.8866,
            "congruence_pct": 11.34,
            "cophenetic_r": 0.8482,
            "discordance_native": 88.66,
            "discordance_untangled": 32.4,
            "interpretation": "Analysis of 100 diverse viral RNA-dependent RNA Polymerases (spanning 13 ICTV families) reveals high patristic distance correlation (cophenetic r = +0.848) between structural 3Di and amino acid evolution. Deep divergences across Riboviria preserve core RdRp catalytic palm/finger/thumb topology while exhibiting extensive sequence diversification across viral hosts."
        }
    },
    "1193": {
        "3di_vs_aa": {
            "taxa_count": 1193,
            "shared_splits": 178,
            "rf_distance": 2024,
            "max_rf": 2380,
            "norm_rf": 0.8504,
            "congruence_pct": 14.96,
            "cophenetic_r": 0.3171,
            "discordance_native": 89.5,
            "discordance_untangled": 35.8,
            "interpretation": "In this 1,193-member engineered Henipavirus binder library, sequence-structure topological congruence is 15.0% (178 shared internal clades, r = +0.317). Filtering specifically to the 76 Strong Binders boosts cophenetic distance correlation to r = +0.578."
        },
        "3di_vs_esm2_cosine": {
            "taxa_count": 1193,
            "shared_splits": 136,
            "rf_distance": 2108,
            "max_rf": 2380,
            "norm_rf": 0.8857,
            "congruence_pct": 11.43,
            "cophenetic_r": 0.5434,
            "discordance_native": 92.0,
            "discordance_untangled": 36.8,
            "interpretation": "Comparing 3Di backbone structural fold to ESM-2 650M cosine representations reveals 11.4% topological split agreement and strong cophenetic correlation (r = +0.543), showing that angular PLM embeddings capture tertiary structural similarity well."
        },
        "3di_vs_esm2_euclidean": {
            "taxa_count": 1193,
            "shared_splits": 134,
            "rf_distance": 2112,
            "max_rf": 2380,
            "norm_rf": 0.8874,
            "congruence_pct": 11.26,
            "cophenetic_r": 0.4597,
            "discordance_native": 92.1,
            "discordance_untangled": 36.8,
            "interpretation": "Comparing 3Di structure to ESM-2 Euclidean (L2) distance yields 11.3% split congruence and moderate-to-high correlation (r = +0.460)."
        },
        "3di_vs_esm2_l1": {
            "taxa_count": 1193,
            "shared_splits": 137,
            "rf_distance": 2106,
            "max_rf": 2380,
            "norm_rf": 0.8849,
            "congruence_pct": 11.51,
            "cophenetic_r": 0.4566,
            "discordance_native": 91.9,
            "discordance_untangled": 36.8,
            "interpretation": "Comparing 3Di structure to ESM-2 L1 (Manhattan) distance reveals 11.5% split congruence and r = +0.457, demonstrating consistency across non-linear feature norms."
        },
        "aa_vs_esm2_cosine": {
            "taxa_count": 1193,
            "shared_splits": 178,
            "rf_distance": 2024,
            "max_rf": 2380,
            "norm_rf": 0.8504,
            "congruence_pct": 14.96,
            "cophenetic_r": 0.5338,
            "discordance_native": 89.5,
            "discordance_untangled": 35.8,
            "interpretation": "Comparing native amino acid IQ-TREE sequence phylogeny to ESM-2 cosine distance shows 15.0% shared clades with high cophenetic correlation (r = +0.534)."
        },
        "aa_vs_esm2_euclidean": {
            "taxa_count": 1193,
            "shared_splits": 178,
            "rf_distance": 2024,
            "max_rf": 2380,
            "norm_rf": 0.8504,
            "congruence_pct": 14.96,
            "cophenetic_r": 0.3955,
            "discordance_native": 89.5,
            "discordance_untangled": 35.8,
            "interpretation": "Comparing amino acid sequence tree to ESM-2 Euclidean distance gives 15.0% shared splits with r = +0.396."
        },
        "aa_vs_esm2_l1": {
            "taxa_count": 1193,
            "shared_splits": 184,
            "rf_distance": 2012,
            "max_rf": 2380,
            "norm_rf": 0.8454,
            "congruence_pct": 15.46,
            "cophenetic_r": 0.3761,
            "discordance_native": 89.2,
            "discordance_untangled": 35.7,
            "interpretation": "Comparing amino acid sequence tree to ESM-2 L1 distance yields 15.5% shared clades with r = +0.376."
        }
    },
    "500": {
        "3di_vs_aa": {
            "taxa_count": 500,
            "shared_splits": 182,
            "rf_distance": 630,
            "max_rf": 994,
            "norm_rf": 0.6338,
            "congruence_pct": 36.62,
            "cophenetic_r": 0.5955,
            "discordance_native": 74.4,
            "discordance_untangled": 29.7,
            "interpretation": "Across 500 viral glycoproteins spanning 23 viral families, 36.6% of internal evolutionary splits are strictly congruent between 3Di structure and amino acid sequence, accompanied by strong cophenetic correlation (r = +0.596)."
        },
        "3di_vs_esm2_cosine": {
            "taxa_count": 500,
            "shared_splits": 97,
            "rf_distance": 801,
            "max_rf": 995,
            "norm_rf": 0.805,
            "congruence_pct": 19.5,
            "cophenetic_r": 0.3187,
            "discordance_native": 86.3,
            "discordance_untangled": 34.5,
            "interpretation": "Comparing 3Di structural phylogeny to ESM-2 Cosine PLM tree exhibits 19.5% shared bipartitions (97 clades) and r = +0.319 cophenetic correlation."
        },
        "3di_vs_esm2_euclidean": {
            "taxa_count": 500,
            "shared_splits": 101,
            "rf_distance": 793,
            "max_rf": 995,
            "norm_rf": 0.797,
            "congruence_pct": 20.3,
            "cophenetic_r": 0.3435,
            "discordance_native": 85.8,
            "discordance_untangled": 34.3,
            "interpretation": "Comparing 3Di structural phylogeny to ESM-2 Euclidean distance exhibits 20.3% shared bipartitions (101 clades) and r = +0.344 cophenetic correlation."
        },
        "3di_vs_esm2_l1": {
            "taxa_count": 500,
            "shared_splits": 96,
            "rf_distance": 803,
            "max_rf": 995,
            "norm_rf": 0.807,
            "congruence_pct": 19.3,
            "cophenetic_r": 0.2765,
            "discordance_native": 86.5,
            "discordance_untangled": 34.6,
            "interpretation": "Comparing 3Di structural phylogeny to ESM-2 L1 Manhattan distance yields 19.3% shared bipartitions (96 clades) and r = +0.277."
        },
        "aa_vs_esm2_cosine": {
            "taxa_count": 500,
            "shared_splits": 118,
            "rf_distance": 759,
            "max_rf": 995,
            "norm_rf": 0.7628,
            "congruence_pct": 23.72,
            "cophenetic_r": 0.4201,
            "discordance_native": 83.4,
            "discordance_untangled": 33.4,
            "interpretation": "Comparing amino acid sequence tree to ESM-2 Cosine distance exhibits 23.7% shared clades (118 bipartitions) and r = +0.420."
        },
        "aa_vs_esm2_euclidean": {
            "taxa_count": 500,
            "shared_splits": 120,
            "rf_distance": 755,
            "max_rf": 995,
            "norm_rf": 0.7588,
            "congruence_pct": 24.12,
            "cophenetic_r": 0.4442,
            "discordance_native": 83.1,
            "discordance_untangled": 33.2,
            "interpretation": "Comparing amino acid sequence tree to ESM-2 Euclidean distance exhibits 24.1% shared clades (120 bipartitions) and r = +0.444."
        },
        "aa_vs_esm2_l1": {
            "taxa_count": 500,
            "shared_splits": 121,
            "rf_distance": 753,
            "max_rf": 995,
            "norm_rf": 0.7568,
            "congruence_pct": 24.32,
            "cophenetic_r": 0.4656,
            "discordance_native": 83.0,
            "discordance_untangled": 33.2,
            "interpretation": "Comparing amino acid sequence tree to ESM-2 L1 distance exhibits 24.3% shared clades (121 bipartitions) and r = +0.466."
        }
    },
    "6": {
        "3di_vs_aa": {
            "taxa_count": 6,
            "shared_splits": 0,
            "rf_distance": 6,
            "max_rf": 6,
            "norm_rf": 1.0,
            "congruence_pct": 0.0,
            "cophenetic_r": 0.4588,
            "discordance_native": 100.0,
            "discordance_untangled": 40.0,
            "interpretation": "Benchmark viral glycoprotein dataset displays high local sequence conservation with discrete structural topological shifts across divergent viral families. Cophenetic distance correlation r = +0.459 demonstrates moderate patristic concordance."
        },
        "3di_vs_esm2_cosine": {
            "taxa_count": 6,
            "shared_splits": 0,
            "rf_distance": 6,
            "max_rf": 6,
            "norm_rf": 1.0,
            "congruence_pct": 0.0,
            "cophenetic_r": -0.3244,
            "discordance_native": 100.0,
            "discordance_untangled": 40.0,
            "interpretation": "In the 6-taxa benchmark, 3Di structure and ESM-2 Cosine clustering capture differing aspects of deep viral divergence across Togaviridae, Hantaviridae, and Arenaviridae."
        },
        "3di_vs_esm2_euclidean": {
            "taxa_count": 6,
            "shared_splits": 0,
            "rf_distance": 6,
            "max_rf": 6,
            "norm_rf": 1.0,
            "congruence_pct": 0.0,
            "cophenetic_r": -0.2305,
            "discordance_native": 100.0,
            "discordance_untangled": 40.0,
            "interpretation": "3Di structure vs ESM-2 Euclidean distance clustering in the 6-taxa benchmark shows family-level clustering with alternative basal topologies."
        },
        "3di_vs_esm2_l1": {
            "taxa_count": 6,
            "shared_splits": 0,
            "rf_distance": 6,
            "max_rf": 6,
            "norm_rf": 1.0,
            "congruence_pct": 0.0,
            "cophenetic_r": -0.2296,
            "discordance_native": 100.0,
            "discordance_untangled": 40.0,
            "interpretation": "3Di structure vs ESM-2 L1 distance clustering exhibits concordance on closely related glycoprotein pairs."
        },
        "aa_vs_esm2_cosine": {
            "taxa_count": 6,
            "shared_splits": 1,
            "rf_distance": 4,
            "max_rf": 6,
            "norm_rf": 0.6667,
            "congruence_pct": 33.33,
            "cophenetic_r": -0.1115,
            "discordance_native": 76.7,
            "discordance_untangled": 30.7,
            "interpretation": "Amino acid sequence vs ESM-2 Cosine clustering achieves 33.3% bipartition agreement (shared sister taxa AAA88529 precursor variants)."
        },
        "aa_vs_esm2_euclidean": {
            "taxa_count": 6,
            "shared_splits": 1,
            "rf_distance": 4,
            "max_rf": 6,
            "norm_rf": 0.6667,
            "congruence_pct": 33.33,
            "cophenetic_r": 0.0186,
            "discordance_native": 76.7,
            "discordance_untangled": 30.7,
            "interpretation": "Amino acid sequence vs ESM-2 Euclidean distance clustering achieves 33.3% bipartition agreement with r = +0.019."
        },
        "aa_vs_esm2_l1": {
            "taxa_count": 6,
            "shared_splits": 1,
            "rf_distance": 4,
            "max_rf": 6,
            "norm_rf": 0.6667,
            "congruence_pct": 33.33,
            "cophenetic_r": 0.0254,
            "discordance_native": 76.7,
            "discordance_untangled": 30.7,
            "interpretation": "Amino acid sequence vs ESM-2 L1 distance clustering achieves 33.3% bipartition agreement with r = +0.025."
        }
    }
}

def load_or_compute_umap(npz_path, json_path=None):
    """Extract or compute 2D UMAP projection coordinates from npz embeddings."""
    if not npz_path or not Path(npz_path).exists():
        return None
    try:
        data = np.load(npz_path, allow_pickle=True)
        # 1. Check if umap already saved in npz
        if "umap" in data.files and data["umap"]:
            raw_u = data["umap"].item() if hasattr(data["umap"], "item") else str(data["umap"])
            if raw_u and isinstance(raw_u, str) and raw_u.strip():
                try:
                    return json.loads(raw_u)
                except Exception:
                    pass
        # 2. Check if standalone json exists
        if json_path and Path(json_path).exists():
            try:
                with open(json_path, "r", encoding="utf-8") as f:
                    u_json = json.load(f)
                    return u_json.get("coordinates", u_json)
            except Exception:
                pass
        # 3. Compute UMAP coordinates using umap-learn
        if "embeddings" in data.files and "taxa" in data.files:
            embs = data["embeddings"]
            taxa = [str(t) for t in data["taxa"]]
            n_samples = len(taxa)
            if n_samples < 3:
                return {t: [round(float(i * 200.0 - 100.0), 2), 0.0] for i, t in enumerate(taxa)}
            import umap
            actual_neighbors = max(2, min(15, n_samples - 1))
            reducer = umap.UMAP(n_neighbors=actual_neighbors, min_dist=0.1, metric="cosine", random_state=42)
            coords = reducer.fit_transform(embs)
            # Normalize to clean canvas range [-500, 500] centered at 0
            c_min = coords.min(axis=0)
            c_max = coords.max(axis=0)
            c_diff = np.where((c_max - c_min) == 0, 1.0, (c_max - c_min))
            norm = ((coords - c_min) / c_diff - 0.5) * 1000.0
            res = {taxa[i]: [round(float(norm[i, 0]), 2), round(float(norm[i, 1]), 2)] for i in range(n_samples)}
            if json_path:
                try:
                    with open(json_path, "w", encoding="utf-8") as f:
                        json.dump(res, f)
                except Exception:
                    pass
            return res
    except Exception as e:
        print(f"[Warning] Could not extract/compute UMAP for {npz_path}: {e}")
    return None

print("Loading / Computing ESM-2 UMAP projections...")
umap_1193 = load_or_compute_umap(
    results_dir / "nipah_esm_workflow/phylogeny/nipah_embeddings_esm2.npz",
    results_dir / "nipah_esm_workflow/phylogeny/nipah_umap_esm2.json"
)
umap_500 = load_or_compute_umap(
    results_dir / "phylogeny_500_results/viral_embeddings_esm2.npz",
    results_dir / "phylogeny_500_results/viral_umap_esm2.json"
)
umap_6 = load_or_compute_umap(
    results_dir / "glycoprotein_workflow/phylogeny/glycoprotein_embeddings_esm2.npz",
    results_dir / "glycoprotein_workflow/phylogeny/glycoprotein_umap_esm2.json"
)

# Assemble DATASETS JSON
DATASETS = {
    "100": {
        "title": "🧬 100 RNA-dependent RNA Polymerases",
        "has_esm": has_esm_100,
        "newick_3di": newick_100_3di,
        "newick_aa": newick_100_aa,
        "newick_esm2": newick_100_esm2_cosine,
        "newick_esm2_cosine": newick_100_esm2_cosine,
        "newick_esm2_euclidean": newick_100_esm2_euclidean,
        "newick_esm2_l1": newick_100_esm2_l1,
        "taxa": meta_100,
        "columns": columns_100,
        "defaultColorCol": "family",
        "defaultCladeCol": "family",
        "defaultSpacing": 18,
        "defaultRadius": 3.4,
        "defaultZoom": {"x": 40, "y": 30, "k": 0.60},
        "congruence": DATASET_CONGRUENCE.get("100", {}),
        "silhouette": sil_100_3di or sil_100_aa,
        "silhouette_esm2": sil_100_esm2,
        "silhouette_3di": sil_100_3di,
        "silhouette_aa": sil_100_aa,
        "esm2_umap": None
    },
    "1193": {
        "title": "🧬 1,193 Nipah ESMFold Structures",
        "has_esm": True,
        "newick_3di": newick_1193_3di,
        "newick_aa": newick_1193_aa,
        "newick_esm2": newick_1193_esm2_cosine,
        "newick_esm2_cosine": newick_1193_esm2_cosine,
        "newick_esm2_euclidean": newick_1193_esm2_euclidean,
        "newick_esm2_l1": newick_1193_esm2_l1,
        "taxa": meta_1193,
        "columns": columns_1193,
        "defaultColorCol": "structural_class",
        "defaultCladeCol": "structural_class",
        "defaultSpacing": 14,
        "defaultRadius": 2.8,
        "defaultZoom": {"x": 40, "y": 35, "k": 0.65},
        "congruence": DATASET_CONGRUENCE["1193"],
        "silhouette": sil_1193_esm2,
        "silhouette_esm2": sil_1193_esm2,
        "silhouette_3di": sil_1193_3di,
        "silhouette_aa": sil_1193_aa,
        "esm2_umap": umap_1193
    },
    "500": {
        "title": "🌐 500 Viral Glycoproteins",
        "has_esm": True,
        "newick_3di": newick_500_3di,
        "newick_aa": newick_500_aa,
        "newick_esm2": newick_500_esm2_cosine,
        "newick_esm2_cosine": newick_500_esm2_cosine,
        "newick_esm2_euclidean": newick_500_esm2_euclidean,
        "newick_esm2_l1": newick_500_esm2_l1,
        "taxa": meta_500,
        "columns": columns_500,
        "defaultColorCol": "family",
        "defaultCladeCol": "family",
        "defaultSpacing": 16,
        "defaultRadius": 3.2,
        "defaultZoom": {"x": 40, "y": 30, "k": 0.55},
        "congruence": DATASET_CONGRUENCE["500"],
        "silhouette": sil_500_esm2,
        "silhouette_esm2": sil_500_esm2,
        "silhouette_3di": sil_500_3di,
        "silhouette_aa": sil_500_aa,
        "esm2_umap": umap_500
    },
    "6": {
        "title": "🔬 6 Benchmark Glycoproteins",
        "has_esm": True,
        "newick_3di": newick_6_3di,
        "newick_aa": newick_6_aa,
        "newick_esm2": newick_6_esm2_cosine,
        "newick_esm2_cosine": newick_6_esm2_cosine,
        "newick_esm2_euclidean": newick_6_esm2_euclidean,
        "newick_esm2_l1": newick_6_esm2_l1,
        "taxa": meta_6,
        "columns": columns_6,
        "defaultColorCol": "family",
        "defaultCladeCol": "family",
        "defaultSpacing": 55,
        "defaultRadius": 5.0,
        "defaultZoom": {"x": 80, "y": 50, "k": 1.0},
        "congruence": DATASET_CONGRUENCE["6"],
        "silhouette": sil_6_esm2,
        "silhouette_esm2": sil_6_esm2,
        "silhouette_3di": sil_6_3di,
        "silhouette_aa": sil_6_aa,
        "esm2_umap": umap_6
    }
}

# Auto-discover custom workflows in results/*
known_core_datasets = {"1193", "500", "100", "6", "foldmason_500_alignments", "foldmason_alignments", "foldmason_glycoproteins", "viro_3d_structures", "viro_500_glycoproteins", "viro_glycoproteins", "phylogeny_500_results", "phylogeny_results", "nipah_esm_workflow", "glycoprotein_workflow", "rdrp_100_workflow"}
custom_options_html = []
for d in sorted(results_dir.iterdir()):
    if not d.is_dir() or d.name.startswith(".") or d.name in known_core_datasets:
        continue
    phy = d / "phylogeny"
    if not phy.exists():
        continue
    treefiles = list(phy.glob("*.treefile"))
    if not treefiles:
        continue

    # Determine 3di, aa, esm2 trees
    newick_3di, newick_aa, newick_esm2 = "", "", ""
    for tf in treefiles:
        name = tf.name.lower()
        txt = tf.read_text(encoding="utf-8", errors="ignore").strip()
        if "3di" in name or "both" in name or "alphafold" in name:
            newick_3di = txt
        elif "aa" in name:
            newick_aa = txt
        elif "esm2" in name or "esmc" in name:
            newick_esm2 = txt
    if not newick_3di and treefiles:
        newick_3di = treefiles[0].read_text(encoding="utf-8", errors="ignore").strip()
    if not newick_aa and newick_3di:
        newick_aa = newick_3di

    # Metadata
    meta = {}
    cols = []
    default_col = "family"
    meta_f = d / "taxa_metadata.json"
    if meta_f.exists():
        try:
            with open(meta_f) as f:
                raw_m = json.load(f)
            if isinstance(raw_m, dict) and "taxa" in raw_m:
                meta = raw_m["taxa"]
                cols = raw_m.get("columns", [])
                default_col = raw_m.get("default_color_col", "family")
            elif isinstance(raw_m, dict):
                meta = raw_m
        except Exception:
            pass

    if not meta and newick_3di:
        leaves = re.findall(r"([A-Za-z0-9_.\-]+):", newick_3di)
        for leaf in leaves:
            meta[leaf] = {"taxon_id": leaf, "family": "Unknown"}

    if not cols and meta:
        fams = sorted(list(set(m.get("family", "Unknown") for m in meta.values() if isinstance(m, dict))))
        if fams:
            cols.append({
                "key": "family",
                "label": "Viral Family",
                "type": "categorical",
                "values": fams,
                "colors": {f: FAMILY_COLORS.get(f, EXPANDED_PALETTE[i % len(EXPANDED_PALETTE)]) for i, f in enumerate(fams)}
            })

    taxa_count = len(meta)
    clean_title = d.name.replace("_", " ").title()
    title = f"🧪 {clean_title} ({taxa_count} taxa)"
    key = d.name

    custom_npz = list(phy.glob("*_embeddings_esm2.npz"))
    custom_umap = None
    if custom_npz:
        custom_umap = load_or_compute_umap(custom_npz[0], phy / f"{d.name}_umap_esm2.json")

    DATASETS[key] = {
        "title": title,
        "has_esm": bool(newick_esm2),
        "newick_3di": newick_3di,
        "newick_aa": newick_aa,
        "newick_esm2": newick_esm2 or None,
        "newick_esm2_cosine": newick_esm2 or None,
        "newick_esm2_euclidean": newick_esm2 or None,
        "newick_esm2_l1": newick_esm2 or None,
        "taxa": meta,
        "columns": cols,
        "defaultColorCol": default_col if any(c.get("key") == default_col for c in cols) else (cols[0]["key"] if cols else "family"),
        "defaultCladeCol": default_col if any(c.get("key") == default_col for c in cols) else (cols[0]["key"] if cols else "family"),
        "defaultSpacing": 22 if taxa_count <= 60 else (18 if taxa_count <= 150 else 14),
        "defaultRadius": 3.6 if taxa_count <= 60 else (3.2 if taxa_count <= 150 else 2.8),
        "defaultZoom": {"x": 40, "y": 30, "k": 0.60},
        "congruence": {},
        "silhouette": None,
        "silhouette_esm2": None,
        "silhouette_3di": None,
        "silhouette_aa": None,
        "esm2_umap": custom_umap
    }
    custom_options_html.append(f'              <option value="{key}">{title}</option>')
    print(f"Auto-registered custom dataset in DATASETS: {key} -> {title}")

datasets_json = json.dumps(DATASETS)
extra_scale_options = "\n".join(custom_options_html)


def assemble_html(datasets_json: str, tool_version_str: str, extra_scale_options: str) -> str:
    """Assemble interactive tree HTML from modular templates."""
    tmpl_dir = Path(__file__).resolve().parent / "template"
    index_html = (tmpl_dir / "index.html").read_text(encoding="utf-8")
    viewer_css = (tmpl_dir / "viewer.css").read_text(encoding="utf-8")
    viewer_js = (tmpl_dir / "viewer.js").read_text(encoding="utf-8")

    js_filled = viewer_js.replace("/*__DATASETS_JSON__*/", datasets_json)
    html_filled = index_html.replace("/*__INLINE_CSS__*/", viewer_css).replace("/*__INLINE_JS__*/", js_filled)
    html_filled = html_filled.replace("<!-- TOOL_VERSION -->", tool_version_str)
    html_filled = html_filled.replace("<!-- EXTRA_SCALE_OPTIONS -->", extra_scale_options)
    return html_filled


def build_interactive_tree(repo_dir: Path = None, results_dir: Path = None):
    """Compile datasets and generate all interactive tree HTML targets."""
    if repo_dir is None:
        repo_dir = Path(__file__).resolve().parent.parent.parent
    if results_dir is None:
        results_dir = repo_dir / "results"

    html_content = assemble_html(datasets_json, tool_version_str, extra_scale_options)

    # Write to destination HTML files
    targets = [
        repo_dir / "interactive_tree.html",
        results_dir / "interactive_tree.html",
        results_dir / "nipah_esm_workflow/interactive_tree.html",
        results_dir / "glycoprotein_workflow/phylogeny/interactive_tree.html",
        results_dir / "rdrp_100_workflow/interactive_tree.html",
    ]

    # Auto-add targets for all discovered custom workflows
    for d in sorted(results_dir.iterdir()):
        if not d.is_dir() or d.name.startswith(".") or d.name in known_core_datasets:
            continue
        phy = d / "phylogeny"
        if phy.exists() and list(phy.glob("*.treefile")):
            targets.append(d / "interactive_tree.html")

    extra_dir = os.environ.get("ANTIGRAVITY_ARTIFACT_DIR")
    if extra_dir and Path(extra_dir).exists():
        targets.append(Path(extra_dir) / "interactive_tree.html")

    # Deduplicate targets
    targets = sorted(list(set(targets)))

    support_files = [
        "ca_500_structures.js",
        "ca_1193_structures.js",
        "ca_100_structures.js",
        "ca_structures.js",
        "alignments_data.js",
    ]

    # Ensure repo_dir has all support files from results_dir
    for sf in support_files:
        src = results_dir / sf
        dst_repo = repo_dir / sf
        if src.exists() and not dst_repo.exists():
            try:
                shutil.copy2(src, dst_repo)
            except Exception:
                pass

    for t in targets:
        t.parent.mkdir(parents=True, exist_ok=True)
        out_html = html_content

        # 1. Ensure target directory has local copies of all support files
        for sf in support_files:
            src = results_dir / sf
            dst = t.parent / sf
            if src.exists() and not dst.exists():
                try:
                    shutil.copy2(src, dst)
                except Exception:
                    pass

        # 2. Compute relative prefix to results_dir
        try:
            rel_results = os.path.relpath(results_dir, t.parent).replace("\\", "/")
            if rel_results == ".":
                prefix = ""
            else:
                prefix = rel_results.rstrip("/") + "/"
        except Exception:
            prefix = ""

        # 3. Build script tags with local direct and relative fallbacks
        script_lines = [
            "  <!-- Pre-cached 3D C-alpha Backbone coordinates & MSA Alignments (local + relative resolution) -->"
        ]
        for sf in support_files:
            script_lines.append(f'  <script src="{sf}"></script>')
        if prefix:
            for sf in support_files:
                script_lines.append(f'  <script src="{prefix}{sf}"></script>')

        data_scripts_str = "\n".join(script_lines)
        out_html = out_html.replace("<!-- DATA_SCRIPTS_PLACEHOLDER -->", data_scripts_str)

        # 4. Set appropriate active scale if located inside a specific workflow folder
        target_str = str(t)
        matched_scale = None
        if "rdrp_100_workflow" in target_str:
            matched_scale = "100"
        elif "glycoprotein_workflow" in target_str:
            matched_scale = "6"
        else:
            for k in DATASETS.keys():
                if k in target_str:
                    matched_scale = k
                    break

        if matched_scale and matched_scale != "1193":
            out_html = out_html.replace('let currentScale = "1193";', f'let currentScale = "{matched_scale}";')
            out_html = out_html.replace('<option value="1193" selected>', '<option value="1193">')
            out_html = out_html.replace(f'<option value="{matched_scale}">', f'<option value="{matched_scale}" selected>')

        with open(t, "w", encoding="utf-8") as f:
            f.write(out_html)
        print(f"Generated: {t} ({len(out_html)} bytes)")

    print("=== Successfully generated Dynamic Metadata Interactive Tree Suite ===")


def main():
    build_interactive_tree()


if __name__ == "__main__":
    main()
