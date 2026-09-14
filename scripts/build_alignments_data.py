#!/usr/bin/env python3
"""Parse FoldMason AA and 3Di FASTA alignments and generate alignments_data.js.
Writes to:
- ./alignments_data.js
- nipah_esm_workflow/alignments_data.js
- glycoprotein_workflow/phylogeny/alignments_data.js
"""
import os
import json
from pathlib import Path

repo_dir = Path(__file__).resolve().parent.parent
results_dir = repo_dir / "results"

def parse_fasta(path):
    res = {}
    curr_id = None
    curr_seq = []
    with open(path, "r", encoding="utf-8") as f:
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                if curr_id:
                    res[curr_id] = "".join(curr_seq)
                curr_id = line[1:].split()[0]
                curr_seq = []
            else:
                curr_seq.append(line)
        if curr_id:
            res[curr_id] = "".join(curr_seq)
    return res

print("Parsing 1,193 Nipah ESMFold alignments...")
nipah_aa = parse_fasta(results_dir / "nipah_esm_workflow/alignment/foldmason.fasta_aa.fa")
nipah_3di = parse_fasta(results_dir / "nipah_esm_workflow/alignment/foldmason.fasta_3di.fa")

print("Parsing 6 Benchmark Glycoprotein alignments...")
bench_aa = parse_fasta(results_dir / "glycoprotein_workflow/alignment/foldmason.fasta_aa.fa")
bench_3di = parse_fasta(results_dir / "glycoprotein_workflow/alignment/foldmason.fasta_3di.fa")

print("Parsing 500 Viral Glycoprotein alignments...")
g500_aa = parse_fasta(results_dir / "foldmason_500_alignments/foldmason.fasta_aa.fa")
g500_3di = parse_fasta(results_dir / "foldmason_500_alignments/foldmason.fasta_3di.fa")

rdrp_100_aa_path = results_dir / "rdrp_100_workflow/alignment/foldmason.fasta_aa.fa"
rdrp_100_3di_path = results_dir / "rdrp_100_workflow/alignment/foldmason.fasta_3di.fa"
rdrp_aa = {}
rdrp_3di = {}
if rdrp_100_aa_path.exists() and rdrp_100_3di_path.exists():
    print("Parsing 100 RdRp alignments...")
    rdrp_aa = parse_fasta(rdrp_100_aa_path)
    rdrp_3di = parse_fasta(rdrp_100_3di_path)

alignments = {
    "1193": {
        "length": len(next(iter(nipah_aa.values()))),
        "taxa_count": len(nipah_aa),
        "aa": nipah_aa,
        "3di": nipah_3di
    },
    "500": {
        "length": len(next(iter(g500_aa.values()))),
        "taxa_count": len(g500_aa),
        "aa": g500_aa,
        "3di": g500_3di
    },
    "6": {
        "length": len(next(iter(bench_aa.values()))),
        "taxa_count": len(bench_aa),
        "aa": bench_aa,
        "3di": bench_3di
    }
}

if rdrp_aa and rdrp_3di:
    alignments["100"] = {
        "length": len(next(iter(rdrp_aa.values()))),
        "taxa_count": len(rdrp_aa),
        "aa": rdrp_aa,
        "3di": rdrp_3di
    }

js_content = "window.ALIGNMENTS_DATA = " + json.dumps(alignments) + ";\n"

targets = [
    repo_dir / "alignments_data.js",
    results_dir / "alignments_data.js",
    results_dir / "nipah_esm_workflow/alignments_data.js",
    results_dir / "glycoprotein_workflow/phylogeny/alignments_data.js",
    results_dir / "rdrp_100_workflow/alignments_data.js",
]
extra_dir = os.environ.get("ANTIGRAVITY_ARTIFACT_DIR")
if extra_dir and Path(extra_dir).exists():
    targets.append(Path(extra_dir) / "alignments_data.js")

for t in targets:
    t.parent.mkdir(parents=True, exist_ok=True)
    with open(t, "w", encoding="utf-8") as f:
        f.write(js_content)
    print(f"Generated {t} ({len(js_content)} bytes)")

print("=== Successfully generated alignments_data.js across all targets ===")
