# Viral Structural Phylogenetics — Testing & Validation Guide

This document defines the automated test suite and a set of manual validation runs for the `viral-phylo` command-line tool. The manual runs use the 6 benchmark structures shipped in [`tests/fixtures/structures/`](fixtures/structures/), so everything except Test 1 works offline.

Run all commands **from the repository root** with the `spt` environment active (`conda activate spt`).

---

## 1. Environment & Prerequisite Check

```bash
python3 --version
foldmason version
iqtree --version
mafft --version
foldseek version
mmseqs version
VeryFastTree -help | head -1
```

### Expected Output:
- Python $\ge 3.10$
- FoldMason, MAFFT, Foldseek and MMseqs2 installed and executable
- IQ-TREE version $\ge 2.3.0$ (or 3.x)
- VeryFastTree (or FastTree) installed

`./install.sh` performs all of these checks, reports whether PyTorch can use the GPU, and runs the automated suite.

---

## 2. Automated Test Suite

```bash
pytest tests/ -v
# or, exactly as CI runs it:
python3 -m unittest discover -s tests -p "test_*.py" -v
```

| File | Covers |
| :--- | :--- |
| `test_alignment.py` | FASTA parsing, 3Di alphabet, gap stripping, matrix files |
| `test_cli.py` | CLI parsing, matrix discovery, embedding interpreter selection |
| `test_clustering.py` | Distance metrics, hierarchical clustering, silhouettes |
| `test_clustering_selection.py` | QC filter, whole-set gate, partition selection and suggestion, two-space silhouette scoring and disagreement, recursion schedule and shatter guard, 3Di/AA split judgement |
| `test_structures.py` | C-alpha extraction from PDB and mmCIF, AA-alignment resolution |
| `test_fasttree_matrix.py` | PAML → FastTree `-trans` matrix conversion; the matrix is really applied |
| `test_start_trees.py` | Starting-tree repair and validation; guided IQ-TREE with identical sequences; fallback when a guided run fails; preflight without FastTree |
| `test_parallel_trees.py` | Tree-job planning within the thread cap, concurrent largest-first scheduling, failure logs, interrupt stops child processes |
| `test_readiness_fixes.py` | Matrix aliasing, manifest never left `running`, viewer support labels, unique cluster names after recursion |
| `test_interactive_tree.py` | Dashboard assembled from the templates; per-alphabet alignment lengths |
| `test_metadata.py`, `test_report.py`, `test_viral_phylo_modules.py` | Metadata parsing, tree reports, module imports |

Integration tests that run IQ-TREE or VeryFastTree/FastTree are skipped when those tools are not installed.

---

## 3. Manual Validation Runs

### Test 1: Viro3D Structure Download (network)
```bash
viral-phylo fetch \
  --qualifier glycoprotein \
  --max-sequences 4 \
  --output-dir test_structures
```

**Validation Checks:**
- [ ] Exit code is `0`.
- [ ] 4 structure files exist in `test_structures/`, each containing `ATOM` coordinate records.
- [ ] `test_structures/taxa_metadata.json` exists.

---

### Test 2: FoldMason Structural Multiple Alignment
```bash
viral-phylo align \
  --folder tests/fixtures/structures \
  --output-dir test_alignment
```

**Validation Checks:**
- [ ] Exit code is `0`.
- [ ] `test_alignment/foldmason.fasta_3di.fa` contains 6 aligned sequences in the 3Di alphabet.
- [ ] `test_alignment/foldmason.fasta_aa.fa` has the same length as the 3Di alignment (FoldMason writes one alignment in two alphabets).
- [ ] `test_alignment/foldmason.fasta.nw` contains a Newick guide tree.

---

### Test 3: FoldMason Guide Tree Extraction
```bash
viral-phylo tree \
  --alignment test_alignment \
  --method foldmason \
  --output-dir test_phylogeny \
  --prefix guide_tree
```

**Validation Checks:**
- [ ] `test_phylogeny/guide_tree_foldmason.nwk` is generated and ends with `;`.

---

### Test 4: IQ-TREE with AlphaFold Matrix (`Q.3Di.AF`)
```bash
viral-phylo tree \
  --alignment test_alignment \
  --method iqtree \
  --matrix alphafold \
  --rate-heterogeneity +G4 \
  --tree-type 3di \
  --bootstrap 1000 \
  --alrt 1000 \
  --threads 2 \
  --output-dir test_phylogeny \
  --prefix af_result
```

**Validation Checks:**
- [ ] `test_phylogeny/af_result_alphafold.treefile` exists with branch lengths.
- [ ] `test_phylogeny/af_result_alphafold.contree` contains bootstrap support values.
- [ ] `test_phylogeny/af_result_alphafold.iqtree` reports the log-likelihood and Gamma shape $\alpha$.

---

### Test 5: IQ-TREE with ESMFold Matrix (`Q.3Di.LLM`)
```bash
viral-phylo tree \
  --alignment test_alignment \
  --method iqtree \
  --matrix esmfold \
  --rate-heterogeneity +G4 \
  --tree-type 3di \
  --bootstrap 1000 \
  --alrt 1000 \
  --threads 2 \
  --output-dir test_phylogeny \
  --prefix esm_result
```

**Validation Checks:**
- [ ] `test_phylogeny/esm_result_esmfold.treefile` is generated under the ESMFold matrix.

---

### Test 6: Automated Model & Rate Heterogeneity Selection
Tests ModelFinder across both matrices and all rate heterogeneity classes. This exercises the `-mset` path, which relies on the upper-case matrix copies in `MATRICES/`.

```bash
viral-phylo tree \
  --alignment test_alignment \
  --method iqtree \
  --matrix both \
  --rate-heterogeneity auto \
  --criterion BIC \
  --tree-type 3di \
  --bootstrap 1000 \
  --alrt 1000 \
  --threads 2 \
  --output-dir test_phylogeny \
  --prefix auto_result
```

**Validation Checks:**
- [ ] `test_phylogeny/auto_result_both.iqtree` contains a model ranking and a line `Best-fit model according to BIC: MATRICES/Q.3DI.AF...` or `...Q.3DI.LLM...`.
- [ ] `MATRICES/Q.3DI.AF` and `MATRICES/Q.3DI.LLM` exist (generated; ignored by git).

---

### Test 7: Dual Phylogeny & Cophylogenetic Tanglegram
```bash
viral-phylo tree \
  --alignment test_alignment/foldmason.fasta_3di.fa \
  --alignment-aa test_alignment/foldmason.fasta_aa.fa \
  --tree-type both \
  --threads 2 \
  --output-dir test_phylogeny \
  --prefix test_dual
```

**Validation Checks:**
- [ ] Both `test_dual_both.treefile` (3Di) and `test_dual_aa.treefile` (AA) are generated.

---

### Test 8: Approximate-ML Tree with the 3Di Matrix (`--method fasttree`)
```bash
viral-phylo tree \
  --alignment test_alignment \
  --method fasttree \
  --matrix alphafold \
  --tree-type both \
  --output-dir test_phylogeny \
  --prefix ft_result
```

**Validation Checks:**
- [ ] `test_phylogeny/ft_result_3di_fasttree.treefile` and `ft_result_aa_fasttree.treefile` are generated.
- [ ] `ft_result_3di_fasttree.treefile.info.json` has `"model": "Q.3Di.AF (converted for -trans) + Gamma"` and a command containing `-trans` (and `-double-precision` for VeryFastTree).

---

### Test 9: End-to-End Pipeline from a Local Folder
```bash
viral-phylo pipeline \
  --input-folder tests/fixtures/structures \
  --tree-type both \
  --matrix alphafold \
  --threads 2 \
  --output-dir results/validation_local_workflow
```

**Validation Checks:**
- [ ] Exit code is `0`, with no network access needed.
- [ ] `results/validation_local_workflow/run_manifest.txt` lists every stage as `ok` or `skipped` (the embed stage is `skipped` because `--embed` was not requested), and the status is `completed_with_skips`.
- [ ] The manifest records the version of every external tool, and the input's file count and content digest.
- [ ] `results/validation_local_workflow/interactive_tree.html` is written.

---

### Test 10: Similarity Clustering (`--multi-alignment`)
```bash
viral-phylo pipeline \
  --input-folder tests/fixtures/structures \
  --multi-alignment \
  --cluster-min-size 2 \
  --tree-type 3di \
  --matrix alphafold \
  --threads 2 \
  --output-dir results/validation_cluster_workflow
```

**Validation Checks:**
- [ ] The console reports the whole-set 3Di gate verdict and `Clustering at evalue=1.0 coverage=0.5 (foldseek structural)`.
- [ ] `results/validation_cluster_workflow/multi_alignment_summary.json` lists the clusters, the QC counts and the operating point.
- [ ] Each `results/validation_cluster_workflow/alignment/cluster_*/` has `cluster_info.json`, and every cluster of 2+ taxa (i.e. every aligned cluster) also has `alignment_info.json`.
- [ ] Clusters with fewer than 4 taxa are recorded in the manifest as skipped, with the reason.

---

### Test 11: Preflight Stops a Misconfigured Run
```bash
viral-phylo pipeline \
  --input-folder tests/fixtures/structures \
  --metadata no_such_file.csv \
  --output-dir results/validation_preflight
```

**Validation Checks:**
- [ ] The run stops within seconds with `[Preflight] BLOCKER: metadata 'no_such_file.csv' does not exist`.
- [ ] `results/validation_preflight/run_manifest.json` has status `failed` and a single `preflight` stage; no alignment was started.

---

### Test 12: Your Own Large Dataset (optional)
The repository does not ship a large cohort. To validate at scale, run the large-dataset command from the README on your own structures and check `run_manifest.txt`: every stage `ok`, the whole-set gate verdict, the cluster count, how many clusters were split or kept above threshold, and the starting-tree outcomes for clusters of 80+ taxa.

---

## 4. End-to-End One-Line Integration Test

```bash
viral-phylo align -i tests/fixtures/structures -o e2e_alignment && \
viral-phylo tree -a e2e_alignment/foldmason.fasta_3di.fa --tree-type both -t 2 -o e2e_phylogeny -p e2e_viral && \
viral-phylo tree -a e2e_alignment --method fasttree --matrix alphafold --tree-type both -o e2e_phylogeny -p e2e_ft && \
echo "=== ALL INTEGRATION TESTS PASSED ==="
```

### Cleanup Test Artifacts
```bash
rm -rf test_structures test_alignment test_phylogeny e2e_alignment e2e_phylogeny \
       results/validation_local_workflow results/validation_cluster_workflow results/validation_preflight
```
