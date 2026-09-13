# Viral Structural Phylogenetics Skill - Testing & Validation Guide

This document defines the test suite and validation benchmarks for the **Viral Structural Phylogenetics** skill (`viro3d_structures.py`, `foldmason_align.py`, and `structural_phylogeny.py`).

---

## 1. Environment & Prerequisite Check

Before executing tests, verify that the environment has `foldmason`, `iqtree`, and `requests` available:

```bash
# Verify Python and binaries:
conda run -n spt python3 --version
conda run -n spt foldmason version
conda run -n spt iqtree --version
```

### Expected Output:
- Python $\ge 3.10$
- FoldMason installed and executable
- IQ-TREE version $\ge 2.3.0$ (or 3.x)

---

## 2. Test Suite Overview

| Test Case | Description | Input | Expected Output |
| :--- | :--- | :--- | :--- |
| **Test 1** | Viro3D Structure Download | Query: `glycoprotein`, 4 structures | 4 valid `.pdb` files in output directory |
| **Test 2** | FoldMason 3Di Structural MSA | 4 PDB files from Test 1 | `foldmason.fasta_3di.fa`, `foldmason.fasta.nw` |
| **Test 3** | FoldMason Guide Tree | `foldmason.fasta_3di.fa` | Formatted `.nwk` cladogram |
| **Test 4** | IQ-TREE with AlphaFold 3Di Matrix | `foldmason.fasta_3di.fa` + `Q.3Di.AF` | `.treefile` with ML branch lengths |
| **Test 5** | IQ-TREE with ESMFold 3Di Matrix | `foldmason.fasta_3di.fa` + `Q.3Di.LLM` | `.treefile` with ESMFold ML tree |
| **Test 6** | Automated ModelFinder Selection | `foldmason.fasta_3di.fa` + both matrices | Optimal model selection (BIC) + bootstrap tree |
| **Test 7** | Full End-to-End Pipeline | Protein search term -> Full ML tree | Clean multi-stage execution |

---

## 3. Detailed Test Cases & Execution

### Test 1: Viro3D Structure Download
Tests API search and structure file downloading from the Viro3D database.

```bash
conda run -n spt python3 viro3d_structures.py \
  --qualifier glycoprotein \
  --max-sequences 4 \
  --output-dir test_structures
```

**Validation Checks:**
- [x] Exit code is `0`.
- [x] Exactly 4 `.pdb` files exist in `test_structures/`.
- [x] Each file contains valid PDB `ATOM` coordinate records and file size $>10\text{ KB}$.

---

### Test 2: FoldMason Structural Multiple Alignment
Tests structural alignment of 3D protein backbones and translation to the 3Di alphabet.

```bash
conda run -n spt python3 foldmason_align.py \
  --folder test_structures \
  --output-dir test_alignment
```

**Validation Checks:**
- [x] Exit code is `0`.
- [x] `test_alignment/foldmason.fasta_3di.fa` exists and contains 4 aligned sequences in the 3Di alphabet (`A, C, D, E, F, G, H, I, K, L, M, N, P, Q, R, S, T, V, W, Y` and `-`).
- [x] `test_alignment/foldmason.fasta.nw` contains a valid Newick guide tree.
- [x] `test_alignment/foldmason.fasta.html` exists for visual verification.

---

### Test 3: FoldMason Guide Tree Extraction
Tests guide tree extraction via `structural_phylogeny.py`.

```bash
conda run -n spt python3 structural_phylogeny.py \
  --alignment test_alignment \
  --method foldmason \
  --output-dir test_phylogeny \
  --prefix guide_tree
```

**Validation Checks:**
- [x] Exit code is `0`.
- [x] `test_phylogeny/guide_tree_foldmason.nwk` is generated.
- [x] The file contains a complete Newick-formatted tree ending with `;`.

---

### Test 4: IQ-TREE with AlphaFold Matrix (`Q.3Di.AF`)
Tests maximum likelihood tree inference using the AlphaFold 3Di substitution matrix with 4 Gamma rate categories (`+G4`) and 1000 ultrafast bootstrap replicates.

```bash
conda run -n spt python3 structural_phylogeny.py \
  --alignment test_alignment \
  --method iqtree \
  --matrix alphafold \
  --rate-heterogeneity +G4 \
  --bootstrap 1000 \
  --alrt 1000 \
  --threads 2 \
  --output-dir test_phylogeny \
  --prefix af_result
```

**Validation Checks:**
- [x] Matrix `matrices/Q.3Di.AF` is verified/loaded.
- [x] `test_phylogeny/af_result_alphafold.treefile` exists and contains substitution branch lengths.
- [x] `test_phylogeny/af_result_alphafold.contree` contains bootstrap support values.
- [x] `test_phylogeny/af_result_alphafold.iqtree` includes log-likelihood and estimated Gamma shape parameter $\alpha$.

---

### Test 5: IQ-TREE with ESMFold Matrix (`Q.3Di.LLM`)
Tests maximum likelihood tree inference using the ESMFold/ProstT5 3Di substitution matrix.

```bash
conda run -n spt python3 structural_phylogeny.py \
  --alignment test_alignment \
  --method iqtree \
  --matrix esmfold \
  --rate-heterogeneity +G4 \
  --bootstrap 1000 \
  --alrt 1000 \
  --threads 2 \
  --output-dir test_phylogeny \
  --prefix esm_result
```

**Validation Checks:**
- [x] Matrix `matrices/Q.3Di.LLM` is verified/loaded.
- [x] Output tree generated under ESMFold substitution rates.
- [x] Output report logs successful convergence of tree search.

---

### Test 6: Automated Model & Rate Heterogeneity Selection
Tests ModelFinder Plus (`-m MFP`) to simultaneously compare both matrices and all rate heterogeneity configurations (`Uniform`, `+I`, `+G4`, `+I+G4`, `+R4`).

```bash
conda run -n spt python3 structural_phylogeny.py \
  --alignment test_alignment \
  --method iqtree \
  --matrix both \
  --rate-heterogeneity auto \
  --criterion BIC \
  --bootstrap 1000 \
  --alrt 1000 \
  --threads 2 \
  --output-dir test_phylogeny \
  --prefix auto_result
```

**Validation Checks:**
- [x] ModelFinder evaluates both `Q.3Di.AF` and `Q.3Di.LLM` across rate classes.
- [x] Winning model is printed to console (e.g. `Best-fit model according to BIC: ...`).
- [x] `test_phylogeny/auto_result_both.iqtree` contains the complete model ranking table sorted by BIC.
- [x] Final tree is inferred using the optimal model.

---

### Test 7: Dual Phylogeny & Cophylogenetic Tanglegram
Tests concurrent inference of structural (3Di) and amino acid sequence phylogenies with IQ-TREE and interactive Tanglegram generation.

```bash
conda run -n spt python3 scripts/viral_phylogenetics.py tree \
  --alignment test_alignment/foldmason.fasta_3di.fa \
  --alignment-aa test_alignment/foldmason.fasta_aa.fa \
  --tree-type both \
  --output-dir test_phylogeny \
  --prefix test_dual
```

**Validation Checks:**
- [x] Both `test_dual_both.treefile` (3Di structural tree) and `test_dual_aa.treefile` (AA sequence tree) are generated.
- [x] IQ-TREE reports model selection for both alignments (`Q.3Di.AF`/`Q.3Di.LLM` for 3Di, and protein models e.g. `PMB`/`LG` for AA).
- [x] Interactive web visualization renders the side-by-side Tanglegram with cophylogenetic connector lines and dynamic node color legend.

---

### Test 8: Local Structure Folder Processing (`--input-folder`)
Tests running the full pipeline from a pre-existing local directory of `.pdb` or `.cif` files without internet access.

```bash
conda run -n spt python3 scripts/viral_phylogenetics.py pipeline \
  --input-folder viro_glycoproteins \
  --tree-type both \
  --threads 2 \
  --output-dir test_local_workflow
```

**Validation Checks:**
- [x] Automatically detects local `.pdb` files in `viro_glycoproteins/`.
- [x] Executes FoldMason alignment without attempting network downloads.
- [x] Successfully builds both 3Di structural tree and amino acid sequence tree.
- [x] Exit code is `0`.

---

### Test 9: Massive Dataset Scaling (500 Structures)
Tests the pipeline on a 500-structure viral glycoprotein cohort with automatic heuristic acceleration.

```bash
conda run -n spt python3 scripts/viral_phylogenetics.py pipeline \
  --input-folder viro_500_glycoproteins \
  --tree-type both \
  --threads 10 \
  --output-dir test_500_workflow
```

**Validation Checks:**
- [x] FoldMason aligns 500 protein structures in $< 35$ seconds.
- [x] Fast heuristic automatically activates (`Q.3Di.AF+G4` and `LG+G4` with `--fast`).
- [x] Generates complete 500-taxa 3Di and AA trees in $< 60$ seconds.
- [x] Output trees contain 500 leaves and 997 branches.

---

## 4. End-to-End One-Line Integration Test

Execute all stages in sequence to confirm pipeline integrity:

```bash
conda run -n spt python3 viro3d_structures.py -q glycoprotein -m 4 -o e2e_structures && \
conda run -n spt python3 foldmason_align.py -i e2e_structures -o e2e_alignment && \
conda run -n spt python3 scripts/viral_phylogenetics.py tree -a e2e_alignment/foldmason.fasta_3di.fa --tree-type both -t 2 -o e2e_phylogeny -p e2e_viral && \
echo "=== ALL INTEGRATION TESTS PASSED ==="
```

### Cleanup Test Artifacts
```bash
rm -rf test_structures test_alignment test_phylogeny e2e_structures e2e_alignment e2e_phylogeny test_local_workflow
```


