---
name: viral-structural-phylogenetics
description: >
  Fetch viral protein structures from Viro3D or AlphaFold DB, align them with FoldMason (or MAFFT) into the 3Di structural alphabet, cluster large or divergent sets by structural similarity (Foldseek) or sequence identity (MMseqs2), and build structural phylogenies with IQ-TREE using AlphaFold and ESMFold 3Di substitution matrices (Q.3Di.AF, Q.3Di.LLM), plus fast whole-set trees with VeryFastTree under the same matrices. Supports ModelFinder rate heterogeneity selection, ESM-2 embedding trees, an interactive dashboard, and a run manifest recording every stage.
---

# Viral Structural Phylogenetics Skill

This skill provides an end-to-end workflow for **structural phylogenetics** of viral proteins. By operating on 3D protein structure coordinates rather than divergent primary amino acid sequences, it enables resolution of deep evolutionary relationships that are otherwise obscured by sequence divergence and substitution saturation.

The pipeline bridges these core technologies:
1. **Viro3D / AlphaFold DB**: Retrieval of AI-predicted viral protein 3D structures (ColabFold/AlphaFold2 and ESMFold).
2. **FoldMason**: Multiple structural alignment (MSTA) encoding protein coordinates into the 20-state **3Di** (3D interaction) structural alphabet.
3. **Foldseek / MMseqs2**: Structural or sequence clustering, so each alignment and tree is built from proteins that genuinely align.
4. **IQ-TREE with 3Di Empirical Matrices**: Maximum likelihood phylogenetic inference using empirical 3Di substitution matrices (`Q.3Di.AF` and `Q.3Di.LLM` from Edmond doi:10.17617/3.1MJJBH) with automated rate heterogeneity selection (Gamma `+G4`, Invariable sites `+I`, FreeRate `+R`) and ultrafast bootstrap support.
5. **VeryFastTree**: Approximate-ML whole-set trees in about a minute for 1,000+ taxa, using the same 3Di matrix, and validated starting trees for IQ-TREE on large clusters.

---

## 1. Prerequisites & Environment

Use the provided environment, which includes every required tool (`foldmason`, `mafft`, `iqtree`, `foldseek`, `mmseqs2`, `veryfasttree`/`fasttree`):

```bash
./install.sh            # creates/updates the `spt` env, installs the package, runs the tests
conda activate spt
```

Or `conda env create -f environment.yml && conda activate spt && pip install -e . --no-deps`.

**Always run from the repository root** — `matrices/` is resolved relative to the working directory. Before doing any work, the pipeline's preflight check verifies every tool, input and interpreter the chosen options need, and stops immediately with a list if something is missing.

### 3Di Substitution Matrices
Empirical 3Di substitution models inferred by Georg Hochberg et al. (*"A general substitution matrix for structural phylogenetics"*, Edmond Dataverse [doi:10.17617/3.1MJJBH](https://doi.org/10.17617/3.1MJJBH)) are located in [`matrices/`](matrices/):
- **`matrices/Q.3Di.AF`** (File ID: `311466`): Inferred from AlphaFold protein structures.
- **`matrices/Q.3Di.LLM`** (File ID: `311467`): Inferred from ESMFold/ProstT5 protein language model translations.

*(If either matrix is missing, the script automatically downloads it from the Edmond repository).*

---

## 2. Fast Track: End-to-End Pipeline

The pipeline can run either from online **Viro3D queries** or from a **local directory of structures** (`.pdb` / `.cif`).

### Option A: Running from Local Folder of Structures
Bypasses the online download step and directly processes existing structures:
```bash
viral-phylo pipeline \
  --input-folder /path/to/my_structures \
  --tree-type both \
  --matrix alphafold \
  --output-dir results/local_workflow
```

### Option B: Running from Online Viro3D Query
```bash
viral-phylo pipeline \
  --qualifier glycoprotein \
  --count 500 \
  --tree-type both \
  --rate-heterogeneity auto \
  --output-dir results/glycoprotein_500_workflow
```

*(For alignments of $\ge 80$ sequences, IQ-TREE skips model selection and uses the first chosen 3Di matrix with `+G4` and `LG+G4` for amino acids. IQ-TREE `--fast` is only switched on automatically when bootstrapping is disabled.)*

### Option C: Large or Divergent Structure Sets (Recommended at Scale)
When the input does not share one fold, aligning everything together yields a mostly-gap alignment. Cluster first:
```bash
viral-phylo pipeline \
  --input-folder /path/to/structures \
  --metadata /path/to/metadata.csv \
  --output-dir results/large_run \
  --multi-alignment --cluster-recurse \
  --tree-type both --matrix alphafold \
  --embed --embed-clustering upgma
```
This removes short fragments, clusters structures with Foldseek (e-value 1.0, coverage 0.5), splits clusters whose own 3Di alignment is still too gappy (only when the split keeps ≥70% of taxa in clusters of 4+), builds an IQ-TREE tree for every cluster of 4+ taxa (clusters of 80+ start from a validated VeryFastTree tree), and a VeryFastTree whole-set tree.

Defaults worth knowing:
- `--threads 8` caps the whole run: every tool (FoldMason, MAFFT, Foldseek/MMseqs2, IQ-TREE, VeryFastTree, PyTorch, BLAS) is limited to it. Do not use `AUTO` on many-core machines: it removes the cap, and IQ-TREE then benchmarks every thread count, which can take longer than the inference.
- Cluster trees run in parallel within that cap, `--tree-threads` (8) each: `--threads 32` builds 4 at a time. Raise `--threads` rather than giving one IQ-TREE run more threads; beyond ~8 IQ-TREE gets slower on short cluster alignments.
- `--master-method fasttree` for the whole-set tree; `--method iqtree` for cluster trees.
- `--embed-clustering nj` produces no silhouette profile; use `upgma` (default) or another linkage method.

---

## 3. Step-by-Step CLI Subcommands

### Step 1: Download Viral Structures (`fetch`)
Downloads AlphaFold/ColabFold PDB coordinate files from the Viro3D database by protein name, product, or virus name:

```bash
viral-phylo fetch \
  --qualifier glycoprotein \
  --max-sequences 6 \
  --output-dir viro_3d_structures
```

### Step 2: Structural Multiple Alignment (`align`)
Aligns the downloaded structures in 3D space using FoldMason, translating backbone conformations into the 20-letter 3Di alphabet:

```bash
viral-phylo align \
  --folder viro_3d_structures \
  --output-dir foldmason_alignments
```

**Key Outputs:**
- `foldmason_alignments/foldmason.fasta_3di.fa`: MSA in the **3Di alphabet** (used for structural phylogenetics).
- `foldmason_alignments/foldmason.fasta_aa.fa`: MSA in standard amino acid alphabet.
- `foldmason_alignments/foldmason.fasta.nw`: Built-in FoldMason progressive guide tree.
- `foldmason_alignments/foldmason.fasta.html`: Interactive structural alignment viewer.

### Step 3: Structural Phylogeny Inference (`tree`)

#### Mode A: Automated Model & Rate Heterogeneity Selection (Recommended)
Uses IQ-TREE's **ModelFinder Plus** (`-m MFP`) to simultaneously test both the AlphaFold and ESMFold matrices across all rate heterogeneity classes (`Uniform`, `+I`, `+G4`, `+I+G4`, `+R`), then computes maximum likelihood branch lengths and bootstrap supports:

```bash
viral-phylo tree \
  --alignment foldmason_alignments/foldmason.fasta_3di.fa \
  --method iqtree \
  --matrix both \
  --rate-heterogeneity auto \
  --criterion BIC \
  --bootstrap 1000 \
  --alrt 1000 \
  --threads 8 \
  --output-dir phylogeny_results \
  --prefix viral_tree
```

#### Mode B: Fast Approximate-ML Tree (`--method fasttree`)
VeryFastTree (or FastTree) with the chosen 3Di matrix applied via `-trans`; about a minute for 1,000+ taxa. Branch support is a single SH-like value, not UFboot:

```bash
viral-phylo tree \
  --alignment foldmason_alignments/foldmason.fasta_3di.fa \
  --method fasttree \
  --matrix alphafold \
  --tree-type 3di \
  --output-dir phylogeny_results \
  --prefix viral_ft
```

#### Mode C: Fast FoldMason Guide Tree Extraction
Extracts and formats FoldMason's built-in progressive guide tree (fast, uncalibrated cladogram):

```bash
viral-phylo tree \
  --alignment foldmason_alignments \
  --method foldmason \
  --output-dir phylogeny_results \
  --prefix fm_guide
```

#### Mode D: Dual Phylogeny & Cophylogenetic Tanglegram (`--tree-type both`)
Infers both the 3Di structural tree (IQ-TREE + `Q.3Di.AF`) and the amino acid sequence tree (IQ-TREE + ModelFinder), enabling direct comparison of tertiary structure conservation versus primary sequence divergence:

```bash
viral-phylo tree \
  --alignment foldmason_alignments/foldmason.fasta_3di.fa \
  --alignment-aa foldmason_alignments/foldmason.fasta_aa.fa \
  --tree-type both \
  --output-dir phylogeny_results \
  --prefix viral_tree
```

---

## 4. Interactive Visualization & Tanglegram (`interactive_tree.html`)

Every pipeline run writes `<output-dir>/interactive_tree.html` (generated output; not part of the repository). It lists every workflow under `results/`; `build-alignments && build-tree-view` rebuilds them all. For a pipeline run the tree shown is the whole-set tree; per-cluster IQ-TREE trees are in `phylogeny/cluster_*/` and can be opened with `viral-phylo report --tree <file>`.

- **5 Visualization Layouts**:
  - **Phylogram**: Cartesian rectangular tree with calibrated branch lengths.
  - **Cladogram**: Uniform-depth branching diagram.
  - **Radial**: Polar projection tree for large viral clades.
  - **Unrooted**: Equal-angle star tree showing topological radiation without ancestral root assumptions.
  - **Tanglegram**: Dual facing trees (3Di Structural Tree on the left, Primary AA Tree on the right) with curved cophylogenetic connector lines highlighting topological congruence and discordance.
- **On-Hover 3D Structure Viewer**:
  - Hovering over any tip node displays the rotating 3D C$\alpha$ backbone trace rendered live in an HTML5 canvas (read from PDB or mmCIF input).
  - Color-coded by AlphaFold pLDDT confidence: Royal Blue ($\ge 90$), Cyan ($70-89$), Yellow ($50-69$), Orange ($< 50$).
  - Full mouse drag-to-rotate interaction and zoom.
- **Dynamic Node Color Legends**:
  - Automatically updates when toggling node color modes: **UFboot Support**, **AlphaFold pLDDT**, or **Viral Taxonomic Family**. Trees with a single support value per branch (VeryFastTree/FastTree) are labelled *single-value support*, not UFboot.
- **Light & Dark Theme Modes**:
  - One-click toggle between high-contrast dark theme (midnight slate) and crisp publication-ready light theme.
  - Theme preference is remembered across sessions via `localStorage`.
- **Interactive Clade Collapsing for Large Trees**:
  - Click any internal node dot or branch to collapse/expand its subtree.
  - Collapsed clades render as stylized triangular wedges colored by dominant viral family, showing leaf counts and average pLDDT.
  - Bulk actions in the sidebar: **"By Family"** (collapses monophyletic family clades into macro-wedges), **"Deep Clades"**, and **"Expand All"**.
  - Real-time visible leaf counter tracks active taxa.

---

## 5. Key Output Files

| File | Description |
| :--- | :--- |
| `run_manifest.json` / `.txt` | The run record: command, resolved settings, tool versions, input fingerprint, and every stage's status (`ok` / `skipped` / `failed`) with reasons. Check this first. Final status: `completed`, `completed_with_skips`, `failed`, or `interrupted`. |
| `multi_alignment_summary.json` | Clusters, QC exclusions, operating point, and recursion record (`--multi-alignment`). |
| `alignment/cluster_*/cluster_info.json` | Per-cluster size, 3Di and AA alignment statistics, recursion depth and parent. |
| `alignment_info.json` | Which aligner really produced an alignment (`foldmason.fasta_*.fa` names are kept even for MAFFT output). |
| `*_3di_fasttree.treefile` (+ `.info.json`) | Whole-set VeryFastTree trees and their program, version, model and command. |
| `*.start_tree.json` | Whether an IQ-TREE run used a starting tree, fell back, or rejected it. |
| `.treefile` | Maximum likelihood tree in Newick format with calibrated substitution branch lengths and branch support values (`SH-aLRT / UFboot`). |
| `.contree` | Ultrafast bootstrap 50% majority-rule consensus tree. |
| `.iqtree` | Comprehensive report: model comparison table, log-likelihood, AIC/BIC scores, estimated rate parameters, and ASCII tree rendering. |
| `.splits.nex` | NEXUS file containing split support values and frequencies. |
| `.mldist` | Pairwise maximum likelihood distance matrix across all viral proteins. |
| `interactive_tree.html` | Interactive web application with 5 layouts, tanglegram, dynamic legends, and on-hover 3D structure viewer. |
| `.log` | Full execution and optimization log. |
