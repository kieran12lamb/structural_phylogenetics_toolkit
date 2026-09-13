---
name: viral-structural-phylogenetics
description: >
  Fetch viral protein structures from Viro3D, perform structural multiple sequence alignment (MSA) using FoldMason to generate 3Di structural alphabet alignments, and construct rigorous structural phylogenetic trees using IQ-TREE with AlphaFold and ESMFold 3Di substitution matrices (Q.3Di.AF, Q.3Di.LLM) or FoldMason guide trees. Supports automated rate heterogeneity and model selection with ModelFinder.
---

# Viral Structural Phylogenetics Skill

This skill provides an end-to-end workflow for **structural phylogenetics** of viral proteins. By operating on 3D protein structure coordinates rather than divergent primary amino acid sequences, it enables resolution of deep evolutionary relationships that are otherwise obscured by sequence divergence and substitution saturation.

The pipeline bridges three core technologies:
1. **Viro3D**: Retrieval of AI-predicted viral protein 3D structures (ColabFold/AlphaFold2 and ESMFold).
2. **FoldMason**: Multiple structural alignment (MSTA) encoding protein coordinates into the 20-state **3Di** (3D interaction) structural alphabet.
3. **IQ-TREE with 3Di Empirical Matrices**: Maximum likelihood phylogenetic inference using empirical 3Di substitution matrices (`Q.3Di.AF` and `Q.3Di.LLM` from Edmond doi:10.17617/3.1MJJBH) with automated rate heterogeneity selection (Gamma `+G4`, Invariable sites `+I`, FreeRate `+R`) and ultrafast bootstrap support.

---

## 1. Prerequisites & Environment

The tools require Python 3.10+ along with `foldmason` and `iqtree` installed via conda/bioconda:

```bash
# Recommended environment setup:
conda create -y -n spt python=3.11
conda install -y -n spt -c bioconda -c conda-forge foldmason iqtree requests
```

### 3Di Substitution Matrices
Empirical 3Di substitution models inferred by Georg Hochberg et al. (*"A general substitution matrix for structural phylogenetics"*, Edmond Dataverse [doi:10.17617/3.1MJJBH](https://doi.org/10.17617/3.1MJJBH)) are located in [`matrices/`](matrices/):
- **`matrices/Q.3Di.AF`** (File ID: `311466`): Inferred from AlphaFold protein structures.
- **`matrices/Q.3Di.LLM`** (File ID: `311467`): Inferred from ESMFold/ProstT5 protein language model translations.

*(If either matrix is missing, the script automatically downloads it from the Edmond repository).*

---

## 2. Fast Track: End-to-End Pipeline

The pipeline can run either from online **Viro3D queries** or from a **local directory of structures** (`.pdb` / `.cif`).

### Option A: Running from Local Folder of Structures (New Feature)
Bypasses the online download step and directly processes existing structures:
```bash
python3 scripts/viral_phylogenetics.py pipeline \
  --input-folder /path/to/my_structures \
  --tree-type both \
  --matrix alphafold \
  --threads 10 \
  --output-dir local_workflow
```

### Option B: Running from Online Viro3D Query
```bash
python3 scripts/viral_phylogenetics.py pipeline \
  --qualifier glycoprotein \
  --count 500 \
  --tree-type both \
  --rate-heterogeneity auto \
  --threads 10 \
  --output-dir glycoprotein_500_workflow
```

*(Note: For large cohorts $\ge 80$ structures such as the 500-structure run, the pipeline automatically activates high-performance heuristic modes: FoldMason lightweight reporting, empirical `Q.3Di.AF+G4` and `LG+G4` models with IQ-TREE `--fast`, executing in under 2 minutes).*


---

## 3. Step-by-Step CLI Subcommands

### Step 1: Download Viral Structures (`fetch`)
Downloads AlphaFold/ColabFold PDB coordinate files from the Viro3D database by protein name, product, or virus name:

```bash
python3 scripts/viral_phylogenetics.py fetch \
  --qualifier glycoprotein \
  --max-sequences 6 \
  --output-dir viro_3d_structures
```

### Step 2: Structural Multiple Alignment (`align`)
Aligns the downloaded structures in 3D space using FoldMason, translating backbone conformations into the 20-letter 3Di alphabet:

```bash
python3 scripts/viral_phylogenetics.py align \
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
python3 scripts/viral_phylogenetics.py tree \
  --alignment foldmason_alignments/foldmason.fasta_3di.fa \
  --method iqtree \
  --matrix both \
  --rate-heterogeneity auto \
  --criterion BIC \
  --bootstrap 1000 \
  --alrt 1000 \
  --threads AUTO \
  --output-dir phylogeny_results \
  --prefix viral_tree
```

#### Mode B: Fast FoldMason Guide Tree Extraction
Extracts and formats FoldMason's built-in progressive guide tree (fast, uncalibrated cladogram):

```bash
python3 scripts/viral_phylogenetics.py tree \
  --alignment foldmason_alignments \
  --method foldmason \
  --output-dir phylogeny_results \
  --prefix fm_guide
```

#### Mode C: Dual Phylogeny & Cophylogenetic Tanglegram (`--tree-type both`)
Infers both the 3Di structural tree (IQ-TREE + `Q.3Di.AF`) and the amino acid sequence tree (IQ-TREE + ModelFinder), enabling direct comparison of tertiary structure conservation versus primary sequence divergence:

```bash
python3 scripts/viral_phylogenetics.py tree \
  --alignment foldmason_alignments/foldmason.fasta_3di.fa \
  --alignment-aa foldmason_alignments/foldmason.fasta_aa.fa \
  --tree-type both \
  --output-dir phylogeny_results \
  --prefix viral_tree
```

---

## 4. Interactive Visualization & Tanglegram (`interactive_tree.html`)

The workflow generates a self-contained, high-performance HTML visualization document:
- **5 Visualization Layouts**:
  - **Phylogram**: Cartesian rectangular tree with calibrated branch lengths.
  - **Cladogram**: Uniform-depth branching diagram.
  - **Radial**: Polar projection tree for large viral clades.
  - **Unrooted**: Equal-angle star tree showing topological radiation without ancestral root assumptions.
  - **Tanglegram**: Dual facing trees (3Di Structural Tree on the left, Primary AA Tree on the right) with curved cophylogenetic connector lines highlighting topological congruence and discordance.
- **On-Hover 3D Structure Viewer**:
  - Hovering over any tip node displays the rotating 3D C$\alpha$ backbone trace rendered live in an HTML5 canvas.
  - Color-coded by AlphaFold pLDDT confidence: Royal Blue ($\ge 90$), Cyan ($70-89$), Yellow ($50-69$), Orange ($< 50$).
  - Full mouse drag-to-rotate interaction and zoom.
- **Dynamic Node Color Legends**:
  - Automatically updates when toggling node color modes: **UFboot Support**, **AlphaFold pLDDT**, or **Viral Taxonomic Family**.
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

| File Extension | Description |
| :--- | :--- |
| `.treefile` | Maximum likelihood tree in Newick format with calibrated substitution branch lengths and branch support values (`SH-aLRT / UFboot`). |
| `.contree` | Ultrafast bootstrap 50% majority-rule consensus tree. |
| `.iqtree` | Comprehensive report: model comparison table, log-likelihood, AIC/BIC scores, estimated rate parameters, and ASCII tree rendering. |
| `.splits.nex` | NEXUS file containing split support values and frequencies. |
| `.mldist` | Pairwise maximum likelihood distance matrix across all viral proteins. |
| `interactive_tree.html` | Interactive web application with 5 layouts, tanglegram, dynamic legends, and on-hover 3D structure viewer. |
| `.log` | Full execution and optimization log. |

