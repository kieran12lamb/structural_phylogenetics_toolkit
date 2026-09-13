# Viral Structural Phylogenetics Suite — Comprehensive Usage Guide

An end-to-end framework for viral protein structural phylogenetics. The suite combines:
1. **Viro3D Automated Retrieval**: Programmatic retrieval of AI-predicted viral 3D structures (ColabFold/AlphaFold2 and ESMFold) and rich clinical/taxonomic metadata.
2. **FoldMason Multiple Structural Alignment (MSTA)**: 3D coordinate superposition and translation into the 20-state **3Di** (3D interaction) structural alphabet.
3. **IQ-TREE Empirical Structural Phylogenetics**: Maximum likelihood phylogenetic inference using empirical 3Di substitution matrices (`Q.3Di.AF`, `Q.3Di.LLM`, `Q.3Di.ESM`) with automated rate heterogeneity selection and ultrafast bootstrap.
4. **Interactive SVG & WebGL Visualization Suite**: High-performance interactive browser interface with dual tanglegrams, dynamic metadata coloring/filtering, topological congruence metrics, synchronized Multiple Sequence Alignment (MSA) viewer with 2D radar overview, and live 3D protein backbone visualization.
5. **Interactive Web Pipeline Studio**: Direct in-browser Viro3D querying, CLI command generation, and background pipeline execution bridge.

---

## Table of Contents
1. [Prerequisites & Installation](#1-prerequisites--installation)
2. [Quick-Start Workflows](#2-quick-start-workflows)
3. [CLI Reference (`scripts/viral_phylogenetics.py`)](#3-cli-reference)
   - [`pipeline` (End-to-End)](#subcommand-pipeline)
   - [`fetch` (Viro3D Download)](#subcommand-fetch)
   - [`align` (FoldMason 3Di Alignment)](#subcommand-align)
   - [`tree` (IQ-TREE 3Di Phylogeny)](#subcommand-tree)
4. [Custom Metadata Integration](#4-custom-metadata-integration)
5. [Interactive Web Visualization Suite (`interactive_tree.html`)](#5-interactive-web-visualization-suite)
   - [Launching the Viewer](#launching-the-viewer)
   - [Tree Layouts & Exploration](#tree-layouts--exploration)
   - [Dynamic Metadata & Clade Grouping](#dynamic-metadata--clade-grouping)
   - [Metadata Filtering (Prune vs Highlight)](#metadata-filtering)
   - [Structural vs Sequence Congruence Analysis](#structural-vs-sequence-congruence-analysis)
   - [Multiple Sequence Alignment (MSA) Viewer & Minimap](#msa-viewer--minimap)
   - [3D C-alpha Backbone Structure Viewer](#3d-structure-viewer)
   - [Zoom-Adaptive Label Scaling](#zoom-adaptive-label-scaling)
6. [Interactive Pipeline Studio Tab](#6-interactive-pipeline-studio-tab)
7. [Troubleshooting & Best Practices](#7-troubleshooting--best-practices)

---

## 1. Prerequisites & Installation

### Environment Setup (Conda / Mamba)
The tool requires Python 3.10+ and standard bioinformatics tools installed via Conda:

```bash
# 1. Create and activate environment
conda create -y -n skills_hackathon python=3.11
conda activate skills_hackathon

# 2. Install core bioinformatics binaries and python libraries
conda install -y -c bioconda -c conda-forge foldmason iqtree
pip install openpyxl pandas biopython
```

### Verification
Ensure `foldmason` and `iqtree` are available on your system `PATH`:
```bash
foldmason version
iqtree --version
```

### Empirical 3Di Substitution Matrices
The pipeline utilizes empirical 3Di substitution models inferred from structural protein databases:
- **`matrices/Q.3Di.AF`**: Inferred from AlphaFold structures (Edmond ID: `311466`).
- **`matrices/Q.3Di.LLM`**: Inferred from ESMFold/ProstT5 translations (Edmond ID: `311467`).

*(If missing, `scripts/viral_phylogenetics.py` automatically downloads them from the Edmond Dataverse repository upon first run).*

---

## 2. Quick-Start Workflows

### Scenario A: Online Viro3D Query to Phylogenetic Trees
Download 50 viral glycoproteins directly from Viro3D, align with FoldMason, and build both 3Di structural and amino acid sequence trees:
```bash
python scripts/viral_phylogenetics.py pipeline \
  --qualifier glycoprotein \
  --count 50 \
  --tree-type both \
  --matrix alphafold \
  --output-dir glycoprotein_50_workflow
```

### Scenario B: Process a Local Directory of PDB Structures
If you already have PDB coordinate files on disk (e.g. from local AlphaFold/ESMFold predictions):
```bash
python scripts/viral_phylogenetics.py pipeline \
  --input-folder path/to/my_structures \
  --tree-type both \
  --matrix both \
  --output-dir my_local_workflow
```

### Scenario C: Launch the Interactive Web Suite
Start the local server bridge and open the visualization suite in your default browser:
```bash
python scripts/serve_interactive.py
```
*(Or simply open `interactive_tree.html` directly in any web browser).*

---

## 3. CLI Reference

The CLI entrypoint is [`scripts/viral_phylogenetics.py`](scripts/viral_phylogenetics.py). It provides four primary subcommands:

```
usage: viral_phylogenetics.py [-h] {fetch,align,tree,pipeline} ...
```

---

### Subcommand: `pipeline`
Runs the complete workflow: fetches/loads structures $\rightarrow$ performs FoldMason MSTA $\rightarrow$ constructs maximum-likelihood trees.

```bash
python scripts/viral_phylogenetics.py pipeline [options]
```

| Argument | Type | Default | Description |
| :--- | :--- | :--- | :--- |
| `-q`, `--qualifier` | str | `glycoprotein` | Protein search term if fetching from Viro3D (e.g. `glycoprotein`, `rdrp`, `spike`, `capsid`). |
| `-c`, `--count` | int | `6` | Maximum number of structures to download from Viro3D. |
| `-i`, `--input-folder` | path | `None` | Path to local folder containing `.pdb` or `.cif` structures (skips Viro3D download). |
| `--aligner` | choice | `foldmason` | Multiple sequence alignment engine: `foldmason` (default, structural MSTA) or `mafft` (sequence & 3Di alignment). |
| `--tree-type` | choice | `both` | Trees to infer: `3di` (structural), `aa` (sequence), `both`, or `tanglegram`. |
| `-m`, `--method` | choice | `iqtree` | Phylogeny method: `iqtree` (ML with 3Di matrices) or `foldmason` (progressive guide tree). |
| `--matrix` | choice | `both` | 3Di substitution matrix: `alphafold`, `esmfold`, or `both` (ModelFinder auto-select). |
| `--rate-heterogeneity`| str | `auto` | Rate heterogeneity model: `auto` (ModelFinder Plus), `+G4`, `+I+G4`, or `+R`. |
| `-b`, `--bootstrap` | int | `1000` | Number of Ultrafast Bootstrap (UFBoot2) replicates (`0` to disable). |
| `--alrt` | int | `1000` | Number of SH-aLRT branch test replicates (`0` to disable). |
| `-t`, `--threads` | str/int | `AUTO` | CPU threads for FoldMason and IQ-TREE (`AUTO` uses all available cores). |
| `-meta`, `--metadata` | path | `None` | Path to custom metadata file (`.xlsx`, `.csv`, `.tsv`, `.json`). |
| `--fast` | flag | `False` | Enables fast heuristic search mode (skips bootstrap for maximum speed). |
| `--embed` | flag | `False` | Extracts PLM embeddings and constructs a hierarchical clustering tree (ESM-2 650M / ESM-C 600M). |
| `--embed-model` | choice | `esm2` | PLM model architecture: `esm2` or `esmc`. |
| `--embed-clustering` | choice | `upgma` | Hierarchical clustering method: `upgma` (average linkage) or `nj` (neighbor-joining). |
| `--embed-metric` | choice | `cosine` | Pairwise distance metric: `cosine` or `euclidean`. |
| `-o`, `--output-dir` | path | `glycoprotein_workflow` | Output folder where structures, alignments, and trees are stored. |

---

### Subcommand: `fetch`
Searches the Viro3D database via REST API, downloads PDB coordinate files concurrently using a thread pool, and saves a JSON metadata index.

```bash
python scripts/viral_phylogenetics.py fetch \
  --qualifier rdrp \
  --max-sequences 100 \
  --output-dir viro_rdrp_structures
```

**Key Outputs:**
- `viro_rdrp_structures/*.pdb`: Relaxed 3D coordinate files.
- `viro_rdrp_structures/taxa_metadata.json`: Index containing protein name, taxonomy family, genus, host, and UniProt IDs.

---

### Subcommand: `align`
Executes multiple structural or sequence alignment on a directory of structures. By default, **FoldMason** is used for rapid multiple structural alignment (MSTA). Alternatively, **MAFFT** can be selected, in which case 3Di structural sequences are aligned using the official Foldseek 3Di scoring substitution matrix ([`matrices/mat3di.out`](matrices/mat3di.out)).

```bash
# Default: FoldMason (fast structural MSTA)
python scripts/viral_phylogenetics.py align \
  --folder viro_rdrp_structures \
  --output-dir rdrp_alignments

# Alternative: MAFFT (with mat3di.out for 3Di structural alignment)
python scripts/viral_phylogenetics.py align \
  --folder viro_rdrp_structures \
  --aligner mafft \
  --output-dir rdrp_alignments
```

| Argument | Type | Default | Description |
| :--- | :--- | :--- | :--- |
| `-i`, `--folder` | path | `viro_3d_structures` | Directory containing `.pdb` or `.cif` coordinate files. |
| `-o`, `--output-dir` | path | `foldmason_alignments` | Output directory where alignments are written. |
| `--aligner` | choice | `foldmason` | Multiple sequence alignment engine: `foldmason` (default, structural MSTA) or `mafft` (sequence & 3Di alignment). |
| `--foldmason-bin` | path | auto-detected | Custom path to the `foldmason` executable. |
| `--mafft-bin` | path | auto-detected | Custom path to the `mafft` executable (searches PATH and Conda envs). |
| `--mafft-matrix` | path | `matrices/mat3di.out` | Substitution matrix for MAFFT 3Di alignment. Automatically fetched from Foldseek if missing. |

**Key Outputs:**
- `rdrp_alignments/foldmason.fasta_3di.fa`: Multiple sequence alignment in the **3Di structural alphabet**.
- `rdrp_alignments/foldmason.fasta_aa.fa`: Multiple sequence alignment in standard **amino acid sequence**.
- `rdrp_alignments/foldmason.fasta.nw`: Guide tree (when using FoldMason).
- `rdrp_alignments/foldmason.fasta.html`: Standalone FoldMason alignment inspection page (when using FoldMason).

---

### Subcommand: `tree`
Constructs maximum likelihood phylogenetic trees from pre-computed FASTA alignments using IQ-TREE with empirical 3Di substitution matrices.

```bash
python scripts/viral_phylogenetics.py tree \
  --alignment rdrp_alignments/foldmason.fasta_3di.fa \
  --tree-type both \
  --matrix alphafold \
  --bootstrap 1000 \
  --output-dir rdrp_phylogeny
```

*(Note: If you pass an alignment directory, e.g. `-s rdrp_alignments`, the pipeline automatically resolves the `foldmason.fasta_3di.fa` file inside it).*

---

### Subcommand: `embed`
Extracts high-dimensional protein representations using Protein Language Models (ESM-2 650M or ESM-C 600M), computes pairwise semantic distance matrices, and constructs hierarchical clustering trees in Newick format.

```bash
python scripts/viral_phylogenetics.py embed \
  --input rdrp_alignments/foldmason.fasta_aa.fa \
  --output-dir rdrp_phylogeny \
  --model esm2 \
  --clustering upgma \
  --metric cosine \
  --prefix rdrp
```

| Argument | Type | Default | Description |
| :--- | :--- | :--- | :--- |
| `-i`, `--input` | path | Required | Path to input FASTA file (e.g. `foldmason.fasta_aa.fa`) or directory of PDB structures. |
| `-o`, `--output-dir` | path | `plm_results` | Directory where the `.treefile` and compressed `.npz` embeddings are saved. |
| `-m`, `--model` | choice | `esm2` | PLM architecture: `esm2` (`facebook/esm2_t33_650M_UR50D`) or `esmc` (`biohub/ESMC-600M`). |
| `-c`, `--clustering` | choice | `upgma` | Hierarchical clustering algorithm: `upgma` (Average Linkage, Ultrametric) or `nj` (Neighbor-Joining, Additive). |
| `--metric` | choice | `cosine` | Pairwise distance metric: `cosine` ($1 - \cos(\theta)$), `euclidean` ($L_2$), or `l1` / `manhattan` (Taxicab norm). |
| `-p`, `--prefix` | str | `viral_plm` | Filename prefix for the output `.treefile` and `.npz` files. |
| `--batch-size` | int | `1` | Inference batch size (optimized for single-sequence unpadded inference on Apple Silicon GPU / CUDA). |
| `--max-length` | int | `1024` | Maximum sequence length for transformer token truncation. |

**Key Outputs:**
- `{output_dir}/{prefix}_tree_{model}_{metric}.treefile`: Hierarchical clustering tree in standard Newick format for the specified metric (or `{prefix}_tree_{model}.treefile` for default cosine).
- `{output_dir}/{prefix}_embeddings_{model}.npz`: Compressed NumPy archive containing the `taxa` identifiers, $1280$-dimensional mean-pooled residue embeddings, model name, distance metric, and clustering metadata. Instant re-clustering across any metric can be computed directly from this archive without re-running the neural model.

---

## 4. Custom Metadata Integration

The pipeline supports arbitrary user-supplied metadata across all common formats (`.xlsx`, `.csv`, `.tsv`, `.json`).

```bash
python scripts/viral_phylogenetics.py pipeline \
  --input-folder custom_pdbs/ \
  --metadata clinical_metadata.xlsx \
  --output-dir custom_annotated_workflow
```

### Flexible Metadata Schema
The metadata parser ([`scripts/metadata_handler.py`](scripts/metadata_handler.py)) automatically detects:
- **Taxon Identifier Column**: Recognizes `Record ID`, `Accession`, `Taxon`, `Name`, `ID`, `Protein ID`, or filename stems.
- **Categorical Attributes**: Discovers categorical columns (e.g. `Viral Family`, `Host Organism`, `Geography`, `Structural Fold`, `Binding Strength`) and assigns high-contrast categorical palettes.
- **Continuous Metrics**: Discovers numerical metrics (e.g. `pLDDT Confidence`, `Sequence Length`, `IC50`, `Resolution`) and computes min/max ranges for continuous gradient coloring.

---

## 5. Interactive Web Visualization Suite

The suite generates an interactive browser application: [`interactive_tree.html`](interactive_tree.html).

### Launching the Viewer
- **Option 1 (Interactive Bridge Server)**:
  ```bash
  python scripts/serve_interactive.py
  ```
  Opens `http://localhost:8000` with direct in-browser pipeline execution and Viro3D querying enabled.
- **Option 2 (Direct Browser Open)**:
  Double-click `interactive_tree.html` or open in Chrome, Firefox, Safari, or Edge.

---

### Tree Layouts & Exploration
Switch between 5 visualization layouts via the **Visualization Layout** buttons in the left sidebar:
1. **Phylo (Rectangular Phylogram)**: Branch lengths represent evolutionary substitutions per structural 3Di site or amino acid substitution.
2. **Clado (Cladogram)**: Equal-length branches highlighting branching order and topological relationships.
3. **Radial Tree**: Circular tree layout projecting outwards from the root, ideal for inspecting broad clades.
4. **Unrooted Star Tree**: Equal-angle unrooted star tree layout showing global evolutionary divergence without rooting bias.
5. **Tangle (Dual Tanglegram)**: Side-by-side comparison of dual phylogenetic topologies with interactive connecting ribbons.

### Tree Datasets & Multi-Metric PLM Embedding Trees
In the **Tree Dataset** dropdown, switch between:
- **🏛️ 3Di Structural Tree**: FoldMason multiple structural alignment + IQ-TREE empirical 3Di substitution models (`Q.3Di.AF` / `Q.3Di.LLM`).
- **🧬 Amino Acid Sequence Tree**: Standard sequence alignment + IQ-TREE maximum likelihood (`LG+G4`).
- **🤖 ESM-2 PLM Tree**: Protein Language Model mean-pooled residue embeddings (`facebook/esm2_t33_650M_UR50D`) clustered with SciPy C-accelerated UPGMA hierarchical clustering.

When **ESM-2 PLM Tree** is selected, a dedicated **Embedding Distance Metric** sub-selector appears directly beneath:
- **📐 Cosine Distance (Directional Semantic Fold)**: Measures angular alignment between normalized embedding vectors ($d = 1 - \frac{u \cdot v}{\|u\|_2 \|v\|_2}$), invariant to sequence length norm scaling.
- **📏 Euclidean Distance ($L_2$ Geometric)**: Measures absolute Pythagorean distance in 1280-dimensional PLM latent space ($d = \|u - v\|_2$).
- **🧱 Manhattan / $L_1$ Distance (Taxicab Norm)**: Sums absolute differences across all embedding dimensions ($d = \sum_k |u_k - v_k|$), offering higher sensitivity to sparse residue activations.

Selecting any metric instantaneously re-renders the phylogenetic tree, adjusts branch lengths, updates phylogenetic metrics, and rescales the viewport.

---

### Dynamic Metadata & Clade Grouping
- **Color Nodes & Clades By Dropdown**: Select any categorical or continuous column from the active dataset. All tree nodes, leaf labels, and branch paths dynamically recolor.
- **Dynamic Legend**: Displays categorical counts or continuous gradient bars. Clicking any categorical badge in the legend collapses or expands all matching leaves in the tree.
- **Group & Collapse Clades By (Clades Tab)**: Automatically aggregates clades sharing $\ge 75\%$ metadata homogeneity into triangular collapsed clade wedges.

---

### Scoped Clade Analysis & Silhouette Cuts (High-Divergence Datasets)
When analyzing massive or highly divergent sequence cohorts (e.g. across disparate viral families or de novo binder folds), viewing the entire cohort simultaneously can obscure fine-grained intra-lineage relationships and create excessive tanglegram ribbon crossings. The **Clades Tab** provides an expert **Scoped Clade Partitioning Engine**:

1. **Universal Partition Sources & Silhouette Metrics**:
   - **🤖 ESM-2 PLM Tree (Silhouette Guided)**: Partitions sequences using hierarchical clustering on 1280-dimensional PLM embeddings.
   - **🏛️ 3Di Structural Tree (Empirical Substitution Distance)**: Partitions sequences using IQ-TREE patristic distances (`Q.3Di.LLM` / `Q.3Di.AF`) with precomputed silhouette evaluation.
   - **🧬 Amino Acid Tree (Sequence Substitution Distance)**: Partitions sequences using IQ-TREE empirical sequence distances (`LG+G4`) with precomputed silhouette evaluation.

2. **Universal Silhouette Curves & Mathematical Cut Suggestions**:
   - For **all three partition sources** (ESM-2, 3Di, and AA), the interface plots an interactive SVG sparkline of the mean Silhouette Score $S(k)$ across cut levels up to $k=35\text{--}40$.
   - **Suggested Cuts**: Local maxima and optimal clustering trade-offs are flagged with ⭐ badges:
     - **3Di Structural Trees**: e.g., $k=10$ ($S=0.464$) and $k=15$ ($S=0.409$) on 500 glycoproteins; $k=3$ ($S=0.505$) on the Nipah library.
     - **Amino Acid Trees**: e.g., $k=35$ ($S=0.396$) and $k=26$ ($S=0.321$) on 500 glycoproteins; $k=35$ ($S=0.220$) on the Nipah library.
     - **ESM-2 PLM Trees**: e.g., $k=2$ ($S=0.720$), $k=5$ ($S=0.539$), $k=16$ ($S=0.477$) on 500 glycoproteins; $k=3$ ($S=0.690$) and $k=5$ ($S=0.549$) on Nipah.
   - Clicking any peak badge instantly jumps the cut slider to that optimal granularity.
   - The cut slider allows fine-grained cuts up to $k=35\text{--}40$ for large cohorts.

3. **Singleton & Minor Subclade Filtering ($\ge 5$ Sequences)**:
   - High-granularity cuts often create solitary taxa or basal singletons that clutter the roster.
   - The clade roster **automatically hides singletons and minor lineages with $<5$ sequences by default**, ensuring only substantial subclades ($\ge 5$ sequences) are immediately prominent.
   - A collapsible drawer at the bottom of the roster (`[ 🔍 X Minor Lineages & Singletons (<5 seqs) ]`) keeps singletons organized and accessible.
   - A header toggle `[ 🛡️ ≥5 Seqs: ON / OFF ]` lets users switch between filtered and exhaustive views at any time.

4. **Scoped Analysis Isolation Mode**:
   - Each resulting clade card displays its taxon count, cohort percentage, majority annotation (e.g., `Mainly Alpha (98%)`), and cohesion score.
   - Clicking **`🔍 Scope Entire Analysis`** isolates the active workspace strictly to that lineage:
     - **Tree Viewports**: All tree layouts (Phylogram, Cladogram, Radial, Unrooted) prune external branches and re-scale to the scoped clade.
     - **Tanglegram**: Side-by-side trees prune to the clade, resolving divergent spaghetti into clean, high-resolution structural vs sequence correspondence.
     - **MSA Drawer**: Dynamically filters to display only sequences from the scoped clade, stripping all-gap columns and recalculating column conservation.
     - **3D Backbone Fold**: Automatically selects and renders the 3D backbone coordinates for a representative structure from the clade.
     - **Top Floating Banner**: A sleek glowing banner (`Scoped Clade View: Clade X`) allows one-click **`✕ Reset to Full Cohort`** at any time.


---

### Structural, Sequence & PLM Congruence Analysis
When in **Tanglegram** mode, select any pairwise comparison from the **Comparison Pair** selector:
- **Structural vs Sequence / PLM**:
  - `🏛️ 3Di Structural vs 🧬 Amino Acid`: Exposes fold conservation despite primary sequence divergence.
  - `🏛️ 3Di Structural vs 📐 ESM-2 (Cosine)`: Compares tertiary backbone structural state against angular semantic PLM representations.
  - `🏛️ 3Di Structural vs 📏 ESM-2 (Euclidean)`: Compares tertiary structure against geometric $L_2$ embedding distances.
  - `🏛️ 3Di Structural vs 🧱 ESM-2 (Manhattan/L1)`: Compares tertiary structure against taxicab feature norm distances.
- **Sequence vs PLM Embeddings**:
  - `🧬 Amino Acid vs 📐 ESM-2 (Cosine)`: Tests how transformer contextual representations align with traditional substitution matrices (`LG+G4`).
  - `🧬 Amino Acid vs 📏 ESM-2 (Euclidean)`: Sequence evolution vs raw Euclidean PLM space.
  - `🧬 Amino Acid vs 🧱 ESM-2 (Manhattan/L1)`: Sequence evolution vs Manhattan PLM norm.
- **PLM Metric Cross-Comparisons**:
  - `📐 ESM-2 (Cosine) vs 📏 ESM-2 (Euclidean)`: Contrasts angular direction with geometric magnitude.
  - `📐 ESM-2 (Cosine) vs 🧱 ESM-2 (Manhattan/L1)`: Contrasts angular direction with $L_1$ sparsity.

Click **`📊 Full Report`** to inspect quantitative discordance metrics for the currently active pair:
- **Robinson-Foulds Congruence ($RF$)**: Strict topological partition distance and shared bipartition count.
- **Cophenetic Correlation ($r$)**: Pearson correlation on pairwise patristic branch substitution distances.
- **Planar Crossing Counter**: Real-time count of intersecting tangle connectors.
- **Tanglegram Alignment Modes**: Toggle between **True Topology** (native branch order), **Min-Crossing Untangled** (optimal barycenter clade rotations), and **Aligned** (horizontal connector matching).

---

### Multiple Sequence Alignment (MSA) Viewer & Minimap
Docked at the bottom of the workspace:
- **Dual Alphabet Toggle**: Switch instantly between **🥩 Amino Acid (AA)** and **🧊 3Di Structural Alphabet** modes.
- **Dynamic All-Gap Column Stripping (`[ ✂️ Strip Gaps: ON ]`)**:
  - Automatically identifies and removes columns that consist entirely of gaps (`-` or `.`) across the currently visible sequence cohort.
  - In filtered subsets or scoped clades, inter-family insertion gaps are instantly eliminated (e.g., in the Nipah library, 81 gap-only columns are removed in the *Mainly Alpha* clade and 151 in the *Mainly Beta* clade).
  - The summary badge highlights the reduction (e.g. `564 seqs • 452 cols • 81 gap-only cols removed • AA`).
  - The column ruler preserves the original global sequence coordinate for each residue, allowing direct reference to wildtype positions.
  - Can be toggled on/off at any time via the toolbar button.
- **Bioinformatics Color Schemes**: ClustalX standard, Zappo physicochemical, Hydrophobicity gradient, Identity consensus, and 3Di Fold Geometry ($\alpha$-helices in blue, $\beta$-strands in amber, turns in green, loops in purple).
- **Docked 2D Radar Minimap**: A whole-alignment density heatmap showing column conservation across the visible dataset. Click or drag anywhere on the radar to navigate across thousands of residues.
- **Decoupled Scrolling**: Alignment navigation is decoupled from tree movement by default (`[ 🔓 Decoupled ]`). Panning and zooming the tree will not interfere with alignment coordinates.

---

### 3D Structure Viewer
- Pinned at the bottom of the left sidebar and in the interactive hover tooltip.
- Renders the interactive **C$\alpha$ backbone tube** in 3D WebGL with per-residue AlphaFold pLDDT confidence coloring ($>90$ very high, $70\text{--}90$ confident, $50\text{--}70$ low, $<50$ disordered).
- Left-click drag to rotate; scroll wheel to zoom.

---

### Zoom-Adaptive Label Scaling
- On large trees ($N \ge 80$), labels on radial and unrooted layouts automatically shrink as you zoom in ($k > 1$), preventing overlapping text from obscuring dense clades and internal nodes.
- Adjust the **Tip Label Size** slider (`0px` to `16px`) in the `⚙️ View` tab to set any preferred font size or set to `0px` (`Hidden`) to view pure branch geometry.

---

## 6. Interactive Pipeline Studio Tab

Located in the **`🚀 Run`** tab in the sidebar of `interactive_tree.html`:

### Features:
1. **Live Viro3D Query Builder**:
   - Enter any protein qualifier (e.g. `glycoprotein`, `rdrp`, `spike`, `capsid`) or click a quick preset.
   - Adjust the target count slider ($6$ to $500$).
   - Click **`🔍 Check Viro3D Available Structures`** to preview total database hits and sample accession records.
2. **Instant Terminal Command Generator**:
   - Generates the exact, syntax-validated CLI command with your chosen matrix, tree type, and thread configuration.
   - Click **`📋 Copy Command`** or **`💾 Download .sh`** to run in any terminal.
3. **In-Browser Execution Bridge**:
   - If `scripts/serve_interactive.py` is running, click **`▶️ Run Pipeline in Background`** to trigger the pipeline directly from the browser!
   - Real-time execution logs stream into the integrated console window.

---

## 7. Troubleshooting & Best Practices

### Issue: `Ultrafast bootstrap does not work with -fast, -te or -n option`
- **Cause**: In IQ-TREE, `--fast` (FastTree-like heuristic) is mutually exclusive with `-B` (Ultrafast bootstrap).
- **Solution**: Do not combine `--fast` with `-b 1000`. In `scripts/viral_phylogenetics.py`, the pipeline automatically resolves this: if bootstrap is requested, `--fast` is omitted; if `--fast` is passed, bootstrap is gracefully suppressed.

### Issue: `Reading alignment file ... ERROR: File not found`
- **Cause**: Passing an invalid or nonexistent relative file path or pointing to a directory instead of a FASTA file.
- **Solution**: The pipeline now features automated path sanitization: it auto-corrects missing leading slashes and automatically resolves directory paths to `foldmason.fasta_3di.fa`.

### Issue: ModelFinder Takes Very Long on Large Cohorts ($N > 100$)
- **Solution**: For datasets with $\ge 80$ structures, avoid `-m MFP` (testing 120+ models). The pipeline automatically pins `-m Q.3Di.AF+G4` for 3Di and `-m LG+G4` for amino acids, accelerating tree inference from hours to under 2 minutes.

---

## Authors & Citation
- **Viro3D**: University of Glasgow Centre for Virus Research (CVR).
- **FoldMason**: Steinegger Lab (*"FoldMason: fast and accurate multiple protein structure alignment"*, Nature Communications).
- **3Di Empirical Matrices**: Georg Hochberg et al. (*"A general substitution matrix for structural phylogenetics"*, Edmond Dataverse [doi:10.17617/3.1MJJBH](https://doi.org/10.17617/3.1MJJBH)).
- **IQ-TREE 2**: Minh et al. (*"IQ-TREE 2: New Models and Efficient Methods for Phylogenetic Inference"*, Molecular Biology and Evolution).
