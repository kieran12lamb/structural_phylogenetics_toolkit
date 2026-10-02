# Viral Structural Phylogenetics Suite — Comprehensive Usage Guide

An end-to-end framework for viral protein structural phylogenetics. The suite combines:
1. **Viro3D / AlphaFold DB Automated Retrieval**: Programmatic retrieval of AI-predicted viral 3D structures (ColabFold/AlphaFold2 and ESMFold) and rich clinical/taxonomic metadata.
2. **FoldMason Multiple Structural Alignment (MSTA)**: 3D coordinate superposition and translation into the 20-state **3Di** (3D interaction) structural alphabet.
3. **Similarity Clustering (Foldseek / MMseqs2)**: Partitioning of large, divergent structure sets into groups that genuinely align, before alignment and tree inference.
4. **IQ-TREE Empirical Structural Phylogenetics**: Maximum likelihood phylogenetic inference using empirical 3Di substitution matrices (`Q.3Di.AF`, `Q.3Di.LLM`) with automated rate heterogeneity selection and ultrafast bootstrap, plus fast whole-set trees with **VeryFastTree** under the same matrices.
5. **Run Manifest & Provenance**: A durable record of every run — command, settings, tool versions, input fingerprint, and the outcome of each stage.
6. **Interactive SVG & WebGL Visualization Suite**: High-performance interactive browser interface with dual tanglegrams, dynamic metadata coloring/filtering, topological congruence metrics, synchronized Multiple Sequence Alignment (MSA) viewer with 2D radar overview, and live 3D protein backbone visualization.
7. **Interactive Web Pipeline Studio**: Direct in-browser Viro3D querying, CLI command generation, and background pipeline execution bridge.

---

## Table of Contents
1. [Prerequisites & Installation](#1-prerequisites--installation)
2. [Quick-Start Workflows](#2-quick-start-workflows)
3. [CLI Reference](#3-cli-reference)
   - [`pipeline` (End-to-End)](#subcommand-pipeline)
   - [`fetch` (Structure Download)](#subcommand-fetch)
   - [`align` (3Di Alignment)](#subcommand-align)
   - [`tree` (Phylogeny)](#subcommand-tree)
   - [`embed` (PLM Embeddings)](#subcommand-embed)
   - [`report` (Tree Report)](#subcommand-report)
4. [Handling Divergent Datasets: Clustering & Recursion](#4-handling-divergent-datasets-clustering--recursion)
5. [Tree Inference Strategy](#5-tree-inference-strategy)
6. [Run Manifest, Preflight & Provenance](#6-run-manifest-preflight--provenance)
7. [Custom Metadata Integration](#7-custom-metadata-integration)
8. [Interactive Web Visualization Suite (`interactive_tree.html`)](#8-interactive-web-visualization-suite)
9. [Interactive Pipeline Studio Tab](#9-interactive-pipeline-studio-tab)
10. [Troubleshooting & Best Practices](#10-troubleshooting--best-practices)

---

## 1. Prerequisites & Installation

### Environment Setup
The recommended route is the installer, which creates the `spt` conda environment from `environment.yml`, installs the package, verifies every tool and package, and runs the tests:

```bash
./install.sh
conda activate spt
```

Or manually:

```bash
conda env create -f environment.yml
conda activate spt
pip install -e . --no-deps
```

The environment provides Python 3.11 and these external tools, all from bioconda:

| Tool | Used for |
| :--- | :--- |
| `foldmason` | Structural alignment, and 3Di extraction when MAFFT aligns |
| `mafft` | Alternative aligner (`--aligner mafft`) |
| `iqtree` | Maximum-likelihood trees for each cluster |
| `foldseek` | Structural clustering (`--cluster-mode structural`, the default) |
| `mmseqs2` | Sequence clustering (`--cluster-mode sequence`) |
| `veryfasttree` / `fasttree` | Whole-set tree and IQ-TREE starting trees (VeryFastTree preferred) |

`pip install -e .` alone installs only the Python packages; the tools above must come from conda. Add `pip install -e ".[ml]"` for `torch` and `transformers`, which `--embed` needs.

### Verification
```bash
foldmason version
iqtree --version
foldseek version
mmseqs version
VeryFastTree -help | head -1
```

You do not need to check these by hand before a run: the pipeline's **preflight** check (see [section 6](#6-run-manifest-preflight--provenance)) verifies every tool the chosen options need and stops with a clear list if anything is missing.

### GPU Acceleration for Embeddings
`--embed` uses PyTorch. If embeddings run on CPU although a GPU is present, the installed torch build most likely targets a newer CUDA version than your GPU driver supports. `install.sh` reports this case and prints the fix — installing a torch build that matches the driver, e.g. `pip install torch --index-url https://download.pytorch.org/whl/cu121`. CPU embedding works, only more slowly (about 20 minutes for ~1,200 ESM-2 sequences).

The pipeline chooses the interpreter for the embedding step by checking which candidate environments can actually import `torch` and `transformers`, preferring one with a working GPU — not simply whichever conda environment is active.

### Empirical 3Di Substitution Matrices
The pipeline utilizes empirical 3Di substitution models inferred from structural protein databases:
- **`matrices/Q.3Di.AF`**: Inferred from AlphaFold structures (Edmond ID: `311466`).
- **`matrices/Q.3Di.LLM`**: Inferred from ESMFold/ProstT5 translations (Edmond ID: `311467`).

*(If missing, they are downloaded automatically from the Edmond Dataverse repository on first run.)*

Run the pipeline **from the repository root**: `matrices/` is resolved relative to the working directory. When IQ-TREE's model selection is used, the pipeline also writes upper-case copies to `MATRICES/`, because IQ-TREE upper-cases `-mset` paths (see [Troubleshooting](#10-troubleshooting--best-practices)). `MATRICES/` is generated and ignored by git.

---

## 2. Quick-Start Workflows

### Scenario A: Online Viro3D Query to Phylogenetic Trees
Download 50 viral glycoproteins directly from Viro3D, align with FoldMason, and build both 3Di structural and amino acid sequence trees:
```bash
viral-phylo pipeline \
  --qualifier glycoprotein \
  --count 50 \
  --tree-type both \
  --matrix alphafold \
  --output-dir results/glycoprotein_50_workflow
```

### Scenario B: Process a Local Directory of Structures
If you already have PDB or mmCIF files on disk (e.g. from local AlphaFold/ESMFold predictions):
```bash
viral-phylo pipeline \
  --input-folder path/to/my_structures \
  --tree-type both \
  --matrix both \
  --output-dir results/my_local_workflow
```

### Scenario C: A Large or Divergent Structure Set
For hundreds to thousands of structures that do not share a single fold, cluster first so every alignment and tree is built from proteins that align:
```bash
viral-phylo pipeline \
  --input-folder path/to/structures \
  --metadata path/to/metadata.csv \
  --output-dir results/large_run \
  --multi-alignment --cluster-recurse \
  --tree-type both --matrix alphafold \
  --embed --embed-clustering upgma
```
See [section 4](#4-handling-divergent-datasets-clustering--recursion) for what each step does and how to tune it.

### Scenario D: Launch the Interactive Web Suite
Every pipeline run writes `<output-dir>/interactive_tree.html`; open it directly in a browser. To use in-browser pipeline execution, start the bridge server instead:
```bash
python scripts/serve_interactive.py
```

---

## 3. CLI Reference

The CLI entrypoint is `viral-phylo` (equivalently [`scripts/viral_phylogenetics.py`](scripts/viral_phylogenetics.py)). It provides these subcommands:

```
usage: viral-phylo [-h] {fetch,align,tree,embed,report,pipeline} ...
```

---

### Subcommand: `pipeline`
Runs the complete workflow: preflight checks $\rightarrow$ fetch or load structures $\rightarrow$ whole-set alignment $\rightarrow$ (optionally) clustering and per-cluster alignment $\rightarrow$ trees $\rightarrow$ metadata $\rightarrow$ (optionally) PLM embeddings $\rightarrow$ interactive dashboard, with a run manifest updated after every stage.

```bash
viral-phylo pipeline [options]
```

#### Input

| Argument | Type | Default | Description |
| :--- | :--- | :--- | :--- |
| `-i`, `--input-folder` | path | `None` | Local folder of `.pdb` / `.cif` structures (skips remote download). |
| `--source` | choice | `viro3d` | Structure repository: `viro3d`, `alphafold` or `afdb` (AlphaFold DB). |
| `-q`, `--qualifier` | str | `glycoprotein` | Protein search term when fetching (e.g. `glycoprotein`, `rdrp`, `spike`, `capsid`). |
| `-u`, `--uniprot` | str | `None` | Comma-separated UniProt ID(s) for AlphaFold DB (e.g. `P00520,P04637`). |
| `--uniprot-file` | path | `None` | Text file with one UniProt ID per line. |
| `-c`, `--count` | int | `6` | Maximum number of structures to download. |
| `--format` | choice | `pdb` | AlphaFold DB coordinate format: `pdb` or `cif`. |
| `--download-pae` | flag | off | Also download the PAE matrix from AlphaFold DB. |
| `-meta`, `--metadata` | path | `None` | Metadata table (`.xlsx`, `.csv`, `.tsv`, `.json`). |
| `-o`, `--output-dir` | path | `glycoprotein_workflow` | Output folder for the run. |

#### Alignment

| Argument | Type | Default | Description |
| :--- | :--- | :--- | :--- |
| `--aligner` | choice | `foldmason` | `foldmason` (structural MSTA) or `mafft`. The choice applies to every alignment in the run — whole-set and per-cluster — and to both alphabets. |
| `--mafft-matrix` | path | `matrices/mat3di.out` | 3Di substitution matrix for MAFFT. |
| `--mafft-bin` | path | auto | Explicit `mafft` path. |
| `--mafft-leavegappyregion` | flag | off | Pass `--leavegappyregion` to MAFFT's amino-acid alignment, which inserts fewer gaps into gap-rich regions. |

#### Clustering (with `--multi-alignment`)

| Argument | Type | Default | Description |
| :--- | :--- | :--- | :--- |
| `--multi-alignment` | flag | off | Partition the structures into clusters, and align and build trees per cluster. |
| `--cluster-mode` | choice | `structural` | `structural` (Foldseek, 3Di+AA), `sequence` (MMseqs2, amino acid), or `coverage` (the legacy greedy partitioner, kept for reproducing older runs). |
| `--cluster-evalue` | float | `1.0` | Foldseek e-value. Depends on database size, so not comparable across datasets of different sizes. |
| `--cluster-coverage` | float | `0.5` | Alignment coverage needed to join a cluster. |
| `--cluster-tmscore` | float | `None` | Use TM-align structural clustering at this TM-score instead of 3Di+AA. Comparable across datasets, but slower. |
| `--cluster-min-seq-id` | float | `0.5` | MMseqs2 identity for `--cluster-mode sequence`. |
| `--cluster-sweep` | flag | off | Try a grid of thresholds and select one, instead of using the operating point above (see [section 4](#4-handling-divergent-datasets-clustering--recursion)). |
| `--cluster-min-size` | int | `10` | Cluster size counted as "usable" when scoring a sweep. |
| `--cluster-min-retained` | float | `0.70` | Fraction of taxa a sweep partition must keep in usable clusters. |
| `--min-seq-length` | int | `50` | Exclude sequences shorter than this before clustering. |
| `--min-length-frac` | float | `0.4` | Also exclude sequences shorter than this fraction of the median length. |
| `--full-set-max-gap` | float | `0.50` | Whole-set alignment is reported unusable at this gap fraction or above. |
| `--full-set-max-expansion` | float | `1.5` | ...or when it exceeds this multiple of the median ungapped length. |
| `--min-coverage` | float | `0.70` | Coverage threshold for the legacy `coverage` mode and `--filter-coverage`. |
| `--filter-coverage` | flag | off | Legacy single-alignment coverage filter. |

#### Recursion (with `--cluster-recurse`)

| Argument | Type | Default | Description |
| :--- | :--- | :--- | :--- |
| `--cluster-recurse` | flag | off | Split any cluster whose own alignment is still too gappy. |
| `--recurse-max-gap` | float | `0.65` | Split a cluster whose alignment reaches this gap fraction. |
| `--recurse-max-expansion` | float | `3.0` | ...or exceeds this multiple of its median ungapped length. |
| `--recurse-max-depth` | int | `2` | How many times a cluster may be split. |
| `--recurse-min-kept` | float | `0.70` | Only accept a split keeping at least this fraction of the cluster's taxa in clusters of 4+. |
| `--recurse-on` | choice | `3di` | Which alignment decides a split: `3di`, or `both` (also split when the AA alignment is over threshold). |

#### Trees

| Argument | Type | Default | Description |
| :--- | :--- | :--- | :--- |
| `--tree-type` | choice | `both` | `3di` (structural), `aa` (sequence), `both`, or `tanglegram`. |
| `--master-method` | choice | `fasttree` | Method for the whole-set tree: `fasttree` (VeryFastTree/FastTree), `iqtree`, or `foldmason` (guide tree). |
| `-m`, `--method` | choice | `iqtree` | Method for the per-cluster trees: `iqtree`, `fasttree`, or `foldmason`. |
| `--matrix` | choice | `both` | 3Di matrix: `alphafold`/`af`, `esmfold`/`llm`, or `both`/`auto` (IQ-TREE ModelFinder chooses between them). |
| `--rate-heterogeneity` | str | `auto` | `auto` (ModelFinder), `+G4`, `+I+G4`, or `+R`. |
| `-b`, `--bootstrap` | int | `1000` | Ultrafast bootstrap replicates (`0` to disable). |
| `--alrt` | int | `1000` | SH-aLRT replicates (`0` to disable). |
| `--guide-trees` / `--no-guide-trees` | flag | on | Start IQ-TREE from a validated VeryFastTree tree on alignments of `--guide-tree-min-taxa` taxa or more. |
| `--guide-tree-min-taxa` | int | `80` | Smallest alignment given a starting tree. |
| `-t`, `--threads` | str/int | `8` | **Maximum CPU threads for the whole run.** Stages run one at a time and every tool is capped at this count: FoldMason, MAFFT, Foldseek/MMseqs2, IQ-TREE, VeryFastTree/FastTree, PyTorch (ESM embeddings) and NumPy/BLAS. The cap is recorded under `resources` in the manifest. `AUTO` removes it: tools then use every core, and IQ-TREE benchmarks every thread count, which is pathologically slow on many-core machines. |
| `--tree-threads` | str | `8` | Threads for each per-cluster tree job. |
| `--tree-jobs` | int | `--threads` ÷ `--tree-threads` | Cluster trees built at once, largest first. Never more than fit within `--threads` (32 ÷ 8 = 4 jobs). Each job logs to `phylogeny/<cluster>/tree_job.log`; the manifest records every job's command, duration and exit status under `trees.cluster_jobs`. |
| `--fast` | flag | off | IQ-TREE fast heuristic search (disables bootstrap). |

#### Embeddings

| Argument | Type | Default | Description |
| :--- | :--- | :--- | :--- |
| `--embed` | flag | off | Compute PLM embeddings and a hierarchical clustering tree. |
| `--embed-model` | choice | `esm2` | `esm2` (ESM-2 650M) or `esmc` (ESM-C 600M). |
| `--embed-clustering` | choice | `upgma` | `upgma`, `average`, `complete`, `single`, `ward`, or `nj`. All but `nj` also produce a silhouette profile. |
| `--embed-metric` | choice | `cosine` | `cosine`, `euclidean`, or `l1` (`cityblock`, `manhattan`). |

#### Fast vs. Thorough Pipeline Profiles:

- **⚡ Fast Mode (Rapid Screening & Prototyping)**:
  ```bash
  viral-phylo pipeline \
    --input-folder tests/fixtures/structures \
    --tree-type both \
    --matrix alphafold \
    --fast \
    --embed \
    --output-dir results/fast_workflow
  ```

- **🔬 Thorough Mode (Publication-Grade Deep Phylogenetics)**:
  ```bash
  viral-phylo pipeline \
    --input-folder tests/fixtures/structures \
    --tree-type both \
    --matrix both \
    --rate-heterogeneity auto \
    --bootstrap 1000 \
    --alrt 1000 \
    --embed \
    --output-dir results/thorough_workflow
  ```

---

### Subcommand: `fetch`
Searches and downloads structures from **AlphaFold Database (AFDB)** or **Viro3D** concurrently using worker threads, saving coordinate files and a standardized `taxa_metadata.json` index.

#### 1. Download from AlphaFold Database (AFDB)
Query by specific UniProt accession ID(s) or protein search term:
```bash
# Fetch specific proteins by UniProt accession with PAE matrices
viral-phylo fetch \
  --source alphafold \
  -u "Q99720, P87666, P0DTC2, P08667, P03452, P03437" \
  --format pdb \
  --download-pae \
  --output-dir afdb_structures

# Query AlphaFold DB via search resolution (UniProt API)
viral-phylo fetch \
  --source alphafold \
  -q "henipavirus glycoprotein" \
  --max-sequences 6 \
  --output-dir afdb_structures
```

#### 2. Download from Viro3D
```bash
viral-phylo fetch \
  --source viro3d \
  --qualifier rdrp \
  --max-sequences 100 \
  --output-dir viro_rdrp_structures
```

**Key Outputs:**
- `afdb_structures/*.pdb` (or `*.cif`): Relaxed 3D coordinate files.
- `afdb_structures/AF-*-predicted_aligned_error.json`: Optional PAE error matrices.
- `afdb_structures/taxa_metadata.json`: Standardized metadata index containing protein name, organism, gene, length, mean pLDDT score, and pLDDT category.

---

### Subcommand: `align`
Executes multiple structural or sequence alignment on a directory of structures. By default, **FoldMason** is used for rapid multiple structural alignment (MSTA). Alternatively, **MAFFT** can be selected, in which case 3Di structural sequences are aligned using the official Foldseek 3Di scoring substitution matrix ([`matrices/mat3di.out`](matrices/mat3di.out)). FoldMason is still used with MAFFT, but only to extract the unaligned 3Di and amino-acid sequences from the structures; MAFFT then aligns each alphabet.

```bash
# Default: FoldMason (fast structural MSTA)
viral-phylo align \
  --folder viro_rdrp_structures \
  --output-dir rdrp_alignments

# Alternative: MAFFT (with mat3di.out for 3Di structural alignment)
viral-phylo align \
  --folder viro_rdrp_structures \
  --aligner mafft \
  --output-dir rdrp_alignments
```

| Argument | Type | Default | Description |
| :--- | :--- | :--- | :--- |
| `-i`, `--folder` | path | `viro_3d_structures` | Directory containing `.pdb` or `.cif` coordinate files. |
| `-o`, `--output-dir` | path | `foldmason_alignments` | Output directory where alignments are written. |
| `--aligner` | choice | `foldmason` | `foldmason` (default, structural MSTA) or `mafft` (sequence & 3Di alignment). |
| `--foldmason-bin` | path | auto-detected | Custom path to the `foldmason` executable. |
| `--mafft-bin` | path | auto-detected | Custom path to the `mafft` executable (searches PATH and Conda envs). |
| `--mafft-matrix` | path | `matrices/mat3di.out` | Substitution matrix for MAFFT 3Di alignment. Automatically fetched from Foldseek if missing. |
| `--multi-alignment` | flag | off | Partition with the **legacy coverage partitioner**. The similarity clustering described in [section 4](#4-handling-divergent-datasets-clustering--recursion) is available through `pipeline --multi-alignment`. |
| `--min-coverage` | float | `0.70` | Coverage threshold for the legacy partitioner and `--filter-coverage`. |
| `--filter-coverage` | flag | off | Keep only sequences meeting `--min-coverage` in a single alignment. |
| `-t`, `--threads` | str | `8` | Maximum CPU threads for FoldMason and MAFFT. `AUTO` = every core. |

**Key Outputs:**
- `rdrp_alignments/foldmason.fasta_3di.fa`: Multiple sequence alignment in the **3Di structural alphabet**.
- `rdrp_alignments/foldmason.fasta_aa.fa`: Multiple sequence alignment in standard **amino acid sequence**.
- `rdrp_alignments/foldmason.fasta.nw`: Guide tree (FoldMason only).
- `rdrp_alignments/foldmason.fasta.html`: Standalone FoldMason alignment inspection page (FoldMason only).
- `rdrp_alignments/mafft.fasta_3di.fa`, `mafft.fasta_aa.fa`: MAFFT output. It is also copied to the `foldmason.fasta_*.fa` names so downstream steps find it; in a pipeline run, `alignment_info.json` records which aligner really produced them.

With FoldMason the AA and 3Di alignments are one alignment written in two alphabets, so their columns correspond. MAFFT aligns each alphabet separately, so their lengths generally differ and columns do not correspond.

---

### Subcommand: `tree`
Constructs phylogenetic trees from pre-computed FASTA alignments: IQ-TREE maximum likelihood with empirical 3Di substitution matrices by default, or VeryFastTree/FastTree approximate maximum likelihood.

```bash
viral-phylo tree \
  --alignment rdrp_alignments/foldmason.fasta_3di.fa \
  --tree-type both \
  --matrix alphafold \
  --bootstrap 1000 \
  --output-dir rdrp_phylogeny
```

*(Note: If you pass an alignment directory, e.g. `--alignment rdrp_alignments`, the `foldmason.fasta_3di.fa` (or `mafft.fasta_3di.fa`) file inside it is used automatically, and the AA alignment beside it is found for `--tree-type both`.)*

| Argument | Type | Default | Description |
| :--- | :--- | :--- | :--- |
| `-a`, `--alignment` | path | `foldmason_alignments/foldmason.fasta_3di.fa` | 3Di alignment file, or a directory containing one. |
| `--alignment-aa` | path | auto-detected | Amino-acid alignment for the AA tree. |
| `--tree-type` | choice | `both` | `3di`, `aa`, `both`, or `tanglegram`. |
| `-m`, `--method` | choice | `iqtree` | `iqtree`, `fasttree` (VeryFastTree/FastTree; the `--matrix` 3Di model is applied via `-trans`), or `foldmason` (guide tree). |
| `--matrix` | choice | `both` | `alphafold`/`af`, `esmfold`/`llm`, or `both`/`auto`. |
| `--rate-heterogeneity` | str | `auto` | `auto` (ModelFinder), `+G4`, `+I+G4`, or `+R`. |
| `--criterion` | choice | `BIC` | Model selection criterion: `BIC`, `AIC`, or `AICc`. |
| `-b`, `--bootstrap` | int | `1000` | Ultrafast bootstrap replicates. |
| `--alrt` | int | `1000` | SH-aLRT replicates. |
| `-t`, `--threads` | str | `8` | Maximum CPU threads for IQ-TREE / VeryFastTree. |
| `--fast` | flag | off | IQ-TREE fast heuristic search. |
| `--guide-trees` / `--no-guide-trees` | flag | on | Start IQ-TREE from a validated VeryFastTree tree on alignments of `--guide-tree-min-taxa`+ taxa. |
| `--guide-tree-min-taxa` | int | `80` | Smallest alignment given a starting tree. |
| `-o`, `--output-dir` | path | `phylogeny_results` | Output directory. |
| `-p`, `--prefix` | str | `viral_tree` | Output filename prefix. |
| `--iqtree-bin` | path | auto-detected | Custom path to `iqtree`. |

---

### Subcommand: `embed`
Extracts high-dimensional protein representations using Protein Language Models (ESM-2 650M or ESM-C 600M), computes pairwise semantic distance matrices, and constructs hierarchical clustering trees in Newick format.

```bash
viral-phylo embed \
  --input rdrp_alignments/foldmason.fasta_aa.fa \
  --output-dir rdrp_phylogeny \
  --model esm2 \
  --clustering upgma \
  --metric cosine \
  --prefix rdrp
```

| Argument | Type | Default | Description |
| :--- | :--- | :--- | :--- |
| `-i`, `--input`, `--fasta` | path | Required | Input FASTA file (gaps are ignored) or directory of structures. |
| `-o`, `--output-dir` | path | `plm_results` | Directory where the `.treefile` and compressed `.npz` embeddings are saved. |
| `-m`, `--model` | choice | `esm2` | PLM architecture: `esm2` (`facebook/esm2_t33_650M_UR50D`) or `esmc` (`biohub/ESMC-600M`). |
| `-c`, `--clustering` | choice | `upgma` | `upgma` (average linkage), `average`, `complete`, `single`, `ward`, or `nj` (neighbour-joining). `nj` builds no linkage matrix, so it produces no silhouette profile; the run says so. |
| `--metric` | choice | `cosine` | Pairwise distance metric: `cosine` ($1 - \cos(\theta)$), `euclidean` ($L_2$), or `l1` / `cityblock` / `manhattan` (Taxicab norm). |
| `-p`, `--prefix` | str | `viral_plm` | Filename prefix for the output `.treefile` and `.npz` files. |
| `--batch-size` | int | `1` | Inference batch size (optimized for single-sequence unpadded inference on Apple Silicon GPU / CUDA). |
| `--max-length` | int | `1024` | Maximum sequence length for transformer token truncation. |
| `-t`, `--threads` | str | `8` | Maximum CPU threads for embedding (PyTorch on CPU) and clustering. `AUTO` = every core. |

**Key Outputs:**
- `{output_dir}/{prefix}_tree_{model}.treefile`: Hierarchical clustering tree in Newick format.
- `{output_dir}/{prefix}_embeddings_{model}.npz`: Compressed NumPy archive containing the `taxa` identifiers, $1280$-dimensional mean-pooled residue embeddings, model name, distance metric, and clustering metadata. Re-clustering across any metric can be computed directly from this archive without re-running the neural model.
- `{output_dir}/{prefix}_silhouette_{model}.json`: Silhouette profile across tree cuts (all clustering methods except `nj`).
- `{output_dir}/{prefix}_umap_{model}.json`: 2D UMAP projection used by the dashboard.

---

### Subcommand: `report`
Builds a standalone interactive HTML report for any existing phylogenetic tree (Newick or NEXUS, e.g. from IQ-TREE, RAxML, BEAST or MrBayes), optionally annotated with a metadata table. The report opens directly in a browser; no server is needed.

```bash
viral-phylo report \
  --tree my_analysis/tree.treefile \
  --metadata my_analysis/samples.csv \
  --output my_analysis/tree_report.html
```

| Argument | Type | Default | Description |
| :--- | :--- | :--- | :--- |
| `-t`, `--tree` | path | Required | Tree file in Newick (`.nwk`, `.treefile`, `.tre`) or NEXUS (`.nex`) format. NEXUS `translate` tables, quoted labels and `[&...]` annotations are handled. |
| `-meta`, `--metadata` | path | None | Optional metadata table (`.xlsx`, `.csv`, `.tsv`, `.json`) with one row per leaf. |
| `-o`, `--output` | path | `<tree>_report.html` | Output HTML path (defaults to next to the tree file). |
| `--title` | str | tree file name | Title shown in the Dataset Cohort selector. |

**Metadata matching:** the ID column is detected automatically as the column whose values best match the tree's leaf names. Matching ignores case, spaces vs. underscores and punctuation, and also accepts IDs that prefix a leaf name (e.g. `AAQ55251.1` matches `AAQ55251.1.2_9251`). Leaves without a metadata row are shown as **Not in metadata**. The sidebar reports how many leaves were matched.

The same thing can be done in any existing report from the browser: under **Dataset Cohort**, click **📂 Load tree file**, then **🏷️ Add metadata** to attach a CSV, TSV or JSON table. Excel files can only be read by the CLI.

---

## 4. Handling Divergent Datasets: Clustering & Recursion

Aligning a large, structurally heterogeneous set as one alignment tends to fail quietly: sequences that share no fold cannot share columns, so the aligner produces a staircase of gaps. On a 1,193-structure Nipah binder set the global alignment was **1,458 columns for a median sequence length of 120 — 91.6% gaps**, and a greedy partition built on that alignment gave 226 clusters, 104 of them singletons. `pipeline --multi-alignment` avoids this in five steps.

### 1. Quality Control
Sequences shorter than `--min-seq-length` (50) or `--min-length-frac` × median (0.4) are excluded before clustering. Every exclusion is recorded in the run manifest with its reason.

### 2. Whole-Set Gate
The full 3Di alignment is still built — it feeds the whole-set tree — and is judged: it is reported **not usable** as a single alignment if its gap fraction is ≥ `--full-set-max-gap` (0.50) or its length exceeds `--full-set-max-expansion` (1.5) × the median ungapped sequence length. Without `--multi-alignment`, a failing gate prints a suggestion to partition first.

### 3. Clustering
The default is **Foldseek structural clustering** with its native 3Di+amino-acid alignment at e-value 1.0 and 50% coverage, with no sequence-identity threshold (requiring sequence identity would defeat the purpose of structural clustering). On the Nipah set this retained 87% of taxa in 85 clusters, 21 of them with 10+ taxa — whereas MMseqs2 amino-acid clustering could not keep even half the taxa in clusters of 10+ at any threshold, because the sequences are too divergent.

Alternatives:
- `--cluster-mode sequence` — MMseqs2 on amino-acid sequences (`--cluster-min-seq-id`).
- `--cluster-tmscore 0.5` — TM-align clustering. Slower and scored worse on the Nipah set, but a TM-score threshold is comparable across datasets whereas an e-value depends on database size.
- `--cluster-mode coverage` — the legacy greedy partitioner, for reproducing older runs.

### 4. Sweep (optional, `--cluster-sweep`)
Instead of one operating point, a grid of thresholds × coverages is clustered and each partition scored. Selection is **constrained**: only partitions keeping at least `--cluster-min-retained` (0.70) of taxa in clusters of `--cluster-min-size` (10)+ are eligible, and among those the best silhouette wins. Maximising silhouette alone would favour the most fragmented partition, because near-singleton clusters are trivially well separated. If no partition meets the constraints, the one retaining the most taxa is chosen and the manifest says the data does not support that granularity. A more inclusive alternative is suggested when one exists.

Each partition is scored in **two independent spaces**:

| Space | Distance | Available |
| :--- | :--- | :--- |
| ESM-2 | cosine distance between ESM-2 embeddings | with `--embed --embed-model esm2`; the embeddings are computed before clustering, in the same run |
| Structural (`--cluster-mode structural`) | 1 − TM-score from an all-vs-all Foldseek 3Di+AA search | always; computed from the structures |
| Identity (`--cluster-mode sequence`) | 1 − fraction identity from an all-vs-all MMseqs2 search | always; computed from the sequences |

ESM-2 is the primary space when present, because it is independent of the clustering and so checks it rather than grading its own work; otherwise the structural (or identity) space decides. When both are available, the run prints both silhouettes for the chosen partition and, if the second space would have chosen a different operating point, a disagreement notice. The notice is informative, not blocking, and is recorded in `multi_alignment_summary.json` and the manifest.

> On the 1,135-structure Nipah binder set with no embeddings, the structural space selected e-value 1.0, coverage 0.5 (85 clusters, 21 with 10+ taxa, 87% retained, silhouette +0.30). A TM-align–based structural space ranked the same partitions with Spearman ρ = 0.99 but took 131 s against 9 s, so the faster 3Di+AA search is used.

### 5. Recursion (`--cluster-recurse`)
After a cluster is aligned, its **3Di** alignment is checked. If it is ≥ `--recurse-max-gap` (0.65) gaps or more than `--recurse-max-expansion` (3.0) × its median length, the cluster is re-clustered — first at the same e-value (which is effectively stricter on the smaller set), then at 0.3, then 0.1. A split is accepted only if it keeps at least `--recurse-min-kept` (0.70) of the cluster's taxa in clusters of 4+; otherwise the cluster is kept whole and reported as not divisible. Up to `--recurse-max-depth` (2) levels.

These thresholds are deliberately looser than the whole-set gate: 40–55% gaps is ordinary for a diverse family, and applying 50%/1.5× to clusters would split almost every one. Each split, each rejected split, and every cluster left above threshold (with the reason — depth limit, not divisible, too small, recursion off) is recorded. Split children record their `parent`.

`--recurse-on both` also splits clusters whose AA alignment is over threshold. By default only 3Di is judged, and tree-eligible clusters with a gappy AA alignment are flagged instead — with MAFFT the AA alignment is built without structural guidance and is typically a few points gappier.

### Cluster Output
Each cluster is written to `alignment/cluster_<n>_<threshold>/` (e.g. `cluster_3_e0.3`) with its structures and `cluster_info.json`, plus its alignments and `alignment_info.json` when it has 2+ taxa. `cluster_info.json` records its size, 3Di and AA alignment statistics, recursion depth, parent, and whether it is tree-eligible. `multi_alignment_summary.json` lists every cluster plus the QC, operating point, selection and recursion records.

---

## 5. Tree Inference Strategy

### Per-Cluster Trees (IQ-TREE)
Every cluster of **4 or more** taxa gets a tree; smaller clusters are recorded as skipped (IQ-TREE cannot bootstrap fewer than 4 sequences). Clusters under 80 taxa use ModelFinder over the chosen 3Di matrices; 80 and over use a fixed model — the first chosen 3Di matrix with `+G4` (e.g. `Q.3Di.AF+G4`) and `LG+G4` for AA — which keeps large trees tractable and involves no model selection.

### Whole-Set Tree (VeryFastTree)
`--master-method fasttree` (the default) builds the whole-set tree with VeryFastTree (or FastTree), which takes about a minute for 1,000+ taxa where IQ-TREE would take days. The 3Di tree uses the requested 3Di matrix through FastTree's `-trans` option: the PAML matrix is converted to FastTree's format, and the conversion was validated by likelihood — on a fixed tree, IQ-TREE `-m Q.3Di.AF` and FastTree `-trans` agree to 0.003 log-likelihood units. On the Nipah set, Q.3Di.AF improved the whole-set log-likelihood from −201,283 (LG) to −151,272. VeryFastTree runs in double precision, which matches IQ-TREE exactly. The AA tree uses LG.

These trees are named for the method — `<prefix>_3di_fasttree.treefile` and `<prefix>_aa_fasttree.treefile` — and each has an `.info.json` recording program, version, model and command. They are an overview: on a gap-rich whole-set alignment most splits are poorly supported.

### Starting Trees for Large Clusters
With `--guide-trees` (on by default), IQ-TREE runs on alignments of 80+ taxa start from a VeryFastTree tree built under the same model. The starting tree is made safe for IQ-TREE first:
- identical sequences are kept (`-keep-ident`) — otherwise IQ-TREE drops duplicates and rejects a tree that contains them;
- multifurcations (FastTree groups identical sequences) are resolved into random binary splits, seeded for reproducibility — IQ-TREE's search requires a bifurcating tree;
- the tree is checked to have exactly the alignment's taxa, be bifurcating, and carry no labels — otherwise it is not used;
- if the guided IQ-TREE run still fails, it is rerun without the starting tree, so no tree is ever lost.

Each guided run writes `<tree>.start_tree.json` (used / fell back / not used, and why). In benchmarks on 175- and 210-taxon clusters, starting trees did **not** reliably shorten IQ-TREE's runtime — its stopping rules dominate — but on the 210-taxon cluster all three guided runs reached higher-likelihood trees than all three unguided runs. Use `--no-guide-trees` to disable.

### Reading Large-Cluster Trees
Replicate IQ-TREE runs on the same large cluster can disagree on many splits at nearly identical likelihoods: the likelihood surface is flat. Treat deep splits in large clusters with caution unless bootstrap support is high.

---

## 6. Run Manifest, Preflight & Provenance

### Preflight
Before any work, the pipeline checks every tool the chosen options need, the input folder and metadata file, the MAFFT matrix when `--aligner mafft` is used (a check that you are in the repository root), and — with `--embed` — that an interpreter can import `torch` and `transformers`. Missing requirements stop the run within seconds with a list. Non-blocking problems are warnings: embeddings on CPU, or no FastTree (starting trees are then disabled and IQ-TREE runs unguided).

### The Manifest
`<output-dir>/run_manifest.json` and a readable `run_manifest.txt` record:
- the exact command and every resolved setting, defaults included;
- host, Python, package version and git commit (and whether the working tree had uncommitted changes);
- the version of every external tool actually used;
- the input directory, file count and a content digest, plus the Hugging Face dataset and revision when the input comes from a Hugging Face cache;
- QC exclusions, the alignment record (whole-set gate, aligner), the clustering record, and tree provenance;
- every stage with its status (`ok`, `skipped`, `failed`), duration, command and outputs;
- warnings, and the final output inventory.

It is rewritten after every stage, so a long run can be followed in `run_manifest.txt`, and an interrupted run still leaves a record. The overall status is `completed`, `completed_with_skips`, `failed`, or `interrupted` (for a crash or `kill`). A failure in the whole-set tree or the embedding stage is recorded and the run continues, so the per-cluster trees and the dashboard are still produced. The console summary at the end matches the manifest.

### Which Aligner Produced a File
The MAFFT path writes `mafft.fasta_*.fa` and copies them to the `foldmason.fasta_*.fa` names that downstream steps read. Every alignment directory therefore has an `alignment_info.json` recording the aligner, its version, the MAFFT matrix and the real producer of each file.

---

## 7. Custom Metadata Integration

The pipeline supports arbitrary user-supplied metadata across all common formats (`.xlsx`, `.csv`, `.tsv`, `.json`).

```bash
viral-phylo pipeline \
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

## 8. Interactive Web Visualization Suite

Each pipeline run generates an interactive browser application at `<output-dir>/interactive_tree.html`. It is generated output, not part of the repository, and it lists every workflow found under `results/` — so a fresh clone has none until a pipeline has run. To rebuild the dashboards for everything already in `results/`, run `build-alignments && build-tree-view` from the repository root.

### Launching the Viewer
- **Option 1 (Interactive Bridge Server)**:
  ```bash
  python scripts/serve_interactive.py
  ```
  Opens `http://localhost:8000` with direct in-browser pipeline execution and Viro3D querying enabled.
- **Option 2 (Direct Browser Open)**:
  Open `<output-dir>/interactive_tree.html` in Chrome, Firefox, Safari, or Edge.

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

For a pipeline run, the 3Di and AA trees shown are the **whole-set** trees (VeryFastTree by default — see [section 5](#5-tree-inference-strategy)). The per-cluster IQ-TREE trees are written to `phylogeny/cluster_*/`; they are not listed in the dashboard, but any of them can be opened with `viral-phylo report --tree <file>`.

**Branch support:** IQ-TREE trees carry an SH-aLRT/UFboot pair per branch; VeryFastTree/FastTree trees carry a single SH-like local support value. The viewer labels the latter as *single-value support* in the legend and tooltip, so it is never presented as UFboot.

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
- **Dual Alphabet Toggle**: Switch instantly between **🥩 Amino Acid (AA)** and **🧊 3Di Structural Alphabet** modes. Each alphabet is shown with its own column count: with MAFFT the AA and 3Di alignments are built separately and usually differ in length, and each is displayed in full.
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

## 9. Interactive Pipeline Studio Tab

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

## 10. Troubleshooting & Best Practices

### First step: read the manifest
For any run that did not do what you expected, open `<output-dir>/run_manifest.txt`. Every stage is listed as `ok`, `skipped` or `failed` with a reason, and warnings are collected at the end.

### Issue: The preflight check stops the run
- **Cause**: A tool, input file, matrix, or embedding package the chosen options need is missing.
- **Solution**: The message lists each problem. Install missing tools from `environment.yml` (or re-run `./install.sh`); check paths; and run from the repository root so `matrices/` is found.

### Issue: `ERROR: File not found MATRICES/Q.3DI.AF`
- **Cause**: IQ-TREE upper-cases the `-mset` argument, so `matrices/Q.3Di.AF` becomes `MATRICES/Q.3DI.AF`, which does not exist on a case-sensitive filesystem. Before this was handled, every tree using model selection failed with exit status 2.
- **Solution**: Handled automatically — an upper-case copy is written to `MATRICES/` and passed to IQ-TREE. If you see this error, check the run was started from the repository root.

### Issue: `It makes no sense to perform bootstrap with less than 4 sequences`
- **Cause**: IQ-TREE cannot bootstrap fewer than 4 sequences.
- **Solution**: Handled automatically: clusters under 4 taxa are recorded as skipped, and a direct `tree` call on such an alignment drops the bootstrap and still builds the tree.

### Issue: The tree stage is extremely slow on a many-core machine
- **Cause**: `--threads AUTO` makes IQ-TREE benchmark every thread count before inferring anything; on a 256-core machine that took over 12 minutes for a 6-taxon tree, repeated for every cluster.
- **Solution**: The default is now `--threads 8`. Leave it, or set an explicit number.

### Issue: The run uses nearly every core on the machine
- **Cause**: Several tools use every core unless told otherwise: FoldMason and Foldseek default to all 256 cores, and PyTorch opens one thread per physical core (128 here) for CPU embedding. Earlier versions passed `--threads` only to IQ-TREE and the clustering tools.
- **Solution**: `--threads` is now a cap for the whole run, passed to every tool and exported as `OMP_NUM_THREADS` and the BLAS equivalents. With `--threads 4` on the test structures, the run peaked at about 4.6 cores (the extra is Foldseek's coordinating thread). A lower `OMP_NUM_THREADS` already set in your shell is kept. Check `resources.thread_cap` in `run_manifest.json`.

### Issue: Cluster trees are slow even with many threads
- **Cause**: IQ-TREE splits its work across threads by alignment site pattern, and cluster alignments are short. On a 167-taxon Nipah binder cluster (303 patterns), 30 search iterations took 88 s at 2 threads, 60 s at 4, **50 s at 8**, 63 s at 16 and 82 s at 32, with identical results. Beyond about 8 threads, the time goes on keeping threads in step.
- **Solution**: Handled by default: each cluster tree gets `--tree-threads` (8) threads, and the rest of `--threads` builds more clusters at once, largest first. For example, `--threads 32` runs 4 cluster trees at a time. Stopping the pipeline stops every running tree job, including its IQ-TREE process.

### Issue: The whole-set IQ-TREE tree takes days
- **Cause**: IQ-TREE on a large, gap-rich whole-set alignment is extremely expensive (a 1,193-taxon run was projected at ~78 hours).
- **Solution**: Use `--master-method fasttree` (the default): VeryFastTree with the same 3Di matrix, in about a minute.

### Issue: Aligning everything produces a mostly-gap alignment
- **Cause**: The input is too structurally diverse to share one alignment. The whole-set gate reports this.
- **Solution**: Add `--multi-alignment` (and usually `--cluster-recurse`) — see [section 4](#4-handling-divergent-datasets-clustering--recursion).

### Issue: Recursion produced many tiny clusters
- **Cause**: Recursion thresholds too strict for the data, or a very large step to the next threshold.
- **Solution**: Loosen `--recurse-max-gap` / `--recurse-max-expansion`. Splits that keep fewer than `--recurse-min-kept` (70%) of a cluster's taxa in clusters of 4+ are already rejected; raise it to be stricter still.

### Issue: No silhouette profile was produced
- **Cause**: `--embed-clustering nj`. Neighbour-joining builds no linkage matrix to cut, so no silhouette can be computed. The run says so in a notice and a manifest warning.
- **Solution**: Use `upgma` (default), `average`, `complete`, `single`, or `ward`.

### Issue: Embeddings run on CPU although a GPU is available
- **Cause**: The torch build targets a newer CUDA than the installed driver supports.
- **Solution**: `./install.sh` detects and reports this. Install a torch build matching the driver, e.g. `pip install torch --index-url https://download.pytorch.org/whl/cu121`. CPU embedding still works, more slowly.

### Issue: The 3D structure viewer is empty
- **Cause**: No C-alpha coordinates could be read from the input structures.
- **Solution**: Both PDB and mmCIF are supported (mmCIF is parsed from its `_atom_site` header, so any column order works); the manifest's `ca-extraction` stage reports how many structures yielded a backbone. Check the input files contain `CA` atoms.

### Issue: `Ultrafast bootstrap does not work with -fast, -te or -n option`
- **Cause**: In IQ-TREE, `--fast` is mutually exclusive with `-B` (Ultrafast bootstrap).
- **Solution**: Do not combine `--fast` with `-b 1000`. The pipeline resolves this automatically: if bootstrap is requested, `--fast` is omitted; if `--fast` is passed, bootstrap is suppressed with a notice.

### Issue: `Reading alignment file ... ERROR: File not found`
- **Cause**: Passing an invalid or nonexistent relative file path, or a directory instead of a FASTA file.
- **Solution**: Paths are sanitized automatically, and alignment directories resolve to the `foldmason.fasta_3di.fa` (or `mafft.fasta_3di.fa`) file inside them.

### Issue: ModelFinder Takes Very Long on Large Cohorts ($N > 100$)
- **Solution**: For alignments with $\ge 80$ sequences the pipeline skips model selection and pins the first chosen 3Di matrix with `+G4` (e.g. `Q.3Di.AF+G4`) and `LG+G4` for amino acids.

### Best practice: keep results out of git
`results/`, the generated dashboards and `.js` data bundles, `MATRICES/` and `nohup.out` are ignored by git. The 6 benchmark structures used by the tests and examples live in `tests/fixtures/structures/`.

---

## Authors & Citation
- **Viro3D**: University of Glasgow Centre for Virus Research (CVR).
- **FoldMason**: Steinegger Lab (*"FoldMason: fast and accurate multiple protein structure alignment"*, Nature Communications).
- **3Di Empirical Matrices**: Georg Hochberg et al. (*"A general substitution matrix for structural phylogenetics"*, Edmond Dataverse [doi:10.17617/3.1MJJBH](https://doi.org/10.17617/3.1MJJBH)).
- **IQ-TREE 2**: Minh et al. (*"IQ-TREE 2: New Models and Efficient Methods for Phylogenetic Inference"*, Molecular Biology and Evolution).
- **Foldseek**: van Kempen et al. (*"Fast and accurate protein structure search with Foldseek"*, Nature Biotechnology, 2024).
- **MMseqs2**: Steinegger & Söding (*"MMseqs2 enables sensitive protein sequence searching for the analysis of massive data sets"*, Nature Biotechnology, 2017).
- **FastTree 2**: Price, Dehal & Arkin (*"FastTree 2 – approximately maximum-likelihood trees for large alignments"*, PLoS ONE, 2010).
- **VeryFastTree**: Piñeiro, Abuín & Pichel (*"Very Fast Tree: speeding up the estimation of phylogenies for large alignments through parallelization and vectorization strategies"*, Bioinformatics, 2020).
- **MAFFT**: Katoh & Standley (*"MAFFT multiple sequence alignment software version 7"*, Molecular Biology and Evolution, 2013).
