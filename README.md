# 🧬 Viral Structural Phylogenetics & Cophylogenetic Tanglegram Suite

[![CI](https://github.com/kieran12lamb/structural_phylogenetics_tool/actions/workflows/ci.yml/badge.svg)](https://github.com/kieran12lamb/structural_phylogenetics_tool/actions/workflows/ci.yml)
[![Python 3.11+](https://img.shields.io/badge/python-3.11+-blue.svg)](https://www.python.org/downloads/)
[![Conda](https://img.shields.io/badge/conda-bioconda%20%7C%20conda--forge-green.svg)](https://anaconda.org)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

An end-to-end framework for **structural phylogenetics** and **protein language model (PLM) clustering** of highly divergent viral proteins.

By computing evolutionary relationships on **tertiary 3D protein structure coordinates** (encoded into the 20-state **3Di** structural alphabet) alongside **primary amino acid sequences** and **PLM semantic embeddings**, this tool resolves deep evolutionary relationships that are otherwise obscured by sequence saturation and extreme viral divergence.

---

## 🌟 Key Capabilities

1. **Automated Structure Fetching (`fetch`)**: Query and download predicted 3D structures from the **AlphaFold Protein Structure Database (AFDB)** via UniProt IDs or search terms, or from [Viro3D](https://viro3d.org) (AlphaFold2, ColabFold, and ESMFold models).
2. **Structural Multiple Sequence Alignment (`align`)**: Fast 3D backbone superposition with **FoldMason** or **MAFFT** using custom 3Di substitution matrices (`mat3di.out`).
3. **Similarity Clustering for Divergent Datasets**: Before aligning, structures are clustered by measured structural similarity (**Foldseek**) or sequence identity (**MMseqs2**), so each alignment and tree is built from proteins that genuinely align. Clusters whose own alignment is still too gappy can be split further.
4. **Maximum Likelihood Structural Phylogenetics (`tree`)**: Rigorous tree inference via **IQ-TREE** with empirical AlphaFold (`Q.3Di.AF`) and ESMFold (`Q.3Di.LLM`) substitution matrices, automated rate heterogeneity model selection (`+G4`, `+I`, `+R`), and ultrafast bootstrap support. Large clusters start from a validated **VeryFastTree** starting tree.
5. **Fast Whole-Set Trees**: An overview tree of the entire input with **VeryFastTree**, using the same 3Di matrix as IQ-TREE — about a minute for 1,000+ taxa rather than days.
6. **PLM Embedding Trees & Clustering (`embed`)**: Extract sequence representations from **ESM-2 (650M)** or **ESM-C (600M)**, compute semantic distance matrices (Cosine, Euclidean, L1/Manhattan), and build hierarchical trees with silhouette profile analysis.
7. **Run Manifest & Provenance**: Every pipeline run records its exact command, resolved settings, tool versions, input fingerprint and the outcome of every stage — including anything skipped or failed, and why.
8. **Interactive Visualization Suite (`interactive_tree.html`)**:
   - **5 Tree Layouts**: Rectangular Phylogram, Cladogram, Radial Tree, Equal-Angle Unrooted Tree, and Cophylogenetic Tanglegram.
   - **Dual Cophylogenetic Tanglegram**: Direct comparison of 3Di Tertiary Structure vs. Primary AA Sequence vs. ESM-2 PLM Distance Trees, complete with Robinson-Foulds and Cophenetic congruence metrics.
   - **Interactive 3D C$\alpha$ Backbone Viewer**: Live rotation and pLDDT confidence coloring directly on hover/click.
   - **Integrated Multiple Sequence Alignment Drawer**: Virtualized MSA rendering with real-time gap-stripping and leaf synchronization.
   - **Multi-Level Clade Scoping & Silhouette Analysis**: Automatic suggestion of optimal tree cut sizes based on silhouette metrics.
9. **📦 Subclade Package Exporter (.zip)**: One-click export of any filtered clade or subset into a standalone, reproducible zip bundle containing alignments, trees, embeddings, structures, an offline interactive HTML viewer, and executable `REPRODUCE.sh` scripts.

---

## 🚀 Quick Start & Installation

### Option 1: One-Click Installer (Recommended)

Clone the repository and run the automated installation script:

```bash
git clone https://github.com/kieran12lamb/structural_phylogenetics_tool.git
cd structural_phylogenetics_tool

# Automated environment setup and test validation:
./install.sh

# Activate the environment:
conda activate spt
```

The installer creates (or updates) the `spt` conda environment, installs the package, checks every external tool and Python package, reports whether PyTorch can use your GPU, and runs the test suite.

### Option 2: Manual Conda Setup

```bash
conda env create -f environment.yml
conda activate spt
pip install -e . --no-deps
```

### Option 3: Editable Python Package Install

```bash
pip install -e .          # core pipeline
pip install -e ".[ml]"    # plus torch/transformers for --embed
```

pip installs the Python packages only. The external tools — `foldmason`, `iqtree`, `mafft`, `mmseqs2`, `foldseek`, and `veryfasttree` (or `fasttree`) — are not on PyPI; install them from bioconda (they are all listed in `environment.yml`).

> **GPU note:** if embeddings run on CPU despite a GPU being present, your torch build probably targets a newer CUDA than your driver supports. `install.sh` detects this and prints the fix.

---

## 💻 CLI Usage

The unified command-line tool `viral-phylo` (equivalently `python3 scripts/viral_phylogenetics.py`) provides modular subcommands. **Run it from the repository root**: the 3Di matrices in `matrices/` are resolved relative to the working directory. The pipeline checks this, and every other requirement, in a preflight step before doing any work.

### 1. Run the Full End-to-End Pipeline

#### 🧭 Recommended for Large or Divergent Datasets

For hundreds to thousands of structures that do not all share a fold — the case where aligning everything together produces a mostly-gap alignment:

```bash
viral-phylo pipeline \
  --input-folder path/to/structures \
  --output-dir results/my_run \
  --metadata path/to/metadata.csv \
  --multi-alignment --cluster-recurse \
  --tree-type both --matrix alphafold \
  --embed --embed-clustering upgma
```

This removes short fragments, clusters the structures with Foldseek, splits any cluster whose alignment is still too gappy, builds an IQ-TREE tree for every cluster of 4+ taxa, a VeryFastTree tree of the whole set, and ESM-2 embeddings, then writes the interactive dashboard and a run manifest. See [How divergent datasets are handled](#-how-divergent-datasets-are-handled) below.

`--threads` (default 8) caps the CPU use of the whole run: stages run one at a time and every tool — FoldMason, MAFFT, Foldseek, IQ-TREE, VeryFastTree, PyTorch — is limited to that many threads. `--threads AUTO` removes the cap. Cluster trees run in parallel within the cap: each gets `--tree-threads` (default 8) threads, so `--threads 32` builds 4 at a time. IQ-TREE gets little from more threads on short cluster alignments.

The pipeline can also be executed in either **Fast Mode** (for rapid exploratory screening) or **Thorough / Slow Mode** (for publication-grade rigorous phylogenetics):

#### ⚡ Fast Analysis Mode (Rapid Screening & Prototyping)
> Uses `--fast` heuristic ML tree search, skips exhaustive model testing by applying verified empirical models (`+G4`), and omits bootstrapping to complete in seconds to a couple minutes.

```bash
# Fast analysis from a local directory of PDB structures:
python3 scripts/viral_phylogenetics.py pipeline \
  --input-folder tests/fixtures/structures \
  --tree-type both \
  --matrix alphafold \
  --fast \
  --embed \
  --output-dir results/fast_workflow

# Fast analysis querying Viro3D directly (e.g. 50 glycoproteins):
python3 scripts/viral_phylogenetics.py pipeline \
  --qualifier glycoprotein \
  --count 50 \
  --tree-type both \
  --matrix auto \
  --fast \
  --embed \
  --output-dir results/fast_glycoprotein_workflow
```

#### 🔬 Thorough / Slow Analysis Mode (Publication-Grade Deep Inference)
> Runs exhaustive ModelFinder rate heterogeneity testing across empirical matrices, infers Maximum Likelihood trees with extensive branch swapping, computes 1,000 Ultrafast Bootstrap replicates (`-b 1000`) and 1,000 SH-aLRT support tests (`--alrt 1000`), evaluates both 3Di matrices (`--matrix both`), and extracts ESM-2 PLM representations (`--embed`).

```bash
# Thorough analysis from a local directory of PDB structures:
python3 scripts/viral_phylogenetics.py pipeline \
  --input-folder tests/fixtures/structures \
  --tree-type both \
  --matrix both \
  --rate-heterogeneity auto \
  --bootstrap 1000 \
  --alrt 1000 \
  --embed \
  --output-dir results/thorough_workflow

# Thorough analysis querying Viro3D directly (e.g. 100 glycoproteins):
python3 scripts/viral_phylogenetics.py pipeline \
  --qualifier glycoprotein \
  --count 100 \
  --tree-type both \
  --matrix auto \
  --rate-heterogeneity auto \
  --bootstrap 1000 \
  --alrt 1000 \
  --embed \
  --output-dir results/thorough_glycoprotein_workflow
```

### 2. Step-by-Step Execution

#### Step A: Download Viral Structures
```bash
python3 scripts/viral_phylogenetics.py fetch \
  --qualifier glycoprotein \
  --max-sequences 50 \
  --output-dir results/structures
```

#### Step B: Multiple Structural Alignment (FoldMason / MAFFT)
```bash
python3 scripts/viral_phylogenetics.py align \
  --folder results/structures \
  --output-dir results/alignments \
  --aligner foldmason
```

#### Step C: Infer Maximum Likelihood Phylogenies (IQ-TREE)
```bash
python3 scripts/viral_phylogenetics.py tree \
  --alignment results/alignments/foldmason.fasta_3di.fa \
  --alignment-aa results/alignments/foldmason.fasta_aa.fa \
  --tree-type both \
  --matrix alphafold \
  --rate-heterogeneity auto \
  --bootstrap 1000 \
  --output-dir results/phylogeny \
  --prefix viral_tree
```

Use `--method fasttree` for a fast approximate-ML tree instead (the 3Di matrix is still applied).

#### Step D: Extract PLM Embeddings & Build Clustering Tree
```bash
viral-phylo embed \
  --input results/alignments/foldmason.fasta_aa.fa \
  --model esm2 \
  --clustering upgma \
  --metric cosine \
  --output-dir results/embeddings \
  --prefix viral_esm2
```

`upgma`, `average`, `complete`, `single` and `ward` also produce a silhouette profile; `nj` (neighbour-joining) does not, because it builds no linkage matrix to cut.

---

## 🧩 How Divergent Datasets Are Handled

With `--multi-alignment`, the pipeline does not force every structure into one alignment:

1. **Quality control** — sequences shorter than 50 residues, or than 40% of the median length, are excluded (and listed in the manifest). Fragments are what stretch a global alignment into a staircase of gaps.
2. **Whole-set gate** — the full 3Di alignment is still built, and judged: if it is ≥50% gaps or more than 1.5× the median sequence length, it is reported as not usable as a single alignment.
3. **Structural clustering** — Foldseek clusters the structures (3Di+AA alignment, e-value 1.0, 50% coverage by default). Use `--cluster-mode sequence` for MMseqs2 amino-acid clustering, or `--cluster-sweep` to try a grid of thresholds and pick the partition with the best silhouette that still keeps enough of the data in usable clusters. Silhouettes are scored in two spaces — ESM-2 embeddings (with `--embed`, computed before clustering in the same run) and a structural 1 − TM-score space computed from the structures — and the run flags when the two disagree.
4. **Recursion** (`--cluster-recurse`) — any cluster whose own 3Di alignment is still too gappy is re-clustered at a stricter threshold, in small steps. A split is only accepted if it keeps most of its taxa in clusters of 4+, so clusters are refined rather than shattered into singletons.
5. **Trees** — IQ-TREE builds a tree for every cluster of 4+ taxa; clusters of 80+ start from a repaired, validated VeryFastTree tree (`--no-guide-trees` to disable). The whole set gets a VeryFastTree overview tree.

On a 1,193-structure Nipah binder set, where the global alignment was 91.6% gaps, structural clustering gave per-cluster 3Di alignments that were typically 40–55% gaps.

---

## 🧾 Outputs & Provenance

Each pipeline run writes to its `--output-dir`:

```
results/my_run/
├── run_manifest.json / .txt     # command, settings, tool versions, input fingerprint, per-stage status
├── interactive_tree.html        # interactive dashboard for this run
├── taxa_metadata.json           # standardized metadata
├── multi_alignment_summary.json # clusters, QC, selection and recursion record (--multi-alignment)
├── alignment/
│   ├── foldmason.fasta_3di.fa   # whole-set alignment (3Di and _aa.fa)
│   ├── alignment_info.json      # which aligner really produced the files
│   └── cluster_<n>_<threshold>/ # one directory per cluster (aligned clusters have their own alignment_info.json)
└── phylogeny/
    ├── *_3di_fasttree.treefile  # whole-set trees (+ .info.json: program, version, model, command)
    └── cluster_<n>_<threshold>/ # IQ-TREE output per cluster (+ .start_tree.json when guided)
```

The manifest is rewritten after every stage, so an interrupted run still leaves a record, and it ends with an explicit status: `completed`, `completed_with_skips`, `failed`, or `interrupted`. Alignment files keep their historical `foldmason.fasta_*.fa` names even when MAFFT produced them — `alignment_info.json` records the real producer.

---

## 🖥️ Interactive Web Visualization Suite

Each pipeline run writes its dashboard to `<output-dir>/interactive_tree.html`; open it in any web browser:

```bash
open results/my_run/interactive_tree.html
```

The dashboard is generated output and is not part of the repository. It lists every workflow found under `results/`, so a fresh clone has none until you run a pipeline. To rebuild the dashboards for everything already in `results/`:

```bash
build-alignments && build-tree-view
```

Or start the optional local bridge server for running pipelines directly from the web interface:

```bash
python3 scripts/serve_interactive.py --port 8000
```

### Visualising Your Own Tree
Any Newick or NEXUS tree, with an optional metadata table (`.xlsx`, `.csv`, `.tsv`, `.json`), can be turned into a standalone report:

```bash
viral-phylo report --tree tree.treefile --metadata samples.csv --output tree_report.html
```

Inside any report you can also use **📂 Load tree file** and **🏷️ Add metadata** (under *Dataset Cohort*) to open a tree and its metadata directly in the browser.

### Branch Support Values
IQ-TREE trees carry an SH-aLRT/UFboot pair per branch; VeryFastTree/FastTree trees carry a single SH-like local support value. The viewer labels the two differently, so a FastTree value is never presented as UFboot.

---

## 📦 Subclade Package Export (.zip)

When exploring a subclade or applying metadata filters in `interactive_tree.html`, click **"📦 Export Subclade ZIP"** to download a self-contained archive:

```
subclade_clade_1_24_taxa_package.zip
├── alignments/
│   ├── subclade_aa.fasta               # Filtered primary sequence alignment
│   └── subclade_3di.fasta              # Filtered structural 3Di alignment
├── trees/
│   ├── subclade_3di.nwk                # Pruned 3Di structural phylogeny
│   ├── subclade_aa.nwk                 # Pruned AA primary sequence phylogeny
│   └── subclade_esm2_cosine.nwk        # Pruned ESM-2 distance tree
├── embeddings/
│   └── subclade_summary.json           # Silhouette and cluster metrics
├── structures/
│   └── subclade_ca_traces.json         # 3D C-alpha coordinates and pLDDT
├── metadata/
│   ├── subclade_metadata.tsv           # TSV metadata table
│   └── subclade_metadata.json          # JSON metadata records
├── subclade_viewer.html                # Offline standalone interactive 3D viewer
├── REPRODUCE.sh                        # Executable reproduction script
└── REPRODUCIBILITY.md                  # Parameters and provenance documentation
```

---

## 📁 Repository Structure

```
.
├── .github/workflows/ci.yml           # GitHub Actions automated CI testing
├── environment.yml                    # Conda environment definition
├── install.sh                         # Automated installation script
├── pyproject.toml                     # Python packaging and CLI entrypoints
├── USAGE_GUIDE.md                     # Comprehensive user and expert guide
├── SKILL.md                           # AI Agent skill specification
├── matrices/                          # 3Di structural substitution matrices
│   ├── mat3di.out                     # MAFFT 3Di substitution matrix
│   ├── Q.3Di.AF                       # IQ-TREE AlphaFold 3Di empirical matrix
│   └── Q.3Di.LLM                      # IQ-TREE ESMFold 3Di empirical matrix
├── viral_phylo/                       # The package
│   ├── cli.py                         # `viral-phylo` subcommands and the pipeline
│   ├── alignment.py                   # FoldMason / MAFFT alignment
│   ├── clustering.py                  # QC, Foldseek/MMseqs2 clustering, recursion
│   ├── tree.py                        # IQ-TREE and VeryFastTree/FastTree trees
│   ├── start_trees.py                 # Validated FastTree starting trees for IQ-TREE
│   ├── matrices.py                    # Matrix download and format conversion
│   ├── manifest.py                    # Run manifest and provenance
│   ├── structures.py                  # C-alpha extraction from PDB and mmCIF
│   ├── embeddings.py                  # ESM-2 / ESM-C embeddings and PLM trees
│   └── web/                           # Dashboard builder and viewer templates
├── scripts/                           # Backward-compatible script entrypoints
├── tests/
│   ├── fixtures/structures/           # 6 benchmark structures used by tests and examples
│   └── test_*.py                      # Automated test suite
├── examples/                          # Example tree reports
└── results/                           # Your pipeline runs (generated; not tracked by git)
```

---

## 🧪 Testing

Run the automated unit test suite:

```bash
# Using pytest:
pytest tests/ -v

# Or using standard python unittest:
python3 -m unittest discover -s tests -p "test_*.py" -v
```

Integration tests that call IQ-TREE or VeryFastTree/FastTree run when those tools are installed and are skipped otherwise.

---

## 📚 References & Citation

- **FoldMason**: Gilchrist, C.L.M., Mirdita, M., & Steinegger, M. (2024). *Simultaneous identification of structural and sequence variation across thousands of proteins*. Bioinformatics.
- **3Di Empirical Substitution Matrices**: Georg Hochberg et al. (2024). *A general substitution matrix for structural phylogenetics*. Edmond Dataverse, doi:10.17617/3.1MJJBH.
- **IQ-TREE 2**: Minh, B.Q. et al. (2020). *IQ-TREE 2: New models and efficient methods for phylogenetic inference in the genomic era*. Molecular Biology and Evolution.
- **Foldseek**: van Kempen, M. et al. (2024). *Fast and accurate protein structure search with Foldseek*. Nature Biotechnology.
- **MMseqs2**: Steinegger, M. & Söding, J. (2017). *MMseqs2 enables sensitive protein sequence searching for the analysis of massive data sets*. Nature Biotechnology.
- **FastTree 2**: Price, M.N., Dehal, P.S. & Arkin, A.P. (2010). *FastTree 2 – approximately maximum-likelihood trees for large alignments*. PLoS ONE.
- **VeryFastTree**: Piñeiro, C., Abuín, J.M. & Pichel, J.C. (2020). *Very Fast Tree: speeding up the estimation of phylogenies for large alignments through parallelization and vectorization strategies*. Bioinformatics.
- **MAFFT**: Katoh, K. & Standley, D.M. (2013). *MAFFT multiple sequence alignment software version 7*. Molecular Biology and Evolution.
- **ESM-2**: Lin, Z. et al. (2023). *Evolutionary-scale prediction of atomic-level protein structure with a language model*. Science.
- **Viro3D**: Steinegger Lab (2024). *AI-predicted viral structural repository*.
