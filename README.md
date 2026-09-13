# 🧬 Viral Structural Phylogenetics & Cophylogenetic Tanglegram Suite

[![CI](https://github.com/kieran12lamb/structural_phylogenetics_tool/actions/workflows/ci.yml/badge.svg)](https://github.com/kieran12lamb/structural_phylogenetics_tool/actions/workflows/ci.yml)
[![Python 3.11+](https://img.shields.io/badge/python-3.11+-blue.svg)](https://www.python.org/downloads/)
[![Conda](https://img.shields.io/badge/conda-bioconda%20%7C%20conda--forge-green.svg)](https://anaconda.org)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

An end-to-end framework for **structural phylogenetics** and **protein language model (PLM) clustering** of highly divergent viral proteins.

By computing evolutionary relationships on **tertiary 3D protein structure coordinates** (encoded into the 20-state **3Di** structural alphabet) alongside **primary amino acid sequences** and **PLM semantic embeddings**, this tool resolves deep evolutionary relationships that are otherwise obscured by sequence saturation and extreme viral divergence.

---

## 🌟 Key Capabilities

1. **Automated Structure Fetching (`fetch`)**: Query and download predicted 3D structures from [Viro3D](https://viro3d.org) (AlphaFold2, ColabFold, and ESMFold models).
2. **Structural Multiple Sequence Alignment (`align`)**: Fast 3D backbone superposition with **FoldMason** or **MAFFT** using custom 3Di substitution matrices (`mat3di.out`).
3. **Maximum Likelihood Structural Phylogenetics (`tree`)**: Rigorous tree inference via **IQ-TREE** with empirical AlphaFold (`Q.3Di.AF`) and ESMFold (`Q.3Di.LLM`) substitution matrices, automated rate heterogeneity model selection (`+G4`, `+I`, `+R`), and ultrafast bootstrap support.
4. **PLM Embedding Trees & Clustering (`embed`)**: Extract sequence representations from **ESM-2 (650M)** or **ESM-C (600M)**, compute semantic distance matrices (Cosine, Euclidean, L1/Manhattan), and build UPGMA hierarchical trees with silhouette profile analysis.
5. **Interactive Visualization Suite (`interactive_tree.html`)**:
   - **5 Tree Layouts**: Rectangular Phylogram, Cladogram, Radial Tree, Equal-Angle Unrooted Tree, and Cophylogenetic Tanglegram.
   - **Dual Cophylogenetic Tanglegram**: Direct comparison of 3Di Tertiary Structure vs. Primary AA Sequence vs. ESM-2 PLM Distance Trees, complete with Robinson-Foulds and Cophenetic congruence metrics.
   - **Interactive 3D C$\alpha$ Backbone Viewer**: Live rotation and pLDDT confidence coloring directly on hover/click.
   - **Integrated Multiple Sequence Alignment Drawer**: Virtualized MSA rendering with real-time gap-stripping and leaf synchronization.
   - **Multi-Level Clade Scoping & Silhouette Analysis**: Automatic suggestion of optimal tree cut sizes based on silhouette metrics.
6. **📦 Subclade Package Exporter (.zip)**: One-click export of any filtered clade or subset into a standalone, reproducible zip bundle containing alignments, trees, embeddings, structures, an offline interactive HTML viewer, and executable `REPRODUCE.sh` scripts.

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

### Option 2: Manual Conda Setup

```bash
conda env create -f environment.yml
conda activate spt
```

### Option 3: Editable Python Package Install

```bash
pip install -e .
```

---

## 💻 CLI Usage

The unified command-line tool `scripts/viral_phylogenetics.py` provides modular subcommands:

### 1. Run the Full End-to-End Pipeline

```bash
# From a local directory of PDB structures (e.g. the 6 benchmark example structures):
python3 scripts/viral_phylogenetics.py pipeline \
  --input-folder results/glycoprotein_workflow/structures \
  --tree-type both \
  --matrix alphafold \
  --embed \
  --output-dir results/my_workflow

# Or query directly from Viro3D:
python3 scripts/viral_phylogenetics.py pipeline \
  --qualifier glycoprotein \
  --count 100 \
  --tree-type both \
  --matrix auto \
  --output-dir results/glycoprotein_100_workflow
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

#### Step D: Extract PLM Embeddings & Build Clustering Tree
```bash
python3 scripts/embed_and_cluster.py \
  --fasta results/alignments/foldmason.fasta_aa.fa \
  --model esm2 \
  --metric cosine \
  --output-dir results/embeddings \
  --prefix viral_esm2
```

---

## 🖥️ Interactive Web Visualization Suite

To explore the precomputed cohorts or your own custom pipeline runs, open `interactive_tree.html` in any web browser:

```bash
open interactive_tree.html
```

Or start the optional local bridge server for running pipelines directly from the web interface:

```bash
python3 scripts/serve_interactive.py --port 8000
```

### Included Cohorts:
- **500 Viral Glycoproteins**: Diverse cross-family viral cohort spanning 22 viral families.
- **1,193 Nipah Virus Glycoproteins**: Comprehensive ESMFold mutant landscape.
- **6 Benchmark Glycoproteins**: Validated structural benchmark set.

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
├── .gitignore                         # Comprehensive Git ignore rules
├── environment.yml                    # Conda environment definition
├── install.sh                         # Automated installation script
├── pyproject.toml                     # Modern Python packaging configuration
├── interactive_tree.html              # Standalone interactive visualization suite
├── USAGE_GUIDE.md                     # Comprehensive user and expert guide
├── SKILL.md                           # AI Agent skill specification
├── matrices/                          # 3Di structural substitution matrices
│   ├── mat3di.out                     # MAFFT 3Di substitution matrix
│   ├── Q.3Di.AF                       # IQ-TREE AlphaFold 3Di empirical matrix
│   └── Q.3Di.LLM                      # IQ-TREE ESMFold 3Di empirical matrix
├── scripts/
│   ├── viral_phylogenetics.py         # Unified CLI pipeline tool
│   ├── embed_and_cluster.py           # ESM-2 / ESM-C PLM tree builder
│   ├── metadata_handler.py            # Universal metadata schema analyzer
│   ├── build_dynamic_interactive_tree.py # Standalone HTML compiler
│   ├── build_alignments_data.py       # Alignment data compiler
│   └── serve_interactive.py           # Local bridge web server
├── tests/                             # Automated unit test suite
│   ├── test_alignment.py              # FASTA, 3Di alphabet, and gap stripping
│   ├── test_clustering.py             # Distance metrics, UPGMA, and silhouettes
│   ├── test_metadata.py               # Metadata parsing and color classification
│   ├── test_cli.py                    # CLI argument parsing and matrix discovery
│   └── test_interactive_tree.py       # HTML compilation and asset verification
└── results/                           # Analysis cohorts, trees, and alignments
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

---

## 📚 References & Citation

- **FoldMason**: Gilchrist, C.L.M., Mirdita, M., & Steinegger, M. (2024). *Simultaneous identification of structural and sequence variation across thousands of proteins*. Bioinformatics.
- **3Di Empirical Substitution Matrices**: Georg Hochberg et al. (2024). *A general substitution matrix for structural phylogenetics*. Edmond Dataverse, doi:10.17617/3.1MJJBH.
- **IQ-TREE 2**: Minh, B.Q. et al. (2020). *IQ-TREE 2: New models and efficient methods for phylogenetic inference in the genomic era*. Molecular Biology and Evolution.
- **ESM-2**: Lin, Z. et al. (2023). *Evolutionary-scale prediction of atomic-level protein structure with a language model*. Science.
- **Viro3D**: Steinegger Lab (2024). *AI-predicted viral structural repository*.
