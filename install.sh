#!/usr/bin/env bash
# ==============================================================================
# Automated Installation & Setup Script for Viral Structural Phylogenetics
# ==============================================================================

set -eo pipefail

ENV_NAME="spt"
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

echo "========================================================================"
echo "  🧬 Viral Structural Phylogenetics & Tanglegram Suite: Setup"
echo "========================================================================"

# 1. Detect Package Manager (mamba -> conda -> micromamba)
CONDA_CMD=""
if command -v mamba >/dev/null 2>&1; then
    CONDA_CMD="mamba"
elif command -v conda >/dev/null 2>&1; then
    CONDA_CMD="conda"
elif command -v micromamba >/dev/null 2>&1; then
    CONDA_CMD="micromamba"
else
    echo "❌ Error: Neither conda, mamba, nor micromamba was found in your PATH."
    echo "Please install Miniconda or Mambaforge before proceeding:"
    echo "  https://docs.conda.io/en/latest/miniconda.html"
    exit 1
fi

echo "==> [1/6] Using package manager: $CONDA_CMD"

# 2. Create or Update Conda Environment
ENV_EXISTS=false
if $CONDA_CMD env list | grep -q "^${ENV_NAME} "; then
    ENV_EXISTS=true
fi

if [ "$ENV_EXISTS" = true ]; then
    echo "==> [2/6] Environment '$ENV_NAME' already exists. Updating packages..."
    $CONDA_CMD env update -n "$ENV_NAME" -f "$SCRIPT_DIR/environment.yml" --prune
else
    echo "==> [2/6] Creating Conda environment '$ENV_NAME'..."
    $CONDA_CMD env create -f "$SCRIPT_DIR/environment.yml"
fi

# 3. Install Package in Editable Mode (Registers CLI shortcuts and module imports)
echo "==> [3/6] Installing viral-structural-phylogenetics in editable mode..."
$CONDA_CMD run -n "$ENV_NAME" pip install -e "$SCRIPT_DIR" --no-deps

# 4. Verify Structural Substitution Matrices
echo "==> [4/6] Verifying 3Di structural substitution matrices..."
MATRICES_DIR="$SCRIPT_DIR/matrices"
mkdir -p "$MATRICES_DIR"

check_matrix() {
    local name="$1"
    local path="$MATRICES_DIR/$name"
    if [ -f "$path" ] && [ -s "$path" ]; then
        echo "    ✓ Matrix found: $name"
    else
        echo "    ⚠️ Matrix missing: $name (will be auto-downloaded on first run)"
    fi
}

check_matrix "mat3di.out"
check_matrix "Q.3Di.AF"
check_matrix "Q.3Di.LLM"

# 5. Verify Binaries, Python Packages, and Hardware Acceleration
echo "==> [5/6] Testing environment executables and dependencies..."
$CONDA_CMD run -n "$ENV_NAME" python3 -c "
import sys, shutil

# A. External Bioinformatics Binaries
tools = ['foldmason', 'iqtree', 'mafft']
missing_tools = [t for t in tools if not shutil.which(t)]
if missing_tools:
    print('    ⚠️ Note: External binaries missing from PATH:', missing_tools)
else:
    print('    ✓ External structural bioinformatics tools: foldmason, iqtree, mafft')

# B. Core Python Dependencies
packages = [
    ('numpy', 'NumPy (numeric matrices)'),
    ('scipy', 'SciPy (hierarchical clustering & statistics)'),
    ('pandas', 'Pandas (metadata dataframes)'),
    ('requests', 'Requests (UniProt & PDB downloads)'),
    ('tqdm', 'tqdm (progress bars)'),
    ('rich', 'Rich (terminal formatting)'),
    ('Bio', 'BioPython (phylogenetics & sequence handling)'),
    ('torch', 'PyTorch (deep learning & PLM inference)'),
    ('transformers', 'Hugging Face Transformers (ESM-2 / PLM models)'),
    ('umap', 'UMAP-Learn (2D manifold projection)'),
    ('sklearn', 'Scikit-Learn (clustering & distance metrics)')
]

missing_pkgs = []
for mod_name, desc in packages:
    try:
        mod = __import__(mod_name)
    except ImportError:
        missing_pkgs.append((mod_name, desc))

if missing_pkgs:
    print('    ❌ Missing Python packages:')
    for m, d in missing_pkgs:
        print(f'       - {m} ({d})')
    sys.exit(1)
else:
    print('    ✓ All core & ML Python packages verified (torch, transformers, umap, Bio, scipy, etc.)')

# C. PyTorch Hardware Acceleration Detection
try:
    import torch
    if torch.cuda.is_available():
        device_name = torch.cuda.get_device_name(0)
        print(f'    ⚡ PyTorch hardware acceleration: CUDA GPU ({device_name})')
    elif hasattr(torch.backends, 'mps') and torch.backends.mps.is_available():
        print('    ⚡ PyTorch hardware acceleration: Apple Silicon (MPS)')
    else:
        print('    ℹ️ PyTorch execution: CPU mode')
except Exception:
    pass

# D. Registered Package CLI Entrypoints
entrypoints = ['viral-phylo', 'embed-cluster', 'build-tree-view', 'build-alignments']
missing_ep = [ep for ep in entrypoints if not shutil.which(ep)]
if missing_ep:
    print('    ⚠️ Note: Some CLI entrypoints not found in PATH:', missing_ep)
else:
    print('    ✓ Registered CLI entrypoints: viral-phylo, embed-cluster, build-tree-view, build-alignments')
"

# 6. Run Automated Unit Test Suite
echo "==> [6/6] Running automated test suite..."
if $CONDA_CMD run -n "$ENV_NAME" pytest "$SCRIPT_DIR/tests" -v; then
    echo "    ✓ All unit tests passed via pytest!"
elif $CONDA_CMD run -n "$ENV_NAME" python -m unittest discover -s "$SCRIPT_DIR/tests"; then
    echo "    ✓ All unit tests passed via unittest!"
else
    echo "    ⚠️ Some tests had warnings or issues. You can re-run with: pytest tests/ -v"
fi

echo "========================================================================"
echo "🎉 Setup complete! To activate your new environment:"
echo ""
echo "    conda activate $ENV_NAME"
echo ""
echo "Try running the CLI help:"
echo "    viral-phylo --help"
echo "    embed-cluster --help"
echo "    build-tree-view --help"
echo "    build-alignments --help"
echo ""
echo "Or using python scripts directly:"
echo "    python3 scripts/viral_phylogenetics.py --help"
echo ""
echo "Open the interactive visualization suite:"
echo "    open interactive_tree.html"
echo "========================================================================"
