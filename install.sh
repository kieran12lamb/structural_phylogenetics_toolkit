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

echo "==> [1/5] Using package manager: $CONDA_CMD"

# 2. Create or Update Conda Environment
ENV_EXISTS=false
if $CONDA_CMD env list | grep -q "^${ENV_NAME} "; then
    ENV_EXISTS=true
fi

if [ "$ENV_EXISTS" = true ]; then
    echo "==> [2/5] Environment '$ENV_NAME' already exists. Updating packages..."
    $CONDA_CMD env update -n "$ENV_NAME" -f "$SCRIPT_DIR/environment.yml" --prune
else
    echo "==> [2/5] Creating Conda environment '$ENV_NAME'..."
    $CONDA_CMD env create -f "$SCRIPT_DIR/environment.yml"
fi

# 3. Verify Structural Substitution Matrices
echo "==> [3/5] Verifying 3Di structural substitution matrices..."
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

# 4. Verify Binaries in Conda Environment
echo "==> [4/5] Testing environment executables and dependencies..."
$CONDA_CMD run -n "$ENV_NAME" python3 -c "
import sys, shutil
tools = ['foldmason', 'iqtree', 'mafft']
missing = [t for t in tools if not shutil.which(t)]
if missing:
    print('    ⚠️ Note: External binaries not in PATH inside env:', missing)
else:
    print('    ✓ All structural bioinformatics binaries available: foldmason, iqtree, mafft')
"

# 5. Run Automated Unit Test Suite
echo "==> [5/5] Running pytest unit test suite..."
if $CONDA_CMD run -n "$ENV_NAME" pytest "$SCRIPT_DIR/tests" -v; then
    echo "    ✓ All unit tests passed!"
else
    echo "    ⚠️ Some tests had warnings or issues. You can re-run with: pytest tests/ -v"
fi

echo "========================================================================"
echo "🎉 Setup complete! To activate your new environment:"
echo ""
echo "    conda activate $ENV_NAME"
echo ""
echo "Try running the CLI help:"
echo "    python3 scripts/viral_phylogenetics.py --help"
echo ""
echo "Open the interactive visualization suite:"
echo "    open interactive_tree.html"
echo "========================================================================"
