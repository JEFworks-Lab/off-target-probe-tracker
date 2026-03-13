#!/usr/bin/env bash
set -euo pipefail

REPO_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

echo "=== OPT + Streamlit GUI Setup ==="
echo ""

# Detect platform
OS="$(uname -s)"

# Step 1: Create conda environment (remove existing if present)
echo "[1/4] Creating conda environment 'opt' from environment.yml..."
if conda env list | grep -q "^opt "; then
    echo "      Removing existing 'opt' environment..."
    conda env remove -n opt -y
fi
conda env create -f "${REPO_DIR}/environment.yml"
echo "      Done."

# Step 2: Activate environment
echo "[2/4] Activating environment..."
eval "$(conda shell.bash hook)"
conda activate opt

# Step 3: Install mummer4 (platform-dependent)
echo "[3/4] Installing mummer4..."
if [ "$OS" = "Linux" ]; then
    conda install -n opt -c bioconda "mummer4>=4.0.1" -y
    echo "      mummer4 installed via conda."
elif [ "$OS" = "Darwin" ]; then
    echo "      macOS detected — mummer4 must be installed via Homebrew (system-wide)."
    echo "      Ensure Homebrew is installed: https://brew.sh"
    brew install autoconf automake libtool md5sha1sum || true
    gem install yaggo || true
    brew install mummer
    echo "      mummer installed via Homebrew."
else
    echo "      WARNING: Unsupported OS '$OS'. Install mummer4 manually."
    echo "      See: https://github.com/mummer4/mummer/blob/master/INSTALL.md"
fi

# Step 4: Install the opt Python package
echo "[4/4] Installing the opt Python package..."
pip install "${REPO_DIR}"
echo "      Done."

# Step 5: Decompress bundled reference data files
echo "[5/5] Decompressing bundled reference data (.gz files)..."
find "${REPO_DIR}/data" -name "*.gz" | while read -r gz_file; do
    out_file="${gz_file%.gz}"
    if [ ! -f "$out_file" ]; then
        echo "      Decompressing $(basename "$gz_file")..."
        file_type="$(file -b "$gz_file")"
        if echo "$file_type" | grep -qi "gzip"; then
            gunzip -k "$gz_file"
        elif echo "$file_type" | grep -qi "zip"; then
            unzip -p "$gz_file" > "$out_file"
        else
            echo "      WARNING: Unknown format for $(basename "$gz_file"), skipping."
        fi
    else
        echo "      Already decompressed: $(basename "$out_file")"
    fi
done
echo "      Done."

echo ""
echo "=== Setup complete ==="
echo ""
echo "To launch the Streamlit GUI:"
echo ""
echo "  conda activate opt"
echo "  streamlit run ${REPO_DIR}/app.py"
echo ""
echo "To host on a shared lab server (accessible from any browser on the network):"
echo ""
echo "  conda activate opt"
echo "  streamlit run ${REPO_DIR}/app.py --server.port 8501 --server.address 0.0.0.0"
echo ""
echo "To raise the file upload size limit (for large reference files):"
echo "  Add the following to ${REPO_DIR}/.streamlit/config.toml"
echo "  [server]"
echo "  maxUploadSize = 1000"
echo ""
