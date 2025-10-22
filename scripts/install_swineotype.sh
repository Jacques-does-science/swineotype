#!/usr/bin/env bash
set -euo pipefail

ENV_NAME="swineotype"

# Check if conda is installed
if ! command -v conda &> /dev/null
then
    echo "[ERROR] conda could not be found"
    echo "Please install Miniconda or Anaconda and add it to your PATH."
    exit 1
fi

echo "==> Creating conda environment: $ENV_NAME"

# Create a new environment with Python 3.11
conda create -n "$ENV_NAME" python=3.11 -y

# Install dependencies
echo "==> Installing dependencies into $ENV_NAME"
conda install -n "$ENV_NAME" -c conda-forge -c bioconda \
    blast \
    samtools \
    bcftools \
    pandas \
    pyyaml \
    snakemake \
    kma \
    click \
    pytest -y

echo "==> Installing swineotype"
conda run -n "$ENV_NAME" pip install .

echo "==> Done."
echo
echo "To use swineotype, run:"
echo "    conda activate $ENV_NAME"
echo
echo "Then you can run the tool, for example:"
echo "    swineotype --help"
