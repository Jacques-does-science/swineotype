#!/usr/bin/env bash
set -euo pipefail

ENV_NAME="swineotype"

echo "==> Creating conda environment: $ENV_NAME"

# Create a new environment with Python 3.11
conda create -n "$ENV_NAME" python=3.11 -y

# Activate environment (note: in scripts, use conda run instead of activate)
echo "==> Installing dependencies into $ENV_NAME"
conda install -n "$ENV_NAME" -c conda-forge -c bioconda \
    blast \
    samtools \
    pandas \
    pyyaml \
    biopython \
    tqdm -y

echo "==> Done."
echo
echo "To use swineotype, run:"
echo "    conda activate $ENV_NAME"
echo
echo "Then test with:"
echo "    swineotype --help"
