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

# Initialize submodule if needed
if [ -z "$(ls -A third_party/serovar_detector)" ]; then
    echo "==> Initializing serovar_detector submodule..."
    git submodule update --init --recursive
fi


# Create a new environment with Python 3.11
conda create -n "$ENV_NAME" python=3.11 -y

# Install dependencies
echo "==> Installing dependencies into $ENV_NAME"
# samtools and bcftools were dropped: Stage 2 used to extract the diagnostic
# base with `samtools faidx`, but it now reads it out of the BLAST alignment
# directly. Neither is referenced by swineotype or by serovar_detector.
# pandas, snakemake and peppy are what serovar_detector (APP) needs to start
# its workflow, which then builds its own BLAST and R environments with conda
# on first use. swineotype itself needs only click and PyYAML, which
# pyproject.toml declares.
conda install -n "$ENV_NAME" -c conda-forge -c bioconda \
    blast \
    pandas \
    pyyaml \
    snakemake \
    click \
    pytest \
    peppy -y

echo "==> Installing swineotype"
conda run -n "$ENV_NAME" pip install -e .

# The bundled submodule is an unmodified serovar_detector release.
echo "==> Installing serovar_detector (APP) from third_party/serovar_detector"
conda run -n "$ENV_NAME" pip install ./third_party/serovar_detector

echo "==> Done."
echo
echo "To use swineotype, run:"
echo "    conda activate $ENV_NAME"
echo
echo "Then you can run the tool, for example:"
echo "    swineotype --help"
