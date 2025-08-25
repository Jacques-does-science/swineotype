#!/usr/bin/env bash
set -euo pipefail

ENV_NAME="swineotype"
THIRDPARTY_DIR="third_party"
SEROVAR_DETECTOR_REPO="https://github.com/KasperThystrup/serovar_detector.git"
SEROVAR_DETECTOR_DIR="${THIRDPARTY_DIR}/serovar_detector"

# Check if conda is installed
if ! command -v conda &> /dev/null
then
    echo "[ERROR] conda could not be found"
    echo "Please install Miniconda or Anaconda and add it to your PATH."
    exit 1
fi

# Check if git is installed
if ! command -v git &> /dev/null
then
    echo "[ERROR] git could not be found"
    echo "Please install git."
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
    kma -y

# Install third-party workflow for APP serotyping
echo "==> Installing third-party workflow for APP serotyping"
mkdir -p "${THIRDPARTY_DIR}"
if [ -d "${SEROVAR_DETECTOR_DIR}" ]; then
    echo "[INFO] ${SEROVAR_DETECTOR_DIR} already exists, skipping clone."
else
    echo "[INFO] Cloning serovar_detector from GitHub..."
    git clone "${SEROVAR_DETECTOR_REPO}" "${SEROVAR_DETECTOR_DIR}"
fi


echo "==> Done."
echo
echo "To use swineotype, run:"
echo "    conda activate $ENV_NAME"
echo
echo "Then you can run the tool, for example:"
echo "    swineotype --help"
