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
conda install -n "$ENV_NAME" -c conda-forge -c bioconda \
    blast \
    samtools \
    bcftools \
    pandas \
    pyyaml \
    snakemake \
    kma \
    click \
    pytest \
    peppy -y

echo "==> Installing swineotype"
conda run -n "$ENV_NAME" pip install -e .

echo "==> Note: If you have system R installed (e.g., via Homebrew),"
echo "    you may need to install R packages manually:"
echo "    R -e 'install.packages(c(\"dplyr\", \"purrr\", \"readr\", \"stringr\", \"tidyr\", \"tibble\", \"yaml\", \"logger\"), repos=\"https://cloud.r-project.org\")'"

echo "==> Done."
echo
echo "To use swineotype, run:"
echo "    conda activate $ENV_NAME"
echo
echo "Then you can run the tool, for example:"
echo "    swineotype --help"
