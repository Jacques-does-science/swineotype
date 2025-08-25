# swineotype

A bioinformatics toolkit for serotyping key porcine bacterial pathogens.

`swineotype` is a lightweight, flat‑install toolkit that:
1.  Natively serotypes **Streptococcus suis** from assemblies.
2.  Wraps a third‑party Snakemake workflow to serotype **Actinobacillus pleuropneumoniae** (APP) via a dedicated adapter, keeping a unified interface and results format.

## Biology & Application Context

### Pathogens
-   ***Streptococcus suis***: Serotyped by the primary `swineotype` tool.
-   ***Actinobacillus pleuropneumoniae*** (APP): Serotyped via the `app_serovar_detector_adapter`.

### Why serotyping matters
For both pathogens, the serotype is critical for vaccine development and epidemiological tracking.
-   *S. suis*: Many closely related serotypes differ by only one or two Single Nucleotide Polymorphisms (SNPs) in the capsule synthesis (cps) genes.
-   *APP*: Detection relies on identifying specific capsule gene clusters.

## Pipeline Architecture & Reasoning

### *S. suis* (Native Pipeline)
The native pipeline uses a two-stage approach to balance speed and accuracy.

-   **Stage 1 (Broad Screen)**: A fast `BLAST` search against a curated `wzx/wzy` gene whitelist assigns a possible serogroup (e.g., the ambiguous groups 1/14 or 2/1/2).
-   **Stage 2 (Fine SNP Resolution)**: For ambiguous cases, a targeted check resolves the serotype. For example, a known SNP at position 483 of the `cpsK` gene distinguishes serotype 1 from 14. This is designed to be resolved by alignment and variant calling tools like `samtools` and `bcftools`.

### APP (Adapter for `serovar_detector`)
The adapter makes the third-party `serovar_detector` tool plug-and-play within the `swineotype` ecosystem. `serovar_detector` uses `kma` (k-mer alignment) to map reads/assemblies to a capsule gene database.

The adapter's role is to automate the setup process by:
-   Generating the required `samples.tsv` and `config.yaml` files.
-   Running the Snakemake workflow with the correct paths and parameters.
-   Optionally merging the APP results back into a `swineotype` summary file for unified reporting.

## Installation

The installation is managed by a shell script that sets up a Conda environment.

**Prerequisites**:
-   A Conda distribution (Miniconda or Anaconda).

To install, run the script from the repository root:
```bash
bash scripts/install_swineotype.sh
```
This script will:
1.  Create a Conda environment named `swineotype`.
2.  Install all necessary dependencies (`blast`, `samtools`, `bcftools`, `kma`, `snakemake`, `click`, `pyyaml`, `pytest`, etc.).

To use the tool, activate the environment:
```bash
conda activate swineotype
```

## Usage

### *S. suis* Serotyping
```bash
swineotype \
  --assembly "/path/to/Ssuis/*.fasta" \
  --out_dir /path/to/results \
  --merged_csv /path/to/results/swineotype_summary.csv \
  --threads 4
```

### APP Serotyping
```bash
python -m swineotype.adapters.app \
  --assembly "/path/to/APP/*.fasta" \
  --out_dir "/path/to/results" \
  --threads 4
```

## Configuration
The tool can be configured via a `config.yaml` file. By default, it will look for a `config.yaml` file in the current directory. You can also specify a path to a custom config file using the `--config` option.

The following options can be configured:
- `data_dir`: The directory containing the data files.
- `wzxwzy_fasta`: The name of the `wzx/wzy` whitelist FASTA file.
- `resolver_refs_fasta`: The name of the resolver references FASTA file.
- `tmp_dir`: The directory to use for temporary files.
- `plurality`: The minimum fraction of the total score that the top hit must have to be considered decisive.
- `delta`: The minimum difference between the top two hits to be considered decisive.
- `require_agreement`: Whether to require that the top hits for `wzx` and `wzy` agree.
- `min_pid`: The minimum percent identity for a BLAST hit to be considered.
- `min_cov`: The minimum coverage for a BLAST hit to be considered.
- `min_res_pid`: The minimum percent identity for a resolver BLAST hit to be considered.
- `min_res_alen`: The minimum alignment length for a resolver BLAST hit to be considered.
- `keep_debug`: Whether to keep the intermediate BLAST results.
- `gzip_debug`: Whether to gzip the intermediate BLAST results.
- `clean_temp`: Whether to clean the temporary directory after running.

## File & Folder Philosophy
The project is organized to separate curated data, custom scripts, external code, and run-specific outputs.
```
swineotype/
├─ bin/                      # Entrypoints
├─ swineotype/               # Core Python package
│  ├─ adapters/              # Glue for third-party workflows
│  ├─ __init__.py
│  ├─ blast.py
│  ├─ config.py
│  ├─ main.py
│  └─ stages.py
├─ data/                     # Curated biological references (stable)
├─ third_party/              # Vendored external tools (unmodified)
└─ results/                  # Ephemeral outputs (run-specific)
```
This design ensures clarity and reproducibility. The `results/` directory can be safely deleted and regenerated on each run.
