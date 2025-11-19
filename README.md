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

Follow these steps to set up `swineotype` on your local machine.

**Prerequisites**:
-   You must have a Conda distribution (Miniconda or Anaconda) installed. For installation instructions, please see the [official Conda documentation](https://conda.io/projects/conda/en/latest/user-guide/install/index.html).
-   You will also need `git` to clone the repository.

**Step-by-step instructions**:

1.  **Clone the repository**:
    Open your terminal and run the following command to clone the repository to your local machine:
    ```bash
    git clone --recursive https://github.com/Jacques-does-science/swineotype.git
    ```

2.  **Navigate to the project directory**:
    ```bash
    cd swineotype
    ```

3.  **Run the installation script**:

    This script will create a dedicated Conda environment with all the necessary dependencies and install the `swineotype` package using the `setup.py` file. This makes the `swineotype` command available in your environment.

    ```bash
    bash scripts/install_swineotype.sh
    ```
    *Note: The script will automatically initialize the `serovar_detector` submodule if it wasn't cloned recursively.*

4.  **Activate the Conda environment**:
    To start using the tools, you must activate the `swineotype` environment:
    ```bash
    conda activate swineotype
    ```

After these steps, you will be ready to use `swineotype`. The `swineotype` command will be available in your `PATH`.

## Usage

`swineotype` can serotype both *S. suis* and *A. pleuropneumoniae* (APP). You can specify the target species using the `--species` flag.

### *S. suis* Serotyping (default)
```bash
swineotype \
  --assembly "/path/to/Ssuis/*.fasta" \
  --out_dir /path/to/results \
  --merged_csv /path/to/results/swineotype_summary.csv \
  --threads 4
```
*Note: `--species suis` is the default and can be omitted.*

### APP Serotyping
```bash
swineotype \
  --species app \
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
├─ setup.py                  # Installation script
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
