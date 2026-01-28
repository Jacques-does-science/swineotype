# swineotype

**A bioinformatics toolkit for serotyping key porcine bacterial pathogens.**

`swineotype` is a user-friendly command-line tool designed to identify the serotypes of two major pig pathogens:
1.  **Streptococcus suis**: Uses a native, fast, and accurate two-stage algorithms built directly into the tool.
2.  **Actinobacillus pleuropneumoniae (APP)**: Uses a built-in **adapter** to automate a sophisticated third-party workflow (`serovar_detector`), handling all the complex configuration for you.

---

## 🚀 Quick Start

### Installation

For users with little command-line experience: you will need **git** (to download the code) and **Conda** (to manage software environments).

1.  **Clone the repository**:
    ```bash
    git clone --recursive https://github.com/Jacques-does-science/swineotype.git
    cd swineotype
    ```
    *Note: `--recursive` is important! It downloads the external APP tools nested inside this project.*

2.  **Install**:
    Run the provided helper script. It creates a secluded "environment" so installing this tool doesn't mess up your other programs.
    ```bash
    bash scripts/install_swineotype.sh
    ```

3.  **Activate**:
    You must run this command every time you open a new terminal to use the tool.
    ```bash
    conda activate swineotype
    ```

### Running the Tool

**For *S. suis*:**
```bash
swineotype \
  --species suis \
  --assembly "data/*.fasta" \
  --out_dir results_suis \
  --merged_csv results_suis/summary.csv
```

**For *APP*:**
```bash
swineotype \
  --species app \
  --assembly "data/*.fasta" \
  --out_dir results_app \
  --merged_csv results_app/final_results.csv
```

---

## 📖 Command Line Manual

When running `swineotype`, you provide **Inputs** (your data) and **Outputs** (where to save results).

| Argument | Required | Description |
| :--- | :---: | :--- |
| `--species` | **Yes** | Choose which pathogen to analyze. Options: `suis` or `app`. |
| `--assembly` | **Yes** | Path to your assembly FASTA files. You can list multiple files or use wildcards (e.g., `"*.fasta"`). **Note: Use quotes around wildcards to let the tool handle them correctly.** |
| `--out_dir` | **Yes** | The folder where `swineotype` will create analysis files. If it doesn't exist, it will be created. |
| `--merged_csv` | No | A convenient summary file collecting results from ALL inputs into one table. Recommended. |
| `--threads` | No | Number of CPU cores to use. Defaults to using half your available cores. |
| `--config` | No | Advanced: Path to a custom `config.yaml` to tweak internal thresholds. |

---

## 📂 Understanding the Output

### *S. suis* Output
When running with `--species suis`, the tool looks for specific marker genes (`wzx`, `wzy`) and sometimes specific mutations (SNPs).

**1. Summary File (`--merged_csv`)**
This is the most important file. Columns include:
- `sample`: The name of your input file.
- `status`:
    - `STAGE1`: Serotype found quickly by gene match.
    - `STAGE2`: Serotype resolved by checking specific mutations (SNPs).
    - `NO_CALL`: Could not determine serotype confidently.
- `final_serotype`: The predicted serotype (e.g., `2`, `1`, `14`, `1/2`).

**2. Output Directory (`--out_dir`)**
Inside `out_dir`, you will find a folder for **each sample**. These contain detailed "debug" info:
- `wzxwzy_vs_asm.tsv`: Raw BLAST results searching for marker genes.
- `resolver_vs_asm.tsv`: (If Stage 2 ran) Raw BLAST results searching for SNP regions.

---

### *APP* Output (and "app_detector")
When running with `--species app`, `swineotype` acts as a **smart wrapper**. It prepares files and runs a complex workflow called `serovar_detector` for you.

**1. Why is there a folder called `app_detector` in my results?**
You will see a structure like `results_app/app_detector`.
- This is **NOT** the software installation.
- This is a **run directory**. The external tool (`serovar_detector`) requires a specific folder structure (configuration files, symbolic links, logs) to run. `swineotype` creates this for you inside your output folder so that the analysis is self-contained.

**2. Key Result Files**
- `results_app/final_results.csv`: (If you used `--merged_csv`) The simple, clean summary of serotypes for your samples.
- `results_app/app_detector/results/serovar.tsv`: The raw output from the underlying `serovar_detector` tool.

---

## 🔬 How It Works (Under the Hood)

### *S. suis* (Native Pipeline)
1.  **Broad Screen**: We blast your assembly against known `wzx` and `wzy` genes. If we find a perfect match that belongs to a unique serotype, we stop there.
2.  **Fine Resolution**: Some serotypes (like 1 vs 14, or 2 vs 1/2) are almost identical. If needed, the tool zooms in on specific DNA letters (SNPs) to tell them apart.

### APP (Wrapper Pipeline)
APP serotyping is complex and uses 3rd party tools (Snakemake, KMA, serovar_detector).


    1.  It creates a temporary "workspace" in your output folder (`app_detector/`).
    2.  It creates the necessary `sample_sheet.csv` and `config.yaml`.
    3.  It runs the workflow for you.
    4.  It reads the result and gives you a clean CSV.

---

## File Structure

```
swineotype/
├─ scripts/                  # Helper scripts (installers)
├─ swineotype/               # The main python code
│  ├─ adapters/              # Code that "wraps" the APP workflow
│  └─ data/                  # Reference databases for S. suis
├─ third_party/              # External tools (where the APP logic lives)
└─ results/                  # Your analysis outputs go here
```
