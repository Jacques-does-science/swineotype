# swineotype

**A unified bioinformatics toolkit for serotyping key porcine bacterial pathogens.**

`swineotype` provides a streamlined command-line interface for the serotyping of **Streptococcus suis** and **Actinobacillus pleuropneumoniae (APP)** from genome assemblies. It integrates a native, high-performance genotyping algorithm for *S. suis* with an automated adapter for the established `serovar_detector` workflow for APP, offering a consistent user experience and unified output format.

---

## Installation

*Prerequisites: git, Conda (Miniconda/Anaconda).*

1.  **Clone the repository** (Recursive clone is required for submodules):
    ```bash
    git clone --recursive https://github.com/Jacques-does-science/swineotype.git
    cd swineotype
    ```
    
2.  **Install Environment**:
    Run the included script to create the `swineotype` Conda environment and install dependencies.
    ```bash
    bash scripts/install_swineotype.sh
    ```

3.  **Activate**:
    ```bash
    conda activate swineotype
    ```
    
---

## Usage

Serotyping can be performed by specifying the target species (`suis` or `app`), input assemblies, and output directory.

### *S. suis* Serotyping
```bash
swineotype \
  --species suis \
  --assembly "data/swine_isolates/*.fasta" \
  --out_dir results_suis \
  --merged_csv results_suis/summary_report.csv \
  --threads 8
```

### APP Serotyping
```bash
swineotype \
  --species app \
  --assembly "data/app_isolates/*.fasta" \
  --out_dir results_app \
  --merged_csv results_app/summary_report.csv \
  --threads 8
```

---

## Molecular Methods & Algorithm

`swineotype` employs distinct molecular strategies optimized for the capsule genetics of each pathogen.

### *Streptococcus suis* (Native Pipeline)

The *S. suis* serotyping module uses a two-stage hierarchical algorithm designed to resolve the significant genetic overlap between specific serotypes (e.g., 2, 1/2, 1, 14).

**Stage 1: *wzx/wzy* Homology Screening**
The tool performs a **BLASTn** search of the input assembly (database) against a curated reference panel of *wzx* (flippase) and *wzy* (polymerase) genes (query).
-   **Scoring**: High-Scoring Pairs (HSPs) are filtered by coverage and percent identity. Valid hits contribute cumulative BitScores to their respective serotypes.
-   **Assignment**: A serotype is assigned if the top-scoring candidate meets the plurality and delta thresholds (score margin > second best).
-   **Ambiguity**: Certain serotypes (e.g., **1 vs 14**, **2 vs 1/2**) are clinically distinct but genetically identical at the *wzx/wzy* loci. These trigger Stage 2.

**Stage 2: SNP-Based Resolution**
For unresolved pairs, the tool targets specific serotype-determining Single Nucleotide Polymorphisms (SNPs).
1.  **Locus Identification**: A targeted BLASTn locates the relevant gene region (e.g., *cpsK*) in the assembly.
2.  **Genotyping**: The specific base at the diagnostic position (e.g., position 483 in *cpsK*) is extracted.
3.  **Resolution**: The base is compared against the reference logic (e.g., `G` = Serotype 14, `C/T` = Serotype 1) to make a definitive call.

### *Actinobacillus pleuropneumoniae* (Adapter Pipeline)

For APP, `swineotype` functions as an automated wrapper for the third-party **serovar_detector** workflow.
-   **KMA Alignment**: Utilizes the K-mer Alignment (KMA) algorithm to map assemblies against a validated APP capsule locus database.
-   **Automation**: `swineotype` handles the complex initialization of the Snakemake workflow, generating the required `samples.tsv` and `config.yaml` manifests dynamically at runtime within the output directory.

---

## 📂 Output Interpretation

### Common Output Files

| File/Directory | Description |
| :--- | :--- |
| `[out_dir]/` | Root directory containing per-sample subdirectories. |
| `--merged_csv` | **Primary Result.** A consolidated CSV table containing results for all input samples. |

### Interpretation of Summary Columns

The summary CSV contains the following key fields:

| Column | Explanation |
| :--- | :--- |
| `sample` | Filename/ID of the input assembly. |
| `final_serotype` | The definitive serotype call (e.g., `2`, `14`, `APP_5`). |
| `status` | Confidence level or method used: <br>• **STAGE1**: Resolved solely by *wzx/wzy* homology. <br>• **STAGE2**: Resolved by SNP analysis (high confidence). <br>• **NO_CALL**: Insufficient evidence for assignment. |
| `stage1_top` | (Debug) The best hit from the initial detailed gene screen. |
| `base` | (Debug) For Stage 2, the specific nucleotide base found at the varying site. |

### Specific Note on APP Results Structure
When running `--species app`, you will observe a subdirectory named `app_detector/`.
*   **Purpose**: This is an encapsulated run-directory required by the external Snakemake workflow.
*   **Contents**: It contains the intermediate `config.yaml`, `sample_sheet.csv`, and symbolic links newly generated for that specific run.
*   **Results**: The raw output from the external tool can be found in `app_detector/results/serovar.tsv`.
