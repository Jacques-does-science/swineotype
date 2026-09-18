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
For unresolved pairs, the tool targets the serotype-determining SNP in *cpsK*, which encodes the glycosyltransferase CpsK. A single residue (161) sets whether the enzyme adds galactose or *N*-acetylgalactosamine to the CPS side chain, and that is the only difference between each pair.

1.  **Locus Identification**: A targeted BLASTn locates *cpsK* in the assembly.
2.  **Genotyping**: The base aligned to the diagnostic position is read **out of the BLAST alignment itself**, walking the gapped alignment column by column. Offset arithmetic on the HSP start would assume an ungapped alignment, and a single upstream indel — the characteristic Oxford Nanopore error — is enough to shift the read-out and flip the call.
3.  **Resolution**: The base is interpreted using the `G_serotype` / `CT_serotype` fields declared in each reference's own FASTA header.

| Locus group | Base | Codon 161 | Residue | Side-chain sugar | Serotype |
| :--- | :--- | :--- | :--- | :--- | :--- |
| 2 / 1&#8239;2 | `G` (pos 483) | `TGG` | Trp | Gal | **2** |
| 2 / 1&#8239;2 | `C`/`T` (pos 483) | `TGY` | Cys | GalNAc | **1/2** |
| 1 / 14 | `G` (pos 492) | `TGG` | Trp | Gal | **14** |
| 1 / 14 | `C`/`T` (pos 492) | `TGY` | Cys | GalNAc | **1** |

The two groups declare different positions for the same residue because the 1/14 references carry 9 extra bases at their 5′ end. The position is therefore declared **per reference** in its header, never hard-coded. See `tests/test_golden_resolver.py`, which asserts this against the shipped data.

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
| `sample` | Input assembly filename without its extension. This is the key the APP results are merged on. |
| `sample_path` | Absolute path of the input assembly, for provenance. |
| `run_dir` | Name of this sample's subdirectory under `[out_dir]`, holding its BLAST debug TSVs. Suffixed with a short hash when two inputs share a filename. |
| `final_serotype` | The definitive serotype call (e.g., `2`, `14`, `APP_5`). Empty when no call was made. |
| `status` | Method used, or why no call was made: <br>• **STAGE1**: Resolved solely by *wzx/wzy* homology. <br>• **STAGE2**: Resolved by SNP analysis. <br>• **NO_CALL_STAGE2**: Insufficient evidence for assignment. <br>• **NO_CALL_PAIR_CONFLICT**: Stage 2 returned a serotype outside the pair Stage 1 pointed at; the call is withheld rather than reported. |
| `warnings` | Semicolon-separated flags. Empty is the normal case. See below. |
| `stage1_top` | (Debug) The best-scoring serotype from the *wzx/wzy* screen. |
| `base` | (Debug) The nucleotide found at the diagnostic site. `-` means the assembly carries a deletion there. |
| `contig`, `contig_pos`, `strand` | (Debug) Where in the assembly the diagnostic site was read. |

### Warnings

| Flag | Meaning |
| :--- | :--- |
| `stage1_pair_ambiguous:<top>/<second>` | The top two *wzx/wzy* hits point at different locus groups. The top hit is used; the call is worth confirming. |
| `resolver_site_deleted` | The assembly has a deletion at the diagnostic site, so no base could be read. |
| `resolver_base_not_gct:<base>` | The diagnostic site is not G, C or T (e.g. an `N`), so the call is withheld. |
| `stage2_outside_stage1_pair:<serotype>` | Internal consistency check failed; the call is withheld. Please report this. |

> **Note on species.** The tool assigns a *cps* type; it does not confirm the species. The reference panel still contains the loci of former serotypes 20, 22, 26 (now *Streptococcus parasuis*), 33 (*S. ruminantium*) and 32, 34 (*S. orisratti*), so a hit to one of those indicates a different organism, not an *S. suis* serotype. Confirm the species independently.

### Specific Note on APP Results Structure
When running `--species app`, you will observe a subdirectory named `app_detector/`.
*   **Purpose**: This is an encapsulated run-directory required by the external Snakemake workflow.
*   **Contents**: It contains the intermediate `config.yaml`, `sample_sheet.csv`, and symbolic links newly generated for that specific run.
*   **Results**: The raw output from the external tool can be found in `app_detector/results/serovar.tsv`.
