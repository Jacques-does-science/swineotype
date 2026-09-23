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

`swineotype` requires **Python ≥ 3.10** and declares its own runtime dependencies (`click`, `pandas`, `PyYAML`) in `setup.py`, so `pip install .` now produces a working install on its own. It previously declared none, and only worked because the Conda script installed them separately. External tools — `blastn` and `makeblastdb` for *S. suis*, plus `snakemake`, `kma` and R for the APP adapter — still come from Conda.

---

## Usage

Serotyping can be performed by specifying the target species (`suis` or `app`), input assemblies, and output directory.

> **Input must be uncompressed FASTA.** Gzipped assemblies are not supported and currently fail with a raw `makeblastdb` error rather than a clear message. Decompress first:
> ```bash
> gunzip -k isolate.fasta.gz
> ```

### *S. suis* Serotyping
```bash
swineotype --species suis --assembly "data/swine_isolates/*.fasta" --out_dir results_suis --threads 8
```

`--merged_csv` is optional: `results_suis/swineotype_summary.csv` and `swineotype_run.json` are written regardless. A pattern that matches no file is now an **error** (exit code 2) rather than a successful run with an empty summary.

`--input_species` records an identification you established elsewhere (ANI, a marker panel). It is stored verbatim and never used to make a call — this tool does not identify species. See [What this does not do](#what-this-does-not-do).

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

-   **Evidence aggregation — "best consistent copy"**: an allele's HSPs are first partitioned into *physical copies* (same contig, same strand, no query overlap, collinear in query and subject). The allele's evidence is the single highest-scoring copy, and **coverage, identity and score are all computed from that one copy**. Previously each of the three came from a different set of HSPs — coverage unioned over all of them, identity length-weighted over all of them including overlapping duplicates, score over a greedy non-overlapping subset — so a full-length 100%-identity match passed on its own but failed once five overlapping 80%-identity copies were added.
-   **Split-contig support**: if the best copy alone falls short of `min_cov`, HSPs on disjoint query intervals are chained in, best-score-first. The result is flagged `split` and records every contig it drew on, because fragments from unrelated locations are not by themselves one intact gene.
-   **Serotype-specific gene required**: only a type whose ***wzy*** was found may be called. *wzx* is conserved across serotypes — one *wzx* routinely matches four different references at 95–98% identity — so it cannot carry a call on its own. *wzy* is the discriminating gene, which is why the published multiplex PCR schemes target it. (The rule is *wzy*-only, not "both genes", because serotype 14 has no *wzx* reference. Configurable via `require_wzy`.)
-   **Families, not individual labels**: confidence is established between **families** before anything is resolved within one. The families are `{1,14}`, `{2,1/2}`, and every other type on its own. Within a family the score is the **best score per marker class** (*wzx*, *wzy*), not the sum: the serotype 1 and serotype 14 references differ at a single site, so summing them would manufacture two pieces of evidence out of one. Across families the scores add up to a denominator, and the leader must clear both `plurality` (share of total family evidence) and `delta` (raw bit-score margin over the runner-up **family**).
-   **Ambiguity inside a family** (**1 vs 14**, **2 vs 1/2**) is expected — those serotypes are clinically distinct but genetically identical outside *cpsK* — and is what Stage 2 exists to settle. **Ambiguity between families is not**, and Stage 2 cannot settle it: reporting `NO_CALL_FAMILY_AMBIGUOUS` is the honest outcome.

> `plurality` and `delta` are **heuristics, not calibrated probabilities**. They are bit-score bookkeeping thresholds that have not been fitted to any labelled collection, and the same is true of `min_pid`, `min_cov`, `min_res_pid` and `min_res_alen`. Treat a passing result as "the evidence was not contradictory", not as a confidence level.

**Stage 0: Species gate**
Before any serotype is assigned, the organism the winning *cps* locus belongs to is checked against the target species. Six of the original 35 serotypes now belong to other taxa; their loci are still in the panel, and a hit to one is reported as that organism rather than as a serotype. See [Species](#species).

**Stage 2: SNP-Based Resolution**
For unresolved pairs, the tool targets the serotype-determining SNP in *cpsK*, which encodes the glycosyltransferase CpsK. A single residue (161) sets whether the enzyme adds galactose or *N*-acetylgalactosamine to the CPS side chain, and that is the only difference between each pair.

1.  **Locus Identification**: A targeted BLASTn locates *cpsK* in the assembly. Qualifying alignments are then collapsed onto **distinct physical loci** — same contig, overlapping subject coordinates — so that the two or three near-identical resolver references that find one copy count as one copy. Different query names are *not* evidence of multiple copies.
2.  **Genotyping**: The whole diagnostic **codon** (query positions `pos-2 … pos`) is read **out of the BLAST alignment itself**, walking the gapped alignment column by column. Offset arithmetic on the HSP start would assume an ungapped alignment, and a single upstream indel — the characteristic Oxford Nanopore error — is enough to shift the read-out and flip the call. Reading only the wobble base was worse still: any triplet ending in G, including `AGG` (Arg), was accepted as `TGG` (Trp).
3.  **Validation**: the three bases must all have been recovered, must be adjacent in the subject in the alignment's own orientation, must be unambiguous, and must form one of the three documented codons. Anything else — `AGG`, `CGG`, `TAG`, an ambiguity code, a deletion at the site, an alignment that stops short of the codon, an insertion inside it — is reported with its own state and **withholds the exact call**.
4.  **Conflict**: if two distinct physical loci imply different serotypes, the tool returns a conflict and leaves `final_serotype` empty. Both pieces of evidence are kept in `resolver_loci`. Bit score and FASTA order do **not** break the tie; previously the highest-scoring HSP silently won.
5.  **Coding integrity**: recovering a coordinate across an indel is not proof that the gene around it still reads. The alignment is checked for a frameshift (unbalanced indels) and for a premature stop in the reference's reading frame. A positively-detected disruption withholds the call; a partial alignment is reported as `UNASSESSED`, which is **not** a synonym for intact.
6.  **Resolution**: the codon is interpreted using the `G_serotype` / `CT_serotype` fields declared in each reference's own FASTA header.

| Locus group | Base | Codon 161 | Residue | Side-chain sugar | Serotype |
| :--- | :--- | :--- | :--- | :--- | :--- |
| 2 / 1&#8239;2 | `G` (pos 483) | `TGG` | Trp | Gal | **2** |
| 2 / 1&#8239;2 | `C`/`T` (pos 483) | `TGY` | Cys | GalNAc | **1/2** |
| 1 / 14 | `G` (pos 492) | `TGG` | Trp | Gal | **14** |
| 1 / 14 | `C`/`T` (pos 492) | `TGY` | Cys | GalNAc | **1** |

The two groups declare different positions for the same residue because the 1/14 references carry 9 extra bases at their 5′ end. The position is therefore declared **per reference** in its header, never hard-coded. The difference between 483 and 492 is intentional and verified; do not unify them. See `tests/test_golden_resolver.py`, which asserts this against the shipped data.

#### Resolver reference provenance

Each of the four resolver references is a verbatim slice of a public record. Accession**.version**, source coordinates, strand, diagnostic position, the SHA-256 of the normalized sequence and the extraction method are pinned in [`swineotype/data/suis_resolver_refs.manifest.json`](swineotype/data/suis_resolver_refs.manifest.json).

| Reference | Source | Coordinates (1-based, inclusive) | Strand | Length | Diagnostic position | Codon |
| :--- | :--- | :--- | :--- | ---: | ---: | :--- |
| `cps1L`   | `AB737817.1` | 13538–14551 | + | 1014 | 492 | `TGT` |
| `cps1/2K` | `AB737816.1` | 15240–16244 | + | 1005 | 483 | `TGT` |
| `cps2K`   | `BR001000.1` | 15237–16241 | + | 1005 | 483 | `TGG` |
| `cps14K`  | `AB737822.1` | 13505–14518 | + | 1014 | 492 | `TGG` |

Verify the shipped file against the manifest — **offline**, no network, and this is what the test suite runs:

```bash
python scripts/resolver_refs.py check
```

Re-extract from ENA and rewrite the FASTA (network; opt-in, never run by the tests). It refuses to write if any slice no longer matches its pinned checksum:

```bash
python scripts/resolver_refs.py regenerate
```

`cps2K` and `cps14K` previously shipped with a **two-base deletion** relative to their source records, which frameshifted the 3′ two-thirds of each gene and introduced five internal stop codons apiece. Because the deletion sat *downstream* of the diagnostic site, every position-level test still passed. The manifest records the checksum of that shipped state as `damaged_sha256`, and `check` rejects it.

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
| `[out_dir]/<sample>/` | Per-sample BLAST debug TSVs (`wzxwzy_vs_asm.tsv`, `resolver_vs_asm.tsv`). Named by the `run_dir` column. |
| `[out_dir]/.swineotype_cache/` | Working directory: the staged (line-ending-normalised) assembly and its BLAST database, both named by a content hash. Safe to delete; it is rebuilt on the next run. |
| `[out_dir]/swineotype_summary.csv` | **Primary Result.** All samples from this run, always written. Overwritten each run. |
| `[out_dir]/swineotype_run.json` | Run record: swineotype version, timestamp, command line, resolved input paths, the **effective configuration** (every threshold actually used), SHA-256 of both reference FASTAs, the resolver manifest version, and per-status counts. |
| `--merged_csv` | Optional additional CSV, **appended** across runs for batch accumulation. |

> **Compatibility change.** The per-run summary and run record are new and are written **by default**. Previously nothing machine-readable was produced unless `--merged_csv` was passed, and nothing recorded which thresholds or reference data produced a call. `--merged_csv` keeps its old append-on-existing behaviour, so use a fresh path per run or delete it first to avoid mixing results.

### Interpretation of Summary Columns

The summary CSV contains the following key fields:

| Column | Explanation |
| :--- | :--- |
| `sample` | Input assembly filename without its extension. This is the key the APP results are merged on. |
| `sample_path` | Absolute path of the input assembly, for provenance. |
| `run_dir` | Name of this sample's subdirectory under `[out_dir]`, holding its BLAST debug TSVs. Suffixed with a short hash when two inputs share a filename. |
| `matched_reference_taxon` | Organism the matched *cps* **reference** is labelled with. This is reference metadata, **not** an identification of the input. Empty when nothing matched. |
| `input_species` | Species of the input as established independently, if `--input_species` was supplied. This tool never measures it. |
| `species_assessment` | `NOT_ASSESSED` (the default — nothing measured the input) or `USER_SUPPLIED`. |
| `species` | **Deprecated alias** of `matched_reference_taxon`, kept so existing readers do not break. |
| `final_serotype` | The exact serotype call (e.g., `2`, `14`, `APP_5`). Empty whenever the exact member was not established. |
| `family_serotype` | The family-level result when the exact member was withheld, e.g. `2 or 1/2`. Empty when an exact call was made or when no family was established. |
| `status` | Method used, or why no call was made: <br>• **STAGE1**: The family contains one type, and it is decisively supported. <br>• **STAGE2**: Family resolved, then the exact member read from the *cpsK* codon. <br>• **FAMILY_ONLY**: The family is established but the exact member is withheld (conflict, invalid codon, detected coding disruption, or no qualifying resolver alignment). `family_serotype` carries the result. <br>• **NO_CALL_FAMILY_AMBIGUOUS**: Competing *cps* families are not separated by the evidence. The *cpsK* site cannot arbitrate between families, so nothing is resolved. <br>• **NO_CALL_PAIR_CONFLICT**: Stage 2 returned a serotype outside the family Stage 1 pointed at; the call is withheld. Please report this. <br>• **NON_TARGET_SPECIES**: The best *cps* match is a reference labelled with another organism; no serotype is reported. See **Species** below. <br>• **NO_WZY_MATCH**: *cps* genes are present but no reference *wzy* matched. See **Novel capsular loci** below. <br>• **NO_CPS_MATCH**: Nothing in the panel matched. No species is inferred. |
| `warnings` | Semicolon-separated flags. Empty is the normal case. See below. |
| `stage1_top` | (Debug) Best-scoring **callable** type from the *wzx/wzy* screen — i.e. one with *wzy* support. |
| `stage1_family`, `stage1_family_label` | (Debug) The winning family (`1_vs_14`, `2_vs_1_2`, or `type:<n>`) and its readable form. |
| `stage1_wzx_only` | (Debug) Best type matched by *wzx* alone. A lead, **not** a serotype: *wzx* is not serotype-specific. Populated on `NO_WZY_MATCH`. |
| `stage2_status` | (Debug) What Stage 2 found: `SKIPPED`, `OK`, `NO_HSP_OR_LOW_QUAL`, `CONFLICTING_COPIES`, `INVALID_TRIPLET`, or `CODING_DISRUPTED`. |
| `triplet` | (Debug) The whole diagnostic codon read from the assembly. `-` marks a deletion, `?` a position the alignment did not reach. |
| `triplet_status` | (Debug) `OK`, `UNEXPECTED_CODON`, `AMBIGUOUS`, `DELETED`, `INCOMPLETE`, or `NON_CONTIGUOUS`. Only `OK` permits an exact call. |
| `coding_status` | (Debug) `INTACT` (alignment spans the whole reference CDS, no frameshift or premature stop), `DISRUPTED` (one was detected), or `UNASSESSED` (partial alignment — **not** a synonym for intact). |
| `resolver_loci` | (Debug) Every distinct physical *cpsK* locus found, with its coordinates, codon, codon status and implied serotype. This is what a conflict looks like. |
| `base` | (Debug) The nucleotide at the diagnostic site — the last base of `triplet`. On its own it is **not** sufficient to make a call. |
| `contig`, `contig_pos`, `strand` | (Debug) Where in the assembly the diagnostic site was read (best-scoring locus). |

> **Compatibility change.** `species` used to be filled from `config["target_species"]` whenever no reference taxon was available, so an assembly that matched nothing was reported as *Streptococcus suis*. It is now empty in that case, and the three columns above separate what the reference says from what (if anything) is known about the input. `NO_CALL_STAGE2` has been replaced by the more specific `FAMILY_ONLY`, `NO_CALL_FAMILY_AMBIGUOUS` and `NO_CPS_MATCH`.

### Warnings

| Flag | Meaning |
| :--- | :--- |
| `competing_cps_families:<a>/<b>` | Two *cps* families are not separated by the evidence. Nothing is resolved; the *cpsK* site cannot arbitrate between families. |
| `family_fraction=<f>;family_delta=<d>` | The two family-level numbers behind that decision, so the margin can be inspected. |
| `family_disagrees_with_top_label:<family>/<pair>` | The best-supported *family* is not the one the single highest-scoring type belongs to. The family decision governs; this records the disagreement. |
| `conflicting_resolver_copies:<a>,<b>` | Distinct physical *cpsK* loci imply different serotypes; a codon outside the scheme (e.g. `AGG`) is listed as itself. Synonymous codons (`TGT`/`TGC`) agree. The exact label is withheld and both loci are kept in `resolver_loci`. |
| `resolver_loci:<n>` | How many distinct physical loci were found. |
| `resolver_triplet_<state>:<codon>` | The diagnostic codon is not interpretable — see `triplet_status`. The call is withheld. |
| `resolver_coding_disrupted:<detail>` | A frameshift or premature stop was detected around the diagnostic site. The call is withheld. |
| `resolver_coding_integrity_unassessed` | The alignment was partial, so coding integrity could not be judged. The call is **not** withheld for this; it is disclosed. |
| `stage2_outside_stage1_pair:<serotype>` | Internal consistency check failed; the call is withheld. Please report this. |
| `no_cps_reference_matched` | Nothing in the panel matched. No species is inferred. |
| `no_reference_wzy_matched` | No serotype-specific gene was found; no serotype can be assigned. |
| `nearest_wzx_relative:cps_type_<n>` | The closest *wzx* relative, as a lead for follow-up. Not a serotype call. |
| `locus_integrity_not_assessed` | A *wzx*-only match says no qualifying reference *wzy* was found. It does **not** say the locus is intact, and it does not say it is novel. |

---

## Novel capsular loci

At least 32 novel *cps* loci (NCLs) have been described in non-serotypeable *S. suis* since 2015, and they are **not** in this reference panel. An isolate carrying one produces `NO_WZY_MATCH`: its *wzx* matches the panel (that gene is conserved), but its *wzy* does not match any of the 29 serotypes at usable identity.

**What `NO_WZY_MATCH` does and does not establish.** It establishes exactly one thing: no qualifying reference *wzy* was found. It does **not** establish that the isolate carries an intact capsular locus, and it does not establish that the locus is novel. A truncated *wzy*, a *wzy* broken across a contig boundary, a deleted *wzy*, a low-quality assembly and a genuinely new polymerase all look identical from here. The tool used to emit a `candidate_novel_capsular_locus` warning that asserted both; it now emits `locus_integrity_not_assessed` instead. An NCL is one hypothesis worth testing, not a finding. To test it:

1. Note `stage1_wzx_only`; it names the closest *wzx* relative.
2. Extract the *cps* locus from the assembly. The `wzxwzy_vs_asm.tsv` in the sample's `run_dir` gives the coordinates of the *wzx* and any partial *wzy* hits; the locus is typically 18–30 kb and its conserved *cpsA/cpsB/cpsD* regulatory block marks one end.
3. BLAST it against NCBI `nt`. Described NCLs were submitted by Zheng, Qiu, Huang and Králová.

Note that the *cps* locus is **not** always flanked by *orfZ–orfX* and *aroA*: seven of the 35 reference loci are flanked by *glf* instead, and five chromosomal arrangements are known ([Okura et al. 2013](https://doi.org/10.1128/AEM.03742-12)). Do not assume a fixed pair of flanking genes when extracting.

---

## Species

*S. suis* was originally described with 35 serotypes. Six have since been moved to other taxa:

| Former serotype | Current species | Reassignment |
| :--- | :--- | :--- |
| 20, 22, 26 | *Streptococcus parasuis* | Nomoto et al. 2015, IJSEM 65:438 |
| 33 | *Streptococcus ruminantium* | Tohya et al. 2017, IJSEM 67:2224 |
| 32, 34 | *Streptococcus orisratti* | Hill et al. 2005, Vet Microbiol 107:63 |

That leaves **29 true *S. suis* serotypes: 1–19, 21, 23–25, 27–31 and 1/2.**

Their *cps* loci are deliberately **kept** in the reference panel — they are how the tool recognises those organisms. Deleting them would turn a wrong answer into a silent `NO_CALL`. Instead, each reference carries a `[species=...]` tag, and a best hit to a non-target species short-circuits the pipeline: the `species` column names the organism, `final_serotype` is left empty, and `status` is `NON_TARGET_SPECIES`.

```
iso_parasuis.fasta  Streptococcus parasuis   (cps type 22)   NON_TARGET_SPECIES
iso_sero9.fasta     Streptococcus suis       9               STAGE1
```

The species a *cps* type belongs to is declared **in the FASTA header, not in the Python**, so correcting or extending it needs no code change:

```
>wzy_AB737828 [locus=SS_CPS] [type_id=22] [species=Streptococcus parasuis] \
  [source_gb=AB737828] [reassigned_from=S. suis serotype 22] \
  [reassignment_ref=Nomoto et al. 2015, IJSEM 65:438]
```

The target species itself is configurable (`target_species`, default `Streptococcus suis`).

### What this does not do

**This is a *cps*-locus inference, not a species identification.** It catches an organism whose capsule locus is in the panel. It does **not** catch:

- an organism outside *S. suis* carrying a *cps* locus that resembles a true *S. suis* serotype — capsule loci transfer horizontally between related streptococci;
- the wider *S. suis* complex. At least a dozen further species and provisional lineages sit alongside *S. suis*, are routinely misidentified by MALDI-TOF, and are **not** reliably separated by the `recN` PCR that most *S. suis* pipelines gate on ([Fittipaldi et al. 2025, *J Clin Microbiol*](https://doi.org/10.1128/jcm.01030-25)).

A genuine species call needs genome-level evidence — ANI against type strains (published thresholds: 93.17% for authentic *S. suis*, 92.33% for cluster delineation within the complex) or a conserved-marker panel. **Confirm the species independently before acting on a serotype.**

## Testing

```bash
python -m pytest tests/ -q
```

The suite is offline by default: reference invariants are checked against the pinned manifest, never against a live service. Tests that shell out to `blastn`/`makeblastdb` (`tests/test_real_blast.py`) **skip with an explicit reason** when those tools are absent — a skipped integration test is not a passing one.

`tests/test_real_blast.py` builds its inputs from the project's own references. That makes it a check on behaviour for sequences of a given shape; it is not a sample of any population and **no accuracy figure may be read off it**.

---

### Specific Note on APP Results Structure
When running `--species app`, you will observe a subdirectory named `app_detector/`.
*   **Purpose**: This is an encapsulated run-directory required by the external Snakemake workflow.
*   **Contents**: It contains the intermediate `config.yaml`, `sample_sheet.csv`, and symbolic links newly generated for that specific run.
*   **Results**: The raw output from the external tool can be found in `app_detector/results/serovar.tsv`.
