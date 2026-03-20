# ATMPI Debugging History

This document records the debugging steps, issues identified, and fixes applied to the ATMPI (ASV-to-MAG Matching Pipeline) during the session on 2026-03-20.

## Overview
The pipeline was initially failing due to multiple issues ranging from input data formatting and script bugs to incompatible tool arguments and unrealistic filtering thresholds.

## Issues and Fixes

### 1. Input Data & Configuration
*   **ASV Table Formatting:** `asv_table.tsv` contained multiple FASTA paths on a single line. 
    *   *Fix:* Reformatted the table to a standard two-column tab-separated format (one ASV per line).
*   **Missing Files:** `config.yaml` referenced a missing `mag_taxonomy` file.
    *   *Fix:* Commented out the optional taxonomy path and updated `Snakefile` to handle its absence gracefully.
*   **Compressed File Support:** MAG files were gzipped (`.gz`), but scripts expected plain text.
    *   *Fix:* Added `gzip` and `shutil` support to `extract_16s.py` and `build_mag_contig_db.py`.

### 2. Snakemake Workflow (`Snakefile`)
*   **Rule Ordering:** Rule `all` was at the bottom, causing Snakemake to only run the first rule it encountered.
    *   *Fix:* Moved `rule all` to the top.
*   **BLAST Database Outputs:** Handled BLAST DB file extensions correctly using `multiext`.
*   **Dependency Tracking:** Updated rule inputs/outputs to ensure scripts rerun when their logic changes.

### 3. Script-Specific Bugs
*   **`extract_16s.py`:**
    *   Removed invalid `barrnap --genomic` option; switched to `--outseq` for sequence extraction.
    *   Fixed `read_fasta_sequences` which was trying to open a sequence string as a file path.
    *   Ensured temporary decompressed files are cleaned up.
*   **`prepare_blast_query.py`:** 
    *   Fixed FASTA parsing to extract only the sequence content, excluding headers.
*   **`build_mag_contig_db.py`:**
    *   Removed `-parse_seqids` from `makeblastdb` to resolve the "local id too long" error (IDs were > 50 characters).
*   **`process_blast_results.py`:**
    *   Improved MAG ID parsing (stripping `_16S_rRNA::` and pipe suffixes) to correctly merge 16S and contig hits.
    *   Relaxed `coverage_threshold` to only apply to the ASV side (80% coverage of a MAG contig is impossible for a ~250bp ASV).
*   **`generate_summary.py`:**
    *   Fixed a pandas aggregation bug where the script attempted to calculate the mean of `asv_id` strings.
*   **`filter_target_species.py`:**
    *   Updated argument parsing and `Snakefile` quoting to handle target species names containing spaces (e.g., "Bifidobacterium longum").

## Final Pipeline State
The pipeline now executes successfully from start to finish.

### Results Summary
*   **Total ASVs:** 1923
*   **Total MAGs:** 1316
*   **Successful Assignments:** 537 ASVs matched to 61 MAGs.
*   **Mean Identity:** 98.86%
*   **Strain Variation:** Identified 50 MAGs with multiple ASV variants.

### Verified Outputs
*   `assignments/asv_to_mag_assignments.tsv`
*   `assignments/asv_strain_resolution.tsv`
*   `summary/asv_mag_summary.tsv`
*   `summary/visualization_report.html`

## How to Rerun
To rerun the pipeline with 16 cores:
```bash
conda run -n env-atmpi snakemake --cores 16 --configfile config.yaml
```
To force a full rerun of all steps:
```bash
conda run -n env-atmpi snakemake --cores 16 --configfile config.yaml --forceall
```
