# RPKM Exon Coverage Curve (DeCovA-style) — Design

**Date:** 2026-06-17
**Project:** KODEX / run 1 NGS reanalysis
**Status:** Approved design, ready for implementation plan

## Purpose

Produce a single R tool that turns the per-exon RPKM matrix from the run 1
run into a coverage curve across all 22 patients. The curve serves two goals at
once:

1. **Coverage QC** — verify every exon is adequately covered in every sample.
2. **Visual CNV detection** — eyeball exon-level deletions/duplications in
   TSC1/TSC2 (and optionally the whole panel) by comparing each patient's RPKM
   to the cohort baseline.

This replicates the visualization of HCL-HUBL's **DeCovA** tool
(`graphAllSample`: one colored depth line per sample, all overlaid; CNV module
centered on the cohort median with del < 0.8x / dup > 1.2x thresholds), but
works directly from the existing RPKM matrix in pure R — no BAMs, no Perl,
no bedtools/GATK.

## Inputs

- **File:** `<run_dir>/KoverGraph/run_1_Blanket.xlsx`
- **Sheet:** `RPKMs`
- **Layout:** columns `Gene, Transcrit, Exon, <sample_1> ... <sample_22>`;
  one row per exon. Notes on the raw layout:
  - `Gene` and `Transcrit` are fill-down (blank/`None` after the first row of
    each gene block) and must be forward-filled.
  - `Exon` rows include padding/UTR sub-rows such as `E01_UTR5_lpad`,
    `E01_lpad`, `E01`, `E01_rpad`. The core exon row is the plain `E0N` form.
  - The sheet spans the whole Twist TSC panel (e.g. `IQGAP2`, `TSC1`, `TSC2`),
    not only the disease genes.

## Package Stack (Approach A)

| Concern | Package |
|---|---|
| Read xlsx sheet | `readxl` |
| Reshape + statistics | `dplyr`, `tidyr` |
| Static figure | `ggplot2` |
| Interactive HTML | `plotly` (`ggplotly()`) |

One ggplot object drives both the static export and the interactive version,
so there is a single source of truth for the plot.

## Data Flow

1. **Load** sheet `RPKMs` with `readxl`.
2. **Forward-fill** `Gene` and `Transcrit` down each gene block.
3. **Filter exon rows** — keep core exons (`E0N`), drop `_lpad`/`_rpad`/`_UTR`
   padding rows by default; a flag (`include_padding`) can keep them.
4. **Subset genes** — default `c("TSC1", "TSC2")`; `genes = NULL` plots the
   whole panel.
5. **Pivot to long form** `(gene, transcrit, exon, patient, rpkm)`; preserve
   exon order along each gene (factor by genomic/exon order, not alphabetical).
6. **Normalize** — per exon, compute the cohort **median** RPKM; derive
   `ratio = rpkm / median`.

## Plot

Two stacked panels sharing the x-axis (exon), faceted by gene:

- **QC panel:** raw RPKM, one line per patient.
- **CNV panel:** `ratio`, one line per patient, with reference lines at
  **1.0** (neutral) and DeCovA thresholds **0.8** (deletion) and **1.2**
  (duplication). Points falling outside the [0.8, 1.2] band are highlighted so
  candidate CNVs stand out.

Rendering:
- **Static:** PNG + PDF (thesis/report quality).
- **Interactive:** `ggplotly()` HTML — hover reveals patient + exon + value;
  zoom into a gene region.

## Structure

A single R script with small, testable functions:

- `load_rpkm(path, sheet = "RPKMs")` → tidy data frame (forward-fill applied).
- `prepare_long(df, genes = c("TSC1","TSC2"), include_padding = FALSE)` →
  long-form data with `ratio` column.
- `plot_coverage_curve(long_df)` → ggplot object (the two-panel figure).
- `main()` / CLI accepting `--xlsx`, `--genes`, `--out`, writing PNG + PDF +
  HTML to an output dir.

Output defaults to a project-local directory (e.g. `output/` under the working
project), never overwriting source data under `<run_dir>/`.

## Out of Scope

- Statistical CNV **calling** with confidence/segmentation — deferred to a
  possible future Bioconductor step (ExomeDepth / cn.mops / CODEX2).
- BAM/coverage computation — already done; we consume the RPKM matrix.
- Re-running DeCovA itself.

## Open Considerations (resolved defaults)

- **Center statistic:** median (matches DeCovA default `med`).
- **Thresholds:** 0.8 / 1.2 (DeCovA defaults `level_del` / `level_dup`).
- **Default gene set:** TSC1 + TSC2, panel available via parameter.
