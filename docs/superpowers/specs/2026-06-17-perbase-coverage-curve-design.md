# Per-base coverage curve (DeCovA-style) — design

**Date:** 2026-06-17
**Status:** approved, pending implementation
**Sibling spec:** `2026-06-17-rpkm-exon-curve-design.md`

## Purpose

A base-resolution coverage curve for the TSC panel across all samples of a
sequencing run, to spot CNVs (deletions / duplications) by eye as coverage
shape. It is the per-base counterpart of `rpkm_exon_curve.R`: overlaid lines,
one color per patient, faceted by gene — but at single-base resolution on a
concatenated-exonic x-axis, plotting **raw depth**.

The intent is a run-level QC / outlier view: the cohort is shown as overlaid
per-sample lines (color-coded by patient), and a sample whose coverage dives
below or jumps above the others over a stretch of exons is a CNV candidate.

## Scope decisions (settled during brainstorming)

- **Source = mosdepth `per-base.bed.gz`**, not `regions.bed.gz`. `regions.bed.gz`
  is mean-depth-per-exon — the same granularity as the RPKM sheet, so it would
  only reproduce the existing exon-level curve. Per-base is what reveals
  intra-exon dropouts and breakpoints.
- **X-axis = concatenated exonic position** (option C): exons laid side by side,
  each exon's width proportional to its length in nucleotides, introns dropped.
- **Y-axis = raw per-base depth.** Per-sample-median normalization is explicitly
  **deferred** to a later session (see TODO). First version uses actual depth.
- **One figure for all samples**, color-coded by patient (legend maps
  color → sample, exactly as in `rpkm_exon_curve.R`). No outlier
  auto-highlighting.

## Data sources

### Panel BED (built once, cached)

The genomic coordinates of the panel exons live in the KoverGraph per-sample
workbooks:

```
<run_dir>/KoverGraph/DATA/<run>/<SAMPLE>/<GENE>/<transcript>.xlsx
  └─ sheet "All_position": Gene | Transcrit | Exon | Chrom | Position | Couverture
```

Coordinates are identical across all samples (same panel, same reference), so a
**single sample** is read to derive the BED. For each `(Gene, Transcrit, Exon)`,
contiguous `Position` runs are collapsed into intervals, producing
`output/panel.bed`:

```
chrom  start  end  gene  transcript  exon  strand
```

Cached after first build; rebuilt only if absent or `--rebuild-bed` is passed.

### Per-base depth

```
<run_dir>/depth/<SAMPLE>.per-base.bed.gz   # mosdepth RLE: chrom start end depth
```

One bgzipped, tabix-indexed file per sample. The panel region is extracted per
sample with `tabix -R output/panel.bed`, so only panel-overlapping records are
loaded — never the genome-wide file.

## Pipeline

1. **Build / load panel BED** from one sample's KoverGraph xlsx (cached).
2. **Extract + annotate per sample:** `tabix -R panel.bed <sample>.per-base.bed.gz`
   yields panel-overlapping RLE intervals; expand each to per-base depth and tag
   every base with its gene / transcript / exon via interval join against the
   panel BED (`data.table::foverlaps`).
3. **Concatenated exonic coordinate:** order exons by genomic position within
   each gene (respecting strand); assign each exon a cumulative offset equal to
   the summed length of preceding exons. `exonic_x = offset + (pos − exon_start)`.
   Exon widths come out proportional to length automatically.
4. **Assemble** one tidy long table:
   `sample | gene | transcript | exon | genomic_pos | exonic_x | depth`.

## Plot

- One figure: `x = exonic_x`, `y = depth`, one **line per sample colored by
  patient**, faceted per gene (`facet_grid(~ gene)` or `facet_wrap`, free x).
- **Gene-model track** beneath each facet: exon segments drawn along x with
  boundaries and exon labels. Because introns are dropped, this is an exon ruler
  rather than the intron-spanning box model of the Lyon reference figure.
- Title per facet carries gene + transcript; tooltip (interactive) carries
  sample, gene, exon, genomic position, depth.
- **Outputs:** PNG + PDF + interactive HTML via `plotly::ggplotly`, reusing the
  `render_outputs` pattern from `rpkm_exon_curve.R`.

## Output location (timestamped, no-overwrite)

Each run writes into its own dated subfolder so a rerun never erases the
previous one:

```
output/<script>/<run_scope>/<DD_HHMM>/
e.g.  output/perbase_coverage/run_1_TSC1-TSC2/17_1024/
```

- **`<run_scope>`** = run name (inferred from the depth dir, e.g. `run 1`) +
  gene scope (`TSC1-TSC2`, or `panel` for the whole panel), so scopes don't mix.
- **`<DD_HHMM>`** = `format(Sys.time(), "%d_%H%M")` — day + hour-minute, enough to
  avoid overwriting (intentionally short, not for long-term archival sorting).
- A `<DD_HHMM>` subfolder is **always** created under the base directory. `--out`
  overrides the base (`output/<script>/<run_scope>`) but the timestamped subfolder
  still applies, so output is never clobbered.

This convention is **also retrofitted to `rpkm_exon_curve.R`**, whose current
default (`output/rpkm_curve`, fixed) overwrites on every run. Its run name is
inferred from the input xlsx filename (e.g. `run_1_Blanket.xlsx` → `run 1`).

## Structure

New script `perbase_coverage_curve.R`, mirroring `rpkm_exon_curve.R`:

- `build_panel_bed(kover_dir, sample = NULL)` — derive/collapse intervals from
  one sample's xlsx; returns + caches the BED.
- `load_depth(sample_path, panel_bed)` — tabix-extract, RLE-expand, annotate.
- `to_exonic_coords(df, panel_bed)` — compute cumulative exonic x per gene.
- `plot_coverage_curve(long_df)` — overlaid colored lines + gene-model track.
- `render_outputs(p, outdir, basename)` — PNG/PDF/HTML (reuse pattern).
- `parse_args` / `main` — CLI: `--depth-dir`, `--kover-dir`, `--genes`
  (default `TSC1,TSC2`; `all`/`panel` keeps whole panel), `--out`,
  `--rebuild-bed`.

## Environment

Combined `kodex` micromamba env (R + Python). Add **`htslib`** (provides
`tabix`) to `environment.yml`. R packages: existing `readxl`, `dplyr`, `tidyr`,
`ggplot2`, `plotly` plus **`data.table`** for fast bgzip reads and interval
joins.

## Constraints

- Read-only on all `<run_dir>/` inputs; all outputs to project-local
  `output/` (never overwrite source data).
- No hardcoded absolute paths beyond documented CLI defaults.
- No external API calls; no micromamba invocation from hooks/automation.

## Out of scope (deferred)

- Per-sample-median normalization of depth (the planned "fix in another
  session").
- Statistical CNV calling (already tracked in `TODO.md`, Approach C).
- Outlier auto-highlighting / CNV-region annotation lines.
