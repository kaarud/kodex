# Per-base Coverage Curve Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Build `perbase_coverage_curve.R`, a base-resolution coverage curve for the TSC panel across all samples of a run (one colored line per patient), and retrofit `rpkm_exon_curve.R` with a timestamped no-overwrite output scheme.

**Architecture:** A panel BED is derived once from one sample's KoverGraph `All_position` xlsx (cached). For each sample, mosdepth `per-base.bed.gz` is region-extracted with `tabix -R`, RLE-expanded to per-base depth, annotated to exons via a non-equi join, and remapped to a concatenated-exonic x-axis (exon width ∝ length, introns dropped). Pure data transforms are unit-tested with synthetic fixtures; I/O wrappers and plotting are verified by a real-data smoke run. Shared output helpers live in `coverage_common.R`, sourced by both scripts.

**Tech Stack:** R 4.3+ (`readxl`, `data.table`, `dplyr`, `ggplot2`, `plotly`, `htmlwidgets`), `testthat` for tests, `tabix` (htslib) for region extraction. Combined `kodex` micromamba env.

**Spec:** `docs/superpowers/specs/2026-06-17-perbase-coverage-curve-design.md`

## Global Constraints

- **R floor:** `r-base>=4.3` (uses vectorized `sequence(nvec, from=)`, R ≥ 4.0).
- **Read-only inputs:** never write into `<run_dir>/`. `.csi` indexes already exist, so no index is created; if one were missing, only an additive sidecar may be written, never a modification.
- **Outputs:** project-local only, under `output/`. Each run writes to `output/<script>/<run_scope>/<DD_HHMM>/` so reruns never overwrite. Timestamp format is exactly `format(now, "%d_%H%M")`.
- **No hardcoded absolute paths** beyond documented CLI defaults (which use `~`).
- **No external API calls; no micromamba invocation from hooks/automation.**
- **Run commands** use the env binary by full path: `Rscript`. `tabix` resolves via the system PATH (homebrew) or the env once `htslib` is added. All commands run from the repo root `~/projets/kodex`.
- **Sample id** = the depth filename with `.per-base.bed.gz` stripped (e.g. `100003_S7`).
- **Coordinates:** KoverGraph `Position` is 1-based inclusive; mosdepth BED is 0-based half-open. Per-base genomic position = `start + 1 .. end`. Panel intervals are stored 1-based inclusive; the on-disk BED for `tabix` is 0-based (`start - 1`).

---

## File Structure

- **Create `coverage_common.R`** — `timestamped_outdir()`, `render_outputs()`. Sourced by both scripts (guarded by `exists()`).
- **Create `perbase_coverage_curve.R`** — the new script (pure transforms + I/O wrappers + CLI).
- **Modify `rpkm_exon_curve.R`** — source `coverage_common.R`, drop its inline `render_outputs`, use timestamped output + run-name inference.
- **Modify `environment.yml`** — add `r-data.table`, `r-testthat`, `htslib`.
- **Create `tests/test_coverage_common.R`** — testthat for the shared helper.
- **Create `tests/test_perbase_coverage.R`** — testthat for the pure transforms.

---

## Task 1: Environment + shared timestamped-output helper

**Files:**
- Modify: `environment.yml`
- Create: `coverage_common.R`
- Create: `tests/test_coverage_common.R`

**Interfaces:**
- Produces: `timestamped_outdir(base, now = Sys.time())` → character path `file.path(base, format(now, "%d_%H%M"))`, directory created, path returned.

- [ ] **Step 1: Add R/htslib deps to `environment.yml`**

Append these three lines under the existing R block (after `- pandoc`):

```yaml
  - r-data.table
  - r-testthat
  - htslib   # provides tabix for per-base.bed.gz region extraction
```

- [ ] **Step 2: Update the env**

Run: `micromamba env update -n kodex -f environment.yml -y`
Expected: completes; `data.table`, `testthat`, and `tabix` become available. (Interactive only — never in a hook.)

- [ ] **Step 3: Write the failing test**

Create `tests/test_coverage_common.R`:

```r
library(testthat)
source("coverage_common.R")

test_that("timestamped_outdir builds DD_HHMM under base and creates it", {
  base <- file.path(tempdir(), "tsout_test")
  on.exit(unlink(base, recursive = TRUE), add = TRUE)
  t <- as.POSIXct("2026-06-17 10:24:00", tz = "UTC")
  out <- timestamped_outdir(base, now = t)
  expect_equal(basename(out), "17_1024")
  expect_true(dir.exists(out))
})
```

- [ ] **Step 4: Run test to verify it fails**

Run: `Rscript -e 'testthat::test_file("tests/test_coverage_common.R")'`
Expected: FAIL — `cannot open file 'coverage_common.R'` (file does not exist yet).

- [ ] **Step 5: Create `coverage_common.R` with both helpers**

```r
# Shared helpers for the coverage-curve scripts (rpkm_exon_curve.R,
# perbase_coverage_curve.R).

suppressPackageStartupMessages({
  library(ggplot2)
  library(plotly)
})

# Timestamped output dir (DD_HHMM) under `base`, so reruns never overwrite.
timestamped_outdir <- function(base, now = Sys.time()) {
  dir <- file.path(base, format(now, "%d_%H%M"))
  if (!dir.exists(dir)) dir.create(dir, recursive = TRUE)
  dir
}

# Render a ggplot to PNG + PDF + interactive HTML. Returns named paths.
render_outputs <- function(p, outdir, basename) {
  if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)
  png_path  <- file.path(outdir, paste0(basename, ".png"))
  pdf_path  <- file.path(outdir, paste0(basename, ".pdf"))
  html_path <- file.path(outdir, paste0(basename, ".html"))

  ggplot2::ggsave(png_path, p, width = 14, height = 8, dpi = 150)
  ggplot2::ggsave(pdf_path, p, width = 14, height = 8)

  ip <- plotly::ggplotly(p, tooltip = "text")
  saved <- tryCatch({
    htmlwidgets::saveWidget(ip, html_path, selfcontained = TRUE); TRUE
  }, error = function(e) {
    htmlwidgets::saveWidget(ip, html_path, selfcontained = FALSE); FALSE
  })
  if (!saved) {
    message("Interactive HTML saved non-self-contained (pandoc unavailable); ",
            "an accompanying *_files/ directory is required to open it.")
  }
  invisible(c(png = png_path, pdf = pdf_path, html = html_path))
}
```

- [ ] **Step 6: Run test to verify it passes**

Run: `Rscript -e 'testthat::test_file("tests/test_coverage_common.R")'`
Expected: PASS (1 passed).

- [ ] **Step 7: Commit**

```bash
git add environment.yml coverage_common.R tests/test_coverage_common.R
git commit -m "feat: shared coverage_common.R + timestamped_outdir helper"
```

---

## Task 2: Retrofit `rpkm_exon_curve.R` to timestamped output

**Files:**
- Modify: `rpkm_exon_curve.R` (remove inline `render_outputs` at lines 128-152; adjust `main` at lines 175-194; add a `script_dir`/`source` block + `infer_run` near the top)

**Interfaces:**
- Consumes: `timestamped_outdir`, `render_outputs` from `coverage_common.R`.
- Produces: `infer_run(path)` → character; extracts a run id like `run 1` from an input filename (`^([A-Za-z]+_[0-9]+)`), else the basename without extension.

- [ ] **Step 1: Add a failing test for `infer_run`**

Create `tests/test_rpkm_infer_run.R`:

```r
library(testthat)
source("coverage_common.R")
source("rpkm_exon_curve.R")

test_that("infer_run pulls the run id from a KoverGraph xlsx name", {
  expect_equal(infer_run("<run_dir>/KoverGraph/run_1_Blanket.xlsx"), "run 1")
  expect_equal(infer_run("/x/run_2_Singularity_Blanket.xlsx"), "run 2")
})

test_that("infer_run falls back to the bare basename", {
  expect_equal(infer_run("/x/cohort.xlsx"), "cohort")
})
```

- [ ] **Step 2: Run test to verify it fails**

Run: `Rscript -e 'testthat::test_file("tests/test_rpkm_infer_run.R")'`
Expected: FAIL — `could not find function "infer_run"`.

- [ ] **Step 3: Add the source block + `infer_run` near the top of `rpkm_exon_curve.R`**

Insert immediately after the `suppressPackageStartupMessages({ ... })` block (currently ending line 18):

```r
script_dir <- function() {
  ca <- commandArgs(FALSE)
  f <- sub("^--file=", "", ca[grepl("^--file=", ca)])
  if (length(f)) dirname(normalizePath(f)) else getwd()
}
if (!exists("timestamped_outdir")) source(file.path(script_dir(), "coverage_common.R"))

# Infer a run id (e.g. run 1) from an input filename; fallback = basename.
infer_run <- function(path) {
  stem <- sub("\\.[^.]+$", "", basename(path))
  m <- regmatches(stem, regexpr("^[A-Za-z]+_[0-9]+", stem))
  if (length(m) == 1 && nzchar(m)) m else stem
}
```

- [ ] **Step 4: Delete the inline `render_outputs` definition**

Remove the entire `render_outputs <- function(p, outdir, basename = "rpkm_exon_curve") { ... }` block (currently lines 128-152). It now comes from `coverage_common.R`.

- [ ] **Step 5: Change the `out` default and timestamped wiring in `main`**

In `parse_args`, change the `out` default (currently `out = "output/rpkm_curve"`) to keep the same base:

```r
    out = "output/rpkm_curve",
```

(unchanged value — it is now treated as the `<script>` base, not the final dir).

In `main`, replace the two output lines (currently `paths <- render_outputs(p, path.expand(opts$out))` and the preceding `p <- plot_coverage_curve(long)`) so the curve is built, then a run-scoped timestamped dir is derived:

```r
  p <- plot_coverage_curve(long)

  scope <- if (is.null(genes)) "panel" else paste(genes, collapse = "-")
  base <- file.path(path.expand(opts$out), paste0(infer_run(path), "_", scope))
  outdir <- timestamped_outdir(base)
  paths <- render_outputs(p, outdir, "rpkm_exon_curve")
```

- [ ] **Step 6: Run the infer_run test to verify it passes**

Run: `Rscript -e 'testthat::test_file("tests/test_rpkm_infer_run.R")'`
Expected: PASS (3 passed).

- [ ] **Step 7: Smoke-test the retrofitted script on real data**

Run: `Rscript rpkm_exon_curve.R --xlsx <run_dir>/KoverGraph/run_1_Blanket.xlsx`
Expected: prints `Wrote:` with three paths under `output/rpkm_curve/run_1_TSC1-TSC2/<DD_HHMM>/`. Confirm the dir exists and the previous fixed `output/rpkm_curve/*.png` is no longer the write target:

```bash
ls -R output/rpkm_curve/
```

- [ ] **Step 8: Commit**

```bash
git add rpkm_exon_curve.R tests/test_rpkm_infer_run.R
git commit -m "feat: timestamped no-overwrite output for rpkm_exon_curve.R"
```

---

## Task 3: Panel BED — `collapse_positions_to_intervals`

**Files:**
- Create: `perbase_coverage_curve.R` (header + library block + source block + this function)
- Create: `tests/test_perbase_coverage.R`

**Interfaces:**
- Produces: `collapse_positions_to_intervals(df)` — `df` has columns `Gene, Transcrit, Exon, Chrom, Position` (1-based int). Returns a `data.frame` with columns `gene, transcript, exon, exon_num, chrom, start, end` (1-based inclusive), keeping only core exons matching `^E[0-9]+$`, ordered by `gene` then `exon_num`.

- [ ] **Step 1: Write the failing test**

Create `tests/test_perbase_coverage.R`:

```r
library(testthat)
library(data.table)
source("coverage_common.R")
source("perbase_coverage_curve.R")

test_that("collapse_positions_to_intervals collapses runs and drops padding", {
  df <- data.frame(
    Gene = "TSC1", Transcrit = "NM_000368",
    Exon = c("E01","E01","E01","E02","E02","E01_lpad"),
    Chrom = "chr9",
    Position = c(100L,101L,102L, 200L,201L, 50L),
    stringsAsFactors = FALSE)
  iv <- collapse_positions_to_intervals(df)
  expect_equal(nrow(iv), 2L)                       # padding row dropped
  expect_equal(iv$exon, c("E01","E02"))            # ordered by exon_num
  expect_equal(iv$start, c(100L, 200L))
  expect_equal(iv$end,   c(102L, 201L))
  expect_equal(iv$exon_num, c(1L, 2L))
  expect_equal(unique(iv$chrom), "chr9")
})
```

- [ ] **Step 2: Run test to verify it fails**

Run: `Rscript -e 'testthat::test_file("tests/test_perbase_coverage.R")'`
Expected: FAIL — `cannot open file 'perbase_coverage_curve.R'`.

- [ ] **Step 3: Create `perbase_coverage_curve.R` header + the function**

```r
#!/usr/bin/env Rscript

# Per-base coverage curve (DeCovA-style) for the TSC panel.
#
# Reads mosdepth per-base depth for every sample of a run, annotates each base
# to its panel exon, remaps onto a concatenated-exonic x-axis (exon width
# proportional to length, introns dropped) and draws one colored line per
# patient, faceted by gene. Raw depth (normalization deferred).
#
# Spec: docs/superpowers/specs/2026-06-17-perbase-coverage-curve-design.md

suppressPackageStartupMessages({
  library(readxl)
  library(data.table)
  library(dplyr)
  library(ggplot2)
  library(plotly)
})

script_dir <- function() {
  ca <- commandArgs(FALSE)
  f <- sub("^--file=", "", ca[grepl("^--file=", ca)])
  if (length(f)) dirname(normalizePath(f)) else getwd()
}
if (!exists("timestamped_outdir")) source(file.path(script_dir(), "coverage_common.R"))

EXON_RE <- "^E[0-9]+$"

# --- pure transforms ---------------------------------------------------------

# Collapse per-base KoverGraph rows into one interval per (gene,transcript,exon).
collapse_positions_to_intervals <- function(df) {
  df <- df[grepl(EXON_RE, df$Exon), , drop = FALSE]
  if (nrow(df) == 0) stop("No core exon rows (E<number>) found in KoverGraph data.")
  dt <- data.table::as.data.table(df)
  out <- dt[, .(chrom = unique(Chrom),
                start = min(Position), end = max(Position)),
            by = .(gene = Gene, transcript = Transcrit, exon = Exon)]
  out[, exon_num := as.integer(sub("^E", "", exon))]
  data.table::setorder(out, gene, exon_num)
  as.data.frame(out)
}
```

- [ ] **Step 4: Run test to verify it passes**

Run: `Rscript -e 'testthat::test_file("tests/test_perbase_coverage.R")'`
Expected: PASS (1 passed).

- [ ] **Step 5: Commit**

```bash
git add perbase_coverage_curve.R tests/test_perbase_coverage.R
git commit -m "feat: collapse KoverGraph positions to panel exon intervals"
```

---

## Task 4: RLE expansion — `expand_rle`

**Files:**
- Modify: `perbase_coverage_curve.R` (add function after `collapse_positions_to_intervals`)
- Modify: `tests/test_perbase_coverage.R` (add test)

**Interfaces:**
- Produces: `expand_rle(dt)` — `dt` has columns `chrom, start, end, depth` (0-based half-open mosdepth intervals). Returns a `data.table` with columns `chrom, pos (1-based int), depth`, one row per covered base. Zero-length intervals contribute no rows.

- [ ] **Step 1: Write the failing test**

Add to `tests/test_perbase_coverage.R`:

```r
test_that("expand_rle expands 0-based half-open intervals to 1-based per-base", {
  dt <- data.table(chrom = c("chr9","chr9"),
                   start = c(100L, 103L),   # 0-based
                   end   = c(103L, 103L),   # half-open; 2nd is zero-length
                   depth = c(7L, 99L))
  out <- expand_rle(dt)
  expect_equal(out$pos, c(101L, 102L, 103L))   # start+1 .. end
  expect_equal(out$depth, c(7L, 7L, 7L))
  expect_equal(unique(out$chrom), "chr9")
})
```

- [ ] **Step 2: Run test to verify it fails**

Run: `Rscript -e 'testthat::test_file("tests/test_perbase_coverage.R")'`
Expected: FAIL — `could not find function "expand_rle"`.

- [ ] **Step 3: Add the function**

```r
# Expand mosdepth RLE intervals (0-based half-open) to per-base depth.
expand_rle <- function(dt) {
  dt <- data.table::as.data.table(dt)
  lens <- dt$end - dt$start
  data.table::data.table(
    chrom = rep(dt$chrom, lens),
    pos   = sequence(lens, from = dt$start + 1L),  # 1-based genomic position
    depth = rep(dt$depth, lens)
  )
}
```

- [ ] **Step 4: Run test to verify it passes**

Run: `Rscript -e 'testthat::test_file("tests/test_perbase_coverage.R")'`
Expected: PASS (2 passed).

- [ ] **Step 5: Commit**

```bash
git add perbase_coverage_curve.R tests/test_perbase_coverage.R
git commit -m "feat: expand mosdepth RLE intervals to per-base depth"
```

---

## Task 5: Exon annotation — `assign_exons`

**Files:**
- Modify: `perbase_coverage_curve.R` (add function after `expand_rle`)
- Modify: `tests/test_perbase_coverage.R` (add test)

**Interfaces:**
- Consumes: per-base long table from `expand_rle` (`chrom, pos, depth`); intervals from `collapse_positions_to_intervals` (`gene, transcript, exon, exon_num, chrom, start, end`).
- Produces: `assign_exons(depth_long, intervals)` → `data.table` with columns `chrom, pos, depth, gene, transcript, exon, exon_num, exon_start, exon_end`, keeping only positions falling inside an exon interval (1-based inclusive).

- [ ] **Step 1: Write the failing test**

Add to `tests/test_perbase_coverage.R`:

```r
test_that("assign_exons keeps only in-exon positions and tags them", {
  depth_long <- data.table(chrom = "chr9", pos = 99:103, depth = 5L)
  intervals <- data.frame(gene = "TSC1", transcript = "NM_000368",
                          exon = "E01", exon_num = 1L, chrom = "chr9",
                          start = 100L, end = 102L, stringsAsFactors = FALSE)
  out <- assign_exons(depth_long, intervals)
  expect_equal(sort(out$pos), c(100L, 101L, 102L))   # 99 and 103 dropped
  expect_true(all(out$exon == "E01"))
  expect_true(all(out$gene == "TSC1"))
  expect_equal(unique(out$exon_start), 100L)
  expect_equal(unique(out$exon_end), 102L)
})
```

- [ ] **Step 2: Run test to verify it fails**

Run: `Rscript -e 'testthat::test_file("tests/test_perbase_coverage.R")'`
Expected: FAIL — `could not find function "assign_exons"`.

- [ ] **Step 3: Add the function (non-equi join)**

```r
# Annotate each per-base position with its exon; drop positions outside the panel.
assign_exons <- function(depth_long, intervals) {
  d  <- data.table::as.data.table(depth_long)
  iv <- data.table::as.data.table(intervals)
  iv[d, on = .(chrom, start <= pos, end >= pos), nomatch = NULL,
     .(chrom, pos = i.pos, depth = i.depth,
       gene, transcript, exon, exon_num,
       exon_start = x.start, exon_end = x.end)]
}
```

- [ ] **Step 4: Run test to verify it passes**

Run: `Rscript -e 'testthat::test_file("tests/test_perbase_coverage.R")'`
Expected: PASS (3 passed).

- [ ] **Step 5: Commit**

```bash
git add perbase_coverage_curve.R tests/test_perbase_coverage.R
git commit -m "feat: annotate per-base depth to panel exons"
```

---

## Task 6: Concatenated-exonic coordinates — `exon_layout` + `to_exonic_coords`

**Files:**
- Modify: `perbase_coverage_curve.R` (add both functions after `assign_exons`)
- Modify: `tests/test_perbase_coverage.R` (add test)

**Interfaces:**
- Consumes: intervals (`gene, transcript, exon, exon_num, chrom, start, end`); annotated per-base table from `assign_exons` (must have `gene, exon, pos`).
- Produces:
  - `exon_layout(intervals)` → `data.frame` `gene, transcript, exon, exon_num, chrom, gstart, gend, len, x_start, x_end, x_mid`. `len = end-start+1`; `x_start` = cumulative summed length of preceding exons within the gene (ordered by `exon_num`); `x_end = x_start + len`; `x_mid = x_start + len/2`.
  - `to_exonic_coords(df, layout)` → input `df` joined to `layout` on `gene, exon`, adds `exonic_x = x_start + (pos - gstart)`.

- [ ] **Step 1: Write the failing test**

Add to `tests/test_perbase_coverage.R`:

```r
test_that("exon_layout gives proportional widths and cumulative offsets", {
  intervals <- data.frame(
    gene = "TSC1", transcript = "NM_000368",
    exon = c("E01","E02"), exon_num = c(1L,2L), chrom = "chr9",
    start = c(100L, 200L), end = c(102L, 201L), stringsAsFactors = FALSE)
  L <- exon_layout(intervals)
  expect_equal(L$len, c(3L, 2L))           # 100-102 -> 3, 200-201 -> 2
  expect_equal(L$x_start, c(0L, 3L))        # E02 starts after E01's 3 bases
  expect_equal(L$x_end, c(3L, 5L))
})

test_that("to_exonic_coords maps genomic pos to concatenated exonic x", {
  intervals <- data.frame(
    gene = "TSC1", transcript = "NM_000368",
    exon = c("E01","E02"), exon_num = c(1L,2L), chrom = "chr9",
    start = c(100L, 200L), end = c(102L, 201L), stringsAsFactors = FALSE)
  L <- exon_layout(intervals)
  df <- data.frame(gene = "TSC1", exon = c("E01","E02"),
                   pos = c(100L, 200L), depth = c(5L, 9L),
                   stringsAsFactors = FALSE)
  out <- to_exonic_coords(df, L)
  x <- out$exonic_x[match(c(100L,200L), out$pos)]
  expect_equal(x, c(0L, 3L))   # first base of E01 -> 0; first base of E02 -> 3
})
```

- [ ] **Step 2: Run test to verify it fails**

Run: `Rscript -e 'testthat::test_file("tests/test_perbase_coverage.R")'`
Expected: FAIL — `could not find function "exon_layout"`.

- [ ] **Step 3: Add both functions**

```r
# Per-exon layout on the concatenated-exonic axis (width proportional to length).
exon_layout <- function(intervals) {
  iv <- data.table::as.data.table(intervals)
  data.table::setorder(iv, gene, exon_num)
  iv[, len := end - start + 1L]
  iv[, x_start := cumsum(c(0L, utils::head(len, -1L))), by = gene]
  iv[, x_end := x_start + len]
  iv[, x_mid := x_start + len / 2]
  as.data.frame(iv[, .(gene, transcript, exon, exon_num, chrom,
                       gstart = start, gend = end, len, x_start, x_end, x_mid)])
}

# Map annotated per-base rows onto the concatenated-exonic x-axis.
to_exonic_coords <- function(df, layout) {
  d <- data.table::as.data.table(df)
  L <- data.table::as.data.table(layout)[, .(gene, exon, x_start, gstart)]
  m <- merge(d, L, by = c("gene", "exon"))
  m[, exonic_x := x_start + (pos - gstart)]
  as.data.frame(m)
}
```

- [ ] **Step 4: Run test to verify it passes**

Run: `Rscript -e 'testthat::test_file("tests/test_perbase_coverage.R")'`
Expected: PASS (5 passed).

- [ ] **Step 5: Commit**

```bash
git add perbase_coverage_curve.R tests/test_perbase_coverage.R
git commit -m "feat: concatenated-exonic coordinate layout and mapping"
```

---

## Task 7: I/O wrappers, plot, CLI, and real-data smoke

**Files:**
- Modify: `perbase_coverage_curve.R` (add `infer_run_scope`, `write_panel_bed`, `build_panel_bed`, `ensure_tabix_index`, `load_depth`, `load_all_samples`, `plot_coverage_curve`, `parse_args`, `main`, and the `if (sys.nframe() == 0) main()` guard)
- Modify: `tests/test_perbase_coverage.R` (add `infer_run_scope` test)

**Interfaces:**
- Consumes: every function from Tasks 3-6; `timestamped_outdir`, `render_outputs` from `coverage_common.R`.
- Produces:
  - `infer_run_scope(depth_dir, genes)` → `"<run>_<scope>"`, run = `basename(dirname(depth_dir))`, scope = `paste(genes, collapse="-")` or `"panel"` when `genes` is `NULL`.
  - `write_panel_bed(intervals, path)` → writes a sorted 0-based BED (`chrom, start-1, end`), returns `path`.
  - `build_panel_bed(kover_dir, bed_path, sample=NULL)` → reads ALL gene xlsx of one sample, returns intervals `data.frame` and writes the BED (always whole-panel; gene filtering happens later).
  - `load_depth(sample_path, bed_path, intervals, sample_id)` → annotated per-base `data.frame` with a `sample` column.
  - `load_all_samples(depth_dir, bed_path, intervals)` → row-bound annotated table for all `*.per-base.bed.gz`.
  - `plot_coverage_curve(long_df, layout)` → ggplot object.

- [ ] **Step 1: Write the failing test for `infer_run_scope`**

Add to `tests/test_perbase_coverage.R`:

```r
test_that("infer_run_scope combines run dir and gene scope", {
  expect_equal(infer_run_scope("/x/run 1/depth", c("TSC1","TSC2")),
               "run_1_TSC1-TSC2")
  expect_equal(infer_run_scope("/x/run 1/depth", NULL), "run_1_panel")
})
```

- [ ] **Step 2: Run test to verify it fails**

Run: `Rscript -e 'testthat::test_file("tests/test_perbase_coverage.R")'`
Expected: FAIL — `could not find function "infer_run_scope"`.

- [ ] **Step 3: Add the I/O wrappers, plot, CLI, and `main`**

Append to `perbase_coverage_curve.R`:

```r
# --- run scope / IO ----------------------------------------------------------

infer_run_scope <- function(depth_dir, genes) {
  run <- basename(dirname(normalizePath(depth_dir, mustWork = FALSE)))
  scope <- if (is.null(genes)) "panel" else paste(genes, collapse = "-")
  paste0(run, "_", scope)
}

write_panel_bed <- function(intervals, path) {
  bed <- data.table::as.data.table(intervals)[, .(chrom, start = start - 1L, end)]
  data.table::setorder(bed, chrom, start)
  data.table::fwrite(bed, path, sep = "\t", col.names = FALSE)
  path
}

# Derive the panel BED from ONE sample's KoverGraph xlsx (coords are identical
# across samples). Always builds the whole panel; plotting filters genes later.
build_panel_bed <- function(kover_dir, bed_path, sample = NULL) {
  samples <- list.dirs(kover_dir, recursive = FALSE)
  if (length(samples) == 0) stop("No sample dirs under ", kover_dir)
  s <- if (is.null(sample)) samples[1] else file.path(kover_dir, sample)
  gene_dirs <- list.dirs(s, recursive = FALSE)
  rows <- lapply(gene_dirs, function(gd) {
    xlsx <- list.files(gd, pattern = "\\.xlsx$", full.names = TRUE)
    if (length(xlsx) == 0) return(NULL)
    data.table::rbindlist(lapply(xlsx, function(f)
      as.data.frame(readxl::read_excel(f, sheet = "All_position"))),
      use.names = TRUE, fill = TRUE)
  })
  df <- data.table::rbindlist(rows, use.names = TRUE, fill = TRUE)
  intervals <- collapse_positions_to_intervals(as.data.frame(df))
  write_panel_bed(intervals, bed_path)
  intervals
}

ensure_tabix_index <- function(gz) {
  if (!file.exists(paste0(gz, ".tbi")) && !file.exists(paste0(gz, ".csi"))) {
    system2("tabix", c("-p", "bed", shQuote(gz)))
  }
  invisible(gz)
}

load_depth <- function(sample_path, bed_path, intervals, sample_id) {
  ensure_tabix_index(sample_path)
  raw <- data.table::fread(
    cmd = paste("tabix", "-R", shQuote(bed_path), shQuote(sample_path)),
    header = FALSE, col.names = c("chrom", "start", "end", "depth"))
  if (nrow(raw) == 0) stop("No depth records for ", sample_id)
  ann <- assign_exons(expand_rle(raw), intervals)
  ann[, sample := sample_id]
  as.data.frame(ann)
}

load_all_samples <- function(depth_dir, bed_path, intervals) {
  files <- list.files(depth_dir, pattern = "\\.per-base\\.bed\\.gz$",
                      full.names = TRUE)
  if (length(files) == 0) stop("No *.per-base.bed.gz under ", depth_dir)
  dfs <- lapply(files, function(f) {
    sid <- sub("\\.per-base\\.bed\\.gz$", "", basename(f))
    load_depth(f, bed_path, intervals, sid)
  })
  as.data.frame(data.table::rbindlist(dfs, use.names = TRUE))
}

# --- plot --------------------------------------------------------------------

plot_coverage_curve <- function(long_df, layout) {
  ld <- to_exonic_coords(long_df, layout)
  L <- data.table::as.data.table(layout)
  track_y <- -0.06 * max(ld$depth, na.rm = TRUE)

  ggplot(ld, aes(x = exonic_x, y = depth, group = sample, color = sample,
                 text = paste0("sample: ", sample, "<br>gene: ", gene,
                               "<br>exon: ", exon, "<br>pos: ", pos,
                               "<br>depth: ", depth))) +
    geom_line(linewidth = 0.25, alpha = 0.8) +
    geom_segment(data = L, inherit.aes = FALSE,
                 aes(x = x_start, xend = x_end, y = track_y, yend = track_y),
                 linewidth = 3, color = "grey70") +
    geom_text(data = L, inherit.aes = FALSE,
              aes(x = x_mid, y = track_y, label = exon),
              size = 1.8, vjust = 1.6) +
    facet_wrap(~ gene, scales = "free_x") +
    labs(x = "Exonic position (concatenated, width proportional to exon length)",
         y = "Depth", color = "Sample",
         title = "Per-base coverage across cohort") +
    theme_bw(base_size = 9) +
    theme(legend.position = "right",
          legend.key.size = unit(0.35, "cm"),
          legend.text = element_text(size = 6))
}

# --- CLI ---------------------------------------------------------------------

parse_args <- function(args) {
  opts <- list(
    depth_dir   = "<run_dir>/depth",
    kover_dir   = "<run_dir>/KoverGraph/DATA/<run>",
    genes       = "TSC1,TSC2",
    out         = "output/perbase_coverage",
    bed         = "output/perbase_coverage/panel.bed",
    bed_sample  = "",
    rebuild_bed = "false"
  )
  i <- 1
  while (i <= length(args)) {
    key <- sub("^--", "", args[i])
    if (key %in% names(opts) && i + 1 <= length(args)) {
      opts[[key]] <- args[i + 1]; i <- i + 2
    } else stop("Unknown or incomplete argument: ", args[i])
  }
  opts
}

main <- function() {
  opts  <- parse_args(commandArgs(trailingOnly = TRUE))
  genes <- if (tolower(opts$genes) %in% c("all", "panel", "")) NULL
           else trimws(strsplit(opts$genes, ",")[[1]])
  kover <- path.expand(opts$kover_dir)
  depth <- path.expand(opts$depth_dir)
  bed   <- path.expand(opts$bed)
  rds   <- sub("\\.bed$", ".intervals.rds", bed)
  if (!dir.exists(dirname(bed))) dir.create(dirname(bed), recursive = TRUE)

  rebuild    <- tolower(opts$rebuild_bed) %in% c("true", "yes", "1")
  bed_sample <- if (nzchar(opts$bed_sample)) opts$bed_sample else NULL
  if (rebuild || !file.exists(bed) || !file.exists(rds)) {
    message("Building panel BED from ", kover)
    intervals <- build_panel_bed(kover, bed, sample = bed_sample)
    saveRDS(intervals, rds)
  } else {
    intervals <- readRDS(rds)
  }
  if (!is.null(genes)) intervals <- intervals[intervals$gene %in% genes, , drop = FALSE]
  if (nrow(intervals) == 0) stop("No exons left after gene filter: ", opts$genes)

  layout <- exon_layout(intervals)
  message("Loading depth from ", depth)
  long <- load_all_samples(depth, bed, intervals)
  message("Samples: ", dplyr::n_distinct(long$sample),
          " | genes: ", paste(unique(intervals$gene), collapse = ", "),
          " | exons: ", nrow(intervals))

  p <- plot_coverage_curve(long, layout)
  base <- file.path(opts$out, infer_run_scope(depth, genes))
  outdir <- timestamped_outdir(base)
  paths <- render_outputs(p, outdir, "perbase_coverage_curve")
  message("Wrote:\n  ", paste(paths, collapse = "\n  "))
}

if (sys.nframe() == 0) main()
```

- [ ] **Step 4: Run the unit suite to verify it passes**

Run: `Rscript -e 'testthat::test_file("tests/test_perbase_coverage.R")'`
Expected: PASS (6 passed).

- [ ] **Step 5: Real-data smoke run (TSC1+TSC2)**

Run: `Rscript perbase_coverage_curve.R`
Expected: builds `output/perbase_coverage/panel.bed` (+ `.intervals.rds`) on first run, then prints `Samples: 22 | genes: TSC1, TSC2 | exons: ...` and `Wrote:` three paths under `output/perbase_coverage/run_1_TSC1-TSC2/<DD_HHMM>/`. Open the PNG and confirm: one colored line per sample, two facets (TSC1, TSC2), exon ruler with proportional widths beneath each facet.

```bash
ls -R output/perbase_coverage/
```

- [ ] **Step 6: Coordinate-concordance spot check**

Confirm the 1-based/0-based mapping is correct by comparing one position's depth against the KoverGraph `Couverture` for the same `Position`. In R:

```r
Rscript -e '
source("coverage_common.R"); source("perbase_coverage_curve.R")
iv <- readRDS("output/perbase_coverage/panel.intervals.rds")
iv <- iv[iv$gene=="TSC1",]
d <- load_depth("<run_dir>/depth/100003_S7.per-base.bed.gz",
                "output/perbase_coverage/panel.bed", iv, "100003_S7")
print(head(d[order(d$pos), c("pos","depth","exon")], 3))
'
```

Then open `<run_dir>/KoverGraph/DATA/<run>/100003_S7/TSC1/NM_000368.xlsx` (`All_position`) and verify the `Couverture` at those exact `Position` values matches `depth`. If off by one, the bug is in `expand_rle` (`from = start + 1L`) — fix and re-run Task 4's test.

- [ ] **Step 7: Commit**

```bash
git add perbase_coverage_curve.R tests/test_perbase_coverage.R
git commit -m "feat: per-base coverage curve script (IO, plot, CLI)"
```

---

## Self-Review (completed during planning)

- **Spec coverage:** source = per-base.bed.gz (Tasks 4-5,7); panel BED from one sample's KoverGraph xlsx, cached (Tasks 3,7); concatenated-exonic x, proportional widths (Task 6); raw depth y (Task 7 plot); one figure, all samples, color per patient (Task 7 plot); gene-model ruler (Task 7 plot); PNG+PDF+HTML (Task 1 `render_outputs`); timestamped no-overwrite output + rpkm retrofit (Tasks 1,2,7); deps htslib/data.table/testthat (Task 1). Deferred items (normalization, CNV calling, outlier highlight) intentionally absent.
- **Coordinate model:** stated once in Global Constraints and exercised by Task 4/Task 6 tests + Task 7 Step 6 spot check.
- **Type consistency:** `intervals` columns (`gene, transcript, exon, exon_num, chrom, start, end`) are produced in Task 3 and consumed unchanged in Tasks 5-7; `layout` columns produced in Task 6 consumed in Task 7 plot; `expand_rle`→`assign_exons`→`to_exonic_coords` column names align across tasks.
- **No tabix index creation** needed for real data (`.csi` present); `ensure_tabix_index` is a guarded no-op there.
