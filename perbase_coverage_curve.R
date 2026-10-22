#!/usr/bin/env Rscript

# Per-base coverage curve (DeCovA-style) for the TSC panel.
#
# Reads mosdepth per-base depth for every sample of a run, annotates each base
# to its panel exon, remaps onto a concatenated-exonic x-axis (exon width
# proportional to length, introns dropped) and draws one colored line per
# patient. Genes are stacked vertically (one row each) by default; --facet none
# aggregates the whole panel onto a single axis (auto for whole-panel runs).
# Default y is depth normalised to each sample's
# median (--normalised); --raw plots absolute depth; --cohort additionally
# divides each position by the cohort median there (DeCovA-style ratio2center),
# flattening shared capture bias so samples cluster and real CNVs stand out.
#
# Three ablation modes take one step out of --cohort, to show what each step
# actually removes: --cohort-only skips the per-sample median; --gc and
# --capture split the cohort reference into its GC-driven and probe-driven
# components (a per-position median cannot separate them on its own, since GC is
# constant at a fixed position) and divide by one component only. --gc and
# --capture need a reference genome (--ref) to compute local GC content.
#
# --smooth trades per-base resolution for readability: an integer is a rolling
# window in bp, confined to each exon; "exon" collapses each exon to its median
# and draws a step function. --smooth-stat picks median (default, preserves the
# front and depth of an abrupt drop) or mean (blurs it into a ramp). Smoothing
# is cosmetic and applied last, it never feeds back into any normalisation.
#
# Panel geometry is a committed, run-independent artifact (panel/panel.tsv).
# Routine runs read it and never touch the KoverGraph dir; pass --rebuild-panel
# (with --kover-dir) only to regenerate it when the sequencing panel changes.
#
# Orientation caveat (v1): exons are laid out in transcript order (by exon_num),
# but WITHIN each exon bases are placed by ascending genomic position. For
# minus-strand genes (e.g. TSC1) this mirrors the within-exon profile relative
# to mRNA 5'->3'. This is harmless for the run-level outlier QC purpose (every
# sample shares the same x-mapping); strand-aware within-exon orientation is
# deferred (see TODO.md).
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

# Normalisation variants. raw/normalised/cohort are the production modes;
# cohortonly/gc/capture are ablations that take one step out of `cohort`.
MODES <- c("normalised", "raw", "cohort", "cohortonly", "gc", "capture")

# Bare CLI flag -> mode. parse_args() maps "-" to "_" first, so --cohort-only
# and --cohortonly both land here.
BARE_MODES <- c(raw = "raw", normalised = "normalised", normalized = "normalised",
                cohort = "cohort", cohort_ref = "cohort",
                cohortonly = "cohortonly", cohort_only = "cohortonly",
                gc = "gc", capture = "capture")

# --- pure transforms ---------------------------------------------------------

# Collapse per-base KoverGraph rows into one interval per (gene,transcript,exon).
# Assumes each exon's core positions form a single contiguous run: the interval
# is min(Position)..max(Position), so a gap within an exon would be bridged.
# Holds for the TSC panel; revisit if a panel has non-contiguous exon blocks.
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

# Annotate each per-base position with its exon; drop positions outside the panel.
assign_exons <- function(depth_long, intervals) {
  d  <- data.table::as.data.table(depth_long)
  iv <- data.table::as.data.table(intervals)
  iv[d, on = .(chrom, start <= pos, end >= pos), nomatch = NULL,
     .(chrom, pos = i.pos, depth = i.depth,
       gene, transcript, exon, exon_num,
       exon_start = x.start, exon_end = x.end)]
}

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

# Lay exons on ONE axis in genomic order, ignoring exon_num.
#
# exon_layout() orders by (gene, exon_num) and restarts the cumsum at each gene,
# which is right for reading a transcript but wrong for a junction view: PKD1 is
# on the minus strand, so its exon numbers run down as the coordinate runs up
# (E46 at chr16:2,089,727, E01 at 2,135,475). Laid out by exon_num next to TSC2
# it would be drawn backwards and the junction would be a fiction. Ordering by
# `start` across every exon of every gene puts adjacent genes where they
# actually adjoin. Same output columns as exon_layout(), so to_exonic_coords()
# and both normalisations apply unchanged.
#
# Introns and intergenic gaps are still dropped: only panel exons are laid down,
# back to back. The 1,117 bp between TSC2 E42 and PKD1 E46 is marked on the
# figure, not plotted (see junction_marker()).
genomic_layout <- function(intervals) {
  iv <- data.table::as.data.table(intervals)
  if (length(unique(iv$chrom)) != 1L)
    stop("genomic_layout needs a single chromosome, got: ",
         paste(sort(unique(iv$chrom)), collapse = ", "))
  data.table::setorder(iv, start)
  iv[, len := end - start + 1L]
  iv[, x_start := cumsum(c(0L, utils::head(len, -1L)))]   # no by = gene
  iv[, x_end := x_start + len]
  iv[, x_mid := x_start + len / 2]
  as.data.frame(iv[, .(gene, transcript, exon, exon_num, chrom,
                       gstart = start, gend = end, len, x_start, x_end, x_mid)])
}

# Where one gene hands over to the next on a genomic layout: the x position of
# the boundary, and the genomic gap it hides. Returns one row per junction.
junction_marker <- function(layout) {
  L <- data.table::as.data.table(layout)
  data.table::setorder(L, x_start)
  g <- L[, .(gstart = min(gstart), gend = max(gend),
             x_lo = min(x_start), x_hi = max(x_end)), by = gene]
  data.table::setorder(g, x_lo)
  if (nrow(g) < 2L) return(data.frame())
  i <- seq_len(nrow(g) - 1L)
  data.frame(x = g$x_hi[i],
             left = g$gene[i], right = g$gene[i + 1L],
             gap = g$gstart[i + 1L] - g$gend[i] - 1L,
             stringsAsFactors = FALSE)
}

# Shift each gene's per-gene x-coords by the cumulative width of the genes
# before it, so all genes lie back-to-back on ONE continuous axis. Used for the
# aggregated (no-facet) whole-panel view; to_exonic_coords then inherits the
# global offset via the layout's x_start, so per-base x becomes global too.
globalize_layout <- function(layout) {
  L <- data.table::as.data.table(layout)
  data.table::setorder(L, gene, exon_num)
  gw <- L[, .(w = max(x_end)), by = gene]
  gw[, off := cumsum(c(0L, utils::head(w, -1L)))]
  L <- merge(L, gw[, .(gene, off)], by = "gene")
  data.table::setorder(L, gene, exon_num)
  L[, `:=`(x_start = x_start + off, x_end = x_end + off, x_mid = x_mid + off)]
  L[, off := NULL]
  as.data.frame(L)
}

# Map annotated per-base rows onto the concatenated-exonic x-axis.
to_exonic_coords <- function(df, layout) {
  d <- data.table::as.data.table(df)
  L <- data.table::as.data.table(layout)[, .(gene, exon, x_start, gstart)]
  m <- merge(d, L, by = c("gene", "exon"))
  m[, exonic_x := x_start + (pos - gstart)]
  as.data.frame(m)
}

# Divide each sample's per-base depth by that sample's own median (over the
# bases being plotted), adding a `y` column. This removes library-size spread
# so curves overlay around 1.0 and real CNVs stand out (del ~0.5, dup ~1.5).
# `depth` is preserved for tooltips.
normalise_by_sample_median <- function(df) {
  dt <- data.table::as.data.table(df)
  dt[, y := depth / stats::median(depth), by = sample]
  as.data.frame(dt)
}

# Cohort-reference normalisation (DeCovA's ratio2center step). After the
# per-sample median divide above has removed library-size spread, divide each
# sample's `y` at a given exonic position by the COHORT median `y` at that same
# position. This cancels the systematic per-base capture/GC signature shared by
# all samples (a hard-to-capture base dips in everyone), so normal samples
# flatten to ~1.0 and only real CNVs deviate (del ~0.5, dup ~1.5). Without it,
# residual per-position bias keeps the curves "wavy" and they don't cluster.
# Expects a `y` column (run normalise_by_sample_median first) and an `exonic_x`
# (run to_exonic_coords first); (gene, exonic_x) is the genome position key.
normalise_by_cohort_reference <- function(df) {
  dt <- data.table::as.data.table(df)
  if (!"y" %in% names(dt))
    stop("normalise_by_cohort_reference needs a `y` column (per-sample median first).")
  if (!"exonic_x" %in% names(dt))
    stop("normalise_by_cohort_reference needs `exonic_x` (to_exonic_coords first).")
  dt[, ref := stats::median(y[y > 0]), by = .(gene, exonic_x)]
  dt[, y := y / ref]
  dt[!is.finite(ref) | ref == 0, y := NA_real_]   # no usable cohort signal here
  dt[, ref := NULL]
  as.data.frame(dt)
}

# --- GC decomposition --------------------------------------------------------

# A GC bin needs at least this many distinct panel positions to contribute a
# point to the fitted GC curve; sparser bins are dropped and interpolated over.
MIN_GC_BIN <- 50L

# Rolling GC fraction over a window centred on each CORE base of `seq`. `seq` is
# the core interval flanked by `half` bases of genomic padding on each side, so
# the returned vector has length nchar(seq) - 2*half. Bases that are not A/C/G/T
# (N, or padding past a contig edge) count in neither numerator nor denominator;
# a window with no defined base yields NA.
rolling_gc <- function(seq, half) {
  b <- strsplit(toupper(seq), "", fixed = TRUE)[[1]]
  n <- length(b)
  if (n < 2L * half + 1L) return(numeric(0))
  cg <- cumsum(c(0L, as.integer(b %in% c("G", "C"))))
  ca <- cumsum(c(0L, as.integer(b %in% c("A", "C", "G", "T"))))
  core <- seq_len(n - 2L * half) + half
  num <- cg[core + half + 1L] - cg[core - half]
  den <- ca[core + half + 1L] - ca[core - half]
  ifelse(den > 0L, num / den, NA_real_)
}

# GC fraction at every panel base, from the reference genome. One `samtools
# faidx -r` call for the whole panel: each exon is fetched with `window %/% 2`
# bases of flanking sequence so the rolling window is defined at exon edges.
# Returns data.frame(chrom, pos, gc), deduplicated (padding of adjacent exons
# can overlap).
gc_by_position <- function(intervals, ref_fasta, window = 101L) {
  half <- as.integer(window) %/% 2L
  iv <- data.table::as.data.table(intervals)
  data.table::setorder(iv, chrom, start)
  lstart <- pmax(1L, iv$start - half)
  regions <- sprintf("%s:%d-%d", iv$chrom, lstart, iv$end + half)

  rf <- tempfile(fileext = ".txt")
  on.exit(unlink(rf), add = TRUE)
  writeLines(regions, rf)
  out <- system2("samtools", c("faidx", shQuote(path.expand(ref_fasta)),
                               "-r", shQuote(rf)), stdout = TRUE)
  hdr <- grepl("^>", out)
  if (sum(hdr) != nrow(iv))
    stop("samtools faidx returned ", sum(hdr), " sequences for ", nrow(iv),
         " exons. Check the reference and the panel chromosome naming.")
  grp <- factor(cumsum(hdr)[!hdr], levels = seq_len(nrow(iv)))
  seqs <- vapply(split(out[!hdr], grp), paste, character(1), collapse = "")

  # Re-pad what samtools clamped at contig edges, so every core base keeps a
  # full window and the returned vector lines up with start:end.
  want <- (iv$end + half) - (iv$start - half) + 1L
  lpad <- lstart - (iv$start - half)
  seqs <- paste0(strrep("N", lpad), seqs,
                 strrep("N", pmax(0L, want - lpad - nchar(seqs))))

  gc <- data.table::rbindlist(lapply(seq_len(nrow(iv)), function(i)
    data.table::data.table(chrom = iv$chrom[i],
                           pos = seq.int(iv$start[i], iv$end[i]),
                           gc = rolling_gc(seqs[i], half))))
  as.data.frame(unique(gc, by = c("chrom", "pos")))
}

# Split the cohort reference into its GC and probe components.
#
# After the per-sample median divide, the cohort median at a position is
# ref(p) ~ G(gc(p)) * P(p): the GC-driven bias and the probe-driven bias, which
# a per-position median cannot tell apart (gc is constant at a fixed position).
# Fitting ref against local GC content estimates G; the residual ref/G is P.
# The fit is run on gene-centred ref, so a tiered panel does not turn the GC
# curve into a between-gene coverage ladder; the gene level stays in P.
# By construction gcurve * probe == ref, so the two components partition the
# cohort reference exactly rather than approximating it twice.
#
# Adds `ref`, `gcurve` and `probe`; does NOT divide `y` — the caller picks the
# divisor according to the mode. Expects `y` (per-sample median first),
# `exonic_x` (to_exonic_coords first) and `chrom`/`pos`.
decompose_cohort_reference <- function(df, gc_tbl) {
  dt <- data.table::as.data.table(df)
  if (!"y" %in% names(dt))
    stop("decompose_cohort_reference needs a `y` column (per-sample median first).")
  if (!"exonic_x" %in% names(dt))
    stop("decompose_cohort_reference needs `exonic_x` (to_exonic_coords first).")
  dt[, ref := stats::median(y[y > 0]), by = .(gene, exonic_x)]

  g <- data.table::as.data.table(gc_tbl)[, .(chrom, pos, gc)]
  dt <- merge(dt, g, by = c("chrom", "pos"), all.x = TRUE)

  # Fit on distinct POSITIONS, not on rows: every position repeats once per
  # sample and would otherwise weight the bin medians by cohort size.
  pos_tbl <- unique(dt[is.finite(ref) & ref > 0 & !is.na(gc),
                       .(gene, exonic_x, ref, gc)])
  # Centre each gene on its own median before fitting. The TSC panel is tiered
  # (TSC1/TSC2/PKD1 sit ~20x above the pathway genes), so pooling raw `ref`
  # across genes makes the bin medians track which tier happens to fall in a GC
  # bin rather than the GC effect itself. The gene level is not discarded, it
  # stays in `probe` where it belongs: it is capture, not GC.
  pos_tbl[, centred := ref / stats::median(ref), by = gene]
  pos_tbl[, bin := round(gc * 100)]
  fit <- pos_tbl[, .(m = stats::median(centred), n = .N), by = bin][n >= MIN_GC_BIN]
  data.table::setorder(fit, bin)
  if (nrow(fit) < 2L)
    stop("Not enough populated GC bins (>= ", MIN_GC_BIN, " positions) to fit a ",
         "GC curve. Widen the gene scope (--genes all) or lower MIN_GC_BIN.")

  dt[, gcurve := stats::approx(fit$bin, fit$m, xout = round(gc * 100), rule = 2)$y]
  dt[, probe := ref / gcurve]
  as.data.frame(dt)
}

# --- run scope / IO ----------------------------------------------------------

infer_run_scope <- function(depth_dir, genes) {
  run <- basename(dirname(normalizePath(depth_dir, mustWork = FALSE)))
  scope <- if (is.null(genes)) "panel" else paste(genes, collapse = "-")
  paste0(run, "_", scope)
}

# Write the 3-column tabix filter BED (0-based start) — a transient intermediate
# derived from the panel geometry, not the source of truth.
write_panel_bed <- function(intervals, path) {
  bed <- data.table::as.data.table(intervals)[, .(chrom, start = start - 1L, end)]
  data.table::setorder(bed, chrom, start)
  data.table::fwrite(bed, path, sep = "\t", col.names = FALSE)
  path
}

# Persist / load the annotated panel geometry (the real panel definition:
# one row per exon with gene/transcript/exon/coords). This is the committed,
# run-independent artifact; routine runs read it and never touch KoverGraph.
PANEL_COLS <- c("chrom", "start", "end", "gene", "transcript", "exon", "exon_num")

write_panel <- function(intervals, path) {
  dt <- data.table::as.data.table(intervals)[, ..PANEL_COLS]
  data.table::setorder(dt, gene, exon_num)
  data.table::fwrite(dt, path, sep = "\t")
  path
}

read_panel <- function(path) {
  dt <- data.table::fread(path, sep = "\t", colClasses = list(
    character = c("chrom", "gene", "transcript", "exon"),
    integer   = c("start", "end", "exon_num")))
  as.data.frame(dt)
}

# Derive the panel geometry from ONE sample's KoverGraph xlsx (coords are
# identical across samples). Always builds the whole panel; gene filtering
# happens at plot time. Only invoked when (re)building the panel file.
build_panel <- function(kover_dir, sample = NULL) {
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
  collapse_positions_to_intervals(as.data.frame(df))
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

# --- shared colour key / display decimation ----------------------------------

# Order samples by their trailing lane number (S1, S2, ... S22).
#
# Plain alphabetical sorting keys on the DEFGEN id that prefixes the lane, so
# the legend comes out ordered by accession number, and S7 lands before S15
# because "1" sorts before "7". Ids carrying no _S<n> suffix keep alphabetical
# order and are listed after the numbered ones.
order_samples <- function(x) {
  s <- sort(unique(as.character(x)))
  has <- grepl("_S[0-9]+$", s)
  n <- rep(NA_integer_, length(s))
  n[has] <- as.integer(sub("^.*_S([0-9]+)$", "\\1", s[has]))
  s[order(!has, n, s)]
}

# One colour per sample, fixed once for the whole report.
#
# Without this each figure builds its own hue scale from the samples it happens
# to contain. That is accidentally consistent while every figure holds the same
# 22 samples, and silently wrong the day one does not: a curve followed from one
# figure to the next would change colour without warning.
sample_palette <- function(samples) {
  s <- order_samples(samples)
  stats::setNames(scales::hue_pal()(length(s)), s)
}

# Thin a per-base series down for interactive display. OFF by default, and it
# should stay off unless the file size is genuinely blocking.
#
# The HTML serialises every plotted point as JSON: six full-resolution per-base
# figures come to ~170 MB. Keeping one base in `step` cuts that roughly tenfold.
#
# The cost is not "small events are missed", it is worse: small events are
# MISDRAWN. A 3 bp heterozygous deletion is a real, readable signal here —
# 100014_S15 carries NM_000548.5:c.4527_4529del and its depth goes 1754 ->
# 1040/1021/1022 -> 1726 at chr16:2,084,981-983, a square notch to ratio 0.64.
# At step 10 exactly one of those three bases survives, and the line joins it to
# neighbours 10 bp away on either side: the notch becomes a V. A figure that
# draws a ramp where the data has a step is more dangerous than one that draws
# nothing. With a one-base offset, the event disappears outright.
#
# So: usable for a coarse overview of large events, never for indel-scale
# reading. Exon first and last bases are always kept, so exon edges stay sharp
# and a per-exon step function stays exact. The PDF and PNGs are always rendered
# from the full data.
decimate_for_display <- function(df, step) {
  step <- as.integer(step)
  if (is.na(step) || step <= 1L) return(df)
  dt <- data.table::as.data.table(df)
  data.table::setorder(dt, sample, gene, exon_num, exonic_x)
  dt[, .keep := seq_len(.N) %% step == 1L | seq_len(.N) == 1L | seq_len(.N) == .N,
     by = .(sample, gene, exon)]
  out <- dt[.keep == TRUE]
  out[, .keep := NULL]
  as.data.frame(out)
}

# --- smoothing ---------------------------------------------------------------

# Rolling smoother for `y` along the exonic axis, `n` bases wide and centred.
#
# After cohort referencing what is left is mostly per-base counting noise, with
# no positional structure: smoothing it away leaves the exon-level signal and
# turns the sawtooth into readable plateaus.
#
# `stat` picks the estimator, and the choice matters more than the window size:
#   "median" (default) is edge-preserving. A drop that starts mid-exon stays a
#     drop, with its front and its depth intact.
#   "mean" blurs an edge into a ramp `n` bases wide and shaves its depth, so an
#     abrupt fall reads as a shallow dip. Kept because it is the smoother most
#     people expect, but it is the wrong tool for spotting breakpoints.
#
# The window is confined to each (sample, gene, exon) block so it never bleeds
# across an exon boundary and blurs a real step; blocks shorter than the window
# use the whole block. Window edges: runmed's Tukey end rule for the median, and
# nearest-value fill for the mean, so each exon's curve stays continuous.
smooth_y <- function(df, n, stat = c("median", "mean")) {
  stat <- match.arg(stat)
  n <- as.integer(n)
  if (is.na(n) || n <= 1L) return(df)
  dt <- data.table::as.data.table(df)
  data.table::setorder(dt, sample, gene, exon_num, exonic_x)
  dt[, y := {
    w <- min(n, .N)
    if (w %% 2L == 0L) w <- w - 1L        # both smoothers need an odd window
    if (w <= 1L) y
    else if (stat == "median") as.numeric(stats::runmed(y, w, endrule = "median"))
    else {
      sm <- data.table::frollmean(y, w, align = "center", na.rm = TRUE)
      data.table::nafill(data.table::nafill(sm, "locf"), "nocb")
    }
  }, by = .(sample, gene, exon)]
  as.data.frame(dt)
}

# Collapse each exon to its median: one level per exon per sample, i.e. a true
# step function. Same per-base data, read at the resolution CNV calling uses.
square_y <- function(df) {
  dt <- data.table::as.data.table(df)
  dt[, y := stats::median(y, na.rm = TRUE), by = .(sample, gene, exon)]
  as.data.frame(dt)
}

# "off"/"0" -> no smoothing; "exon" -> per-exon median; an integer -> window bp.
parse_smooth <- function(x) {
  x <- tolower(trimws(as.character(x)))
  if (x %in% c("", "0", "off", "none", "false")) return(list(kind = "none", n = 0L))
  if (x == "exon") return(list(kind = "exon", n = 0L))
  n <- suppressWarnings(as.integer(x))
  if (is.na(n) || n < 0L) stop("--smooth expects off | exon | a window size in bp, got: ", x)
  if (n <= 1L) list(kind = "none", n = 0L) else list(kind = "roll", n = n)
}

# --- plot --------------------------------------------------------------------

plot_coverage_curve <- function(long_df, layout, mode = MODES,
                                facet = c("gene", "none"), gc_tbl = NULL,
                                smooth = "off", smooth_stat = "median",
                                axis = c("exonic", "genomic"),
                                palette = NULL, show_legend = TRUE,
                                decimate = 0L) {
  mode  <- match.arg(mode, MODES)
  facet <- match.arg(facet)
  axis  <- match.arg(axis)
  # One continuous axis. Only for an exonic layout: globalize_layout re-offsets
  # genes in alphabetical order, which would undo a genomic layout and put PKD1
  # before TSC2.
  if (facet == "none" && axis == "exonic") layout <- globalize_layout(layout)
  ld <- to_exonic_coords(long_df, layout)
  # Titles and axis labels are in French: the figures circulate to the lab and
  # the thesis jury. Code and comments stay in English.
  if (mode == "raw") {
    ld$y  <- ld$depth
    ylab  <- "Profondeur (nombre de lectures par base)"
    title <- "Profondeur de séquençage brute"
  } else if (mode == "cohortonly") {
    ld$y  <- ld$depth
    ld    <- normalise_by_cohort_reference(ld)
    ylab  <- "Profondeur / médiane de la cohorte"
    title <- "Couverture rapportée à la cohorte, sans normalisation par échantillon"
  } else {
    ld <- normalise_by_sample_median(ld)
    if (mode == "cohort") {
      ld <- normalise_by_cohort_reference(ld)
      ylab  <- "Couverture relative (1 = couverture attendue pour deux copies)"
      title <- "Couverture rapportée à la cohorte"
    } else if (mode %in% c("gc", "capture")) {
      if (is.null(gc_tbl))
        stop("mode '", mode, "' needs gc_tbl (run gc_by_position() first).")
      ld  <- decompose_cohort_reference(ld, gc_tbl)
      div <- if (mode == "gc") ld$gcurve else ld$probe
      ld$y <- ld$y / div
      ld$y[!is.finite(div) | div == 0] <- NA_real_
      if (mode == "gc") {
        ylab  <- "Couverture relative, composante GC retirée"
        title <- "Couverture rapportée à la cohorte, correction du GC seule"
      } else {
        ylab  <- "Couverture relative, composante de capture retirée"
        title <- "Couverture rapportée à la cohorte, correction des sondes seule"
      }
    } else {
      ylab  <- "Profondeur / médiane de l'échantillon"
      title <- "Couverture normalisée par échantillon"
    }
  }
  sm <- parse_smooth(smooth)
  if (sm$kind == "roll") {
    ld <- smooth_y(ld, sm$n, smooth_stat)
    stat_fr <- if (smooth_stat == "median") "médiane" else "moyenne"
    title <- paste0(title, ", ", stat_fr, " glissante sur ", sm$n, " pb")
  } else if (sm$kind == "exon") {
    ld <- square_y(ld)
    title <- paste0(title, ", médiane par exon")
  }

  ld <- decimate_for_display(ld, decimate)   # display only; never before a stat
  # Legend order follows the factor levels, so set them by lane number here,
  # after every statistic has been computed.
  ld$sample <- factor(ld$sample, levels = order_samples(ld$sample))

  L <- data.table::as.data.table(layout)
  track_y <- -0.06 * max(ld$y, na.rm = TRUE)
  # 87 exons on one axis cannot all be labelled; thin them out.
  Ltxt <- if (axis == "genomic") L[seq(1L, .N, by = 3L)] else L

  xlab <- if (axis == "genomic")
    paste("Ordre génomique (exons du panel mis bout à bout ; introns et régions",
          "intergéniques non représentés)")
  else
    paste("Position exonique (exons mis bout à bout, largeur proportionnelle à",
          "leur longueur ; introns non représentés)")

  p <- ggplot(ld, aes(x = exonic_x, y = y, group = sample, color = sample,
                 text = paste0("échantillon : ", sample, "<br>gène : ", gene,
                               "<br>exon : ", exon, "<br>position : ", pos,
                               "<br>profondeur : ", depth))) +
    # DeCovA alert lines, only where y is a ratio to the cohort (centred on 1).
    (if (mode %in% c("cohort", "cohortonly", "gc", "capture"))
       geom_hline(yintercept = c(DEL_THRESHOLD, DUP_THRESHOLD),
                  linetype = "dashed", colour = "grey40", linewidth = 0.3)) +
    geom_line(linewidth = 0.25, alpha = 0.8) +
    geom_segment(data = L, inherit.aes = FALSE,
                 aes(x = x_start, xend = x_end, y = track_y, yend = track_y),
                 linewidth = 3, color = "grey70") +
    geom_text(data = Ltxt, inherit.aes = FALSE,
              aes(x = x_mid, y = track_y, label = exon),
              size = 1.8, vjust = 1.6) +
    labs(x = xlab, y = ylab, color = "Échantillon", title = title) +
    theme_bw(base_size = 9) +
    theme(legend.position = if (show_legend) "right" else "none",
          legend.key.size = unit(0.35, "cm"),
          legend.text = element_text(size = 6))
  if (!is.null(palette))
    p <- p + scale_colour_manual(values = palette,
                                 guide = if (show_legend) "legend" else "none")

  # On a genomic axis, name each gene and mark where one hands over to the next:
  # a deletion crossing the junction must read as one event, not a coincidence.
  if (axis == "genomic") {
    gb <- L[, .(x = (min(x_start) + max(x_end)) / 2), by = gene]
    top <- max(ld$y, na.rm = TRUE)
    jm <- junction_marker(layout)
    if (nrow(jm)) {
      p <- p +
        geom_vline(data = jm, inherit.aes = FALSE, aes(xintercept = x),
                   linetype = "dashed", colour = "grey30") +
        geom_text(data = jm, inherit.aes = FALSE,
                  aes(x = x, y = top, label = paste0(left, " | ", right,
                                                     "  (", gap, " pb)")),
                  size = 2.2, angle = 90, vjust = -0.4, hjust = 1,
                  colour = "grey30")
    }
    p <- p + geom_text(data = gb, inherit.aes = FALSE,
                       aes(x = x, y = top, label = gene),
                       size = 3.2, fontface = "bold", vjust = 1.2,
                       colour = "grey20")
  }

  # Stack genes vertically (one row each) so neither is cramped; "none"
  # aggregates the whole panel onto a single axis with no per-gene split.
  if (facet == "gene") p <- p + facet_wrap(~ gene, ncol = 1, scales = "free_x")
  p
}

# --- CLI ---------------------------------------------------------------------

parse_args <- function(args) {
  opts <- list(
    depth_dir     = "",            # required: dir of mosdepth *.per-base.bed.gz
    kover_dir     = "",            # KoverGraph DATA/<run> dir, --rebuild-panel only
    genes         = "TSC1,TSC2",
    out           = "output/perbase_coverage",
    panel         = file.path(script_dir(), "panel", "panel.tsv"),
    panel_sample  = "",
    rebuild_panel = "false",
    mode          = "normalised",  # see BARE_MODES below
    facet         = "",            # gene | none ; default auto (none for whole panel)
    ref           = "",            # reference FASTA, required by gc/capture only
    gc_window     = "101",         # bp, centred on each base; gc/capture only
    smooth        = "off",         # off | exon | <bp> rolling window
    smooth_stat   = "median",      # median (edge-preserving) | mean
    html_step     = "1"            # keep 1 base in N for the HTML only; 1 = off
  )
  given <- character()   # keys the caller actually passed; selects report vs single
  i <- 1
  while (i <= length(args)) {
    key <- gsub("-", "_", sub("^--", "", args[i]))  # accept --depth-dir or --depth_dir
    if (key %in% names(BARE_MODES)) {               # bare flags, e.g. --cohort-only
      opts$mode <- unname(BARE_MODES[key])
      given <- c(given, "mode")
      i <- i + 1
    } else if (key %in% names(opts) && i + 1 <= length(args)) {
      opts[[key]] <- args[i + 1]; given <- c(given, key); i <- i + 2
    } else stop("Unknown or incomplete argument: ", args[i])
  }
  opts$given <- given
  opts
}

# The report is the default. Asking for a specific view (a mode flag, a
# smoothing, a facet) means the caller wants that one figure, as before.
SINGLE_FIGURE_KEYS <- c("mode", "smooth", "smooth_stat", "facet")
wants_single_figure <- function(opts) any(SINGLE_FIGURE_KEYS %in% opts$given)

# Genes whose junction gets its own figures. TSC2 and PKD1 are adjacent on
# 16p13.3 and their joint deletion is the contiguous gene deletion syndrome
# (TSC + polycystic kidney disease), which only reads as one event on a shared
# genomic axis.
JUNCTION_GENES <- c("TSC2", "PKD1")

# Load the panel geometry, rebuilding it from KoverGraph only when asked.
load_panel <- function(opts) {
  panel <- path.expand(opts$panel)
  rebuild <- tolower(opts$rebuild_panel) %in% c("true", "yes", "1")
  if (rebuild || !file.exists(panel)) {
    if (!nzchar(opts$kover_dir)) stop("--rebuild-panel needs --kover-dir")
    kover <- path.expand(opts$kover_dir)
    message("Building panel geometry from ", kover)
    iv <- build_panel(kover, sample = if (nzchar(opts$panel_sample)) opts$panel_sample else NULL)
    if (!dir.exists(dirname(panel))) dir.create(dirname(panel), recursive = TRUE)
    write_panel(iv, panel)
    message("Wrote panel geometry: ", panel)
    iv
  } else {
    read_panel(panel)
  }
}

load_depth_for <- function(intervals, depth_dir) {
  bed <- tempfile(fileext = ".bed")
  on.exit(unlink(bed), add = TRUE)
  write_panel_bed(intervals, bed)
  load_all_samples(depth_dir, bed, intervals)
}

# The three views applied to one scope: the readable per-exon step function, the
# per-base view the pipeline reports, and absolute depth for "is there enough
# depth to read a ratio at all".
scope_views <- function(long, layout, facet, axis, prefix,
                        palette = NULL, decimate = 0L) {
  mk <- function(mode, smooth) plot_coverage_curve(
    long, layout, mode = mode, facet = facet, smooth = smooth, axis = axis,
    palette = palette, decimate = decimate)
  stats::setNames(
    list(mk("cohort", "exon"), mk("cohort", "off"), mk("raw", "off")),
    paste0(prefix, c("exon", "cohort", "raw")))
}

# Ratio series behind a view, for the contrast table. Mirrors what
# plot_coverage_curve computes for mode "cohort", at the same smoothing.
ratio_series <- function(long, layout, smooth) {
  ld <- normalise_by_cohort_reference(
    normalise_by_sample_median(to_exonic_coords(long, layout)))
  if (smooth == "exon") ld <- square_y(ld)
  ld
}

main <- function() {
  opts  <- parse_args(commandArgs(trailingOnly = TRUE))
  genes <- if (tolower(opts$genes) %in% c("all", "panel", "")) NULL
           else trimws(strsplit(opts$genes, ",")[[1]])
  if (!nzchar(opts$depth_dir)) stop("--depth-dir is required (dir of *.per-base.bed.gz)")
  depth <- path.expand(opts$depth_dir)
  panel <- load_panel(opts)

  iv_scope <- if (is.null(genes)) panel else panel[panel$gene %in% genes, , drop = FALSE]
  if (nrow(iv_scope) == 0) stop("No exons left after gene filter: ", opts$genes)
  facet <- if (nzchar(opts$facet)) match.arg(tolower(opts$facet), c("gene", "none"))
           else if (is.null(genes)) "none" else "gene"

  # ---- single figure: an explicit view was requested -----------------------
  if (wants_single_figure(opts)) {
    layout <- exon_layout(iv_scope)
    message("Loading depth from ", depth)
    long <- load_depth_for(iv_scope, depth)
    message("Samples: ", dplyr::n_distinct(long$sample),
            " | genes: ", paste(unique(iv_scope$gene), collapse = ", "),
            " | exons: ", nrow(iv_scope))
    mode <- match.arg(tolower(opts$mode), MODES)
    gc_tbl <- NULL
    if (mode %in% c("gc", "capture")) {
      win <- as.integer(opts$gc_window)
      if (!nzchar(opts$ref)) stop("--gc and --capture need --ref <reference FASTA>")
      message("Computing local GC (", win, " bp window) from ", opts$ref)
      gc_tbl <- gc_by_position(iv_scope, opts$ref, window = win)
    }
    sm <- parse_smooth(opts$smooth)
    message("Mode: ", mode, " | facet: ", facet,
            if (sm$kind != "none") paste0(" | smoothing: ", opts$smooth) else "")
    p <- plot_coverage_curve(long, layout, mode = mode, facet = facet,
                             gc_tbl = gc_tbl, smooth = opts$smooth,
                             smooth_stat = opts$smooth_stat)
    tag <- switch(sm$kind, none = "", exon = "_exon",
                  roll = paste0("_s", sm$n, substr(opts$smooth_stat, 1, 3)))
    outdir <- timestamped_outdir(
      file.path(opts$out, paste0(infer_run_scope(depth, genes), "_", mode, tag)))
    paths <- render_outputs(p, outdir, "perbase_coverage_curve")
    message("Wrote:\n  ", paste(paths, collapse = "\n  "))
    return(invisible(paths))
  }

  # ---- report: the default -------------------------------------------------
  junction <- all(JUNCTION_GENES %in% panel$gene)
  iv_junc  <- panel[panel$gene %in% JUNCTION_GENES, , drop = FALSE]
  # One depth pass over the union; each figure then subsets BEFORE normalising,
  # since normalise_by_sample_median takes the median over the bases plotted.
  iv_union <- panel[panel$gene %in% unique(c(iv_scope$gene,
                                             if (junction) JUNCTION_GENES)), , drop = FALSE]
  message("Loading depth from ", depth)
  long <- load_depth_for(iv_union, depth)
  message("Samples: ", dplyr::n_distinct(long$sample),
          " | scope: ", paste(unique(iv_scope$gene), collapse = ", "),
          if (junction) " | junction: TSC2-PKD1" else "",
          " | exons: ", nrow(iv_union))

  pal  <- sample_palette(long$sample)
  step <- as.integer(opts$html_step)

  lay_scope  <- exon_layout(iv_scope)
  long_scope <- long[long$gene %in% iv_scope$gene, , drop = FALSE]
  plots  <- scope_views(long_scope, lay_scope, facet, "exonic", "", pal)
  light  <- scope_views(long_scope, lay_scope, facet, "exonic", "", pal, step)
  sizes <- list(exon = c(14, 10), cohort = c(14, 10), raw = c(14, 10))
  contrast <- rbind(
    exon_contrast(ratio_series(long_scope, lay_scope, "exon"), "1_exon"),
    exon_contrast(ratio_series(long_scope, lay_scope, "off"),  "2_cohort"))

  if (junction) {
    lay_j  <- genomic_layout(iv_junc)
    long_j <- long[long$gene %in% JUNCTION_GENES, , drop = FALSE]
    plots <- c(plots, scope_views(long_j, lay_j, "none", "genomic", "junction_", pal))
    light <- c(light, scope_views(long_j, lay_j, "none", "genomic", "junction_",
                                  pal, step))
    for (n in c("junction_exon", "junction_cohort", "junction_raw"))
      sizes[[n]] <- c(18, 7)
    contrast <- rbind(contrast,
      exon_contrast(ratio_series(long_j, lay_j, "exon"), "4_junction_exon"),
      exon_contrast(ratio_series(long_j, lay_j, "off"),  "5_junction_cohort"))
  }

  outdir <- timestamped_outdir(
    file.path(opts$out, paste0(infer_run_scope(depth, genes), "_report")))
  if (step > 1L)
    message("WARNING: HTML thinned to 1 base in ", step, ". Events a few bases ",
            "wide (small indels) will be misdrawn as spikes or lost entirely; ",
            "read those on the PDF or PNG, which keep full resolution.")
  paths <- render_report(plots, outdir, "coverage_report", sizes = sizes,
                         html_plots = light)

  # Top contrasts per view: the same event scored under each rendering, so the
  # views can be ranked on evidence rather than on impression.
  top <- do.call(rbind, lapply(split(contrast, contrast$view),
                               function(d) utils::head(d, 20)))
  tsv <- file.path(outdir, "discernibility.tsv")
  utils::write.table(top, tsv, sep = "\t", quote = FALSE, row.names = FALSE)

  message("Wrote:\n  ", paste(c(paths, contrast = tsv), collapse = "\n  "))
  invisible(c(paths, contrast = tsv))
}

if (sys.nframe() == 0) main()
