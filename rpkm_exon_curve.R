#!/usr/bin/env Rscript

# RPKM exon coverage curve (DeCovA-style) for the TSC panel.
#
# Reads the per-exon RPKM matrix (Gene/Transcrit/Exon + one column per sample),
# normalizes each exon to the cohort median, and renders a single vertical
# column of panels — one per gene — as PNG + PDF + interactive HTML. One
# colored line per patient, overlaid; scroll to compare genes.
#
# Note: KoverGraph Blanket ships the RPKMs sheet already cohort-normalised
# (median == 1 on every exon), so the ratio equals the sheet value there. The
# ratio is still computed so non-normalised inputs plot correctly.
#
# Design: docs/superpowers/specs/2026-06-17-rpkm-exon-curve-design.md

suppressPackageStartupMessages({
  library(readxl)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(plotly)
})

script_dir <- function() {
  ca <- commandArgs(FALSE)
  f <- sub("^--file=", "", ca[grepl("^--file=", ca)])
  if (length(f)) dirname(normalizePath(f)) else getwd()
}
if (!exists("timestamped_outdir")) source(file.path(script_dir(), "coverage_common.R"))

# Infer a run id (e.g. run_1) from an input filename; fallback = basename.
infer_run <- function(path) {
  stem <- sub("\\.[^.]+$", "", basename(path))
  m <- regmatches(stem, regexpr("^[A-Za-z]+_[0-9]+", stem))
  if (length(m) == 1 && nzchar(m)) m else stem
}

META_COLS <- c("Gene", "Transcrit", "Exon")

# Read the RPKMs sheet and forward-fill the fill-down meta columns.
# "None" and blank cells in Gene/Transcrit mark continuation rows.
load_rpkm <- function(path, sheet = "RPKMs") {
  df <- readxl::read_excel(path, sheet = sheet)
  df[META_COLS] <- lapply(df[META_COLS], function(col) {
    col <- as.character(col)
    col[col == "None"] <- NA
    col
  })
  df <- tidyr::fill(df, dplyr::all_of(c("Gene", "Transcrit")), .direction = "down")

  sample_cols <- setdiff(names(df), META_COLS)
  df[sample_cols] <- lapply(df[sample_cols], function(x) suppressWarnings(as.numeric(x)))
  df
}

# Filter to core exon rows, subset genes, pivot to long form, add median ratio.
#   genes = NULL keeps the whole panel; include_padding keeps _lpad/_rpad/_UTR.
prepare_long <- function(df, genes = c("TSC1", "TSC2"), include_padding = FALSE) {
  if (!include_padding) {
    df <- dplyr::filter(df, grepl("^E[0-9]+$", .data$Exon))
  }
  if (!is.null(genes)) {
    df <- dplyr::filter(df, .data$Gene %in% genes)
  }
  if (nrow(df) == 0) {
    stop("No exon rows left after filtering. Check the `genes` argument and sheet contents.")
  }

  # Preserve genomic exon order (order of appearance), not alphabetical.
  exon_levels <- unique(df$Exon)
  df$Exon <- factor(df$Exon, levels = exon_levels)
  df$Gene <- factor(df$Gene, levels = unique(df$Gene))

  sample_cols <- setdiff(names(df), META_COLS)
  long <- tidyr::pivot_longer(
    df,
    cols = dplyr::all_of(sample_cols),
    names_to = "patient",
    values_to = "rpkm"
  )

  long %>%
    dplyr::group_by(.data$Gene, .data$Exon) %>%
    dplyr::mutate(
      cohort_median = stats::median(.data$rpkm, na.rm = TRUE),
      ratio = dplyr::if_else(.data$cohort_median > 0,
                             .data$rpkm / .data$cohort_median,
                             NA_real_)
    ) %>%
    dplyr::ungroup()
}

# One panel per gene, stacked in a single column. Plots the cohort ratio only:
# on KoverGraph Blanket input the RPKMs sheet is already cohort-normalised
# (median == 1 on every exon), so a separate raw-RPKM panel was bit-for-bit
# identical to the ratio panel. Plotting `ratio` keeps the figure correct if a
# future input ships non-normalised RPKMs — but note that in that case the
# absolute-depth QC view is no longer on the figure.
plot_coverage_curve <- function(long_df) {
  genes <- levels(long_df$Gene)

  # `tip` is built once as a column so every layer — including the red CNV
  # markers, which draw on top and would otherwise swallow the hover with no
  # tooltip of their own — carries the same sample-identifying hover text.
  plot_df <- long_df %>%
    dplyr::mutate(
      value = .data$ratio,
      tip = paste0("patient: ", .data$patient,
                   "<br>gene: ", .data$Gene,
                   "<br>exon: ", .data$Exon,
                   "<br>ratio: ", round(.data$ratio, 3),
                   "<br>rpkm: ", round(.data$rpkm, 3))
    )

  # Threshold reference lines, repeated on every gene panel.
  thresh_df <- do.call(rbind, lapply(genes, function(g) data.frame(
    Gene = factor(g, levels = genes),
    y = c(DEL_THRESHOLD, 1.0, DUP_THRESHOLD),
    kind = c("del", "neutral", "dup")
  )))

  # CNV candidates: points outside the [del, dup] band. Tag the tooltip so the
  # red markers say why they are red, on top of naming the sample.
  flagged <- plot_df %>%
    dplyr::filter(!is.na(.data$value),
                  .data$value < DEL_THRESHOLD | .data$value > DUP_THRESHOLD) %>%
    dplyr::mutate(tip = paste0(.data$tip, "<br><b>CNV candidate (",
                               dplyr::if_else(.data$value > 1, "DUP", "DEL"), ")</b>"))

  ggplot(plot_df, aes(x = .data$Exon, y = .data$value, group = .data$patient,
                      color = .data$patient, text = .data$tip)) +
    geom_hline(data = thresh_df, aes(yintercept = .data$y, linetype = .data$kind),
               inherit.aes = FALSE, color = "grey40") +
    geom_line(linewidth = 0.35, alpha = 0.8) +
    geom_point(size = 0.6, alpha = 0.7) +
    geom_point(data = flagged,
               aes(x = .data$Exon, y = .data$value, text = .data$tip),
               inherit.aes = FALSE, color = "red", size = 1.4) +
    facet_wrap(~ Gene, ncol = 1, scales = "free") +
    scale_linetype_manual(values = c(del = "dashed", neutral = "solid", dup = "dashed"),
                          guide = "none") +
    labs(x = "Exon", y = "RPKM ratio to cohort median", color = "Patient",
         title = "RPKM exon coverage across cohort (ratio to cohort median)") +
    theme_bw(base_size = 9) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 5),
          legend.position = "right",
          legend.key.size = unit(0.35, "cm"),
          legend.text = element_text(size = 6))
}

# Minimal `--key value` CLI parser (avoids an optparse dependency).
parse_args <- function(args) {
  opts <- list(
    xlsx = "",   # required: KoverGraph <run>_Blanket.xlsx
    genes = "TSC1,TSC2",
    out = "output/rpkm_curve",
    include_padding = "false"
  )
  i <- 1
  while (i <= length(args)) {
    key <- sub("^--", "", args[i])
    if (key %in% names(opts) && i + 1 <= length(args)) {
      opts[[key]] <- args[i + 1]
      i <- i + 2
    } else {
      stop("Unknown or incomplete argument: ", args[i])
    }
  }
  opts
}

main <- function() {
  opts <- parse_args(commandArgs(trailingOnly = TRUE))
  genes <- if (tolower(opts$genes) %in% c("all", "panel", "")) {
    NULL
  } else {
    trimws(strsplit(opts$genes, ",")[[1]])
  }
  include_padding <- tolower(opts$include_padding) %in% c("true", "yes", "1")

  if (!nzchar(opts$xlsx)) stop("--xlsx is required (KoverGraph <run>_Blanket.xlsx)")
  path <- path.expand(opts$xlsx)
  message("Reading: ", path)
  long <- prepare_long(load_rpkm(path), genes = genes, include_padding = include_padding)
  message("Exons: ", dplyr::n_distinct(long$Exon),
          " | patients: ", dplyr::n_distinct(long$patient),
          " | genes: ", paste(levels(long$Gene), collapse = ", "))

  p <- plot_coverage_curve(long)

  scope <- if (is.null(genes)) "panel" else paste(genes, collapse = "-")
  base <- file.path(path.expand(opts$out), paste0(infer_run(path), "_", scope))
  outdir <- timestamped_outdir(base)
  # One panel per gene stacked in a single column: scale height with the panel
  # count so each row keeps room (≈8 in/panel).
  n_panels <- dplyr::n_distinct(long$Gene)
  paths <- render_outputs(p, outdir, "rpkm_exon_curve",
                          width = 12, height = 8 * n_panels)
  message("Wrote:\n  ", paste(paths, collapse = "\n  "))
}

if (sys.nframe() == 0) {
  main()
}
