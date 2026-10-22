#!/usr/bin/env Rscript

# Memory-safe, full-resolution coverage curves.
#
# Unlike the original implementation, this keeps the panel annotation once and
# stores depth as a compact integer matrix (bases x samples).  It never creates
# a data frame containing one copy of chrom/gene/exon for every sample.  PNG and
# PDF are drawn one patient at a time, so every base is kept without retaining a
# second, plot-sized table in memory.  HTML is deliberately not emitted: a
# full-resolution interactive graph serialises every point again as JSON.

suppressPackageStartupMessages(library(data.table))

script_dir <- function() {
  x <- commandArgs(FALSE)
  f <- sub("^--file=", "", x[grepl("^--file=", x)])
  if (length(f)) dirname(normalizePath(f)) else getwd()
}

usage <- function() {
  cat("Usage:\n  Rscript perbase_coverage_curve_memory_safe.R --depth-dir DIR --out DIR [--panel FILE] [--genes TSC1,TSC2|all] [--reindex true]\n\n",
      "The depth directory must contain *.per-base.bed.gz files.\n",
      "The default panel is panel/panel.tsv next to this script.\n", sep = "")
}

parse_args <- function(args) {
  out <- list(depth_dir = NULL, out = "output/perbase_coverage",
              panel = file.path(script_dir(), "panel", "panel.tsv"),
              genes = "TSC1,TSC2", format = "both", mode = "report",
              reindex = "false")
  bare <- c(raw = "raw", normalised = "normalised", normalized = "normalised",
            cohort = "cohort")
  i <- 1L
  while (i <= length(args)) {
    key <- gsub("-", "_", sub("^--", "", args[i]))
    if (key %in% names(bare)) {
      out$mode <- unname(bare[key]); i <- i + 1L
    } else if (key %in% names(out) && i < length(args)) {
      out[[key]] <- args[i + 1L]; i <- i + 2L
    } else if (key %in% c("help", "h")) {
      usage(); quit(status = 0L)
    } else stop("Unknown or incomplete argument: ", args[i])
  }
  if (is.null(out$depth_dir)) stop("--depth-dir is required")
  out$format <- match.arg(tolower(out$format), c("png", "pdf", "both"))
  out
}

read_panel <- function(path) {
  if (!file.exists(path)) stop("Panel not found: ", path,
                              "\nPass --panel /path/to/panel.tsv")
  p <- fread(path, sep = "\t")
  needed <- c("chrom", "start", "end", "gene", "transcript", "exon", "exon_num")
  if (!all(needed %in% names(p))) stop("panel.tsv lacks required columns")
  p[, start := as.integer(start)][, end := as.integer(end)]
  p
}

# One small annotation row per panel base.  This is intentionally stored once.
make_base_map <- function(intervals) {
  iv <- copy(as.data.table(intervals))
  setorder(iv, gene, exon_num, start)
  iv[, len := end - start + 1L]
  iv[, x_start := cumsum(c(0L, head(len, -1L))), by = gene]
  pieces <- lapply(seq_len(nrow(iv)), function(i) {
    n <- iv$len[i]
    data.table(chrom = iv$chrom[i], pos = seq.int(iv$start[i], iv$end[i]),
               gene = iv$gene[i], exon = iv$exon[i], exon_num = iv$exon_num[i],
               x = iv$x_start[i] + seq_len(n) - 1L)
  })
  ans <- rbindlist(pieces, use.names = TRUE)
  ans[, row_id := .I]
  setkey(ans, chrom, pos)
  ans
}

write_filter_bed <- function(intervals) {
  f <- tempfile(fileext = ".bed")
  b <- as.data.table(intervals)[, .(chrom, start = start - 1L, end)]
  setorder(b, chrom, start)
  fwrite(b, f, sep = "\t", col.names = FALSE)
  f
}

ensure_tabix_index <- function(gz, force = FALSE) {
  if (force || (!file.exists(paste0(gz, ".tbi")) && !file.exists(paste0(gz, ".csi")))) {
    # mosdepth per-base files are 0-based, half-open BED intervals.
    status <- system2("tabix", c("-f", "-0", "-s", "1", "-b", "2", "-e", "3", gz))
    if (status != 0L) stop("Could not index: ", gz)
  }
}

# Read one sample only. The temporary expanded table is discarded immediately;
# the retained result is one integer column of the compact depth matrix.
read_depth_column <- function(gz, bed, map, force_reindex = FALSE) {
  ensure_tabix_index(gz, force = force_reindex)
  # -0 applies to the BED region file passed with -R as well as to indexing.
  # Without it tabix warns on start=0 and may interpret the filter incorrectly.
  raw <- fread(cmd = paste("tabix -0 -R", shQuote(bed), shQuote(gz)), header = FALSE,
               col.names = c("chrom", "start", "end", "depth"), showProgress = FALSE)
  if (!nrow(raw)) stop("No panel depth records in ", gz)
  lengths <- raw$end - raw$start
  expanded <- data.table(chrom = rep(raw$chrom, lengths),
                         pos = sequence(lengths, from = raw$start + 1L),
                         depth = rep.int(as.integer(raw$depth), lengths))
  setkey(expanded, chrom, pos)
  aligned <- expanded[map, on = .(chrom, pos)]
  if (anyNA(aligned$depth))
    stop("Missing depth at ", sum(is.na(aligned$depth)), " panel bases in ", basename(gz))
  as.integer(aligned$depth)
}

load_depth_matrix <- function(depth_dir, intervals, force_reindex = FALSE) {
  files <- list.files(depth_dir, pattern = "\\.per-base\\.bed\\.gz$", full.names = TRUE)
  if (!length(files)) stop("No *.per-base.bed.gz files under ", depth_dir)
  map <- make_base_map(intervals)
  bed <- write_filter_bed(intervals)
  on.exit(unlink(bed), add = TRUE)
  d <- matrix(NA_integer_, nrow = nrow(map), ncol = length(files),
              dimnames = list(NULL, sub("\\.per-base\\.bed\\.gz$", "", basename(files))))
  message("Loading ", ncol(d), " samples × ", nrow(d), " bases into a compact matrix")
  for (i in seq_along(files)) {
    message("  [", i, "/", length(files), "] ", basename(files[i]))
    d[, i] <- read_depth_column(files[i], bed, map, force_reindex)
  }
  list(map = map, depth = d)
}

column_medians <- function(d) vapply(seq_len(ncol(d)), function(i) median(d[, i]), numeric(1))
row_medians_positive <- function(x) {
  # Run in chunks to bound transient allocations even on large panels.
  ans <- numeric(nrow(x)); block <- 10000L
  for (lo in seq.int(1L, nrow(x), by = block)) {
    hi <- min(nrow(x), lo + block - 1L)
    ans[lo:hi] <- apply(x[lo:hi, , drop = FALSE], 1L, function(z) median(z[z > 0]))
  }
  ans
}

sample_colours <- function(n) grDevices::hcl.colors(n, palette = "Dynamic")

draw_view <- function(map, depth, yfun, title, outfile, type) {
  genes <- unique(map$gene)
  samples <- colnames(depth)
  cols <- sample_colours(ncol(depth))
  if (type == "png") {
    grDevices::png(outfile, width = 4200, height = max(1100, 950 * length(genes)), res = 300)
  } else grDevices::pdf(outfile, width = 14, height = max(4, 3.2 * length(genes)))
  on.exit(grDevices::dev.off(), add = TRUE)
  op <- par(mfrow = c(length(genes), 1L), mar = c(3.4, 4.1, 3.0, 1.0), oma = c(0, 0, 2, 0))
  on.exit(par(op), add = TRUE)
  for (g in genes) {
    idx <- which(map$gene == g)
    # Determine limits from a small number of vectors, never from a long table.
    ymax <- 0
    for (i in seq_len(ncol(depth))) {
      z <- yfun(i)[idx]
      ymax <- max(ymax, stats::quantile(z[is.finite(z)], 0.995, names = FALSE, na.rm = TRUE))
    }
    if (!is.finite(ymax) || ymax <= 0) ymax <- 1
    plot(map$x[idx], rep(NA_real_, length(idx)), type = "n", xlab = "Position exonique (pb)",
         ylab = if (grepl("brute", title)) "Profondeur" else "Couverture relative",
         ylim = c(0, ymax * 1.05), main = g)
    abline(h = if (grepl("brute", title)) 0 else 1, col = "grey75", lty = 2)
    # DeCovA alert levels (coverage_common.R), on the cohort-ratio view only.
    if (grepl("cohorte", title)) abline(h = c(0.8, 1.2), col = "grey40", lty = 2)
    for (i in seq_len(ncol(depth))) lines(map$x[idx], yfun(i)[idx], col = cols[i], lwd = 0.45)
  }
  mtext(title, outer = TRUE, line = 0.3, cex = 1.2, font = 2)
}

main <- function() {
  o <- parse_args(commandArgs(trailingOnly = TRUE))
  panel <- read_panel(path.expand(o$panel))
  genes <- if (tolower(o$genes) %in% c("all", "panel", "")) unique(panel$gene)
           else trimws(strsplit(o$genes, ",", fixed = TRUE)[[1]])
  iv <- panel[gene %in% genes]
  if (!nrow(iv)) stop("No exons left after --genes ", o$genes)
  reindex <- tolower(o$reindex) %in% c("true", "yes", "1")
  obj <- load_depth_matrix(path.expand(o$depth_dir), iv, force_reindex = reindex)
  map <- obj$map; depth <- obj$depth; rm(obj); gc()
  med <- column_medians(depth)
  no_usable_coverage <- !is.finite(med) | med <= 0
  if (any(no_usable_coverage)) {
    message("WARNING: zero/absent panel coverage; omitted from normalised and cohort ",
            "views: ", paste(colnames(depth)[no_usable_coverage], collapse = ", "))
    med[no_usable_coverage] <- NA_real_
  }
  # This numeric matrix has one value per base/sample; it replaces the much
  # larger repeated-character long table used in the original implementation.
  normalised <- sweep(depth, 2L, med, "/")
  cohort_ref <- NULL
  if (o$mode %in% c("report", "cohort")) {
    message("Computing the per-base cohort reference")
    cohort_ref <- row_medians_positive(normalised)
  }
  modes <- if (o$mode == "report") c("normalised", "cohort", "raw") else o$mode
  stamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
  outdir <- file.path(path.expand(o$out), paste0("perbase_full_resolution_", stamp))
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  for (mode in modes) {
    if (mode == "raw") {
      fn <- function(i) as.numeric(depth[, i]); title <- "Profondeur de séquençage brute"
    } else if (mode == "normalised") {
      fn <- function(i) normalised[, i]; title <- "Couverture normalisée par échantillon"
    } else {
      fn <- function(i) normalised[, i] / cohort_ref; title <- "Couverture rapportée à la cohorte"
    }
    for (type in if (o$format == "both") c("png", "pdf") else o$format) {
      f <- file.path(outdir, paste0("perbase_", mode, ".", type))
      message("Writing ", f)
      draw_view(map, depth, fn, title, f, type)
    }
    gc()
  }
  message("Wrote full-resolution files to: ", outdir)
}

if (sys.nframe() == 0L) main()
