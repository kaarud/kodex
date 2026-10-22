# Shared helpers for the coverage-curve scripts (rpkm_exon_curve.R,
# perbase_coverage_curve.R).

suppressPackageStartupMessages({
  library(ggplot2)
  library(plotly)
})

# DeCovA alert levels (level_del / level_dup) on a ratio to the cohort median.
# Shared by the exon (RPKM) and per-base views so both figures read alike.
DEL_THRESHOLD <- 0.8
DUP_THRESHOLD <- 1.2

# Timestamped output dir (DD_HHMM) under `base`, so reruns never overwrite.
timestamped_outdir <- function(base, now = Sys.time()) {
  dir <- file.path(base, format(now, "%d_%H%M"))
  if (!dir.exists(dir)) dir.create(dir, recursive = TRUE)
  dir
}

# Render a ggplot to PNG + PDF + interactive HTML. Returns named paths.
# width/height size the static PNG/PDF (inches); tall single-column layouts
# need a larger height so panels aren't squished.
render_outputs <- function(p, outdir, basename, width = 14, height = 8) {
  if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)
  png_path  <- file.path(outdir, paste0(basename, ".png"))
  pdf_path  <- file.path(outdir, paste0(basename, ".pdf"))
  html_path <- file.path(outdir, paste0(basename, ".html"))

  ggplot2::ggsave(png_path, p, width = width, height = height, dpi = 150)
  ggplot2::ggsave(pdf_path, p, width = width, height = height)

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

# One ggplot -> one interactive block of the report page.
as_spaced_widget <- function(p, wrap = 110L, gap = "48px", bottom = 95) {
  w <- plotly::plotly_build(plotly::ggplotly(p, tooltip = "text"))
  m <- w$x$layout$margin          # plotly's computed margins: keep, never drop

  # ggplotly drops labs(subtitle = ) silently. Ours carry interpretation caveats
  # (the PKD1 pseudogene warning) and the figure travels on its own, so losing it
  # in the HTML while keeping it in the PDF is not acceptable. Re-attached under
  # the title as a wrapped <sup> line.
  title <- p$labels$title
  t_mar <- m$t
  sub   <- p$labels$subtitle
  if (!is.null(sub) && nzchar(sub)) {
    lines <- strwrap(sub, width = wrap)
    title <- paste0(title, "<br><sup>", paste(lines, collapse = "<br>"), "</sup>")
    t_mar <- 40 + 16 * (length(lines) + 1)
  }

  # Every margin is restated. `layout(margin = list(t = ...))` REPLACES the
  # margin object rather than merging into it, so naming only `t` silently reset
  # bottom/left/right to plotly's defaults — which is what clipped the x-axis
  # title. The bottom is additionally floored: on a faceted plot ggplotly places
  # the axis title as an annotation that the default bottom margin cuts off.
  #
  # plotly_build is not optional either: layout() on a ggplotly object only
  # queues the change in `layoutAttrs`, applied at print time. Building bakes it
  # into `layout`, which is what save_html serialises.
  w <- plotly::plotly_build(plotly::layout(
    w,
    title  = list(text = title, x = 0, xanchor = "left"),
    margin = list(t = t_mar, b = max(m$b, bottom), l = m$l, r = m$r)))

  # Stacked widgets otherwise sit flush, dropping one figure's title onto the
  # previous figure's axis label.
  htmltools::div(style = paste0("margin-bottom:", gap, ";"), w)
}

# Render several ggplots as ONE document: a multi-page PDF, one interactive HTML
# holding every figure, and one PNG per figure for previewing.
#
# `plots` is a named list; the names order the pages and suffix the PNGs.
# `sizes` is an optional list of c(width, height) per plot, recycled from the
# `width`/`height` defaults when absent — a single-axis junction figure wants a
# different shape from a two-panel faceted one.
#
# The HTML is deliberately NOT self-contained: per-base figures run to ~20 MB of
# point data each and embedding them fails. It needs the `lib/` directory beside
# it, so move the two together.
render_report <- function(plots, outdir, basename, sizes = list(),
                          width = 14, height = 10, html_plots = NULL) {
  if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)
  nms <- names(plots)
  if (is.null(nms) || any(!nzchar(nms)))
    stop("render_report needs a fully named list of plots.")
  size_of <- function(n) if (!is.null(sizes[[n]])) sizes[[n]] else c(width, height)

  pdf_path  <- file.path(outdir, paste0(basename, ".pdf"))
  html_path <- file.path(outdir, paste0(basename, ".html"))

  # One page per figure. grDevices::pdf takes a single page size, so use the
  # largest requested; smaller figures simply get more margin.
  page <- do.call(pmax, lapply(nms, size_of))
  grDevices::pdf(pdf_path, width = page[1], height = page[2], onefile = TRUE)
  for (n in nms) print(plots[[n]])
  grDevices::dev.off()

  png_paths <- vapply(seq_along(nms), function(i) {
    n <- nms[i]; sz <- size_of(n)
    f <- file.path(outdir, sprintf("%s_%d_%s.png", basename, i, n))
    ggplot2::ggsave(f, plots[[n]], width = sz[1], height = sz[2], dpi = 150)
    f
  }, character(1))

  # The HTML is served the lighter plots when they exist: full per-base data
  # for six figures runs to ~170 MB of JSON, which no browser opens comfortably.
  src <- if (is.null(html_plots)) plots else html_plots
  widgets <- lapply(nms, function(n) as_spaced_widget(src[[n]]))
  htmltools::save_html(htmltools::tagList(widgets), file = html_path,
                       libdir = "lib")
  message("Interactive HTML needs the adjacent 'lib/' directory; move both together.")

  invisible(c(pdf = pdf_path, html = html_path,
              stats::setNames(png_paths, paste0("png_", nms))))
}

# How far each (sample, exon) stands out from the cohort cloud, ON THE VALUES A
# GIVEN VIEW ACTUALLY DRAWS:
#
#   contrast = |median y of this sample over the exon - 1| / MAD(all points the
#              other samples draw over that exon)
#
# The denominator is the visual thickness of the cohort band: it carries both
# the spread between samples and the per-base scatter within each of them. That
# is the point of measuring on the displayed values rather than on exon medians
# — the exon medians are identical whether or not the curve was smoothed, so a
# median-based contrast would score every view the same and say nothing about
# discernibility. Smoothing thins the band without moving the medians, and this
# ratio registers that.
#
# Only meaningful for ratio views centred on 1; callers skip it for raw depth.
# The MAD is floored so a degenerate exon ranks high rather than yielding Inf.
exon_contrast <- function(df, view, mad_floor = 1e-3) {
  dt <- data.table::as.data.table(df)[!is.na(y)]
  em <- dt[, .(m = stats::median(y)), by = .(sample, gene, exon)]

  # Band thickness per exon, measured over every sample including the deviant
  # one. Leaving it out would be more correct in principle and pointless in
  # practice: the MAD has a 50% breakdown point, so one outlier among twenty-odd
  # samples does not move it — and the leave-one-out join costs millions of rows.
  band <- dt[, .(spread = pmax(stats::mad(y), mad_floor)), by = .(gene, exon)]
  em <- merge(em, band, by = c("gene", "exon"), all.x = TRUE)

  em[, `:=`(view = view,
            contrast = abs(m - 1) / spread,
            direction = data.table::fifelse(m > 1, "DUP", "DEL"))]
  data.table::setorder(em, -contrast)
  as.data.frame(em[, .(view, sample, gene, exon,
                       exon_median = round(m, 3),
                       cohort_mad = round(spread, 4),
                       contrast = round(contrast, 2), direction)])
}
