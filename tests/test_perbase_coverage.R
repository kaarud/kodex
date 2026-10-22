library(testthat)
library(data.table)
source("../coverage_common.R")
source("../perbase_coverage_curve.R")

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

test_that("write_panel / read_panel round-trip preserves geometry and types", {
  iv <- data.frame(
    chrom = "chr9", start = c(100L, 200L), end = c(102L, 201L),
    gene = "TSC1", transcript = "NM_000368",
    exon = c("E01","E02"), exon_num = c(1L, 2L), stringsAsFactors = FALSE)
  f <- tempfile(fileext = ".tsv")
  on.exit(unlink(f))
  write_panel(iv, f)
  back <- read_panel(f)
  expect_equal(back[PANEL_COLS], iv[PANEL_COLS])
  expect_true(is.integer(back$start))
  expect_true(is.integer(back$exon_num))
})

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

test_that("globalize_layout shifts later genes onto one continuous axis", {
  intervals <- data.frame(
    gene = c("TSC1","TSC1","TSC2"), transcript = "NM",
    exon = c("E01","E02","E01"), exon_num = c(1L,2L,1L), chrom = "chr9",
    start = c(100L,200L, 50L), end = c(102L,201L, 54L), stringsAsFactors = FALSE)
  L  <- exon_layout(intervals)                 # per-gene: each gene starts at 0
  G  <- globalize_layout(L)
  t1 <- G[G$gene == "TSC1", ]
  t2 <- G[G$gene == "TSC2", ]
  expect_equal(t1$x_start, c(0L, 3L))          # TSC1 unchanged (first gene)
  expect_equal(max(t1$x_end), 5L)              # TSC1 total width = 3 + 2
  expect_equal(t2$x_start, 5)                  # TSC2 shifted past TSC1's width
  expect_equal(t2$x_end, 10)                   # 5 + 5 (TSC2 exon length)
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

test_that("normalise_by_sample_median divides each sample by its own median", {
  df <- data.frame(
    sample = c("A","A","A", "B","B","B"),
    depth  = c(10L,20L,30L, 100L,200L,300L),
    stringsAsFactors = FALSE)
  out <- normalise_by_sample_median(df)
  # median(A)=20, median(B)=200 -> both collapse onto 0.5/1/1.5
  expect_equal(out$y, c(0.5, 1, 1.5, 0.5, 1, 1.5))
  expect_equal(out$depth, df$depth)   # raw depth preserved for tooltips
})

test_that("normalise_by_cohort_reference divides each exonic position by the cohort median", {
  # 3 samples x 2 positions. A is a flat 0.5x (deletion), B normal, C 1.5x (dup).
  df <- data.frame(
    sample   = c("A","A", "B","B", "C","C"),
    gene     = "TSC1",
    exonic_x = c(0L,1L, 0L,1L, 0L,1L),
    y        = c(2.0, 0.5,  4.0, 1.0,  6.0, 1.5),
    stringsAsFactors = FALSE)
  out <- normalise_by_cohort_reference(df)
  o <- out[order(out$sample, out$exonic_x), ]
  # pos0 cohort median = median(2,4,6) = 4 ; pos1 = median(0.5,1,1.5) = 1
  expect_equal(o$y, c(0.5, 0.5,  1.0, 1.0,  1.5, 1.5))
})

test_that("normalise_by_cohort_reference yields NA where the cohort has no signal", {
  df <- data.frame(
    sample   = c("A","B"),
    gene     = "TSC1",
    exonic_x = c(0L, 0L),
    y        = c(0, 0),          # all samples zero at this position
    stringsAsFactors = FALSE)
  out <- normalise_by_cohort_reference(df)
  expect_true(all(is.na(out$y)))
})

test_that("infer_run_scope combines run dir and gene scope", {
  expect_equal(infer_run_scope("/x/run_1/depth", c("TSC1","TSC2")),
               "run_1_TSC1-TSC2")
  expect_equal(infer_run_scope("/x/run_1/depth", NULL), "run_1_panel")
})

# --- GC decomposition --------------------------------------------------------

test_that("rolling_gc computes a centred window over core bases only", {
  # half = 1 -> 3 bp window; first and last base are flanking padding.
  expect_equal(rolling_gc("AGCGT", 1L), c(2/3, 3/3, 2/3))
  # half = 0 -> each base is its own window.
  expect_equal(rolling_gc("ACGT", 0L), c(0, 1, 1, 0))
  # soft-masked reference sequence must count the same as upper case.
  expect_equal(rolling_gc("agcgt", 1L), rolling_gc("AGCGT", 1L))
})

test_that("rolling_gc excludes N from both numerator and denominator", {
  # core bases N,G,N -> windows ANG (1 GC / 2 defined), NGN (1/1), GNA (1/2)
  expect_equal(rolling_gc("ANGNA", 1L), c(0.5, 1, 0.5))
  expect_true(is.na(rolling_gc("ANNNA", 1L)[2]))      # window NNN: no signal
})

test_that("decompose_cohort_reference partitions the cohort reference exactly", {
  set.seed(1)
  n_pos <- 2400L   # 40 GC bins x 60 positions, above MIN_GC_BIN
  gc_tbl <- data.frame(chrom = "chr9", pos = seq_len(n_pos),
                       gc = rep(seq(0.30, 0.69, by = 0.01), each = 60L))
  base <- 1 + 0.5 * (gc_tbl$gc - 0.5) + rnorm(n_pos, sd = 0.02)  # GC + probe noise
  df <- do.call(rbind, lapply(c("A", "B", "C"), function(s) data.frame(
    sample = s, chrom = "chr9", pos = gc_tbl$pos, gene = "TSC1",
    exonic_x = gc_tbl$pos - 1L, y = base, stringsAsFactors = FALSE)))
  out <- decompose_cohort_reference(df, gc_tbl)
  expect_equal(out$gcurve * out$probe, out$ref)   # the partition identity
  expect_false(any(is.na(out$gcurve)))
})

test_that("decompose_cohort_reference attributes a pure GC trend to the GC component", {
  # The GC curve is fitted on gene-centred ref, so the residual `probe` is the
  # gene's own coverage level: constant, which is what "nothing left over" means.
  n_pos <- 2400L
  gc_tbl <- data.frame(chrom = "chr9", pos = seq_len(n_pos),
                       gc = rep(seq(0.30, 0.69, by = 0.01), each = 60L))
  base <- 1 + 0.5 * (gc_tbl$gc - 0.5)             # ref depends on GC only
  df <- do.call(rbind, lapply(c("A", "B", "C"), function(s) data.frame(
    sample = s, chrom = "chr9", pos = gc_tbl$pos, gene = "TSC1",
    exonic_x = gc_tbl$pos - 1L, y = base, stringsAsFactors = FALSE)))
  out <- decompose_cohort_reference(df, gc_tbl)
  expect_equal(diff(range(out$probe)), 0, tolerance = 1e-8)  # nothing left over
})

test_that("decompose_cohort_reference fits on positions, not on sample-rows", {
  # 40 positions x 2 samples = 80 rows; with MIN_GC_BIN = 50 a row-weighted fit
  # would keep the bin, a position-weighted one must not.
  gc_tbl <- data.frame(chrom = "chr9", pos = 1:40, gc = 0.50)
  df <- do.call(rbind, lapply(c("A", "B"), function(s) data.frame(
    sample = s, chrom = "chr9", pos = 1:40, gene = "TSC1",
    exonic_x = 0:39, y = 1, stringsAsFactors = FALSE)))
  expect_error(decompose_cohort_reference(df, gc_tbl), "GC bins")
})

test_that("cohortonly keeps library-size spread that cohort absorbs", {
  # B is sequenced 2x deeper than A and C; both are otherwise identical.
  df <- data.frame(
    sample   = rep(c("A", "B", "C"), each = 2L),
    gene     = "TSC1",
    exonic_x = rep(c(0L, 1L), 3L),
    depth    = c(10, 20,  20, 40,  10, 20),
    stringsAsFactors = FALSE)
  cohort_only <- normalise_by_cohort_reference(transform(df, y = df$depth))
  expect_equal(cohort_only$y[cohort_only$sample == "B"], c(2, 2))
  full <- normalise_by_cohort_reference(normalise_by_sample_median(df))
  expect_equal(full$y[full$sample == "B"], c(1, 1))
})

test_that("parse_args maps the ablation flags and the gc options", {
  expect_equal(parse_args(c("--cohort-only"))$mode, "cohortonly")
  expect_equal(parse_args(c("--cohortonly"))$mode, "cohortonly")
  expect_equal(parse_args(c("--gc"))$mode, "gc")
  expect_equal(parse_args(c("--capture"))$mode, "capture")
  expect_equal(parse_args(c("--cohort"))$mode, "cohort")
  expect_equal(parse_args(c("--gc-window", "201"))$gc_window, "201")
  expect_equal(parse_args(c("--ref", "/x/hg38.fa"))$ref, "/x/hg38.fa")
})

test_that("decompose_cohort_reference fits GC on gene-centred ref, not on the tier ladder", {
  # Two genes, same GC-to-coverage relation, but B sits 20x above A (tiered
  # panel). Pooling raw ref would make the curve track the tier; centring first
  # must give the same GC curve as either gene alone.
  gcv <- rep(seq(0.30, 0.69, by = 0.01), each = 60L)
  shape <- 1 + 0.5 * (gcv - 0.5)
  mk <- function(g, level) do.call(rbind, lapply(c("A", "B", "C"), function(s)
    data.frame(sample = s, chrom = "chr9", pos = seq_along(gcv), gene = g,
               exonic_x = seq_along(gcv) - 1L, y = shape * level,
               stringsAsFactors = FALSE)))
  df <- rbind(mk("TSC1", 1), mk("PDPK1", 20))
  gc_tbl <- data.frame(chrom = "chr9", pos = seq_along(gcv), gc = gcv)
  out <- decompose_cohort_reference(df, gc_tbl)
  expect_equal(out$gcurve * out$probe, out$ref)          # partition still exact
  # one probe level per gene, and the 20x ladder lands there, not in gcurve
  lv <- tapply(out$probe, out$gene, function(x) diff(range(x)))
  expect_true(all(lv < 1e-8))
  expect_equal(unname(tapply(out$probe, out$gene, mean)["PDPK1"] /
                      tapply(out$probe, out$gene, mean)["TSC1"]), 20,
               tolerance = 1e-6)
  expect_lt(max(out$gcurve), 2)                           # no tier ladder in the curve
})

# --- smoothing ---------------------------------------------------------------

test_that("parse_smooth reads off, exon and a window size", {
  expect_equal(parse_smooth("off")$kind, "none")
  expect_equal(parse_smooth("0")$kind, "none")
  expect_equal(parse_smooth("")$kind, "none")
  expect_equal(parse_smooth("exon")$kind, "exon")
  expect_equal(parse_smooth("101"), list(kind = "roll", n = 101L))
  expect_equal(parse_smooth("1")$kind, "none")       # a 1 bp window is a no-op
  expect_error(parse_smooth("wide"), "--smooth expects")
})

test_that("smooth_y averages within an exon and never across exon boundaries", {
  # E01 is a flat 1; E02 is a flat 5. A window spanning the join would bleed.
  df <- data.frame(
    sample = "A", gene = "TSC1",
    exon = rep(c("E01", "E02"), each = 5L), exon_num = rep(c(1L, 2L), each = 5L),
    exonic_x = 0:9, y = c(rep(1, 5), rep(5, 5)), stringsAsFactors = FALSE)
  out <- smooth_y(df, 3L, "mean")
  expect_equal(out$y[out$exon == "E01"], rep(1, 5))
  expect_equal(out$y[out$exon == "E02"], rep(5, 5))
})

test_that("smooth_y mean spreads a single-base spike, median deletes it", {
  df <- data.frame(
    sample = "A", gene = "TSC1", exon = "E01", exon_num = 1L,
    exonic_x = 0:8, y = c(1, 1, 1, 1, 10, 1, 1, 1, 1), stringsAsFactors = FALSE)
  avg <- smooth_y(df, 3L, "mean")
  expect_false(any(is.na(avg$y)))       # centred window leaves no NA at the ends
  expect_equal(avg$y[4:6], c(4, 4, 4))  # the spike, spread over 3 bp
  # the median removes a one-base outlier outright rather than spreading it
  expect_equal(smooth_y(df, 3L, "median")$y, rep(1, 9))
})

test_that("smooth_y with a window wider than the exon collapses it to its mean", {
  # mean only: the median of a monotone ramp is the ramp itself, by design.
  df <- data.frame(
    sample = "A", gene = "TSC1", exon = "E01", exon_num = 1L,
    exonic_x = 0:4, y = c(1, 2, 3, 4, 5), stringsAsFactors = FALSE)
  expect_equal(smooth_y(df, 101L, "mean")$y, rep(3, 5))
  expect_equal(smooth_y(df, 101L, "median")$y, c(1, 2, 3, 4, 5))
})

test_that("square_y gives one level per exon and per sample", {
  df <- data.frame(
    sample = rep(c("A", "B"), each = 4L), gene = "TSC1",
    exon = rep(c("E01", "E01", "E02", "E02"), 2L),
    exon_num = rep(c(1L, 1L, 2L, 2L), 2L),
    exonic_x = rep(0:3, 2L),
    y = c(1, 3, 10, 12,  0.4, 0.6, 1, 1), stringsAsFactors = FALSE)
  out <- square_y(df)
  expect_equal(out$y[out$sample == "A"], c(2, 2, 11, 11))
  expect_equal(out$y[out$sample == "B"], c(0.5, 0.5, 1, 1))
})

test_that("parse_args reads --smooth", {
  expect_equal(parse_args(c("--smooth", "exon"))$smooth, "exon")
  expect_equal(parse_args(c("--smooth", "201"))$smooth, "201")
  expect_equal(parse_args(character(0))$smooth, "off")
})

test_that("smooth_y median keeps an abrupt drop, mean blurs it into a ramp", {
  # A clean 1 -> 0.5 step in the middle of one exon, plus per-base noise.
  set.seed(42)
  n <- 401L
  step <- c(rep(1, 200), rep(0.5, 201))
  df <- data.frame(
    sample = "A", gene = "TSC1", exon = "E01", exon_num = 1L,
    exonic_x = seq_len(n) - 1L, y = step + rnorm(n, sd = 0.05),
    stringsAsFactors = FALSE)

  med <- smooth_y(df, 31L, "median")$y
  avg <- smooth_y(df, 31L, "mean")$y

  # depth of the step, measured well away from the transition
  expect_equal(mean(med[1:150]) - mean(med[251:401]), 0.5, tolerance = 0.02)
  # the median keeps the front sharp: the transition spans a few bases, the
  # mean smears it over the whole window
  front <- function(v) sum(v > 0.6 & v < 0.9)
  expect_lt(front(med), front(avg))
  # and both actually removed noise
  expect_lt(sd(med[1:150]), sd(df$y[1:150]))
  expect_lt(sd(avg[1:150]), sd(df$y[1:150]))
})

test_that("smooth_y median leaves no NA and respects exon blocks", {
  df <- data.frame(
    sample = "A", gene = "TSC1",
    exon = rep(c("E01", "E02"), each = 5L), exon_num = rep(c(1L, 2L), each = 5L),
    exonic_x = 0:9, y = c(rep(1, 5), rep(5, 5)), stringsAsFactors = FALSE)
  out <- smooth_y(df, 3L, "median")
  expect_false(any(is.na(out$y)))
  expect_equal(out$y[out$exon == "E01"], rep(1, 5))
  expect_equal(out$y[out$exon == "E02"], rep(5, 5))
})

test_that("parse_args reads --smooth-stat", {
  expect_equal(parse_args(c("--smooth-stat", "mean"))$smooth_stat, "mean")
  expect_equal(parse_args(character(0))$smooth_stat, "median")
})

# --- genomic layout ----------------------------------------------------------

# TSC2 (plus strand, exon_num up with coordinate) then PKD1 (minus strand,
# exon_num DOWN with coordinate) — the case exon_layout gets wrong.
two_gene_panel <- function() data.frame(
  chrom = "chr16",
  start = c(2048616L, 2088446L, 2089727L, 2135475L),
  end   = c(2048753L, 2088610L, 2090194L, 2135689L),
  gene  = c("TSC2", "TSC2", "PKD1", "PKD1"),
  transcript = c("NM_000548", "NM_000548", "NM_001009944", "NM_001009944"),
  exon = c("E02", "E42", "E46", "E01"), exon_num = c(2L, 42L, 46L, 1L),
  stringsAsFactors = FALSE)

test_that("genomic_layout orders by coordinate, not by exon number", {
  L <- genomic_layout(two_gene_panel())
  expect_equal(L$exon, c("E02", "E42", "E46", "E01"))   # genomic order
  expect_equal(L$gene, c("TSC2", "TSC2", "PKD1", "PKD1"))
  expect_true(all(diff(L$x_start) > 0))
  # exon_layout would put PKD1 E01 before E46 and restart TSC2 at 0
  expect_equal(exon_layout(two_gene_panel())$exon, c("E01", "E46", "E02", "E42"))
})

test_that("genomic_layout lays exons back to back with proportional widths", {
  L <- genomic_layout(two_gene_panel())
  expect_equal(L$len, L$gend - L$gstart + 1L)
  expect_equal(L$x_end, L$x_start + L$len)
  expect_equal(L$x_start[-1], utils::head(L$x_end, -1))   # no gaps, no overlap
  expect_equal(L$x_start[1], 0L)
})

test_that("genomic_layout refuses more than one chromosome", {
  iv <- two_gene_panel(); iv$chrom[1] <- "chr9"
  expect_error(genomic_layout(iv), "single chromosome")
})

test_that("to_exonic_coords places adjacent genes adjacently on a genomic layout", {
  L <- genomic_layout(two_gene_panel())
  # last base of TSC2 E42 and first base of PKD1 E46 must be neighbours on x
  df <- data.frame(
    sample = "A", chrom = "chr16", gene = c("TSC2", "PKD1"),
    exon = c("E42", "E46"), pos = c(2088610L, 2089727L), depth = 100,
    stringsAsFactors = FALSE)
  out <- to_exonic_coords(df, L)
  x <- out$exonic_x[match(c(2088610L, 2089727L), out$pos)]
  expect_equal(diff(x), 1)   # the 1117 bp intergenic gap is not drawn
})

test_that("junction_marker reports the boundary and the gap it hides", {
  j <- junction_marker(genomic_layout(two_gene_panel()))
  expect_equal(nrow(j), 1L)
  expect_equal(j$left, "TSC2")
  expect_equal(j$right, "PKD1")
  expect_equal(j$gap, 1116L)    # 2089727 - 2088610 - 1
  expect_equal(j$x, 138L + 165L)  # widths of TSC2 E02 and E42
})

test_that("junction_marker is empty for a single gene", {
  iv <- two_gene_panel()[1:2, ]
  expect_equal(nrow(junction_marker(genomic_layout(iv))), 0L)
})

# --- report assembly ---------------------------------------------------------

test_that("wants_single_figure distinguishes a view request from a plain run", {
  expect_false(wants_single_figure(parse_args(character(0))))
  expect_false(wants_single_figure(parse_args(c("--genes", "TSC1"))))
  expect_true(wants_single_figure(parse_args(c("--cohort"))))
  expect_true(wants_single_figure(parse_args(c("--smooth", "exon"))))
  expect_true(wants_single_figure(parse_args(c("--facet", "none"))))
})

test_that("render_report writes one PDF, one HTML and a PNG per figure", {
  skip_if_not_installed("htmltools")
  d <- file.path(tempdir(), "report_test")
  on.exit(unlink(d, recursive = TRUE), add = TRUE)
  toy <- function(t) ggplot2::ggplot(data.frame(x = 1:3, y = 1:3),
                                     ggplot2::aes(x, y)) +
    ggplot2::geom_line() + ggplot2::labs(title = t)
  paths <- render_report(list(alpha = toy("A"), beta = toy("B")), d, "rep",
                         sizes = list(alpha = c(6, 4)), width = 5, height = 3)
  expect_true(file.exists(paths[["pdf"]]))
  expect_true(file.exists(paths[["html"]]))
  expect_true(dir.exists(file.path(d, "lib")))
  expect_true(file.exists(file.path(d, "rep_1_alpha.png")))
  expect_true(file.exists(file.path(d, "rep_2_beta.png")))
})

test_that("render_report refuses an unnamed plot list", {
  toy <- ggplot2::ggplot(data.frame(x = 1, y = 1), ggplot2::aes(x, y))
  expect_error(render_report(list(toy), tempdir(), "x"), "fully named")
})

# --- contrast ----------------------------------------------------------------

# n_other samples flat at 1, one sample at `level`, `jitter` of per-base noise.
contrast_fixture <- function(level, jitter, n_other = 5L, n_base = 40L) {
  set.seed(7)
  mk <- function(s, lv) data.frame(
    sample = s, gene = "TSC2", exon = "E10",
    y = lv + rnorm(n_base, sd = jitter), stringsAsFactors = FALSE)
  rbind(mk("case", level),
        do.call(rbind, lapply(seq_len(n_other), function(i) mk(paste0("ctrl", i), 1))))
}

test_that("exon_contrast scores a lone deviant sample high", {
  out <- exon_contrast(contrast_fixture(0.5, 0.02), "v")
  expect_equal(out$sample[1], "case")
  expect_equal(out$direction[1], "DEL")
  expect_gt(out$contrast[1], 5)
})

test_that("exon_contrast scores a cohort-wide shift low", {
  # every sample at 0.5: the cohort centre moves, nothing stands out
  df <- contrast_fixture(0.5, 0.02)
  df$y[df$sample != "case"] <- df$y[df$sample != "case"] - 0.5
  out <- exon_contrast(df, "v")
  expect_lt(max(out$contrast) / min(out$contrast), 3)
})

test_that("exon_contrast rewards a thinner band for the same event", {
  # identical event, different per-base scatter: smoothing thins the band
  noisy  <- exon_contrast(contrast_fixture(0.5, 0.10), "noisy")
  smooth <- exon_contrast(contrast_fixture(0.5, 0.01), "smooth")
  expect_gt(smooth$contrast[1], noisy$contrast[1])
})

test_that("exon_contrast survives a perfectly flat cohort", {
  df <- contrast_fixture(0.5, 0)
  out <- exon_contrast(df, "v")
  expect_true(all(is.finite(out$contrast)))
})

test_that("a genomic layout survives facet = none (globalize_layout must not run)", {
  # globalize_layout re-offsets genes alphabetically; on a genomic layout that
  # would place PKD1 (higher coordinates) before TSC2 and invert the junction.
  iv <- two_gene_panel()
  L <- genomic_layout(iv)
  df <- data.frame(
    sample = "A", chrom = "chr16", gene = c("TSC2", "PKD1"),
    exon = c("E42", "E46"), pos = c(2088610L, 2089727L), depth = 100,
    exon_num = c(42L, 46L), stringsAsFactors = FALSE)
  p <- plot_coverage_curve(df, L, mode = "raw", facet = "none", axis = "genomic")
  x <- p$data$exonic_x[match(c("TSC2", "PKD1"), p$data$gene)]
  expect_lt(x[1], x[2])                       # TSC2 stays left of PKD1
  expect_equal(diff(x), 1)                    # and they stay adjacent
  expect_equal(junction_marker(L)$left, "TSC2")
})

test_that("as_spaced_widget re-attaches the subtitle ggplotly drops", {
  p <- ggplot2::ggplot(data.frame(x = 1:3, y = 1:3), ggplot2::aes(x, y)) +
    ggplot2::geom_line() +
    ggplot2::labs(title = "Titre", subtitle = paste(rep("avertissement", 20), collapse = " "))
  # ggplotly alone loses it
  bare <- plotly::ggplotly(p)
  expect_false(grepl("avertissement", bare$x$layout$title$text))
  # the wrapper keeps it, wrapped onto several lines, with room in the margin
  d <- as_spaced_widget(p, wrap = 40L)
  w <- d$children[[1]]
  expect_true(grepl("avertissement", w$x$layout$title$text))
  expect_gt(length(gregexpr("<br>", w$x$layout$title$text)[[1]]), 1)
  expect_gt(w$x$layout$margin$t, 60)
})

test_that("as_spaced_widget separates stacked figures and leaves plain plots alone", {
  p <- ggplot2::ggplot(data.frame(x = 1:3, y = 1:3), ggplot2::aes(x, y)) +
    ggplot2::geom_line() + ggplot2::labs(title = "Sans sous-titre")
  d <- as_spaced_widget(p)
  expect_true(grepl("margin-bottom", d$attribs$style))
  expect_equal(d$children[[1]]$x$layout$title$text, "Sans sous-titre")
})

# --- shared colour key / display decimation ----------------------------------

test_that("sample_palette is stable and named, whatever the input order", {
  a <- sample_palette(c("B", "A", "C"))
  b <- sample_palette(c("C", "C", "A", "B"))
  expect_equal(a, b)                       # same colour per sample either way
  expect_equal(names(a), c("A", "B", "C"))
  expect_length(unique(a), 3L)
})

test_that("a shared palette gives one sample the same colour in every figure", {
  # a sample must keep its colour from figure to figure, so that a curve can be
  # followed across the report; a figure missing one sample would otherwise
  # reassign every colour after it
  pal <- sample_palette(c("S1", "S2", "S3"))
  expect_equal(pal[["S3"]], sample_palette(c("S1", "S2", "S3"))[["S3"]])
  expect_false(pal[["S3"]] == sample_palette(c("S1", "S3"))[["S3"]])  # why it is needed
})

test_that("decimate_for_display thins the middle but keeps both exon edges", {
  df <- data.frame(
    sample = "A", gene = "TSC1", exon = "E01", exon_num = 1L,
    exonic_x = 0:99, y = 1, stringsAsFactors = FALSE)
  out <- decimate_for_display(df, 10L)
  expect_lt(nrow(out), 15)
  expect_true(0 %in% out$exonic_x)     # first base kept
  expect_true(99 %in% out$exonic_x)    # last base kept: exon edges stay sharp
  expect_equal(decimate_for_display(df, 1L), df)   # step 1 is a no-op
})

test_that("decimate_for_display is lossless on a per-exon step function", {
  # square_y makes y constant within an exon; thinning must not change what is
  # drawn, only how many points carry it
  df <- data.frame(
    sample = rep(c("A", "B"), each = 20L), gene = "TSC1",
    exon = rep(rep(c("E01", "E02"), each = 10L), 2L),
    exon_num = rep(rep(c(1L, 2L), each = 10L), 2L),
    exonic_x = rep(0:19, 2L),
    y = rep(rep(c(0.5, 1.0), each = 10L), 2L), stringsAsFactors = FALSE)
  out <- decimate_for_display(df, 5L)
  lv <- unique(out[, c("sample", "exon", "y")])
  expect_equal(nrow(lv), 4L)                       # 2 samples x 2 exons
  expect_setequal(lv$y, c(0.5, 1.0))
})

test_that("decimate_for_display never mixes exons or samples", {
  df <- data.frame(
    sample = rep(c("A", "B"), each = 6L), gene = "TSC1",
    exon = rep(rep(c("E01", "E02"), each = 3L), 2L),
    exon_num = rep(rep(c(1L, 2L), each = 3L), 2L),
    exonic_x = rep(0:5, 2L), y = 1, stringsAsFactors = FALSE)
  out <- decimate_for_display(df, 2L)
  # every (sample, exon) block keeps at least its two edges
  n <- table(out$sample, out$exon)
  expect_true(all(n >= 2))
})

test_that("as_spaced_widget keeps every margin and floors the bottom one", {
  # layout(margin = list(t = ...)) replaces the margin object instead of merging,
  # so naming only `t` used to reset bottom/left/right and clip the axis title.
  p <- ggplot2::ggplot(data.frame(x = 1:3, y = 1:3, g = c("a", "a", "b")),
                       ggplot2::aes(x, y)) +
    ggplot2::geom_line() + ggplot2::facet_wrap(~ g, ncol = 1) +
    ggplot2::labs(title = "T", subtitle = "S", x = "Exonic position")
  m <- as_spaced_widget(p)$children[[1]]$x$layout$margin
  expect_gte(m$b, 95)
  expect_true(all(vapply(m[c("t", "b", "l", "r")], is.numeric, logical(1))))
  expect_gt(m$l, 0)
  expect_gt(m$r, 0)
})

test_that("display decimation is off by default", {
  expect_equal(parse_args(character(0))$html_step, "1")
})

test_that("decimation misdraws a 3 bp deletion — why it is off by default", {
  # Real event: 100014_S15, NM_000548.5:c.4527_4529del, chr16:2,084,981-983.
  # Depth 1754 -> 1040/1021/1022 -> 1726, i.e. a square notch to ratio ~0.64.
  pos <- 2084961:2085000
  y <- rep(1.08, length(pos))
  y[match(2084981:2084983, pos)] <- c(0.64, 0.63, 0.64)
  df <- data.frame(sample = "100014_S15", gene = "TSC2", exon = "E35",
                   exon_num = 35L, exonic_x = seq_along(pos), pos = pos, y = y,
                   stringsAsFactors = FALSE)

  full <- decimate_for_display(df, 1L)
  expect_equal(sum(full$y < 0.7), 3L)          # the notch, three bases wide

  thin <- decimate_for_display(df, 10L)
  expect_lt(sum(thin$y < 0.7), 3L)             # at best one base of it survives
  # and the survivor is joined to neighbours ten bases away: a V, not a notch
  kept <- thin$pos[thin$y < 0.7]
  if (length(kept)) expect_gt(min(diff(sort(c(kept, thin$pos[thin$y > 0.7])))), 1L)
})

test_that("order_samples sorts by lane number, not by accession", {
  x <- c("100014_S15", "100003_S7", "100016_S22", "100008_S2")
  expect_equal(order_samples(x),
               c("100008_S2", "100003_S7", "100014_S15", "100016_S22"))
  # what plain sorting would have given, and why this exists
  expect_equal(sort(x),
               c("100003_S7", "100008_S2", "100014_S15", "100016_S22"))
})

test_that("order_samples puts S7 before S15 and keeps unsuffixed ids last", {
  expect_equal(order_samples(c("a_S15", "a_S7", "a_S2")),
               c("a_S2", "a_S7", "a_S15"))
  expect_equal(order_samples(c("b_S3", "control", "a_S1")),
               c("a_S1", "b_S3", "control"))
  expect_equal(order_samples(c("z", "a")), c("a", "z"))
})

test_that("sample_palette follows the lane order", {
  pal <- sample_palette(c("100014_S15", "100003_S7", "100008_S2"))
  expect_equal(names(pal), c("100008_S2", "100003_S7", "100014_S15"))
})

test_that("plot_coverage_curve orders its legend by lane number", {
  L <- exon_layout(data.frame(
    chrom = "chr9", start = 100L, end = 104L, gene = "TSC1",
    transcript = "NM", exon = "E01", exon_num = 1L, stringsAsFactors = FALSE))
  df <- data.frame(
    sample = rep(c("100014_S15", "100003_S7", "100008_S2"), each = 5L),
    chrom = "chr9", gene = "TSC1", exon = "E01", exon_num = 1L,
    pos = rep(100:104, 3L), depth = 100, stringsAsFactors = FALSE)
  p <- plot_coverage_curve(df, L, mode = "raw", facet = "none")
  expect_equal(levels(p$data$sample),
               c("100008_S2", "100003_S7", "100014_S15"))
})

test_that("cohort-ratio views draw the DeCovA alert lines, raw and per-sample views do not", {
  intervals <- data.frame(
    gene = "TSC1", transcript = "NM_000368",
    exon = c("E01","E02"), exon_num = c(1L,2L), chrom = "chr9",
    start = c(100L, 200L), end = c(102L, 201L), stringsAsFactors = FALSE)
  L <- exon_layout(intervals)
  df <- data.frame(sample = rep(c("A", "B", "C"), each = 5L), gene = "TSC1",
                   exon = rep(c("E01","E01","E01","E02","E02"), 3L),
                   pos = rep(c(100L,101L,102L,200L,201L), 3L),
                   depth = c(10,10,10,10,10, 20,20,20,20,20, 10,10,5,10,10),
                   stringsAsFactors = FALSE)
  hlines <- function(p) unlist(unname(lapply(p$layers, function(l)
    if (inherits(l$geom, "GeomHline")) l$data$yintercept)))
  for (m in c("cohort", "cohortonly"))
    expect_equal(sort(hlines(plot_coverage_curve(df, L, mode = m))),
                 c(DEL_THRESHOLD, DUP_THRESHOLD))
  for (m in c("raw", "normalised"))
    expect_null(hlines(plot_coverage_curve(df, L, mode = m)))
})
