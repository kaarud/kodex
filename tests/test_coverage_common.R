library(testthat)
source("../coverage_common.R")

test_that("timestamped_outdir builds DD_HHMM under base and creates it", {
  base <- file.path(tempdir(), "tsout_test")
  on.exit(unlink(base, recursive = TRUE), add = TRUE)
  t <- as.POSIXct("2026-06-17 10:24:00", tz = "UTC")
  out <- timestamped_outdir(base, now = t)
  expect_equal(basename(out), "17_1024")
  expect_true(dir.exists(out))
})
