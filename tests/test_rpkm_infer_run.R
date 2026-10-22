library(testthat)
source("../coverage_common.R")
source("../rpkm_exon_curve.R")

test_that("infer_run pulls the run id from a KoverGraph xlsx name", {
  expect_equal(infer_run("~/data/run_1/KoverGraph/run_1_Blanket.xlsx"), "run_1")
  expect_equal(infer_run("/x/run_2_Singularity_Blanket.xlsx"), "run_2")
})

test_that("infer_run falls back to the bare basename", {
  expect_equal(infer_run("/x/cohort.xlsx"), "cohort")
})
