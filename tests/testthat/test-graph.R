library(testthat)
library(gfpop)
context("graph")

test_that("char states are ok", {
  ecg.graph <- gfpop::graph(
    gfpop::edge("Q", "R", "up",   0),
    gfpop::edge("R", "S", "down", 0),
    gfpop::edge("S", "Q", "up",   0),
    all.null=TRUE)
  expect_identical(sum(ecg.graph$type=="null"), 3L)
  data(ECG, package="gfpop")
  fit <- gfpop::gfpop(vectData = ECG$data$millivolts, mygraph = ecg.graph, type = "gauss")
  expect_is(fit$states.chr, "character")
})
