library(testthat)
library(gfpop)
context("graph")

data(ECG, package="gfpop")
test_that("gfpop returns character states", {
  myGraph <- graph(
    Edge("beforeQ", "Q",      "down", 80000000, gap=0),
    Edge("Q",       "R",      "up",   0,        gap=2000),
    Edge("R",       "S",      "down", 0,        gap=5000),
    Edge("S",       "S1",     "up",   0,        gap=2000),
    Edge("S1",      "S2",     "up",   0,        gap=1000),
    Edge("S2",      "peak",   "up",   0,        gap=0),
    Edge("peak",    "after",  "down", 0,        gap=0),
    Edge("after",   "last",   "down", 0,        gap=0),
    Edge("last",    "beforeQ","up",   0,        gap=0), all.null.edges = TRUE)
  fit <- gfpop(ECG$data$millivolts, mygraph = myGraph, type = "mean")
  expect_true(all(fit$states %in% myGraph$state1))
})
