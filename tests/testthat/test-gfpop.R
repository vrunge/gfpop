library(testthat)
library(gfpop)
library(devtools)
devtools::install_github("vrunge/gfpop.data")
library(gfpop.data)
context("graph")

data(ECG, package="gfpop.data")
test_that("gfpop returns character states", {
  myGraph <- graph(
    Edge("beforeQ", "Q",      "down", penalty = 80000000, gap=0),
    Edge("Q",       "R",      "up",   penalty = 0,        gap=2000),
    Edge("R",       "S",      "down", penalty = 0,        gap=5000),
    Edge("S",       "S1",     "up",   penalty = 0,        gap=2000),
    Edge("S1",      "S2",     "up",   penalty = 0,        gap=1000),
    Edge("S2",      "peak",   "up",   penalty = 0,        gap=0),
    Edge("peak",    "after",  "down", penalty = 0,        gap=0),
    Edge("after",   "last",   "down", penalty = 0,        gap=0),
    Edge("last",    "beforeQ","up",   penalty = 0,        gap=0), all.null.edges = TRUE)
  fit <- gfpop(ECG$data$millivolts, mygraph = myGraph, type = "mean")
  expect_true(all(fit$states %in% myGraph$state1))
})
