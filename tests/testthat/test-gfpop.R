library(testthat)
library(gfpop)
context("graph")

data(ECG, package="gfpop")
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

test_that("gfpop runs on a graph imported from csv with stringsAsFactors=T", {
  # The following structure was produced from creating a a standard 'isotonic' graph,
  # exporting to csv, re-importing, and using dput to make portable.
  test_graph <- structure(list(
    state1 = structure(c(1L, 1L), .Label = "Iso", class = "factor"),
    state2 = structure(c(1L, 1L), .Label = "Iso", class = "factor"),
    type = structure(1:2, .Label = c("null", "up"), class = "factor"),
    parameter = 1:0, penalty = c(0L, 15L), K = c(Inf, Inf), a = c(
      0L,
      0L
    ), min = c(NA, NA), max = c(NA, NA)
  ), row.names = c(
    NA,
    -2L
  ), class = c("data.frame", "graph"))

  # Test that the imported graph can be used to run gfpop
  test_data <- gfpop::dataGenerator(50, changepoints = c(1), parameters = c(1))
  expect_equal(gfpop(data = test_data, mygraph = test_graph)$changepoints, 50)
})
