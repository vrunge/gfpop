library(testthat)
library(gfpop)
context("edge")

test_that("edge is vectorized", {
  expect_silent({
    edges <- gfpop::Edge(0:1, 0:1, "null")
  })
  expect_identical(edges$state1, 0:1)
  expect_identical(edges$state2, 0:1)
  expect_identical(edges$type, c("null", "null"))
  expect_identical(edges$penalty, c(0, 0))
})
