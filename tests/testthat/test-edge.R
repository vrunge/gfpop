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

test_that("StartEnd same col names as Edge", {
  se <- StartEnd("down", "down")
  e <- Edge("down", "down", "std")
  expect_identical(names(se), names(e))
})

test_that("StartEnd is vectorized", {
  start.vec <- "A"
  end.vec <- c("A","B")
  se <- StartEnd(start.vec, end.vec)
  start.rows <- subset(se, type=="start")
  end.rows <- subset(se, type=="end")
  expect_identical(start.rows$state1, start.vec)
  expect_identical(end.rows$state1, end.vec)
})

