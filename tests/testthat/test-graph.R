library(testthat)
library(gfpop)
context("graph")

test_that("graph ok with no args", {
  g <- graph()
  expect_is(g, "graph")
  expect_is(g, "data.frame")
  expect_identical(nrow(g), 0L)
})

test_that("graph type=updown ok", {
  pen.val <- 1.1
  g <- graph(type="updown", penalty=pen.val)
  expect_is(g, "graph")
  expect_is(g, "data.frame")
  expect_identical(g$state1, c("D", "U", "D", "U"))
  expect_identical(g$state2, c("U", "D", "D", "U"))
  expect_identical(g$type, c("up", "down", "null", "null"))
  expect_identical(g$penalty, c(pen.val, pen.val, 0, 0))
})

test_that("custom 3 state graph ok, default all null edges", {
  g <- graph(
    Edge("Q", "R", "up",   1.5,  gap=1000),
    Edge("R", "S", "down", 0,    gap=5000),
    Edge("S", "Q", "up",   0,    gap=2000))
  expect_is(g, "graph")
  expect_is(g, "data.frame")
  expect_identical(g$state1, c("Q", "R", "S", "Q", "R", "S"))
  expect_identical(g$state2, c("R", "S", "Q", "Q", "R", "S"))
  expect_identical(g$type, c("up", "down", "up", "null", "null", "null"))
  expect_identical(g$penalty, c(1.5, 0, 0, 0, 0, 0))
  expect_identical(g$parameter, c(1000, 5000, 2000, 1, 1, 1))
})

test_that("custom 3 state graph with no null edges", {
  g <- graph(
    Edge("Q", "R", "up",   1.5,  gap=1000),
    Edge("R", "S", "down", 0,    gap=5000),
    Edge("S", "Q", "up",   0,    gap=2000),
    all.null.edges=FALSE)
  expect_is(g, "graph")
  expect_is(g, "data.frame")
  expect_identical(g$state1, c("Q", "R", "S"))
  expect_identical(g$state2, c("R", "S", "Q"))
  expect_identical(g$type, c("up", "down", "up"))
  expect_identical(g$penalty, c(1.5, 0, 0))
  expect_identical(g$parameter, c(1000, 5000, 2000))
})
