library(testthat)
library(gfpop)
context("graph")

test_that("graph ok with no args", {
  g <- graph()
  expect_is(g, "graph")
  expect_is(g, "data.frame")
  expect_identical(nrow(g), 0L)
})

test_that("error for non-graph args", {
  expect_error({
    graph(data.frame(state1="foo", state2="bar"))
  }, error="please use gfpop::Edge to specify graph, problem arguments: 1")
})

D <- "Dw"
U <- "Up"
pen.val <- 1.1
test_that("graph type=updown ok", {
  g <- graph(type="updown", penalty=pen.val)
  expect_is(g, "graph")
  expect_is(g, "data.frame")
  expect_identical(g$state1, c(D, U, D, U))
  expect_identical(g$state2, c(D, U, U, D))
  expect_identical(g$type, c("null", "null", "up", "down"))
  expect_identical(g$penalty, c(0, 0, pen.val, pen.val))
})

test_that("graph type=updown with StartEnd ok", {
  g <- rbind(
    graph(type="updown", penalty=pen.val),
    StartEnd(D, D))
  expect_is(g, "graph")
  expect_is(g, "data.frame")
  expect_identical(g$state1, c(D, U, D, U, D, D))
  expect_identical(g$state2, c(D, U, U, D, NA, NA))
  expect_identical(g$type, c("null", "null", "up", "down", "start", "end"))
  expect_identical(g$penalty, c(0, 0, pen.val, pen.val, NA, NA))
})

test_that("custom 3 state graph ok, default all null edges", {
  g <- graph(
    Edge("Q", "R", "up",   penalty=1.5,  gap=1000),
    Edge("R", "S", "down", penalty=0,    gap=5000),
    Edge("S", "Q", "up",   penalty=0,    gap=2000), all.null.edges = TRUE)
  expect_is(g, "graph")
  expect_is(g, "data.frame")
  expect_identical(g$state1, c( "S", "R", "Q", "Q", "R", "S"))
  expect_identical(g$state2, c("S", "R", "Q", "R", "S", "Q"))
  expect_identical(g$type, c("null", "null", "null", "up", "down", "up"))
  expect_identical(g$penalty, c(0, 0, 0, 1.5, 0, 0))
  expect_identical(g$parameter, c(1, 1, 1, 1000, 5000, 2000))
})

test_that("custom 3 state graph with no null edges", {
  g <- graph(
    Edge("Q", "R", "up",   penalty=1.5,  gap=1000),
    Edge("R", "S", "down", penalty=0,    gap=5000),
    Edge("S", "Q", "up",   penalty=0,    gap=2000),
    all.null.edges=FALSE)
  expect_is(g, "graph")
  expect_is(g, "data.frame")
  expect_identical(g$state1, c("Q", "R", "S"))
  expect_identical(g$state2, c("R", "S", "Q"))
  expect_identical(g$type, c("up", "down", "up"))
  expect_identical(g$penalty, c(1.5, 0, 0))
  expect_identical(g$parameter, c(1000, 5000, 2000))
})
