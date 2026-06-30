library(testthat)
library(gfpop)
context("rule")

test_that("Edge has an integer rule column, NA by default", {
  e <- Edge("a", "b", "std")
  expect_true("rule" %in% names(e))
  expect_identical(names(e)[length(names(e))], "rule")
  expect_type(e$rule, "integer")
  expect_true(is.na(e$rule))
})

test_that("Edge(rule=) sets the rule id", {
  e <- Edge("a", "b", "std", rule = 2)
  expect_identical(e$rule, 2L)
})

test_that("Edge is vectorized over rule", {
  e <- Edge(c("a", "c"), c("b", "d"), "std", rule = 1)
  expect_identical(e$rule, c(1L, 1L))
})

test_that("Edge rejects non-positive rule ids", {
  expect_error(Edge("a", "b", "std", rule = 0), "rule must be a positive integer")
  expect_error(Edge("a", "b", "std", rule = -1), "rule must be a positive integer")
})

test_that("StartEnd and Node carry a rule column matching Edge", {
  se <- StartEnd("a", "a")
  nd <- Node("a", 0, 1)
  e <- Edge("a", "a", "std")
  expect_identical(names(se), names(e))
  expect_identical(names(nd), names(e))
  expect_true(is.na(se$rule[1]))
  expect_true(is.na(nd$rule[1]))
})

test_that("predefined graphs have an all-NA rule column", {
  g <- graph(type = "updown", penalty = 5)
  expect_true("rule" %in% names(g))
  expect_true(all(is.na(g$rule)))
})

test_that("custom graph preserves per-edge rule values", {
  g <- graph(
    Edge("a", "a", "null", rule = 1),
    Edge("a", "a", "std", penalty = 5, rule = 2))
  expect_identical(g$rule, c(1L, 2L))
})

test_that("graphReorder preserves rule values with their edges", {
  g <- graph(
    Edge("a", "a", "null", rule = 1),
    Edge("a", "a", "std", penalty = 5, rule = 2),
    StartEnd("a", "a"))
  reordered <- gfpop:::graphReorder(g)$graph
  expect_true("rule" %in% names(reordered))
  expect_identical(reordered$rule[reordered$type == "null"], 1L)
  expect_identical(reordered$rule[reordered$type == "std"], 2L)
})

test_that("gfpop still solves a standard problem with the rule column present", {
  set.seed(42)
  x <- c(rnorm(40, 0), rnorm(40, 10))
  g <- graph(type = "std", penalty = 2 * log(80))
  fit <- gfpop(x, g, type = "mean")
  expect_s3_class(fit, "gfpop")
  ## last changepoint is the data length
  expect_identical(fit$changepoints[length(fit$changepoints)], 80L)
  ## one internal change near the true location (index 40)
  internal <- fit$changepoints[-length(fit$changepoints)]
  expect_true(any(internal >= 35 & internal <= 45))
  ## parameters and changepoints are aligned in length
  expect_identical(length(fit$parameters), length(fit$changepoints))
})
