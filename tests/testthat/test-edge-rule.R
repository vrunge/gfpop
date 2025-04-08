library(testthat)
library(gfpop)
context("edge-rule")

test_that("edge with default rule works", {
  e <- Edge("A", "B")
  expect_equal(e$rule, 1)
})

test_that("edge with custom rule works", {
  e <- Edge("A", "B", rule = 2)
  expect_equal(e$rule, 2)
})

test_that("edge rule is vectorized", {
  edges <- Edge(c("A", "B"), c("B", "C"), rule = c(1, 2))
  expect_equal(edges$rule, c(1, 2))
})

test_that("invalid rule values cause error", {
  expect_error(Edge("A", "B", rule = -1), "rule must be positive")
  expect_error(Edge("A", "B", rule = "1"), "rule must be numeric")
  expect_error(Edge("A", "B", rule = 1.5), "rule must be a positive integer")
})

test_that("StartEnd includes rule column", {
  se <- StartEnd("A", "B")
  expect_true("rule" %in% names(se))
})

test_that("Node includes rule column", {
  n <- Node("A", 0, 1)
  expect_true("rule" %in% names(n))
})

test_that("graph construction works with rule", {
  g <- graph(Edge("A", "B", rule = 2))
  expect_true("rule" %in% names(g))
  expect_equal(g$rule, 2)
})

test_that("graph combines edges with different rules", {
  g <- graph(Edge("A", "B", rule = 2), Edge("B", "C", rule = 3))
  expect_equal(g$rule, c(2, 3))
})

# Additional robust test cases

test_that("large rule values work correctly", {
  e <- Edge("A", "B", rule = 1000000)
  expect_equal(e$rule, 1000000)
})

test_that("rules work with all edge types", {
  edge_types <- c("null", "std", "up", "down", "abs")
  for (type in edge_types) {
    e <- Edge("A", "B", type = type, rule = 5)
    expect_equal(e$rule, 5)
  }
})

test_that("rules are preserved in predefined graphs", {
  g <- graph(type = "updown")
  expect_true("rule" %in% names(g))
  expect_equal(length(unique(g$rule)), 1) # All should have default rule=1
  
  # Test isotonic graph
  g <- graph(type = "isotonic")
  expect_true("rule" %in% names(g))
  expect_equal(length(unique(g$rule)), 1)
})

test_that("all.null.edges preserves rule values", {
  g <- graph(Edge("A", "B", rule = 3), Edge("B", "C", rule = 4), all.null.edges = TRUE)
  expect_true("rule" %in% names(g))
  # The null edges should have rule=1 (default), others should keep their values
  expect_true(1 %in% g$rule)
  expect_true(3 %in% g$rule)
  expect_true(4 %in% g$rule)
})

test_that("complex graph with multiple rules works", {
  g <- graph(
    Edge("A", "B", "up", rule = 2),
    Edge("B", "C", "down", rule = 3),
    Edge("C", "A", "std", rule = 4),
    StartEnd(start = "A", end = "C")
  )
  expect_equal(sort(unique(g[g$type %in% c("up", "down", "std"), "rule"])), c(2, 3, 4))
})
