library(testthat)
library(gfpop)
context("missing")

set.seed(1)
N <- 10
no.missing <- c(rnorm(N, -2), rnorm(N, 2))
g <- gfpop::graph(type="std", penalty=0)
test_that("gfpop works fine without missing data", {
  fit <- gfpop::gfpop(no.missing, g)
  expect_equal(fit$parameters, no.missing)
})

one.missing <- c(NA, no.missing)
test_that("gfpop error for missing data", {
  expect_error({
    gfpop::gfpop(one.missing, g)
  }, "data has missing values, please remove them")
})
