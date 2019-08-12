library(gfpop)
library(testthat)
context("sn")

library(data.table)
data(profile614chr2, package="gfpop")

sngraph <- function(n.segs, type, gap){
  stopifnot(is.integer(n.segs), length(n.segs)==1, n.segs >= 1)
  s <- n.segs-1
  null.vec <- 0:s
  seg.vec <- 1:s
  gfpop::graph(
    gfpop::StartEnd(start=0, end=s),
    gfpop::Edge(null.vec, null.vec, "null"),
    if(1 < n.segs)gfpop::Edge(seg.vec-1, seg.vec, type, gap=gap),
    all.null.edges=FALSE)
}

test_that("absSup model with 1 segment returned", {
  g1 <- sngraph(1L, "absSup", 1)
  fit1 <- gfpop::gfpop(
    profile614chr2$probes$logratio,
    mygraph = g1, type = "mean")
  expect_identical(fit1$changepoints, nrow(profile614chr2$probes))
  expect_equal(fit1$parameters, mean(profile614chr2$probes$logratio))
})

test_that("absSup model with 2 segments returned", {
  g2 <- sngraph(2L, "abs", 1)
  fit2 <- gfpop::gfpop(
    profile614chr2$probes$logratio,
    mygraph = g2, type = "mean")
  expect_identical(length(fit2$changepoints), 2L)
})

test_that("absSup model with 3 segments returned", {
  g3 <- sngraph(3L, "abs", 1)
  fit3 <- gfpop::gfpop(
    profile614chr2$probes$logratio,
    mygraph = g3, type = "mean")
  expect_identical(length(fit3$changepoints), 3L)
})

test_that("std model with 1 segment returned", {
  g1 <- sngraph(1L, "std", 0)
  fit1 <- gfpop::gfpop(
    profile614chr2$probes$logratio,
    mygraph = g1, type = "mean")
  expect_identical(fit1$changepoints, nrow(profile614chr2$probes))
  expect_equal(fit1$parameters, mean(profile614chr2$probes$logratio))
})

test_that("std model with 2 segments returned", {
  g2 <- sngraph(2L, "std", 0)
  fit2 <- gfpop::gfpop(
    profile614chr2$probes$logratio,
    mygraph = g2, type = "mean")
  expect_identical(length(fit2$changepoints), 2L)
})

test_that("std model with 3 segments returned", {
  g3 <- sngraph(3L, "std", 0)
  fit3 <- gfpop::gfpop(
    profile614chr2$probes$logratio,
    mygraph = g3, type = "mean")
  expect_identical(length(fit3$changepoints), 3L)
})

sngraph2 <- function(n.segs, type, gap){
  stopifnot(is.integer(n.segs), length(n.segs)==1, n.segs >= 1)
  s <- n.segs-1
  seg.vec <- 1:s
  gfpop::graph(
    gfpop::StartEnd(start=0, end=s),
    if(1 < n.segs)gfpop::Edge(seg.vec-1, seg.vec, type, gap=gap),
    all.null.edges=TRUE)
}

test_that("all.null=TRUE std model with 2 segments returned", {
  g2 <- sngraph2(2L, "std", 0)
  fit2 <- gfpop::gfpop(
    profile614chr2$probes$logratio,
    mygraph = g2, type = "mean")
  expect_identical(length(fit2$changepoints), 2L)
})

test_that("all.null=TRUE std model with 3 segments returned", {
  g3 <- sngraph2(3L, "std", 0)
  fit3 <- gfpop::gfpop(
    profile614chr2$probes$logratio,
    mygraph = g3, type = "mean")
  expect_identical(length(fit3$changepoints), 3L)
})
