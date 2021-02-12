library(gfpop)
library(gfpop.data)
library(testthat)
context("sn")

library(data.table)
data(profile614chr2, package="gfpop.data")

### reduce data size
profile614chr2$probes <- profile614chr2$probes[1:10000,]

g <- gfpop::graph(type="std")
x <- rnorm(10)
gfpop::gfpop(x, g)

######
###### sngraph function
######

sngraph <- function(n.segs, type, gap)
{
  stopifnot(is.integer(n.segs), length(n.segs)==1, n.segs >= 1)
  s <- n.segs-1
  null.vec <- 0:s
  seg.vec <- 1:s
  edge.df <- gfpop::Edge(paste0("seg", seg.vec-1), paste0("seg", seg.vec), type, gap=gap)
  gfpop::graph(
    gfpop::StartEnd(start=paste0("seg", 0), end=paste0("seg", s)),
    gfpop::Edge(paste0("seg", null.vec), paste0("seg", null.vec), "null"),
    if(1 < n.segs)edge.df else edge.df[0,],
    all.null.edges=FALSE)
}

######
###### tests with abs, 1,2,3 segments
######

test_that("abs model with 1 segment returned", {
  g1 <- sngraph(1L, "abs", 1)
  fit1 <- gfpop::gfpop(
    profile614chr2$probes$logratio,
    mygraph = g1, type = "mean")
  expect_identical(fit1$changepoints, nrow(profile614chr2$probes))
  expect_equal(fit1$parameters, mean(profile614chr2$probes$logratio))
  expect_identical(fit1$states, "seg0")
})

test_that("abs model with 2 segments returned", {
  g2 <- sngraph(2L, "abs", 1)
  fit2 <- gfpop::gfpop(
    profile614chr2$probes$logratio,
    mygraph = g2, type = "mean")
  expect_identical(length(fit2$changepoints), 2L)
  expect_identical(fit2$changepoints[2], nrow(profile614chr2$probes))
  expect_identical(length(fit2$parameters), 2L)
  expect_identical(fit2$states, c("seg0", "seg1"))
})

test_that("abs model with 3 segments returned", {
  g3 <- sngraph(3L, "abs", 1)
  fit3 <- gfpop::gfpop(
    profile614chr2$probes$logratio,
    mygraph = g3, type = "mean")
  expect_identical(length(fit3$changepoints), 3L)
  expect_identical(fit3$changepoints[3], nrow(profile614chr2$probes))
  expect_identical(length(fit3$parameters), 3L)
  expect_identical(fit3$states, c("seg0", "seg1", "seg2"))
})

######
###### tests with std, 1,2,3 segments
######

test_that("std model with 1 segment returned", {
  g1 <- sngraph(1L, "std", 0)
  fit1 <- gfpop::gfpop(
    profile614chr2$probes$logratio,
    mygraph = g1, type = "mean")
  expect_identical(fit1$changepoints, nrow(profile614chr2$probes))
  expect_equal(fit1$parameters, mean(profile614chr2$probes$logratio))
  expect_identical(fit1$states, "seg0")
})

test_that("std model with 2 segments returned", {
  g2 <- sngraph(2L, "std", 0)
  fit2 <- gfpop::gfpop(
    profile614chr2$probes$logratio,
    mygraph = g2, type = "mean")
  expect_identical(length(fit2$changepoints), 2L)
  expect_identical(fit2$changepoints[2], nrow(profile614chr2$probes))
  expect_identical(length(fit2$parameters), 2L)
  expect_identical(fit2$states, c("seg0", "seg1"))
})

test_that("std model with 3 segments returned", {
  g3 <- sngraph(3L, "std", 0)
  fit3 <- gfpop::gfpop(
    profile614chr2$probes$logratio,
    mygraph = g3, type = "mean")
  expect_identical(length(fit3$changepoints), 3L)
  expect_identical(fit3$changepoints[3], nrow(profile614chr2$probes))
  expect_identical(length(fit3$parameters), 3L)
  expect_identical(fit3$states, c("seg0", "seg1", "seg2"))
})
