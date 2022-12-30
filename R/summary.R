##  GPL-3 License
## Copyright (c) 2022 Vincent Runge

#' summary.gfpop
#' @description Returns all information contains in gfpop object
#' @param object a gfpop class object
#' @param ... Other parameters
#' @return summary of gfpop object
#' @examples
#' n <- 1000 #data length
#' data <- dataGenerator(n, c(0.3, 0.4, 0.7, 0.95, 1), c(1, 3, 1, -1, 4), "mean", sigma = 3)
#' myGraph <- graph(type = "relevant", gap = 0.5, penalty = 2 * sdDiff(data) ^ 2 * log(n))
#' g <- gfpop(data, myGraph, type = "mean")
#' summary(g)
summary.gfpop <- function(object,...)
{
  cat("change-point positions:", object$changepoints, "\n")
  cat("segment parameter values:", object$parameters, "\n")
  cat("number of segments:", length(object$changepoints), "\n")
  cat("global cost:", object$globalCost, "\n")
  cat("segment states:", object$states, "\n")
  cat("state transition status (constrained or not T/F):", object$forced, "\n")
}

