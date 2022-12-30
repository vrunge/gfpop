##  GPL-3 License
## Copyright (c) 2022 Vincent Runge

#' print.gfpop
#' @description Printing the resulting change-point model found by gfpop function
#' @param x a gfpop class object
#' @param ... Other parameters
#' @return print the gfpop object
#' @examples
#' n <- 1000 #data length
#' data <- dataGenerator(n, c(0.3, 0.4, 0.7, 0.95, 1), c(1, 3, 1, -1, 4), "mean", sigma = 3)
#' myGraph <- graph(type = "relevant", gap = 0.5, penalty = 2 * sdDiff(data) ^ 2 * log(n))
#' g <- gfpop(data, myGraph, type = "mean")
#' print(g)
print.gfpop <- function(x,...)
{
  cat("$changepoints \n")
  print(x$changepoints)
  cat("$states \n")
  print(x$states)
  cat("$forced \n")
  print(x$forced)
  cat("$globalCost \n")
  print(x$globalCost)
}






