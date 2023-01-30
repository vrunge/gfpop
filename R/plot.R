##  GPL-3 License
## Copyright (c) 2022 Vincent Runge

#' plot.gfpop
#' @description Plotting inferred segment parameters (the result of gfpop) and data.
#' @param x a gfpop class object
#' @param ... Other parameters
#' @param data the data from which we get the gfpop result
#' @param multiple if \code{TRUE} we plot data and the model on different graphs.
#' Only with \code{"mean"} and \code{"poisson"} cost functions (as in that case the parameter
#' values represent the data mean value over each segment) we allow the User
#' to plot signal and data on a single graph.
#' @return plot data and the inferred gfpop segments
#' @examples
#' n <- 1000 #data length
#' data <- dataGenerator(n, c(0.3, 0.4, 0.7, 0.95, 1), c(1, 3, 1, -1, 4), "mean", sigma = 3)
#' myGraph <- graph(type = "relevant", gap = 0.5, penalty = 2 * sdDiff(data) ^ 2 * log(n))
#' g <- gfpop(data, myGraph, type = "mean")
#' plot(x = g, data = data, multiple = FALSE)
#'
#' data <- dataGenerator(n, c(0.4, 0.8, 1), c(1, 1.7, 2.3), "exp")
#' g <- gfpop(data,graph(type = "isotonic", penalty = 2 * sdDiff(data) ^ 2 * log(n)), type = "exp")
#' plot(x = g, data = data, multiple = TRUE)
#'
#' data <- dataGenerator(n, c(0.22, 0.75, 1), c(1.4,1,0.8), "poisson")
#' g <- gfpop(data, paperGraph(8), type = "poisson")
#' plot(x = g, data = data, multiple = TRUE)
plot.gfpop <- function(x, ..., data, multiple = TRUE)
{
  oparops <- par(no.readonly = TRUE)
  on.exit(par(oparops))
  par(mar = c(3, 3, 2, 2))

  n <- 1:length(data)
  p <- length(x$changepoints)
  xbis <- c(1, x$changepoints)
  y <- x$parameters
  limParam <- c(min(x$parameters), max(x$parameters))
  limData <- c(min(data), max(data))

  type <- attributes(x)$type

  if(!any(c("mean","poisson") %in% type) && (multiple == FALSE))
    {stop('Only with "mean" and "poisson" cost type functions can the User plot signal and data on a unique graph.')}

  if(!any(c("mean","poisson") %in% type) || (multiple == TRUE))
  {
    par(mfrow = c(2,1))
    plot(1:length(data), data, pch = '+', ylim = limData, xlab = "", ylab = "", main = "data & changepoints")
    for(i in 1:p){abline(v = xbis[i+1], col= "#3514B3", lty = 1, lwd = 2)}
    plot(1:length(data), ylim =limParam, xlab = "", ylab = "", col = "white", main = "parameters & changepoints")
    for(i in 1:p){segments(xbis[i], y[i], xbis[i+1], y[i], col= "#C73617", lty = 1, lwd = 4)}
    for(i in 1:p){abline(v = xbis[i+1], col= "#3514B3", lty = 1, lwd = 2)}
   }
  else
  {
    plot(1:length(data), data, pch = '+', ylim = c(min(limParam[1], limData[1]), max(limParam[2], limData[2])), xlab = "", ylab = "", main = "data, parameters & changepoints")
    for(i in 1:p){segments(xbis[i], y[i], xbis[i+1], y[i], col= "#C73617", lty = 1, lwd = 4)}
    for(i in 1:p){abline(v = xbis[i+1], col= "#3514B3", lty = 1, lwd = 2)}
  }

}
