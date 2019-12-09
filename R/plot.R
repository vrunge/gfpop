
#' plot.gfpop
#' @description Plot the result of the gfpop function with the data-points
#' @param x a gfpop class object
#' @param ... Other parameters
#' @param data the data from which we get the gfpop result
#' @return plot data and the inferred gfpop segments
plot.gfpop <- function(x, ..., data)
{
  n <- 1:length(data)
  p <- length(x$changepoints)
  xbis <- c(1, x$changepoints)
  y <- x$parameters
  lim <- c(min(min(x$parameters,data)), max(max(x$parameters,data)))

  plot(1:length(data), data, pch = '+', ylim = lim)
  for(i in 1:p)
  {
    segments(xbis[i], y[i], xbis[i+1], y[i], col= 2, lty = 1, lwd = 5)
  }
}

