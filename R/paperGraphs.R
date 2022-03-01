##  GPL-3 License
## Copyright (c) 2022 Vincent Runge

#' Graphs of our paper in JSS
#'
#' @description Graphs of our paper in JSS (Journal of Statistical Software)
#' @param nb the number of the Figure in paper
#' @param penalty the penalty to use for change-point
#' @param decay a nonnegative number to give the strength of the exponential decay into the segment
#' @param gap a nonnegative number to constrain the size of the gap in the change of state
#' @param oneValue the value for parameter when we consider the collective anomalies problem
#' @param K a positive number. Threshold for the Biweight robust loss
#' @param a a positive number. Slope for the Huber robust loss
#' @return a dataframe with 9 variables (columns are named \code{"state1", "state2", "type", "parameter", "penalty", "K", "a", "min", "max"}) with additional \code{"graph"} class.
#'
paperGraph <- function(nb, penalty = 0, decay = 1, gap = 0, oneValue = 0, K = Inf, a = 0)
{
  mygraph <- NULL
  if(nb == 4)
  {
    mygraph <- graph(type = "std", penalty = penalty)
  }
  if(nb == 5)
  {
    mygraph <- graph(type = "isotonic", penalty = penalty)
  }
  if(nb == 6)
  {
    mygraph <- graph(
      Edge("Dw", "Dw", type = "null", decay = decay),
      Edge("Up", "Up", type = "null", decay = decay),
      Edge("Dw", "Up", type = "up", penalty = penalty),
      Edge("Up", "Dw", type = "down", penalty = penalty),
      StartEnd(start = "Dw", end = "Dw")
    )
  }
  if(nb == 7)
  {
    mygraph <- graph(type = "isotonic", decay = decay, gap = gap, penalty = penalty, K = K, a = a)
  }
  if(nb == 8)
  {
    mygraph <- graph(
      Edge("1", "2", type = "std", penalty = 0),
      Edge("2", "3", type = "std", penalty = 0),
      Edge("1", "1", type = "null", penalty = 0),
      Edge("2", "2", type = "null", penalty = 0),
      Edge("3", "3", type = "null", penalty = 0),
      StartEnd(start = "1", end = "3")
    )
  }
  if(nb == 9)
  {
    mygraph <- graph(
      Edge("Wait", "Seg", type = "null", penalty = 0),
      Edge("Seg", "Seg", type = "null", penalty = 0),
      Edge("Seg", "Wait", type = "std", penalty = penalty),
      StartEnd(start = "Wait", end = "Seg")
    )
  }
  if(nb == 11 || nb == 19)
  {
    mygraph <- graph(type = "relevant", decay = decay, gap = gap, penalty = penalty)
  }
  if(nb == 12)
  {
    mygraph <- graph(
      Edge("Dw", "Dw", type = "down", penalty = 0),
      Edge("Up", "Up", type = "up",  penalty =  0),
      Edge("Dw", "Up", type = "up",  penalty = penalty),
      Edge("Up", "Dw", type = "down",  penalty = 0)
      )
  }
  if(nb == 13)
  {
    mygraph <- graph(
      Edge(0, 1, "down", penalty = 8e7, gap=0),
      Edge(1, 2, "up", penalty = 0, gap = 2000),
      Edge(2, 3, "down", penalty = 0, gap = 5000),
      Edge(3, 4, "up", penalty = 0, gap = 2000),
      Edge(4, 5, "up", penalty = 0, gap = 1000),
      Edge(5, 6, "up", penalty = 0, gap = 0),
      Edge(6, 7, "down", penalty = 0, gap = 0),
      Edge(7, 8, "down", penalty = 0, gap = 0),
      Edge(8, 0, "up", penalty =0, gap = 0),
      all.null.edges = TRUE
      )
  }
  ###
  if(nb == 17)
  {
    mygraph <- graph(
      Edge("mu0", "mu0", type = "null", penalty = 0, K = K),
      Edge("mu0", "Coll", type = "std", penalty = penalty),
      Edge("Coll", "Coll", type = "null"),
      Edge("Coll", "mu0", type = "std",  K = K),
      StartEnd(start = "mu0", end = c("mu0", "Coll")),
      Node("mu0", min = oneValue, max = oneValue)
    )
  }
  if(nb == 18)
  {
    mygraph <- graph(
      Edge("Wait1", "Wait2", type = "null", penalty = 0),
      Edge("Wait2", "Seg", type = "null", penalty = 0),
      Edge("Seg", "Seg", type = "null", penalty = 0),
      Edge("Seg", "Wait1", type = "std", penalty = penalty),
      StartEnd(start = "Wait1", end = "Seg")
    )
  }

  if(nb == 20)
  {
    mygraph <- graph(
      Edge("Dw", "Dw", type = "null", penalty = 0),
      Edge("Dw", "Up", type = "up", penalty = penalty, gap = gap),
      Edge("Up", "Up", type = "null", penalty = 0),
      Edge("Up", "Dw", type = "down", penalty = penalty, gap = gap),
      StartEnd(start = "Dw", end = "Dw")
    )
  }
  if(nb == 21)
  {
    mygraph <- graph(
      Edge("Dw1", "Dw", type = "null", penalty = 0),
      Edge("Dw", "Dw", type = "null", penalty = 0),
      Edge("Dw", "Up1", type = "up", penalty = penalty),
      Edge("Up1", "Up", type = "null", penalty = 0),
      Edge("Up", "Up", type = "null", penalty = 0),
      Edge("Up", "Dw1", type = "down", penalty = penalty),
      StartEnd(start = "Dw1", end = "Dw")
    )
  }
  if(nb == 22)
  {
    mygraph <- graph(
      Edge("Dw", "Dw", type = "down", penalty = 0),
      Edge("Dw", "Up", type = "up", penalty = penalty),
      Edge("Up", "Dw", type = "down", penalty = 0),
      StartEnd(start = "Dw", end = "Dw")
    )
  }
  if(nb == 23)
  {
    mygraph <- graph(
      Edge("Dw", "Dw", type = "down", penalty = 0),
      Edge("Dw", "Up", type = "up", penalty = penalty),
      Edge("Up", "Up", type = "up", penalty = 0),
      Edge("Up", "Dw", type = "down", penalty = penalty),
      StartEnd(start = "Dw", end = "Dw")
    )
  }

  return(mygraph)
}

