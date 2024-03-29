% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/graph.R
\name{graph}
\alias{graph}
\title{Graph generation}
\usage{
graph(
  ...,
  type = "empty",
  decay = 1,
  gap = 0,
  penalty = 0,
  K = Inf,
  a = 0,
  all.null.edges = FALSE
)
}
\arguments{
\item{...}{This is a list of edges definied by functions \code{Edge}, \code{StartEnd} and \code{Node}. See gfpop functions \code{\link[gfpop:Edge]{gfpop::Edge()}}, \code{\link[gfpop:StartEnd]{gfpop::StartEnd()}} and  \code{\link[gfpop:Node]{gfpop::Node()}}}

\item{type}{a string equal to \code{"std"}, \code{"isotonic"}, \code{"updown"} or \code{"relevant"} to build a predefined classic graph}

\item{decay}{a nonnegative number to give the strength of the exponential decay into the segment}

\item{gap}{a nonnegative number to constrain the size of the gap in the change of state}

\item{penalty}{a nonnegative number equals to the common penalty to use for all edges}

\item{K}{a positive number. Threshold for the Biweight robust loss}

\item{a}{a positive number. Slope for the Huber robust loss}

\item{all.null.edges}{a boolean. Add null edges to all nodes automatically}
}
\value{
a dataframe with 9 variables :
columns are named \code{"state1"}, \code{"state2"}, \code{"type"}, \code{"parameter"}, \code{"penalty"}, \code{"K"}, \code{"a"}, \code{"min"}, \code{"max"} with additional \code{"graph"} class.
}
\description{
Graph creation using component functions \code{Edge}, \code{StartEnd} and \code{Node}
}
\examples{
graph(type = "updown", gap = 1.3, penalty = 5)

graph(Edge("Dw","Dw"),
      Edge("Up","Up"),
      Edge("Dw","Up","up", gap = 0.5, penalty = 10),
      Edge("Up","Dw","down"),
      StartEnd("Dw","Dw"),
      Node("Dw",0,1),
      Node("Up",0,1))

graph(Edge("1", "2", type = "std"),
      Edge("2", "3", type = "std"),
      Edge("3", "4", type = "std"),
      StartEnd(start = "1", end = "4"),
      all.null.edges = TRUE)
}
