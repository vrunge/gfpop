#' Edge generation
#'
#' @description Edge creation
#' @param state1 a nonnegative integer defining the starting state of the edge
#' @param state2 a nonnegative integer defining the ending state of the edge
#' @param type a string equals to "std", "up", "down", "absInf" or "absSup"
#' @param penalty a nonnegative number. The penalization of this change of state
#' @param parameter a nonnegative number to constraint the size of the gap in the change of state
#' @return a list (with the additional "edge" class) with five components equal to the five parameters
#' @examples
#' edge(0, 1, "up", 10, 1)
edge <- function(state1, state2, type, penalty, parameter = 0)
{
  ###STOP###
  if(state1%%1 != 0){stop('state1 is not an integer.')}
  if(state2%%1 != 0){stop('state1 is not an integer.')}
  if(state1 < 0){stop('state1 must be nonnegative')}
  if(state2 < 0){stop('state2 must be nonnegative')}

  if(type != "std" && type != "up" && type != "down" && type != "absInf" && type != "absSup")
    {stop('Argument not appropriate. Choose a type among the following: "std", "up", "down", "absInf", "absSup".')}

  if(!is.double(penalty)){stop('penalty is not a double.')}
  if(!is.double(parameter)){stop('parameter is not an double.')}

  if(penalty < 0){stop('penalty must be nonnegative')}
  if(parameter < 0){stop('parameter must be nonnegative')}

  ###myNewEdge is a list of class edge
  myNewEdge <- list(state1, state2, type, penalty, parameter)

  class(myNewEdge) <- c(class(myNewEdge), "edge")
  return(myNewEdge)
}

###############################################

#' Graph generation
#'
#' @description Graph creation
#' @param penalty a nonnegative number equals to the common penalty to use for all edges
#' @param type a string equal to "std", "isotonic", "updown", "infsup". to build a predefined classic graph
#' @return a dataframe with edges in rows (columns are named "state1", "state2", "type", "penalty", "parameter") with additional "graph" class.
#' @examples
#' myGraph <- graph(penalty = 10, "updown")
#' myEmptyGraph <- graph()
graph <- function(penalty = 0, type = "empty")
{
  ###STOP###
  if(!is.double(penalty)){stop('penalty is not a double.')}
  if(penalty < 0){stop('penalty must be nonnegative')}

  if(type != "empty" && type != "std" && type != "isotonic" && type != "updown" && type != "infsup")
    {stop('Arugment "type" not appropriate. Choose among "std", "isotonic", "updown", "infsup"')}

  ###Empty Graph###
  myNewGraph <- data.frame(numeric(), numeric(0), character(), numeric(0), numeric(0), stringsAsFactors = FALSE)
  names(myNewGraph)<- c("state1", "state2", "type", "penalty", "parameter")
  class(myNewGraph) <- c(class(myNewGraph), "graph")

  ###Usual graphs###
  if(type == "std") {myNewGraph[1, ] <- edge(0, 0, "std", penalty, 0)}

  if(type == "isotonic"){myNewGraph[1, ] <- edge(0, 0, "up", penalty, 0)}

  if(type == "updown")
  {
    myNewGraph[1, ] <- edge(0, 1, "up", penalty, 0)
    myNewGraph[2, ] <- edge(1, 0, "down", penalty, 0)
  }

  if(type == "infsup")
  {
    myNewGraph[1, ] <- edge(0, 1, "absInf", penalty, 0)
    myNewGraph[2, ] <- edge(1, 0, "absSup", penalty, 0)
  }
  return(myNewGraph)
}

###############################################

#' Adding edge to graph
#'
#' @description Adding edge to a graph
#' @param graph a dataframe of class graph
#' @param edge a vector of class edge
#' @return the graph with the additional edge "edge"
#' @examples
#' myGraph <- graph()
#' myGraph <- addEdge(myGraph, edge(0, 1, "up", 10))
#' myGraph <- addEdge(myGraph, edge(1, 0, "down", 0))
addEdge <- function(graph, edge)
{
  if(!any(class(graph) == "graph")){stop('Your graph is not a graph...')}
  if(!any(class(edge) == "edge")){stop('Your edge is not an edge...')}

  if(nrow(graph) > 0)
    {for(i in 1:nrow(graph)){if(all(graph[i,1:3] == edge[1:3])){stop('This edge is already the graph.')}}}

  graph[nrow(graph) + 1, ] <- edge

  return(graph)
}

###############################################

#' Constraint starting and ending states to a graph
#'
#' @description Adding constraints on the starting and ending states to a graph
#' @param graph a dataframe of class graph
#' @param start a nonnegative integer. The first state contrained in the changepoint inference
#' @param end a nonnegative integer. The end state contrained in the changepoint inference
#' @return the graph with these new constraints
#' @examples
#' myGraph <- graph()
#' myGraph <- addEdge(myGraph, edge(0, 1, "up", 10))
#' myGraph <- addEdge(myGraph, edge(1, 0, "down", 0))
#' myGraph <- addStartEnd(myGraph, 0, 0)
addStartEnd <- function(graph, start = - 1, end = - 1)
{
  ###STOP###
  if(!any(class(graph) == "graph")){stop('Your graph is not a graph...')}

  if(start != -1){if(any(graph[,3] == "start")){stop('You already have a start edge in your graph.')}}
  if(end != -1){if(any(graph[,3] == "end")){stop('You already have an end edge in your graph.')}}

  if(start != -1){if(!any(graph[,1:2] == start)){stop('Your start edge is not in your graph.')}}
  if(end != -1){if(!any(graph[,1:2] == end)){stop('Your end edge is not in your graph.')}}

  ###Add the start/end edges###
  if(start != -1){graph[nrow(graph) + 1, ] <- list(start, NA, "start", NA, NA)}
  if(end != -1){graph[nrow(graph) + 1, ] <- list(end, NA, "end", NA, NA)}

  return(graph)
}

