#' Edge generation
#'
#' @description Edge creation for graph
#' @param state1 a nonnegative integer defining the starting state of the edge
#' @param state2 a nonnegative integer defining the ending state of the edge
#' @param type a string equal to "null", "std", "up", "down", "absInf" or "absSup"
#' @param penalty a nonnegative number. The penality associated to this state transition
#' @param decay a nonnegative number to give the strength of the exponential decay into the segment
#' @param gap a nonnegative number to constraint the size of the gap in the change of state
#' @return a dataframe with five components equal to the five parameters
#' @examples
#' edge(0, 1, "up", 10, gap = 1)
edge <- function(state1, state2, type = "null", penalty = 0, decay = 1, gap = 0)
{
  ###STOP###
  if(state1%%1 != 0){stop('state1 is not an integer.')}
  if(state2%%1 != 0){stop('state1 is not an integer.')}
  if(state1 < 0){stop('state1 must be a nonnegative integer')}
  if(state2 < 0){stop('state2 must be a nonnegative integer')}

  if(type == "null" && (state1 != state2))
    {stop('You can not build a null edge between two different states".')}

  if(type != "null" && type != "std" && type != "up" && type != "down" && type != "absInf" && type != "absSup")
    {stop('Argument not appropriate. Choose a type among the following: "null", "std", "up", "down", "absInf", "absSup".')}



  if(!is.double(penalty)){stop('penalty is not a double.')}
  if(!is.double(decay)){stop('decay is not a double.')}
  if(!is.double(gap)){stop('gap is not a double.')}

  if(penalty < 0){stop('penalty must be nonnegative')}
  if(decay < 0){stop('decay must be nonnegative')}
  if(gap < 0){stop('gap must be nonnegative')}

  if(type == "null"){parameter <- decay}else{parameter <- gap}

  ###response = a dataframe with a unique row
  df <- data.frame(state1, state2, type, penalty, parameter, stringsAsFactors = FALSE)
  return(df)
}


###############################################

#' Constraint the starting and ending states of a graph
#'
#' @description Adding constraints on the starting and ending states of a graph
#' @param start a vector of nonnegative integers. The starting vertices in the changepoint inference
#' @param end a vector of nonnegative integers. The ending vertices in the changepoint inference
#' @return a dataframe with five components (as for edge) with only `state1`` and `type` = start or end defined.
#' @examples
#' StartEnd(start = 0, end = c(1,2))

StartEnd <- function(start = NULL, end = NULL)
{
  df <- data.frame(numeric(), numeric(0), character(), numeric(0), numeric(0), stringsAsFactors = FALSE)
  colnames(df) <- c("state1", "state2", "type", "penalty", "parameter")
  if(length(start) != 0)
  {
    for(i in 1:length(start))
    {
      if(start[i]%%1 != 0){stop('The vector start ontains a non-integer element')}
      if(start[i] < 0){stop('The vector start contains a negative integer')}
      df[i,] <- list(start[i], NA, "start", NA, NA)
    }
  }
  if(length(end) != 0)
  {
    for(i in 1:length(end))
    {
      if(end[i]%%1 != 0){stop('The vector end contains a non-integer element')}
      if(end[i] < 0){stop('The vector end contains a negative integer')}
      df[i + length(start),] <- list(end[i], NA, "end", NA, NA)
    }
  }
  return(df)
}


###############################################

#' Graph generation
#'
#' @description Graph creation
#' @param ... This is a list of edges definied by functions edge and StartEnd
#' @param penalty a nonnegative number equals to the common penalty to use for all edges
#' @param type a string equal to "std", "isotonic", "updown", "infsup". to build a predefined classic graph
#' @return a dataframe with edges in rows (columns are named "state1", "state2", "type", "penalty", "parameter") with additional "graph" class.
#' @examples
#' UpDownGraph <- graph(penalty = 10, type = "updown")
#' MyGraph <- graph(edge(0,0), edge(1,1), edge(0,1,"up",10,gap=0.5), edge(1,0,"down"), StartEnd(0,0))

graph <- function(..., penalty = 0, type = "empty")
{
  myNewGraph <- rbind(...)

  if(is.null(myNewGraph) == TRUE)
  {
    ###STOP###
    if(!is.double(penalty)){stop('penalty is not a double.')}
    if(penalty < 0){stop('penalty must be nonnegative')}

    if(type != "empty" && type != "std" && type != "isotonic" && type != "updown" && type != "infsup")
      {stop('Arugment "type" not appropriate. Choose among "std", "isotonic", "updown", "infsup"')}

    myNewGraph <- data.frame(numeric(), numeric(0), character(), numeric(0), numeric(0), stringsAsFactors = FALSE)
    names(myNewGraph) <- c("state1", "state2", "type", "penalty", "parameter")

    if(type == "std")
    {
      myNewGraph[1, ] <- edge(0, 0, "null")
      myNewGraph[2, ] <- edge(0, 0, "std", penalty, 0)
    }

    if(type == "isotonic")
    {
      myNewGraph[1, ] <- edge(0, 0, "null")
      myNewGraph[2, ] <- edge(0, 0, "up", penalty, 0)
    }

    if(type == "updown")
    {
      myNewGraph[1, ] <- edge(0, 0, "null")
      myNewGraph[2, ] <- edge(1, 1, "null")
      myNewGraph[3, ] <- edge(0, 1, "up", penalty, 0)
      myNewGraph[4, ] <- edge(1, 0, "down", penalty, 0)
    }

    if(type == "infsup")
    {
      myNewGraph[1, ] <- edge(0, 0, "null")
      myNewGraph[2, ] <- edge(1, 1, "null")
      myNewGraph[3, ] <- edge(0, 1, "absSup", penalty, 0)
      myNewGraph[4, ] <- edge(1, 0, "absInf", penalty, 0)
    }
  }

  class(myNewGraph) <- c("data.frame", "graph")

  return(myNewGraph)
}
