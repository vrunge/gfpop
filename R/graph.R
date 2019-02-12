#' Edge generation
#'
#' @description Edge creation for graph
#' @param state1 a nonnegative integer defining the starting state of the edge
#' @param state2 a nonnegative integer defining the ending state of the edge
#' @param type a string equal to "null", "std", "up", "down", "absInf" or "absSup"
#' @param penalty a nonnegative number. The penality associated to this state transition
#' @param parameter a nonnegative number to constraint the size of the gap in the change of state
#' @return a dataframe with five components equal to the five parameters
#' @examples
#' edge(0, 1, "up", 10, 1)
edge <- function(state1, state2, type = "null", penalty = 0, parameter = 0)
{
  ###STOP###
  if(state1%%1 != 0){stop('state1 is not an integer.')}
  if(state2%%1 != 0){stop('state1 is not an integer.')}
  if(state1 < 0){stop('state1 must be a nonnegative integer')}
  if(state2 < 0){stop('state2 must be a nonnegative integer')}

  if(type != "null" && type != "std" && type != "up" && type != "down" && type != "absInf" && type != "absSup")
    {stop('Argument not appropriate. Choose a type among the following: "null", "std", "up", "down", "absInf", "absSup".')}

  if(!is.double(penalty)){stop('penalty is not a double.')}
  if(!is.double(parameter)){stop('parameter is not a double.')}

  if(penalty < 0){stop('penalty must be nonnegative')}
  if(parameter < 0){stop('parameter must be nonnegative')}

  ###response = a dataframe with a unique row
  df <- data.frame(state1, state2, type, penalty, parameter, stringsAsFactors = FALSE)
  return(df)
}


###############################################

#' Constraint starting and ending states to a graph
#'
#' @description Adding constraints on the starting and ending states to a graph
#' @param start a nonnegative integer. The first vertex in the changepoint inference
#' @param end a nonnegative integer. The end vertex in the changepoint inference
#' @return a dataframe with five components (as for edge) with only state1 and type = start or end defined.
#' @examples
#' StartEnd(start = 0, end = 1)

StartEnd <- function(start = -1, end = -1)
{
  if(start%%1 != 0){stop('start is not an integer.')}
  if(end%%1 != 0){stop('end is not an integer.')}
  if(start < -1){stop('start must be a nonnegative integer')}
  if(end < -1){stop('end must be a nonnegative integer')}

  df1 <- data.frame(start, NA, "start", NA, NA, stringsAsFactors = FALSE)
  df2 <- data.frame(end, NA, "end", NA, NA, stringsAsFactors = FALSE)
  colnames(df1) <- c("state1", "state2", "type", "penalty", "parameter")
  colnames(df2) <- c("state1", "state2", "type", "penalty", "parameter")
  df <- NULL
  if(start != -1){df <- df1}
  if(end != -1){df <- rbind(df,df2)}
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
#' MyGraph <- graph(edge(0,0), edge(1,1), edge(0,1,"up",10), edge(1,0,"down",0), StartEnd(0,0))

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
      myNewGraph[1, ] <- edge(0, 0, "null", 0, 0)
      myNewGraph[2, ] <- edge(0, 0, "std", penalty, 0)
    }

    if(type == "isotonic")
    {
      myNewGraph[1, ] <- edge(0, 0, "null", 0, 0)
      myNewGraph[2, ] <- edge(0, 0, "up", penalty, 0)
    }

    if(type == "updown")
    {
      myNewGraph[1, ] <- edge(0, 0, "null", 0, 0)
      myNewGraph[2, ] <- edge(1, 1, "null", 0, 0)
      myNewGraph[3, ] <- edge(0, 1, "up", penalty, 0)
      myNewGraph[4, ] <- edge(1, 0, "down", penalty, 0)
    }

    if(type == "infsup")
    {
      myNewGraph[1, ] <- edge(0, 0, "null", 0, 0)
      myNewGraph[2, ] <- edge(1, 1, "null", 0, 0)
      myNewGraph[3, ] <- edge(0, 1, "absSup", penalty, 0)
      myNewGraph[4, ] <- edge(1, 0, "absInf", penalty, 0)
    }
  }

  class(myNewGraph) <- c("data.frame", "graph")

  return(myNewGraph)
}
