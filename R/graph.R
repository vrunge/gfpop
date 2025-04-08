##  GPL-3 License
## Copyright (c) 2022 Vincent Runge

#' Edge generation
#'
#' @description Edge creation for gfpop-R-Package graph
#' @param state1 a string defining the starting state of the edge
#' @param state2 a string defining the ending state of the edge
#' @param type a string equal to \code{"null"}, \code{"std"}, \code{"up"}, \code{"down"} or \code{"abs"}. Default type is \code{"null"},
#' the transition to stay on the same segment.
#' @param decay a nonnegative number to give the strength of the exponential decay into the segment
#' @param gap a nonnegative number to constrain the size of the gap in the change of state
#' @param penalty a nonnegative number. The penality associated to this state transition
#' @param K a positive number. Threshold for the Biweight robust loss
#' @param a a positive number. Slope for the Huber robust loss
#' @param rule a positive integer identifying the rule ID for time-dependent constraints
#' @return a one-row dataframe with 10 variables
#' @examples
#' Edge("Dw", "Up", "up", gap = 1, penalty = 10, K = 3)
#'
#' Edge(0, 1, "abs", penalty = 2, gap = 1)
#'
#' Edge(0, 0, "null", penalty = 0, K = 2, a = 1)
#'
#' Edge("Dw", "Dw", type = "null", decay = 0.997)
#' 
#' Edge("Dw", "Up", "up", rule = 2) # Use with time-dependent constraints
Edge <- function(state1, state2, type = "null", decay = 1, gap = 0, penalty = 0, K = Inf, a = 0, rule = 1)
{
  allowed.types <- c("null", "std", "up", "down", "abs")
  if(!type %in% allowed.types){stop('type must be one of: ', paste(allowed.types, collapse=", "))}
  if(!is.double(decay)){stop('decay is not a double.')}
  if(!is.double(gap)){stop('gap is not a double.')}
  if(!is.double(penalty)){stop('penalty is not a double.')}
  if(!is.double(K)){stop('K is not a double.')}
  if(!is.double(a)){stop('a is not a double.')}
  
  # Add validation for rule parameter
  if(any(!is.numeric(rule))){stop('rule must be numeric')}
  if(any(rule <= 0)){stop('rule must be positive')}
  if(any(rule != floor(rule))){stop('rule must be a positive integer')}

  if(any(type == "null" && decay == 0, na.rm=TRUE))stop('decay must be non-zero')
  if(any(decay < 0, na.rm=TRUE)){stop('decay must be nonnegative')}
  if(any(gap < 0, na.rm=TRUE)){stop('gap must be nonnegative')}
  if(any(penalty < 0, na.rm=TRUE)){stop('penalty must be nonnegative')}
  if(any(K <= 0, na.rm=TRUE)){stop('K must be positive')}
  if(any(a < 0, na.rm=TRUE)){stop('a must be nonnegative')}

  #fill parameter variable
  if(type == "null"){parameter <- decay}else{parameter <- gap}
  
  data.frame(state1, state2, type, parameter, penalty, K, a, min=NA, max=NA, rule=rule, stringsAsFactors = FALSE)
}

#' Start and End nodes for the graph
#'
#' @description Defining the beginning and ending states of a graph
#' @param start a vector of states. The beginning nodes for the changepoint inference
#' @param end a vector of states. The ending nodes for the changepoint inference
#' @return dataframe with 10 variables with only \code{state1} and \code{type = "start"} or \code{"end"} defined (not \code{NA}).
#' @examples
#' StartEnd(start = "A", end = c("A","B"))
#'
#' StartEnd(start = 0)
#'
#' StartEnd(start = 1, end = 1)
#'
#' StartEnd(start = "v0", end = "v3")
#'
#' StartEnd(end = "s0")
StartEnd <- function(start = NULL, end = NULL)
{
  start <- unique(start)
  end <- unique(end)

  df <- data.frame(character(), character(), character(), numeric(0), numeric(0), numeric(0), numeric(0), numeric(0), numeric(0), numeric(0), stringsAsFactors = FALSE)
  colnames(df) <- c("state1", "state2", "type", "parameter", "penalty", "K", "a", "min", "max", "rule")
  if(length(start) != 0)
  {
    for(i in 1:length(start))
      {df[i,] <- list(start[i], NA, "start", NA, NA, NA, NA, NA, NA, NA)}
  }
  if(length(end) != 0)
  {
    for(i in 1:length(end))
      {df[i + length(start),] <- list(end[i], NA, "end", NA, NA, NA, NA, NA, NA, NA)}
  }
  return(df)
}

#' Node Values
#'
#' @description Constrain the range of values to consider at a node
#' @param state a string defining the state to constrain
#' @param min minimal value for the inferred parameter
#' @param max maximal value for the inferred parameter
#' @return a dataframe with 10 variables with only \code{state1}, \code{min} and \code{max} defined (not \code{NA}).
#' @examples
#' Node(state = "s0", min = 0, max = 2)
#'
#' Node(state = 0, min = -1, max = 1)
#'
#' Node(state = "positive", min = 0)
#'
#' Node(state = "mu0", min = 0.5, max = 0.5)
Node <- function(state = NULL, min = -Inf, max = Inf)
{
  if(!is.double(min)){stop('min is not a double.')}
  if(!is.double(max)){stop('max is not a double.')}
  if(min > max){stop('min is greater than max')}

  df <- data.frame(character(), character(), character(), numeric(0), numeric(0), numeric(0), numeric(0), numeric(0), numeric(0), numeric(0), stringsAsFactors = FALSE)
  colnames(df) <- c("state1", "state2", "type", "parameter", "penalty", "K", "a", "min", "max", "rule")
  df [1,] <- data.frame(state, state, "node", NA, NA, NA, NA, min, max, NA, stringsAsFactors = FALSE)
  return(df)
}

#' Graph generation
#'
#' @description Graph creation using component functions \code{Edge}, \code{StartEnd} and \code{Node}
#' @param ... This is a list of edges defined by functions \code{Edge}, \code{StartEnd} and \code{Node}
#' @param type a string equal to \code{"std"}, \code{"isotonic"}, \code{"updown"} or \code{"relevant"} 
#' @param decay a nonnegative number to give the strength of the exponential decay
#' @param gap a nonnegative number to constrain the size of the gap in the change of state
#' @param penalty a nonnegative number equals to the common penalty for all edges
#' @param K a positive number. Threshold for the Biweight robust loss
#' @param a a positive number. Slope for the Huber robust loss
#' @param all.null.edges a boolean. Add null edges to all nodes automatically
#' @return a dataframe with 10 variables with additional \code{"graph"} class.
#' @examples
#' graph(type = "updown", gap = 1.3, penalty = 5)
#'
#' graph(Edge("Dw","Dw"),
#'       Edge("Up","Up"),
#'       Edge("Dw","Up","up", gap = 0.5, penalty = 10),
#'       Edge("Up","Dw","down"),
#'       StartEnd("Dw","Dw"),
#'       Node("Dw",0,1),
#'       Node("Up",0,1))
graph <- function(..., type = "empty", decay = 1, gap = 0, penalty = 0, K = Inf, a = 0, all.null.edges = FALSE)
{
  myNewGraph <- rbind(...)
  
  if(!is.null(myNewGraph) && all.null.edges) {
    # Extract all unique state names from both columns and remove NAs
    all_states <- unique(c(as.character(myNewGraph$state1), as.character(myNewGraph$state2)))
    all_states <- all_states[!is.na(all_states)]
    
    # Save original edges unchanged (preserve their rule values)
    original_edges <- myNewGraph
    
    # Create a data frame of auto-generated null edges.
    # Use reverse sorted order to match expected ordering.
    null_edges <- data.frame()
    for(state in rev(sort(all_states))) {
      null_edge <- Edge(state, state, "null", decay = decay, rule = 1)
      null_edges <- rbind(null_edges, null_edge)
    }
    
    # Combine auto null edges (first) with original edges (later)
    myNewGraph <- rbind(null_edges, original_edges)
  }
  
  if(is.null(myNewGraph) == TRUE)
  {
    allowed.graphs <- c("empty", "std", "isotonic", "updown", "relevant")
    if(!type %in% allowed.graphs){stop('type must be one of: ', paste(allowed.graphs, collapse=", "))}

    if(!is.double(decay)){stop('decay is not a double.')}
    if(!is.double(gap)){stop('gap is not a double.')}
    if(decay < 0){stop('decay must be nonnegative')}
    if(gap < 0){stop('gap must be nonnegative')}
    if(!is.double(penalty)){stop('penalty is not a double.')}
    if(penalty < 0){stop('penalty must be nonnegative')}

    myNewGraph <- data.frame(character(), character(), character(), numeric(0), numeric(0), numeric(0), numeric(0), numeric(0), numeric(0), numeric(0), stringsAsFactors = FALSE)
    names(myNewGraph) <- c("state1", "state2", "type", "parameter", "penalty", "K", "a", "min", "max", "rule")

    if(type == "std")
    {
      myNewGraph[1, ] <- Edge("Std", "Std", "null", decay = decay, K = K, a = a)
      myNewGraph[2, ] <- Edge("Std", "Std", "std", penalty = penalty, K = K, a = a)
    }
    else if(type == "isotonic")
    {
      myNewGraph[1, ] <- Edge("Iso", "Iso", "null", decay = decay, K = K, a = a)
      myNewGraph[2, ] <- Edge("Iso", "Iso", "up", gap = gap, penalty = penalty, K = K, a = a)
    }
    else if(type == "updown")
    {
      myNewGraph[1, ] <- Edge("Dw", "Dw", "null", decay = decay, K = K, a = a)
      myNewGraph[2, ] <- Edge("Up", "Up", "null", decay = decay, K = K, a = a)
      myNewGraph[3, ] <- Edge("Dw", "Up", "up", gap = gap, penalty = penalty, K = K, a = a)
      myNewGraph[4, ] <- Edge("Up", "Dw", "down", gap = gap, penalty = penalty, K = K, a = a)
    }
    else if(type == "relevant")
    {
      myNewGraph[1, ] <- Edge("Abs", "Abs", "null", decay = decay, K = K, a = a)
      myNewGraph[2, ] <- Edge("Abs", "Abs", "abs", gap = gap, penalty = penalty, K = K, a = a)
    }
  }
  class(myNewGraph) <- c("graph", "data.frame")
  return(myNewGraph)
}

graphReorder <- function(mygraph)
{
  graphNA <- mygraph[is.na(mygraph[,5]),]
  graphVtemp <-  mygraph[!is.na(mygraph[,5]),]
  myVertices <- unique(c(graphVtemp[,1], graphVtemp[,2]))

  if(!all(is.element(mygraph[is.na(mygraph[,5]), 1], myVertices))){stop("Some start-end-node names not related to edges")}

  absEdge <- graphVtemp[,3] == "abs"

  if(!all(absEdge == FALSE))
  {
    graphVtemp[absEdge,3] <- "down"
    addToGraphVV <- graphVtemp[absEdge,]
    addToGraphVV[,3] <- "up"
    graphV <- rbind(graphVtemp, addToGraphVV)
  }else
  {
    graphV <- graphVtemp
  }

  myNewGraph <- graph()
  selectNull <- graphV[, 3] == "null"
  graphV[selectNull, 5] <- -1

  for(vertex in myVertices)
  {
    selectRaw <- graphV[graphV[,2]==vertex, ]
    ordre <- order(selectRaw[,5])
    selectRaw <- selectRaw[ordre,]
    myNewGraph <- rbind(myNewGraph, selectRaw)
  }

  myNewGraph <- rbind(myNewGraph, graphNA)
  selectNull <- myNewGraph[, 3] == "null"
  myNewGraph[selectNull, 5] <- 0

  for(i in 1:dim(myNewGraph)[1])
  {
    myNewGraph[i,1] <- which(myNewGraph[i,1] == myVertices) - 1
    if(!is.na(myNewGraph[i,2])){myNewGraph[i,2] <- which(myNewGraph[i,2] == myVertices) - 1}
  }

  class(myNewGraph$state1) <- "numeric"
  class(myNewGraph$state2) <- "numeric"

  response <- list(graph = myNewGraph, vertices = myVertices)
  return(response)
}

explore <- function(mygraph)
{
  graph <- mygraph$graph
  graph[,1] <- graph[,1] + 1
  graph[,2] <- graph[,2] + 1
  len <- length(mygraph$vertices)

  theStart <- graph[graph[,3] == "start", 1]
  theEnd <- graph[graph[,3] == "end", 1]
  if(length(theStart) == 0){theStart <- 1:len}
  if(length(theEnd) == 0){theEnd <- 1:len}

  recNodes <- graph[which(graph$state1 == graph$state2),]$state1
  seenNodes <- NULL

  for(i in 1:len)
  {
    Vi <- visit(graph, i)
    if(length(intersect(Vi,theEnd)) == 0){stop('Not all nodes lead to an end node')}

    if(i %in% theStart && length(intersect(Vi,recNodes)) == 0){stop('Not all path have a recursive edge')}
    if(i %in% theStart){seenNodes <- c(seenNodes, Vi)}
  }
  seenNodes <- sort(unique(seenNodes))
  if(length(seenNodes) != len){stop('One or more nodes is/are not seen by the algorithm')}
}

visit <- function(graph, startNode)
{
  visited <- NULL
  toVisit <- c(startNode)
  while(length(toVisit) > 0)
  {
    visited <- c(visited, toVisit[1])
    newToVisit <- graph[graph[,1] == toVisit[1], 2]
    newToVisit <- newToVisit[!is.na(newToVisit)]
    newToVisit <- setdiff(newToVisit, visited)
    newToVisit <- setdiff(newToVisit, toVisit)
    toVisit <- toVisit[-1]
    toVisit <- c(newToVisit, toVisit)
  }
  return(visited)
}
