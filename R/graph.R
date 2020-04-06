##  GPL-3 License
## Copyright (c) 2019 Vincent Runge

#' Edge generation
#'
#' @description Edge creation for graph
#' @param state1 a string defining the starting state of the edge
#' @param state2 a string defining the ending state of the edge
#' @param type a string equal to "null", "std", "up", "down" or "abs"
#' @param decay a nonnegative number to give the strength of the exponential decay into the segment
#' @param gap a nonnegative number to constrain the size of the gap in the change of state
#' @param penalty a nonnegative number. The penality associated to this state transition
#' @param K a positive number. Threshold for the Biweight robust loss
#' @param a a positive number. Slope for the Huber robust loss
#' @return a one-row dataframe with 9 variables
#' @examples
#' Edge("Dw", "Up", "up", gap = 1, penalty = 10, K = 3)
Edge <- function(state1, state2, type = "null", decay = 1, gap = 0, penalty = 0, K = Inf, a = 0)
{
  allowed.types <- c("null", "std", "up", "down", "abs")
  if(!type %in% allowed.types){stop('type must be one of: ', paste(allowed.types, collapse=", "))}
  if(!is.double(decay)){stop('decay is not a double.')}
  if(!is.double(gap)){stop('gap is not a double.')}
  if(!is.double(penalty)){stop('penalty is not a double.')}
  if(!is.double(K)){stop('K is not a double.')}
  if(!is.double(a)){stop('a is not a double.')}

  if(any(type == "null" && decay == 0, na.rm=TRUE))stop('decay must be non-zero')
  if(any(decay < 0, na.rm=TRUE)){stop('decay must be nonnegative')}
  if(any(gap < 0, na.rm=TRUE)){stop('gap must be nonnegative')}
  if(any(penalty < 0, na.rm=TRUE)){stop('penalty must be nonnegative')}
  if(any(K <= 0, na.rm=TRUE)){stop('K must be positive')}
  if(any(a < 0, na.rm=TRUE)){stop('a must be nonnegative')}

  #fill parameter variable

  if(type == "null"){parameter <- decay}else{parameter <- gap}

  data.frame(state1, state2, type, parameter, penalty, K, a, min=NA, max=NA, stringsAsFactors = FALSE)

}


###############################################

#' Start and End nodes for the graph
#'
#' @description Defining the beginning and ending states of a graph
#' @param start a vector of states. The beginning nodes for the changepoint inference
#' @param end a vector of states. The ending nodes for the changepoint inference
#' @return dataframe with 9 variables with only `state1` and `type` = start or end defined.
#' @examples
#' StartEnd(start = "A", end = c("A","B"))

StartEnd <- function(start = NULL, end = NULL)
{
  ### delete repetitions if any
  start <- unique(start)
  end <- unique(end)

  df <- data.frame(character(), character(), character(), numeric(0), numeric(0), numeric(0), numeric(0), numeric(0), numeric(0), stringsAsFactors = FALSE)
  colnames(df) <- c("state1", "state2", "type", "parameter", "penalty", "K", "a", "min", "max")
  if(length(start) != 0)
  {
    for(i in 1:length(start))
      {df[i,] <- list(start[i], NA, "start", NA, NA, NA, NA, NA, NA)}
  }
  if(length(end) != 0)
  {
    for(i in 1:length(end))
      {df[i + length(start),] <- list(end[i], NA, "end", NA, NA, NA, NA, NA, NA)}
  }
  return(df)
}


###############################################

#' Node Values
#'
#' @description Constrain the range of values to consider at a node
#' @param state a string defining the state to constrain
#' @param min minimal value for the inferred parameter
#' @param max maximal value for the inferred parameter
#' @return a dataframe with 9 variables with only `state1`, `min` and `max` defined.
#' @examples
#' Node(state = "s0", min = 0, max = 2)

Node <- function(state = NULL, min = -Inf, max = Inf)
{
  if(!is.double(min)){stop('min is not a double.')}
  if(!is.double(max)){stop('max is not a double.')}
  if(min > max){stop('min is greater than max')}

  df <- data.frame(character(), character(), character(), numeric(0), numeric(0), numeric(0), numeric(0), numeric(0), numeric(0), stringsAsFactors = FALSE)
  colnames(df) <- c("state1", "state2", "type", "parameter", "penalty", "K", "a", "min", "max")
  df [1,] <- data.frame(state, state, "node", NA, NA, NA, NA, min, max, stringsAsFactors = FALSE)
  return(df)
}


###############################################

#' Graph generation
#'
#' @description Graph creation using component functions "Edge", "StartEnd" and "Node"
#' @param ... This is a list of edges definied by functions "Edge", "StartEnd" and "Node"
#' @param type a string equal to "std", "isotonic", "updown", "relevant". to build a predefined classic graph
#' @param decay a nonnegative number to give the strength of the exponential decay into the segment
#' @param gap a nonnegative number to constrain the size of the gap in the change of state
#' @param penalty a nonnegative number equals to the common penalty to use for all edges
#' @param K a positive number. Threshold for the Biweight robust loss
#' @param a a positive number. Slope for the Huber robust loss
#' @param all.null.edges a boolean. Add null edges to all nodes automatically
#' @return a dataframe with 9 variables (columns are named "state1", "state2", "type", "parameter", "penalty", "K", "a", "min", "max") with additional "graph" class.
#' @examples
#' UpDownGraph <- graph(type = "updown", gap = 1.3, penalty = 10)
#' MyGraph <- graph(Edge("Dw","Dw"), Edge("Up","Up"), Edge("Dw","Up","up", gap = 0.5, penalty = 10),
#' Edge("Up","Dw","down"), StartEnd("Dw","Dw"), Node("Dw",0,1), Node("Up",0,1))

graph <- function(..., type = "empty", decay = 1, gap = 0, penalty = 0, K = Inf, a = 0, all.null.edges = FALSE)
{
  #### build the graph with the collection ... of edges
  myNewGraph <- rbind(...)
  if(is.null(myNewGraph) == FALSE)
  {
    not.special <- subset(myNewGraph, ! type %in% c("node", "null", "start", "end"))
    if(isTRUE(all.null.edges) && nrow(not.special > 0))
    {
      myNewGraph <- subset(myNewGraph, type != "null")
      u.states <- with(not.special, unique(c(state1, state2)))
      for(i in u.states){myNewGraph <- rbind(Edge(i, i, "null", decay = decay), myNewGraph)}
    }
  }

  #### user specified graph
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

    myNewGraph <- data.frame(character(), character(), character(), numeric(0), numeric(0), numeric(0), numeric(0), numeric(0), numeric(0), stringsAsFactors = FALSE)
    names(myNewGraph) <- c("state1", "state2", "type", "parameter", "penalty", "K", "a", "min", "max")

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
  class(myNewGraph) <- c("data.frame", "graph")
  return(myNewGraph)
}


###############################################
###############################################
###############################################
# invisible function for the user
#Order the graph and create integer state values


graphReorder <- function(mygraph)
{
  ### BUILD an ordered Graph : myOrderedGraph ###
  ##separate start, end, node from vertices
  graphNA <- mygraph[is.na(mygraph[,5]),] ## Start End nodes and range values nodes
  graphVtemp <-  mygraph[!is.na(mygraph[,5]),] ## Edges of the graph
  myVertices <- unique(c(graphVtemp[,1], graphVtemp[,2]))

  if(!all(is.element(mygraph[is.na(mygraph[,5]), 1], myVertices))){stop("Some start-end-node names not related to edges")}

  ###transform the abs edge into two edges (up and down)
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

  ##create a new graph
  myNewGraph <- graph()
  selectNull <- graphV[, 3] == "null" ### => penalty = 0
  graphV[selectNull, 5] <- -1 #set penalty to -1

  for(vertex in myVertices)
  {
    selectRaw <- graphV[graphV[,2]==vertex, ]
    ordre <- order(selectRaw[,5])
    selectRaw <- selectRaw[ordre,]
    myNewGraph <- rbind(myNewGraph, selectRaw)
  }

  myNewGraph <- rbind(myNewGraph, graphNA)
  selectNull <- myNewGraph[, 3] == "null"
  myNewGraph[selectNull, 5] <- 0 #for ordering

  ###Label the vertices with integers from 0 to nbVertices
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



###############################################
# invisible function for the user
# to use after graphAnalysis and test whether the graph can be used in the gfpop function

explore <- function(mygraph)
{
  graph <- mygraph$graph
  graph[,1] <- graph[,1] + 1
  graph[,2] <- graph[,2] + 1 ### states from 1 to nbState
  len <- length(mygraph$vertices)

  theStart <- graph[graph[,3] == "start", 1]
  theEnd <- graph[graph[,3] == "end", 1]
  if(length(theStart) == 0){theStart <- 1:len}
  if(length(theEnd) == 0){theEnd <- 1:len}

  recNodes <- graph[which(graph$state1 == graph$state2),]$state1 ### recursive nodes
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



###############################################

visit <- function(graph, startNode)
{
  visited <- NULL
  toVisit <- c(startNode)
  while(length(toVisit) > 0)
  {
    visited <- c(visited, toVisit[1]) #we visit toVisit[1]
    newToVisit <- graph[graph[,1] == toVisit[1], 2] #to visit from toVisit[1]
    newToVisit <- newToVisit[!is.na(newToVisit)] #remove NA
    newToVisit <- setdiff(newToVisit, visited)
    newToVisit <- setdiff(newToVisit, toVisit)
    toVisit <- toVisit[-1] #we remove toVisit[1] in node to visit
    toVisit <- c(newToVisit, toVisit)
  }
  return(visited)
}
