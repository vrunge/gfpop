#' Edge generation
#'
#' @description Edge creation for graph
#' @param state1 a string defining the starting state of the edge
#' @param state2 a string defining the ending state of the edge
#' @param type a string equal to "null", "std", "up", "down", "absInf" or "absSup"
#' @param penalty a nonnegative number. The penality associated to this state transition
#' @param decay a nonnegative number to give the strength of the exponential decay into the segment
#' @param gap a nonnegative number to constraint the size of the gap in the change of state
#' @param oneValue the value of the oneValue edge
#' @param K a positive number. Threshold for the Biweight robust loss
#' @param a a positive number. Slope for the Huber robust loss
#' @return a one-row dataframe with 9 variables
#' @examples
#' Edge("D", "U", "up", penalty = 10, gap = 1)
Edge <- function(state1, state2, type = "null", penalty = 0, decay = 1, gap = 0, K = Inf, a = Inf)
{
  ###STOP###
  if(type != "null" && type != "std" && type != "up" && type != "down" && type != "abs")
  {stop('Argument not appropriate. Choose a type among the following: "null", "std", "up", "down", "abs".')}

  if(!is.double(penalty)){stop('penalty is not a double.')}
  if(!is.double(decay)){stop('decay is not a double.')}
  if(!is.double(gap)){stop('gap is not a double.')}
  if(!is.double(K)){stop('K is not a double.')}
  if(!is.double(a)){stop('a is not a double.')}

  if(penalty < 0){stop('penalty must be nonnegative')}
  if(decay < 0){stop('decay must be nonnegative')}
  if(type == "null"){if(decay == 0){stop('decay must be non-zero')}}
  if(gap < 0){stop('gap must be nonnegative')}
  if(K <= 0){stop('K must be positive')}
  if(a < 0){stop('a must be nonnegative')}

  if(type == "null"){parameter <- decay}else{parameter <- gap}

  ###response = a dataframe with a unique row
  df <- data.frame(character(), character(), character(), numeric(0), numeric(0), numeric(0), numeric(0), numeric(0), numeric(0), stringsAsFactors = FALSE)
  colnames(df) <- c("state1", "state2", "type", "penalty", "parameter", "K", "a", "min", "max")
  df[1,] <- data.frame(state1, state2, type, penalty, parameter, K, a, NA, NA, stringsAsFactors = FALSE)
  return(df)
}


###############################################

#' Start and End nodes for the graph
#'
#' @description Adding constraints on the beginning and ending states of a graph
#' @param start a vector of states. The beginning vertices in the changepoint inference
#' @param end a vector of states. The ending vertices in the changepoint inference
#' @return dataframe with 9 variables with only `state1` and `type` = start or end defined.
#' @examples
#' StartEnd(start = "A", end = c("A","B"))

StartEnd <- function(start = NULL, end = NULL)
{
  df <- data.frame(character(), character(), character(), numeric(0), numeric(0), numeric(0), numeric(0), numeric(0), numeric(0), stringsAsFactors = FALSE)
  colnames(df) <- c("state1", "state2", "type", "penalty", "parameter", "K", "a", "min", "max")
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
#' @description Constraint the rabge of values to consider for a node
#' @param state a string defining the state to constrain
#' @param min minimal value for the infered parameter
#' @param max maximal value for the infered parameter
#' @return a dataframe with 9 variables with only `state1`, `min` and `max` defined.
#' @examples
#' Node(state = "S", min = 0, max = 2)

Node <- function(state = NULL, min = -Inf, max = Inf)
{
  if(!is.double(min)){stop('max is not a double.')}
  if(!is.double(max)){stop('min is not a double.')}
  if(min > max){stop('min is greater than max')}
  df <- data.frame(character(), character(), character(), numeric(0), numeric(0), numeric(0), numeric(0), numeric(0), numeric(0), stringsAsFactors = FALSE)
  colnames(df) <- c("state1", "state2", "type", "penalty", "parameter", "K", "a", "min", "max")
  df [1,] <- data.frame(state, NA, "node", NA, NA, NA, NA, min, max, stringsAsFactors = FALSE)
  return(df)
}




###############################################

#' Graph generation
#'
#' @description Graph creation using component function Edge, StartEnd, Node
#' @param ... This is a list of edges definied by functions edge and StartEnd
#' @param penalty a nonnegative number equals to the common penalty to use for all edges
#' @param type a string equal to "std", "isotonic", "updown", "infsup". to build a predefined classic graph
#' @return a dataframe with edges in rows (columns are named "state1", "state2", "type", "penalty", "parameter") with additional "graph" class.
#' @examples
#' UpDownGraph <- graph(penalty = 10, type = "updown")
#' MyGraph <- graph(Edge(0,0), Edge(1,1), Edge(0,1,"up",10,gap=0.5), Edge(1,0,"down"), StartEnd(0,0), Node(0,0,1), Node(1,0,1))

graph <- function(..., penalty = 0, type = "empty")
{
  myNewGraph <- rbind(...)

  if(is.null(myNewGraph) == TRUE)
  {
    ###STOP###
    if(!is.double(penalty)){stop('penalty is not a double.')}
    if(penalty < 0){stop('penalty must be nonnegative')}

    if(type != "empty" && type != "std" && type != "isotonic" && type != "updown")
      {stop('Arugment "type" not appropriate. Choose among "std", "isotonic", "updown"')}

    myNewGraph <- data.frame(character(), character(), character(), numeric(0), numeric(0), numeric(0), numeric(0), numeric(0), numeric(0), stringsAsFactors = FALSE)
    names(myNewGraph) <- c("state1", "state2", "type", "penalty", "parameter", "K", "a", "min", "max")

    if(type == "std")
    {
      myNewGraph[1, ] <- Edge("S", "S", "null")
      myNewGraph[2, ] <- Edge("S", "S", "std", penalty)
    }

    if(type == "isotonic")
    {
      myNewGraph[1, ] <- Edge("S", "S", "null")
      myNewGraph[2, ] <- Edge("S", "S", "up", penalty)
    }

    if(type == "updown")
    {
      myNewGraph[1, ] <- Edge("D", "D", "null")
      myNewGraph[2, ] <- Edge("U", "U", "null")
      myNewGraph[3, ] <- Edge("D", "U", "up", penalty)
      myNewGraph[4, ] <- Edge("U", "D", "down", penalty)
    }
  }
  class(myNewGraph) <- c("data.frame", "graph")
  return(myNewGraph)
}




###############################################
# invisible function for the user

graphAnalysis <- function(mygraph)
{
  ### BUILD an ordered Graph : myOrderedGraph ###
  ##separate startend from vertices
  graphNA <- mygraph[is.na(mygraph[,2]),] ## Start End vertices
  graphV <-  mygraph[!is.na(mygraph[,2]),] ## Edges of the graph

  myVertices <- unique(c(graphV[,1], graphV[,2]))

  ##create a new graph
  myNewGraph <- graph()

  selectNull <- graphV[, 3] == "null" ### => penalty = 0
  graphV[selectNull, 4] <- -1 #set penalty to -1


  for(vertex in myVertices)
  {
    selectRaw <- graphV[graphV[,2]==vertex, ]
    ordre <- order(selectRaw[,4])
    selectRaw <- selectRaw[ordre,]
    myNewGraph <- rbind(myNewGraph, selectRaw)
  }

  myNewGraph <- rbind(myNewGraph, graphNA)
  selectNull <- myNewGraph[, 3] == "null"
  myNewGraph[selectNull, 4] <- 0 #for ordering


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
# to use after graphAnalysis

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
