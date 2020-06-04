##  GPL-3 License
## Copyright (c) 2019 Vincent Runge

#' Graph-Constrained Functional Pruning Optimal Partitioning (gfpop)
#'
#' @description Functional pruning optimal partitioning with a graph structure to take into account constraints on consecutive segment parameters. The user has to specify the graph he wants to use (see the graph function) and a type of cost function. This is the main function of the gfpop package.
#' @param data vector of data to segment
#' @param mygraph dataframe of class "graph" to constrain the changepoint inference
#' @param type a string defining the cost model to use: "mean", "variance", "poisson", "exp", "negbin"
#' @param weights vector of weights (positive numbers), same size as data
#' @param testMode boolean. False by default. Used to debug the code
#' @return a gfpop object = (changepoints, states, forced, parameters, globalCost)
#' \describe{
#' \item{\code{changepoints}}{is the vector of changepoints (we give the last element of each segment)}
#' \item{\code{states}}{is the vector giving the state of each segment}
#' \item{\code{forced}}{is the vector specifying whether the constraints of the graph are active (=1) or not (=0)}
#' \item{\code{parameters}}{is the vector of successive parameters of each segment}
#' \item{\code{globalCost}}{is a number equal to the total loss: the minimal cost for the optimization problem with all penalty values excluded}
#'  }
gfpop <- function(data, mygraph, type = "mean", weights = NULL, testMode = FALSE)
{
  #enforce factor to string if necessary
  mygraph$state1 <- as.character(mygraph$state1)
  mygraph$state2 <- as.character(mygraph$state2)
  mygraph$type <- as.character(mygraph$type)

  ############
  ### STOP ###
  ############

  if(!any(class(mygraph) == "graph")){stop('Your graph is not a graph created with the graph function of the gfpop package.')}

  allowed.types <- c("mean", "variance", "poisson", "exp", "negbin")
  if(!type %in% allowed.types){stop('type must be one of: ', paste(allowed.types, collapse=", "))}

  ### if we have weights
  if(!is.null(weights))
  {
    if(length(data) != length(weights)){stop('data vector and weights vector have different sizes')}
    if(!all(weights > 0)){stop('weights vector has non strictly positive components')}
  }
  else{weights <- 0} #to send a double in gfpopTransfer
  if(length(data) < 2){stop('data vector length is less than 2...')}
  if(any(is.na(data)))stop("data has missing values, please remove them")

  ######################
  ### GRAPH ANALYSIS ###
  ######################
  mynewgraph <- graphReorder(mygraph) ### reorder the edges
  explore(mynewgraph) ### test if the graph can be used

  newGraph <- mynewgraph$graph
  vertices <- mynewgraph$vertices

  ###########################
  ### CALL Rcpp functions ###
  ###########################

  res <- gfpopTransfer(data, newGraph, type, weights, testMode)

  ############################
  ### Response class gfpop ###
  ############################

  if(length(res$changepoints) == 1) ##### best output state
  {
    response <- list(changepoints = c(rev(res$changepoints[[1]][-1]), length(data)), states = vertices[rev(res$states[[1]])+1], forced = rev(res$forced[[1]]), parameters = rev(res$param[[1]]), globalCost = res$cost[[1]])
  }
  else  ##### multiple output state
  {
    p <- length(res$changepoints)
    lastStates <- NULL
    for(i in 1:p)
    {
      res$changepoints[[i]] <- c(rev(res$changepoints[[i]][-1]), length(data))
      res$states[[i]] <- vertices[rev(res$states[[i]])+1]
      lastStates <- c(lastStates, rev(res$states[[i]])[1])
      res$forced[[i]] <- rev(res$forced[[i]])
      res$param[[i]] <- rev(res$param[[i]])
    }
    names(res$changepoints) <- lastStates
    names(res$states) <- lastStates
    names(res$forced) <- lastStates
    names(res$param) <- lastStates
    names(res$cost) <- lastStates
    response <- list(changepoints = res$changepoints, states = res$states, forced = res$forced, parameters = res$param, globalCost = res$cost)
  }

  attr(response, "class") <- "gfpop"
  return(response)
}


########################################################################################
# mygraph has penalties of type = sigma^2 or const * sigma^2

#' Graph-constrained functional pruning optimal partitioning iterated
#'
#' @description Functional pruning optimal partitioning with a graph structure to take into account constraints on consecutive segment parameters. This is an iterated version of the main gfpop function using a Birgé Massart like penalty
#' @param data vector of data to segment
#' @param mygraph dataframe of class graph to constrain the changepoint inference
#' @param type a string defining the cost model to use: "gauss", "variance", "poisson", "exp", "negbin"
#' @param weights vector of weights (positive numbers), same size as data
#' @param iter.max maximal number of iteration of the gfpop function
#' @param D.init initialisation of the number of segments
#' @return a gfpop object = (changepoints, states, forced, parameters, globalCost, Dvect)
#' \describe{
#' \item{\code{changepoints}}{is the vector of changepoints (we give the last element of each segment)}
#' \item{\code{states}}{is the vector giving the state of each segment}
#' \item{\code{forced}}{is the vector specifying whether the constraints of the graph are active (=1) or not (=0)}
#' \item{\code{parameters}}{is the vector of successive parameters of each segment}
#' \item{\code{globalCost}}{is a number equal to the total loss: the minimal cost for the optimization problem with all penalty values excluded}
#' \item{\code{Dvect}}{is a vector of integers. The successive tested D in the Birgé Massart penalty until convergence}
#'  }
itergfpop <- function(data, mygraph, type = "mean", weights = NULL, iter.max = 100, D.init = 1)
{
  ############
  ### STOP ###
  ############
  if(!any(class(mygraph) == "graph")){stop('Your graph is not a graph created with the graph function in gfpop package...')}

  allowed.types <- c("mean", "variance", "poisson", "exp", "negbin")
  if(!type %in% allowed.types){stop('type must be one of: ', paste(allowed.types, collapse=", "))}

  ### if we have weights
  if(!is.null(weights))
  {
    if(length(data) != length(weights)){stop('data vector and weights vector have different sizes')}
    if(!all(weights>0)){stop('weights vector has non strictly positive components')}
  }
  if(length(data) < 2){stop('data vector length is less than 2...')}

  ######################
  ### GRAPH ANALYSIS ###
  ######################
  mynewgraph <- graphReorder(mygraph) ### reorder the edges
  explore(mynewgraph) ### test if the graph can be used

  newGraph <- mynewgraph$graph
  vertices <- mynewgraph$vertices

  ######################
  ### ITERATED GFPOP ###
  ######################
  iter <- 0
  beta_old <- 0
  n <- length(data)
  beta <- getDerivativePenalty(D.init, n)

  #update newGraph
  newGraph[,4] <- newGraph[,4] * beta
  Dvect <- D.init

  while(beta_old != beta & iter <= iter.max)
  {
    ###########################
    ### CALL Rcpp functions ###
    ###########################
    res <- gfpopTransfer(data, newGraph, type, weights)
    beta_old <- beta

    newGraph[,4] <- newGraph[,4] / beta
    D <- length(res$changepoints)
    Dvect <- c(Dvect, D)
    beta <- getDerivativePenalty(D, n)
    newGraph[,4] <- newGraph[,4] * beta

    iter <- iter + 1
  }

  ############################
  ### Response class gfpop ###
  ############################
  response <- list(changepoints = c(rev(res$changepoints[-1]), length(data)), states = vertices[rev(res$states)+1], forced = rev(res$forced), parameters = rev(res$param), globalCost = res$cost, Dvect = rev(rev(Dvect)[-1]))
  attr(response, "class") <- "gfpop"

  return(response)
}
