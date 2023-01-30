##  GPL-3 License
## Copyright (c) 2022 Vincent Runge

#' Graph-Constrained Functional Pruning Optimal Partitioning (gfpop)
#'
#' @description Functional pruning optimal partitioning with a graph structure to take
#' into account constraints on consecutive segment parameters. The user has to specify
#' the graph he wants to use (see the graph function) and a type of cost function.
#' This is the main function of the gfpop package. Its result can be plotted using
#' the S3 gfpop function \code{\link[gfpop:plot.gfpop]{gfpop::plot()}}
#' @param data vector of data to segment. For simulation studies, Data can be generated using gfpop package function \code{\link[gfpop:dataGenerator]{gfpop::dataGenerator()}}
#' @param mygraph dataframe of class "graph" to constrain the changepoint inference, see \code{\link[gfpop:graph]{gfpop::graph()}}
#' @param type a string defining the cost model to use: \code{"mean"}, \code{"variance"}, \code{"poisson"}, \code{"exp"}, \code{"negbin"}
#' @param weights vector of weights (positive numbers), same size as data
#' @param testMode boolean. \code{FALSE} by default. Used to debug the code
#' @return a gfpop object = (\code{changepoints, states, forced, parameters, globalCost})
#' \describe{
#' \item{\code{changepoints}}{is the vector of changepoints (we give the last element of each segment)}
#' \item{\code{states}}{is the vector giving the state of each segment}
#' \item{\code{forced}}{is the vector specifying whether the constraints of the graph are active (\code{= TRUE}) or not (\code{= FALSE})}
#' \item{\code{parameters}}{is the vector of successive parameters of each segment}
#' \item{\code{globalCost}}{is a number equal to the total loss: the minimal cost for the optimization problem with all penalty values excluded}
#'  }
#' @details
#' The constrained optimization problem for n data points takes the following general form:
#' \deqn{Q_n = min (with constraints) (\sum_{t=1}^n (\gamma(e[t])(y[t], \mu[t]) + \beta(e[t]))}
#' with data points \eqn{y[t]}, edges \eqn{e[t]}, edge-dependent penalties \eqn{\beta(e[t])} and cost functions \eqn{\gamma}.
#' The cost function can take three different forms for parameter x and constants (A, B, C):
#' \itemize{
#'   \item \emph{quadratic}, with representation \eqn{Ax^2 + Bx +C} with  \eqn{x} in R
#'   \item \emph{log-linear}, with representation \eqn{Ax  - B log(x) +C} with  \eqn{x \ge 0}
#'   \item \emph{log-log}, with representation \eqn{- A log(x) - B log(1-x) +C} with  \eqn{0 \le x \le 1}
#' }
#' For each optimization problem, we consider a unique cost representation.
#' However, the User can define robustness values (K and a) specific to each edge, making the cost function edge-dependent.
#' We give the atomic form of each of the five available types (for one data point of value y with weight w)
#' \itemize{
#'   \item \code{"mean"} :  \eqn{A = w}, \eqn{B = -2wy}, \eqn{C = wy^2}
#'   \item \code{"variance"} : \eqn{A = wy^2}, \eqn{B = w}, \eqn{C = 0}
#'   \item \code{"poisson"} : \eqn{A = w}, \eqn{B = wy}, \eqn{C = 0}
#'   \item \code{"exp"} : \eqn{A = wy}, \eqn{B = w}, \eqn{C = 0}
#'   \item \code{"negbin"} : \eqn{A = w}, \eqn{B = wy}, \eqn{C = 0}
#' }
#' @seealso
#' \itemize{
#'   \item \code{\link[gfpop:dataGenerator]{gfpop::dataGenerator()}} to generate data for multiple change-point simulations
#'   \item \code{\link[gfpop:graph]{gfpop::graph()}} to create graphs complying with the gfpop function
#'   \item \code{\link[gfpop:plot.gfpop]{gfpop::plot()}} to plot the gfpop object and visualize inferred changepoints and parameters
#' }
#' @examples
#' n <- 1000 #data length
#' ### EXAMPLE 1 ### updown graph + poisson loss
#'  myData <- dataGenerator(n, c(0.1, 0.3, 0.5, 0.8, 1), c(1, 2, 1, 3, 1), type = "poisson")
#'  myGraph <- graph(penalty = 2 * sdDiff(myData)^2 * log(n), type = "updown")
#'  gfpop(data = myData, mygraph = myGraph, type = "poisson")
#'
#' ### EXAMPLE 2 ### relevant graph with min gap = 2 + poisson loss
#'  myData <- dataGenerator(n, c(0.1, 0.3, 0.5, 0.8, 1), c(1, 2, 3, 5, 3), type = "poisson")
#'  myGraph <- graph(type = "relevant", penalty = 2 * log(n), gap = 2)
#'  gfpop(data =  myData, mygraph = myGraph, type = "poisson")
#'
#' ### EXAMPLE 3 ### std graph with robust loss + variance loss
#'  myData <- dataGenerator(n, c(0.1, 0.3, 0.5, 0.8, 1), c(1, 5, 1, 5, 1), type = "variance")
#'  outliers <- 5 * rbinom(n, 1, 0.05) - 5 * rbinom(n, 1, 0.05)
#' ### with robust parameter K
#'  myGraph <- graph(type = "std", penalty = 2 * log(n), K = 10)
#'  gfpop(data =  myData + outliers, mygraph = myGraph, type = "variance")
#' ### no K
#'  myGraph <- graph(type = "std", penalty = 2 * log(n))
#'  gfpop(data =  myData, mygraph = myGraph, type = "variance")
#'
#' ### EXAMPLE 4 ###  3-segment graph with mean (Gaussian) loss
#'  myData <- dataGenerator(n, c(0.12, 0.31, 0.53, 0.88, 1), c(1, 2, 0, 1, 2), type = "mean")
#'  outliers <- 5 * rbinom(n, 1, 0.05) - 5 * rbinom(n, 1, 0.05)
#'  gfpop(data =  myData + outliers, mygraph = paperGraph(8, penalty = 2 * log(n)), type = "mean")
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
  attr(response, "type") <- type

  return(response)
}


########################################################################################
# mygraph has penalties of type = sigma^2 or const * sigma^2

#' Graph-constrained functional pruning optimal partitioning iterated
#'
#' @description Functional pruning optimal partitioning with a graph structure to take into account constraints on consecutive segment parameters.
#' This is an iterated version of the main gfpop function using a Birgé Massart like penalty
#' @param data vector of data to segment. For simulation studies, Data can be generated using gfpop package function \code{\link[gfpop:dataGenerator]{gfpop::dataGenerator()}}
#' @param mygraph dataframe of class "graph" to constrain the changepoint inference, see \code{\link[gfpop:graph]{gfpop::graph()}}
#' @param type a string defining the cost model to use: \code{"mean"}, \code{"variance"}, \code{"poisson"}, \code{"exp"}, \code{"negbin"}
#' @param weights vector of weights (positive numbers), same size as data
#' @param iter.max maximal number of iteration of the gfpop function
#' @param D.init initialisation of the number of segments
#' @return a gfpop object = (\code{changepoints, states, forced, parameters, globalCost})
#' \describe{
#' \item{\code{changepoints}}{is the vector of changepoints (we give the last element of each segment)}
#' \item{\code{states}}{is the vector giving the state of each segment}
#' \item{\code{forced}}{is the vector specifying whether the constraints of the graph are active (\code{= TRUE}) or not (\code{= FALSE})}
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
  attr(response, "type") <- type

  return(response)
}
