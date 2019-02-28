
###############################################

#' Graph-constrained functional pruning optimal partitioning
#'
#' @description Graph-contstrained functional pruning optimal partitionning
#' @param vectData vector of data to segment
#' @param vectWeight vector of weights (positive numbers) same size as vectData
#' @param mygraph dataframe of class graph to constraint the changepoint dynamic programming algorithm
#' @param type a string defining the type of cost to use. "gauss", "poisson" or "binomial"
#' @param K a positive number. Threshold for the Biweight robust loss
#' @param a a positive number. Slope for the Huber robust loss
#' @param min minimal bound for the infered means
#' @param max maximal bound for the infered means
#' @return a gfpop object = (changepoints, states, forced, means).
#' 'changepoints' is the vector of changepoints (we give the last element of each segment).
#' 'states' is the vector giving the state of each segment
#' 'forced' is the vector specifying whether the constraints of the graph are active (=1) or not (=0)
#' 'means' is the vector of successive means of each segment
#' 'cost' is a number equal to the global cost of the graph-constrained segmentation

gfpop <- function(vectData = c(0), vectWeight = c(0), mygraph, type = "gauss", K = Inf, a = 0, min = -Inf, max = Inf)
{
  ### STOP ###
    if(!any(class(mygraph) == "graph")){stop('Your graph is not a graph...')}
    if(type != "gauss" && type != "poisson" && type != "binomial")
    {stop('Arugment "type" not appropriate. Choose among "gauss", "poisson", "binomial"')}

    if(!is.double(K)){stop('K is not a double.')}
    if(K <= 0){stop('K must be positive (= Inf if not robust loss)')}

    if(!is.double(K)){stop('a is not a double.')}
    if(a < 0){stop('a must be nonnegative')}

    if(!is.double(min)){stop('min is not an double.')}
    if(!is.double(max)){stop('max is not an double.')}
    if(max <= min){stop('max is less than min...')}

    if(length(vectWeight) > 1)
    {
      if(length(vectData) != length(vectWeight)){stop('vectData and vectWeight have different size')}
      if(!all(vectWeight>0)){stop('vectWeight has non strictly positive components')}
    }

  if(type == "poisson"){stop('poisson loss not yet available')}
  if(type == "binomial"){stop('binomial loss not yet available')}


  ### BUILD an ordered Graph : myOrderedGraph ###
  startend <- mygraph[is.na(mygraph[,2]),]
  mygraph <-  mygraph[!is.na(mygraph[,2]),]
  maxVertex <- max(mygraph[,c(1,2)])

  myOrderedGraph <- graph()
  selectNullDecay <- mygraph[, 3] == "null"

  mygraph[selectNullDecay, 4] <- -1 #for ordering
  for(i in 0:maxVertex)
  {
    selectRaw <- mygraph[mygraph[,2]==i, ]
    ordre <- order(selectRaw[,4])
    selectRaw <- selectRaw[ordre,]
    myOrderedGraph <- rbind(myOrderedGraph, selectRaw)
  }

  myOrderedGraph <- rbind(myOrderedGraph, startend)
  selectNullDecay <- myOrderedGraph[, 3] == "null"
  myOrderedGraph[selectNullDecay, 4] <- 0 #for ordering

  ###CALL Rcpp functions###
  res <- gfpopTransfer(vectData, vectWeight, myOrderedGraph, type, K, a, min, max)

  ###Response class gfpop###
  response <- list(changepoints = c(rev(res$changepoints[-1]), length(vectData)), states = rev(res$states), forced = rev(res$forced), means = rev(res$means), cost = res$cost)
  attr(response, "class") <- "gfpop"

  return(response)
}



########################################################################################
# mygraph as penalties of type = sigma^2 or const * sigma^2

#' Graph-constrained functional pruning optimal partitioning iterated
#'
#' @description Graph-contstrained functional pruning optimal partitionning iterated with a Birgé Massart like penalty
#' @param vectData vector of data to segment
#' @param vectWeight vector of weights (positive numbers) same size as vectData
#' @param mygraph dataframe of class graph to constraint the changepoint dynamic programming algorithm
#' @param type a string defining the type of cost to use. "gauss", "poisson" or "binomial"
#' @param K a positive number. Threshold for the Biweight robust loss
#' @param a a positive number. Slope for the Huber robust loss
#' @param min minimal bound for the infered means
#' @param max maximal bound for the infered means
#' @param iter.max maximal number of iteration of the gfpop function
#' @param D.init initialisation of the number of segments
#' @return a gfpop object = (changepoints, states, forced, means).
#' 'changepoints' is the vector of changepoints (we give the last element of each segment).
#' 'states' is the vector giving the state of each segment
#' 'forced' is the vector specifying whether the constraints of the graph are active (=1) or not (=0)
#' 'means' is the vector of successive means of each segment
#' 'cost' is a number equal to the global cost of the graph-constrained segmentation
#' 'Dvect' is a vector of integers. The successive tested D in the Birgé Massart penalty until convergence

itergfpop <- function(vectData = c(0), vectWeight = c(0), mygraph, type = "gauss", K = Inf, a = 0, min = -Inf, max = Inf, iter.max = 100, D.init = 1)
{
  ### STOP ###
  if(!any(class(mygraph) == "graph")){stop('Your graph is not a graph...')}
  if(type != "gauss" && type != "poisson" && type != "binomial")
  {stop('Arugment "type" not appropriate. Choose among "gauss", "poisson", "binomial"')}

  if(!is.double(K)){stop('K is not a double.')}
  if(K <= 0){stop('K must be positive (= Inf if not robust loss)')}

  if(!is.double(K)){stop('a is not a double.')}
  if(a < 0){stop('a must be nonnegative')}

  if(!is.double(min)){stop('min is not an double.')}
  if(!is.double(max)){stop('max is not an double.')}
  if(max <= min){stop('max is less than min...')}

  if(length(vectWeight) > 1)
  {
    if(length(vectData) != length(vectWeight)){stop('vectData and vectWeight have different size')}
    if(!all(vectWeight>0)){stop('vectWeight has non strictly positive components')}
  }

  if(type == "poisson"){stop('poisson loss not yet available')}
  if(type == "binomial"){stop('binomial loss not yet available')}


  ### BUILD an ordered Graph : myOrderedGraph ###
  startend <- mygraph[is.na(mygraph[,2]),]
  mygraph <-  mygraph[!is.na(mygraph[,2]),]
  maxVertex <- max(mygraph[,c(1,2)])

  myOrderedGraph <- graph()
  for(i in 0:maxVertex)
  {
    selectRaw <- mygraph[mygraph[,2]==i, ]
    ordre <- order(selectRaw[,4])
    selectRaw <- selectRaw[ordre,]
    myOrderedGraph <- rbind(myOrderedGraph, selectRaw)
  }
  myOrderedGraph <- rbind(myOrderedGraph, startend)

  ####################################
  iter <- 0
  beta_old <- 0
  n <- length(vectData)
  beta <- getDerivativePenalty(D.init, n)

  #update myOrderedGraph
  myOrderedGraph[,4] <- myOrderedGraph[,4] * beta
  Dvect <- D.init

  while(beta_old != beta & iter <= iter.max)
  {
    #UPDATE GRAPH !!!

    ###CALL Rcpp functions###
    res <- gfpopTransfer(vectData, vectWeight, myOrderedGraph, K, a, min, max)

    beta_old <- beta

    myOrderedGraph[,4] <- myOrderedGraph[,4] / beta
    D <- length(res$changepoints)
    Dvect <- c(Dvect, D)
    beta <- getDerivativePenalty(D, n)
    myOrderedGraph[,4] <- myOrderedGraph[,4] * beta

    iter <- iter + 1
  }

  ####################################

  ###Response class gfpop###
  response <- list(changepoints = c(rev(res$changepoints[-1]), length(vectData)), states = rev(res$states), forced = rev(res$forced), means = rev(res$means), cost = res$cost, Dvect = rev(rev(Dvect)[-1]))
  attr(response, "class") <- "gfpop"

  return(response)
}





