
###############################################

#' Graph-contstrained functional pruning optimal partitioning
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

gfpop <- function(vectData = c(0), vectWeight = c(0), mygraph, type = "gauss", K = Inf, a = 0, min = -Inf, max = Inf)
{
  ###STOP###
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


  ###CALL Rcpp functions###
  if(type == "gauss"){res <- gfpopTransfert_Gauss(vectData, vectWeight, mygraph, K, a, min, max)}
  if(type == "poisson"){res <- gfpopTransfert_Poisson(vectData, vectWeight, mygraph, K, a, min, max)}
  if(type == "binomial"){res <- gfpopTransfert_Binomial(vectData, vectWeight, mygraph, K, a, min, max)}

  ###Response class gfpop###
  response <- list(changepoints = c(rev(res$changepoints[-1]), length(vectData)), states = rev(res$states), forced = rev(res$forced), means = rev(res$means), cost = res$cost)
  attr(response, "class") <- "gfpop"

  return(response)
}


