#ISOTONIC REGRESSION

library(isotone) ## slow but more general
library(UniIsoRegression)
library(stats)
#https://www.sciencedirect.com/science/article/pii/S0167947308003964?via%3Dihub
## noise = corrupted simular to figure in https://arxiv.org/pdf/1707.09157.pdf
#reg_1d -> package UniIsoRegression
#isoreg -> package stats

#devtools::install_github("vrunge/gfpop", force = TRUE)
library(gfpop)

####################################################################################################
#########################
# DATA SIMULATION FUNCTION

simuData <- function(n = 10^3, D = 10, type = "isostep", noise = "student", dfree = 3, pcorr = 0.3, sigma = 10, length = "regular", min.length = 10)
{
  # n = nb points
  # k = nb steps
  # type = isostep OR isolinear
  # noise = student OR corrupted
  # dfree = for student distribution : degree of freedom
  # pcorr = percent of corrupted data
  ###
  ### type of function. DATA FROM 0 TO 100 !
  ###
  if(length == "regular")
  {
    if(type == "isostep")
    {
      data <- list()
      data$signal <- (floor(seq(from = 0, to = D-0.000001, length.out = n)))*100/D
      data$signal <-  data$signal - mean(data$signal)
      data$y <- data$signal ## modified later
      data$w <- rep(1, length(data$y))
      data$D <- D
    }
    if(type == "isolinear")
    {
      data <- list()
      data$signal <- (0:(n-1))*(100/n)
      data$signal <-  data$signal - mean(data$signal)
      data$y <- data$signal ## modified later
      data$w <- rep(1, length(data$y))
      data$D <- D
    }
  }

  if(length == "nonregular") #only case type == "isolinear"
  {
    data <- list()
    LengthSegment <- diff(c(0, sort(sample(n - D*min.length, D-1, replace = FALSE)), n - D*min.length)) + min.length
    data$signal <- rep(0:(D-1)*100/D, LengthSegment)
    data$signal <-  data$signal - mean(data$signal)
    data$y <- data$signal ## modified later
    data$w <- rep(1, length(data$y))
    data$D <- D
  }

  ###
  ### noise
  ###
  if(noise == "student")
  {
    data$y <- data$y + sigma*sqrt((dfree-2)/dfree)*rt(length(data$y), df = dfree) #Student t distribution (!= gauss; heavy tails)
  }

  if(noise == "corrupted")
  {
    data$y <- data$y + sigma*rnorm(length(data$y)) #gaussian distribution
    isp <- sample(1:n, trunc(pcorr*n))
    data$y[isp] <- -data$y[isp] #opposition for data indices in isp
  }

  if(noise == "gauss")
  {
    data$y <- data$y + sigma*rnorm(length(data$y)) #gaussian distribution
  }


  return(data)
}


################################################################################################
#########################
# ESTIMATE TIME AND FIT

estimateAllSeg <- function(simuData, type = "mean", K = 2, pen = 2)
{
  time <- c()
  fit <- list()

  time[1]  <- system.time(fit[[1]] <- isoreg(simuData$y)$yf)[3]
  time[2]  <- system.time(fit[[2]] <- reg_1d(simuData$y, simuData$w, metric = 2, unimodal = FALSE, decreasing = FALSE))[3]
  time[3]  <- system.time(fit[[3]] <- reg_1d(simuData$y, simuData$w, metric = 1, unimodal = FALSE, decreasing = FALSE))[3]

  time[4]  <- system.time(isoFpop4 <- gfpop(data = simuData$y, weights = simuData$w, mygraph = graph(penalty = 0, type = "isotonic"), type = type))[3]
  time[5]  <- system.time(isoFpop5 <- gfpop(data = simuData$y, weights = simuData$w, mygraph = graph(penalty = 0, type = "isotonic", K = K), type = type))[3]
  #time[6]  <- time[4]
  #time[7]  <- time[5]
  #isoFpop6 <- isoFpop4
  #isoFpop7 <- isoFpop5

  time[6]  <- system.time(isoFpop6 <- gfpop(data = simuData$y, weights = simuData$w, mygraph = graph(penalty = pen, type = "isotonic"), type = type))[3]
  time[7]  <- system.time(isoFpop7 <- gfpop(data = simuData$y, weights = simuData$w, mygraph = graph(penalty = pen, type = "isotonic", K = K), type = type))[3]

  fit[[4]] <- rep(isoFpop4$parameters, diff(c(0, isoFpop4$changepoints)))
  fit[[5]] <- rep(isoFpop5$parameters, diff(c(0, isoFpop5$changepoints)))
  fit[[6]] <- rep(isoFpop6$parameters, diff(c(0, isoFpop6$changepoints)))
  fit[[7]] <- rep(isoFpop7$parameters, diff(c(0, isoFpop7$changepoints)))

  response <- list(time =  time, fit = fit)
  return(response)
}




MSEall <- function(data, estimate)
{
  n <- length(data$signal)
  res <- data.frame(matrix(0, nrow = 1, ncol = 7))

  for(i in 1:7)
  {
    u <- data$signal - estimate$fit[[i]]
    res[1,i] <- (1/n)*(t(u)%*%u)
  }
  return(res)
}



MSEdataSignal <- function(data, signal, method = "NONE")
{
  if(method == "NONE")
  {
    n <- length(data)
    u <- data - signal
    res <- (1/n)*(t(u)%*%u)
  }
  if(method == "linear")
  {
    n <- length(data)
    x <- 1:n
    linearMod <- lm(data ~ x)
    line <- linearMod$coefficients[1] + linearMod$coefficients[2]*x
    u <- line - signal
    res <- (1/n)*(t(u)%*%u)
  }
  return(res)
}

####################################################################################################
#########################
# getK nb of segments for each estimate

getD <- function(estimate)
{
  res <- sapply(estimate$fit, FUN = function(x) sum(diff(x) != 0))
  return(res + 1)
}

################################################################################################
################################################################################################
################################################################################################

varDiff <- function(x)
{
  res = mad(diff(x)/sqrt(2))
  return(res)
}

#varDiff sous-estime la variance dans le cas student


varDiff <- function(x, method = 'HALL'){
  if(method == "SD")
  {
    return(sd(diff(x)/sqrt(2)))
  }
  if(method == "MAD")
  {
    return(mad(diff(x)/sqrt(2)))
  }
  if(method == "HALL")
  {
    n = length(x)
    wei <- c(0.1942, 0.2809, 0.3832, -0.8582)
    mat <- wei %*% t(x)

    mat[2, -n] = mat[2, -1]
    mat[3, -c(n-1, n)] = mat[3, -c(1, 2)]
    mat[4, -c(n-2, n-1, n)] = mat[4, -c(1, 2, 3)]
    return(sqrt(sum(apply(mat[, -c(n-2, n-1, n)], 2, sum)^2) / (n-3)))
  }
}



plotAllSeg <- function(data, estimate)
{
  plot(data$y, pch = "+", cex = 1.8) #cex = 1.2
  lty <- c(1, 2, 3, 4, 5, 6, 7)
  col <- c("black", "blue", "red", "yellow", "orange", "magenta", "green")
  for(i in 1:7) lines(estimate$fit[[i]], lty = lty[i], lwd = 2, col = col[i])
}




plotSeg <- function(data, estimate)
{
  plot(data$y, pch = ".", cex = 1) #cex = 1.2
  lty <- c(1, 3, 7)
  col <- c("red", "magenta", "blue")
  ttype <- c(4, 1, 5)
  lines(data$signal, lwd = 1.9, col = "black")
  for(i in 1:3) lines(estimate$fit[[lty[i]]], lwd = 2, lty = ttype[i], col = col[i])
  legend("bottomright", legend=c("true signal", "isoreg", "reg_1d L1", "gfpop4"),
         col=c("black",col), lty=c(1,ttype), lwd = 2, cex=1.5)
}


