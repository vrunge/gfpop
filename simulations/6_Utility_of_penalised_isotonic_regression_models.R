
###
### How to install the gfpop PACKAGE
###

#(method 1) ON GITHUB 
#devtools::install_github("vrunge/gfpop", force = TRUE)

#(method 2) ON CRAN
install.packages("gfpop")
library(gfpop) 

############################################################################
## (SECTION 6. JSS PAPER) Utility of penalized isotonic regression models ## 
############################################################################

#-----------------------------------------------------#
## DATA SIMULATION FUNCTION for ISOTONIC SIMULATIONS ##
#-----------------------------------------------------#

simuData <- function(n = 10^3, 
                     D = 10, 
                     type = "isostep", 
                     noise = "student", 
                     dfree = 3, 
                     pcorr = 0.3, 
                     sigma = 10, 
                     min.length = 10)
{
  # n = nb points
  # k = nb steps
  # type = isostep OR isolinear
  # noise = student OR corrupted
  # dfree = for student distribution : degree of freedom
  # pcorr = percent of corrupted data
  ###
  ### ALL signals between 0 and 100 !
  ###

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
  
  ### noise
  if(noise == "gauss")
  {
    data$y <- data$y + sigma*rnorm(length(data$y)) #gaussian distribution
  }
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

  return(data)
}


#-----------------#
## ESTIMATE FIT  ##
#-----------------#


estimateAllSeg <- function(simuData, type = "mean", K = 2, pen = 2)
{
  response <- list()
  
  response[[1]] <- isoreg(simuData$y)$yf
  response[[2]] <- reg_1d(simuData$y, simuData$w, metric = 2, unimodal = FALSE, decreasing = FALSE)
  response[[3]] <- reg_1d(simuData$y, simuData$w, metric = 1, unimodal = FALSE, decreasing = FALSE)
  
  isoFpop4 <- gfpop(data = simuData$y, weights = simuData$w, mygraph = graph(penalty = 0, type = "isotonic"), type = type)
  isoFpop5 <- gfpop(data = simuData$y, weights = simuData$w, mygraph = graph(penalty = 0, type = "isotonic", K = K), type = type)
  isoFpop6 <- gfpop(data = simuData$y, weights = simuData$w, mygraph = graph(penalty = pen, type = "isotonic"), type = type)
  isoFpop7 <- gfpop(data = simuData$y, weights = simuData$w, mygraph = graph(penalty = pen, type = "isotonic", K = K), type = type)
  
  response[[4]] <- rep(isoFpop4$parameters, diff(c(0, isoFpop4$changepoints)))
  response[[5]] <- rep(isoFpop5$parameters, diff(c(0, isoFpop5$changepoints)))
  response[[6]] <- rep(isoFpop6$parameters, diff(c(0, isoFpop6$changepoints)))
  response[[7]] <- rep(isoFpop7$parameters, diff(c(0, isoFpop7$changepoints)))
  
  return(response)
}

#----------------------------#
## Compute MSE of 7 methods ##
#----------------------------#

MSEall <- function(data, estimate)
{
  n <- length(data$signal)
  res <- data.frame(matrix(0, nrow = 1, ncol = 7))
  
  for(i in 1:7)
  {
    u <- data$signal - estimate[[i]]
    res[1,i] <- (1/n)*(t(u)%*%u)
  }
  return(res)
}

#-----------------#
## MSEdataSignal ##
#-----------------#

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

#--------------------------#
## get number of segments ##
#--------------------------#

getD <- function(estimate)
{
  res <- sapply(estimate, FUN = function(x) sum(diff(x) != 0))
  return(res + 1)
}

#-----------#
## varDiff ##
#-----------#

varDiff <- function(x, method = 'HALL')
{
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


#--------------#
## plotAllSeg ##
#--------------#

plotAllSeg <- function(data, estimate)
{
  plot(data$y, pch = "+", cex = 1.8) #cex = 1.2
  lty <- c(1, 2, 3, 4, 5, 6, 7)
  col <- c("black", "blue", "red", "yellow", "orange", "magenta", "green")
  for(i in 1:7) lines(estimate[[i]], lty = lty[i], lwd = 2, col = col[i])
}


#--------------#
##  plotSeg   ##
#--------------#

plotSeg <- function(data, estimate)
{
  plot(data$y, pch = ".", cex = 1) #cex = 1.2
  lty <- c(1, 3, 7)
  col <- c("red", "magenta", "blue")
  ttype <- c(4, 1, 5)
  lines(data$signal, lwd = 1.9, col = "black")
  for(i in 1:3) lines(estimate[[lty[i]]], lwd = 2, lty = ttype[i], col = col[i])
  legend("bottomright", legend=c("true signal", "isoreg", "reg_1d L1", "gfpop4"),
         col=c("black",col), lty=c(1,ttype), lwd = 2, cex=1.5)
}




#############
############# one.simu
#############

one.simu <- function(i,
                     n = 10^2,
                     D = 10,
                     type = "isolinear",
                     noise = "gauss",
                     sigma = 10,
                     df = 3,
                     pcorr = 0.3,
                     estimating = "MSE",
                     sigma.estimate = FALSE)
{
  # PENALTY + K
  myK <- 3 * sigma
  myPen <- 2 * sigma * sigma * log(n)
  
  #DATA
  data <- simuData(n, D = D, type = type, noise = noise, 
                   pcorr = pcorr, df = df, sigma = sigma)
  if(sigma.estimate == TRUE){sigma <- sdDiff(data$y)}#SIGMA ESTIMATION
  
  #Estimation signal
  estim <- estimateAllSeg(data, K = myK, pen = myPen)
 
  #----------------------------------#
  if(estimating == "MSE")
  {
    #all MSE (x9)
    ref <- MSEdataSignal(data$y, data$signal)
    ref_lin <- MSEdataSignal(data$y, data$signal, method = "linear")
    mseAll <- MSEall(data, estim)
    #ind, method, mse, beta, K
    df <- data.frame(numeric(0), character(), numeric(0), numeric(0), numeric(0), stringsAsFactors = FALSE)
    colnames(df) <- c("index", "method", "MSE", "penalty", "K")
    df[1,] <- c(i,"MSE", ref, 0, Inf)
    df[2,] <- c(i,"MSE_linearFit", ref_lin, 0, Inf)
    df[3,] <- c(i,"isoreg", mseAll[1], 0, Inf)
    df[4,] <- c(i,"reg_1d_L2", mseAll[2], 0, Inf)
    df[5,] <- c(i,"reg_1d_L1", mseAll[3], 0, Inf)
    df[6,] <- c(i,"gfpop", mseAll[4], 0, Inf)
    df[7,] <- c(i,"gfpop", mseAll[5], 0, myK)
    df[8,] <- c(i,"gfpop", mseAll[6], myPen, Inf)
    df[9,] <- c(i,"gfpop", mseAll[7], myPen, myK)
    return(df)
  }
  #----------------------------------#
  if(estimating == "D")
  {
    #all D (x7)
    myD <- getD(estim)
    #ind, method, mse, beta, K
    df <- data.frame(numeric(0), character(), numeric(0), numeric(0), numeric(0), stringsAsFactors = FALSE)
    colnames(df) <- c("index", "method", "D", "penalty", "K")
    df[1,] <- c(i,"isoreg", myD[1], 0, Inf)
    df[2,] <- c(i,"reg_1d_L2", myD[2], 0, Inf)
    df[3,] <- c(i,"reg_1d_L1", myD[3], 0, Inf)
    df[4,] <- c(i,"gfpop", myD[4], 0, Inf)
    df[5,] <- c(i,"gfpop", myD[5], 0, myK)
    df[6,] <- c(i,"gfpop", myD[6], myPen, Inf)
    df[7,] <- c(i,"gfpop", myD[7], myPen, myK)
    return(df)
  }
}



#------------------------------------------------------
##### ISOTONIC REGRESSION ##### 

library(isotone) ## slow but more general
library(UniIsoRegression)
library(stats)
#https://www.sciencedirect.com/science/article/pii/S0167947308003964?via%3Dihub
## noise = corrupted simular to figure in https://arxiv.org/pdf/1707.09157.pdf
#reg_1d -> package UniIsoRegression
#isoreg -> package stats



#----------------------------------#
##            Figure 16           ##
#----------------------------------#
n = 10^4
sigma = 10
myK <- (2 * sigma) ^ 2
myPen <- 2 * sigma * sigma * log(n)
#DATA
data <- simuData(n, 
                 D= 10, 
                 type = "isostep", 
                 noise = "corrupted", 
                 pcorr = 0.3, 
                 sigma = sigma)

estim <- estimateAllSeg(data, K = myK, pen = myPen)

plot(data$y, pch = ".", cex = 1) #cex = 1.2
lty <- c(1, 3, 7)
col <- c("red", "magenta", "blue")
ttype <- c(4, 1, 5)
lines(data$signal, lwd = 1.9, col = "black")
for(i in 1:3){lines(estim[[lty[i]]], lwd = 2, lty = ttype[i], col = col[i])}

legend("bottomright", 
       legend=c("true signal", "isoreg", "reg_1d L1", "gfpop4"),
       col=c("black",col), 
       lty=c(1,ttype), 
       lwd = 2, 
       cex=1.5)

######### ######### ######### ######### ######### 
######### ######### ######### ######### ######### 
######### ######### ######### ######### ######### 


library(parallel)
cores <- detectCores()
cores <- 1 ### for windows user only
#cores <- 50


#########################################################################
nbSimu <- 100 #number of simulations
my_n = 10^4  #time-series data length
########################################################################

#############################
### SIMULATION FOR TABLE 1
lres1_lin <- mclapply(1:nbSimu, FUN = one.simu,
                      n = my_n,
                      type = "isolinear",
                      noise = "gauss",
                      sigma = 1,
                      estimating = "MSE",
                      mc.cores = cores)

lres2_lin <- mclapply(1:nbSimu, FUN = one.simu,
                      n = my_n,
                      type = "isolinear",
                      noise = "student",
                      sigma = 1,
                      df = 3,
                      estimating = "MSE",
                      mc.cores = cores)

lres3_lin <- mclapply(1:nbSimu, FUN = one.simu,
                      n = my_n,
                      type = "isolinear",
                      noise = "corrupted",
                      sigma = 1,
                      pcorr = 0.3,
                      estimating = "MSE",
                      mc.cores = cores)

df_linear_gauss <- do.call(rbind, lres1_lin)
df_linear_student <- do.call(rbind, lres2_lin)
df_linear_corrupted <- do.call(rbind, lres3_lin)

#############################
### SIMULATION FOR TABLE 2
lres1_step <- mclapply(1:nbSimu, FUN = one.simu,
                       n = my_n,
                       type = "isostep",
                       noise = "gauss",
                       sigma = 1,
                       estimating = "MSE",
                       mc.cores = cores)

lres2_step <- mclapply(1:nbSimu, FUN = one.simu,
                       n = my_n,
                       type = "isostep",
                       noise = "student",
                       sigma = 1,
                       df = 3,
                       estimating = "MSE",
                       mc.cores = cores)

lres3_step <- mclapply(1:nbSimu, FUN = one.simu,
                       n = my_n,
                       type = "isostep",
                       noise = "corrupted",
                       sigma = 1,
                       pcorr = 0.3,
                       estimating = "MSE",
                       mc.cores = cores)


df_step_gauss <- do.call(rbind, lres1_step)
df_step_student <- do.call(rbind, lres2_step)
df_step_corrupted <- do.call(rbind, lres3_step)


#############################
### SIMULATION FOR TABLE 3
lres1_stepD <- mclapply(1:nbSimu, FUN = one.simu,
                        n = my_n,
                        D = 10,
                        type = "isostep",
                        noise = "gauss",
                        sigma = 1,
                        estimating = "D",
                        mc.cores = cores)

lres2_stepD <- mclapply(1:nbSimu, FUN = one.simu.getD,
                        n = my_n,
                        D = 10,
                        type = "isostep",
                        noise = "student",
                        sigma = 1,
                        df = 3,
                        estimating = "D",
                        mc.cores = cores)

lres3_stepD <- mclapply(1:nbSimu, FUN = one.simu.getD,
                        n = my_n,
                        D = 10,
                        type = "isostep",
                        noise = "corrupted",
                        sigma = 1,
                        pcorr = 0.3,
                        estimating = "D",
                        mc.cores = cores)

df_D_gauss <- do.call(rbind, lres1_stepD)
df_D_student <- do.call(rbind, lres2_stepD)
df_D_corrupted <- do.call(rbind, lres3_stepD)


##########################################
#----------------------#
## TABLE 2: isolinear ##
#----------------------#

meth <- c(df_linear_gauss[1:5,2],"gfpop1","gfpop2","gfpop3","gfpop4")
df_linear_gauss[,2] <- rep(meth, nbSimu)
df_linear_student[,2] <-  rep(meth, nbSimu)
df_linear_corrupted[,2] <-  rep(meth, nbSimu)

df_linear_gauss$MSE <- as.numeric(df_linear_gauss$MSE)
df_linear_student$MSE <- as.numeric(df_linear_student$MSE)
df_linear_corrupted$MSE <- as.numeric(df_linear_corrupted$MSE)

with(df_linear_gauss, tapply(MSE, method, mean))
with(df_linear_student, tapply(MSE, method, mean))
with(df_linear_corrupted, tapply(MSE, method, mean))
with(df_linear_gauss, tapply(MSE, method, sd))
with(df_linear_student, tapply(MSE, method, sd))
with(df_linear_corrupted, tapply(MSE, method, sd))



#--------------------#
## TABLE 2: isostep ##
#--------------------#

meth <- c(df_step_gauss[1:5,2],"gfpop1","gfpop2","gfpop3","gfpop4")
df_step_gauss[,2] <- rep(meth, nbSimu)
df_step_student[,2] <-  rep(meth, nbSimu)
df_step_corrupted[,2] <-  rep(meth, nbSimu)

df_step_gauss$MSE <- as.numeric(df_step_gauss$MSE)
df_step_student$MSE <- as.numeric(df_step_student$MSE)
df_step_corrupted$MSE <- as.numeric(df_step_corrupted$MSE)

with(df_step_gauss, tapply(MSE, method, mean))
with(df_step_student, tapply(MSE, method, mean))
with(df_step_corrupted, tapply(MSE, method, mean))
with(df_step_gauss, tapply(MSE, method, sd))
with(df_step_student, tapply(MSE, method, sd))
with(df_step_corrupted, tapply(MSE, method, sd))

#------------------#
## TABLE 3: D HAT ##
#------------------#

meth <- c(df_D_gauss[1:3,2],"gfpop1","gfpop2","gfpop3","gfpop4")
df_D_gauss[,2] <- rep(meth, nbSimu)
df_D_student[,2] <-  rep(meth, nbSimu)
df_D_corrupted[,2] <-  rep(meth, nbSimu)

df_D_gauss$D <- as.numeric(df_D_gauss$D)
df_D_student$D <- as.numeric(df_D_student$D)
df_D_corrupted$D <- as.numeric(df_D_corrupted$D)

with(df_D_gauss, tapply(D, method, mean))
with(df_D_student, tapply(D, method, mean))
with(df_D_corrupted, tapply(D, method, mean))
with(df_D_gauss, tapply(D, method, sd))
with(df_D_student, tapply(D, method, sd))
with(df_D_corrupted, tapply(D, method, sd))



#---------------------------------------------------------------#
## Figure 17: Violin plots of the MSE for iso-step simulations ##
#---------------------------------------------------------------#

library(ggplot2)

#create factor
df_step_student$method <- as.factor(df_step_student$method)

#force numeric values
df_step_student$MSE <- as.numeric(df_step_student$MSE)

#partial plot
me <- c("isoreg", "reg_1d_L1", "gfpop2", "gfpop3", "gfpop4")
df <- df_step_student[which(df_step_student$method %in% me),]
p <- ggplot(df, aes(x = method, y = MSE, color=method)) + geom_violin()
p