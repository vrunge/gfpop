rm(list=ls())

source("simulations/iso_functions.R")
library(parallel)
cores <- detectCores()
cores <- 40
#############
############# one.simu
#############

one.simu.getD <- function(i,
                     n = 10^2,
                     D = 10,
                     type = "isolinear",
                     noise = "gauss",
                     sigma = 10,
                     df = 3,
                     pcorr = 0.3)
{
  # PENALTY + K
  myK <- 3 * sigma
  myPen <- 2 * sigma * sigma * log(n)
  
  #DATA
  data <- simuData(n, D = D, type = type, noise = noise, pcorr = pcorr, df = df, sigma = sigma)
  sigma <- sdDiff(data$y)
  
  #Estimation signal
  estim <- estimateAllSeg(data, K = myK, pen = myPen)
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


#########################################################################
nbSimu <- 10
lres_1 <- mclapply(1:nbSimu, FUN = one.simu.getD,
                   n = 10^3,
                   D = 10,
                   type = "isostep",
                   noise = "gauss",
                   sigma = 1,
                   mc.cores = cores)

lres_2 <- mclapply(1:nbSimu, FUN = one.simu.getD,
                   n = 10^3,
                   D = 10,
                   type = "isostep",
                   noise = "student",
                   sigma = 1,
                   df = 3,
                   mc.cores = cores)

lres_3 <- mclapply(1:nbSimu, FUN = one.simu.getD,
                   n = 10^3,
                   D = 10,
                   type = "isostep",
                   noise = "corrupted",
                   sigma = 1,
                   pcorr = 0.3,
                   mc.cores = cores)

df_D_gauss <- do.call(rbind, lres_1)
df_D_student <- do.call(rbind, lres_2)
df_D_corrupted <- do.call(rbind, lres_3)

save(df_D_gauss, file="df_D_gauss.RData")
save(df_D_student, file="df_D_student.RData")
save(df_D_corrupted, file="df_D_corrupted.RData")

