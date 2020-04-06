rm(list=ls())


source("simulations/iso_functions.R")
library(parallel)
cores <- detectCores()
cores <- 40
#############
############# one.simu
#############

one.simu <- function(i,
                     n = 10^2,
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
  data <- simuData(n, type = type, noise = noise, pcorr = pcorr, df = df, sigma = sigma)
  ## sigma <- sdDiff(data$y) # SIGMA ESTIMATION
  
  #Estimation signal
  estim <- estimateAllSeg(data, K = myK, pen = myPen)
  
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


#########################################################################
nbSimu <- 100
lres_1 <- mclapply(1:nbSimu, FUN = one.simu,
                   n = 10^4,
                   type = "isostep",
                   noise = "gauss",
                   sigma = 1,
                   mc.cores = cores)

lres_2 <- mclapply(1:nbSimu, FUN = one.simu,
                   n = 10^4,
                   type = "isostep",
                   noise = "student",
                   sigma = 1,
                   df = 3,
                   mc.cores = cores)

lres_3 <- mclapply(1:nbSimu, FUN = one.simu,
                   n = 10^4,
                   type = "isostep",
                   noise = "corrupted",
                   sigma = 1,
                   pcorr = 0.3,
                   mc.cores = cores)


df_step_gauss <- do.call(rbind, lres_1)
df_step_student <- do.call(rbind, lres_2)
df_step_corrupted <- do.call(rbind, lres_3)

save(df_step_gauss, file="df_step_gauss.RData")
save(df_step_student, file="df_step_student.RData")
save(df_step_corrupted, file="df_step_corrupted.RData")
