
##### ISOTONIC REGRESSION ##### 

library(isotone) ## slow but more general
library(UniIsoRegression)
library(stats)
#https://www.sciencedirect.com/science/article/pii/S0167947308003964?via%3Dihub
## noise = corrupted simular to figure in https://arxiv.org/pdf/1707.09157.pdf
#reg_1d -> package UniIsoRegression
#isoreg -> package stats

##### gfpop ##### 

#devtools::install_github("vrunge/gfpop", force = TRUE)
library(gfpop)

##### Functions for isotonic simulations

#set working directory by automatic method
setwd(dirname(rstudioapi::getSourceEditorContext()$path))
source("6_Functions_for_isotonic_simulations.R")


######### ######### ######### ######### ######### 
######### ######### Figure 16 ######### ######### 
######### ######### ######### ######### ######### 
n = 10^4
sigma = 10
myK <- (2* sigma)^2
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
for(i in 1:3) lines(estimate$fit[[lty[i]]], lwd = 2, lty = ttype[i], col = col[i])

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
#cores <- 40


#########################################################################
nbSimu <- 100
my_n = 10^4
#########################################################################
### SIMULATION FOR TABLE 1
lres1_lin <- mclapply(1:nbSimu, FUN = one.simu,
                   n = my_n,
                   type = "isolinear",
                   noise = "gauss",
                   sigma = 1,
                   mc.cores = cores)

lres2_lin <- mclapply(1:nbSimu, FUN = one.simu,
                   n = my_n,
                   type = "isolinear",
                   noise = "student",
                   sigma = 1,
                   df = 3,
                   mc.cores = cores)

lres3_lin <- mclapply(1:nbSimu, FUN = one.simu,
                   n = my_n,
                   type = "isolinear",
                   noise = "corrupted",
                   sigma = 1,
                   pcorr = 0.3,
                   mc.cores = cores)

df_linear_gauss <- do.call(rbind, lres1_lin)
df_linear_student <- do.call(rbind, lres2_lin)
df_linear_corrupted <- do.call(rbind, lres3_lin)

#########################################################################
### SIMULATION FOR TABLE 2
lres1_step <- mclapply(1:nbSimu, FUN = one.simu,
                   n = my_n,
                   type = "isostep",
                   noise = "gauss",
                   sigma = 1,
                   mc.cores = cores)

lres2_step <- mclapply(1:nbSimu, FUN = one.simu,
                   n = my_n,
                   type = "isostep",
                   noise = "student",
                   sigma = 1,
                   df = 3,
                   mc.cores = cores)

lres3_step <- mclapply(1:nbSimu, FUN = one.simu,
                   n = my_n,
                   type = "isostep",
                   noise = "corrupted",
                   sigma = 1,
                   pcorr = 0.3,
                   mc.cores = cores)


df_step_gauss <- do.call(rbind, lres1_step)
df_step_student <- do.call(rbind, lres2_step)
df_step_corrupted <- do.call(rbind, lres3_step)


#########################################################################
### SIMULATION FOR TABLE 3
lres1_stepD <- mclapply(1:nbSimu, FUN = one.simu.getD,
                   n = my_n,
                   D = 10,
                   type = "isostep",
                   noise = "gauss",
                   sigma = 1,
                   mc.cores = cores)

lres2_stepD <- mclapply(1:nbSimu, FUN = one.simu.getD,
                   n = my_n,
                   D = 10,
                   type = "isostep",
                   noise = "student",
                   sigma = 1,
                   df = 3,
                   mc.cores = cores)

lres3_stepD <- mclapply(1:nbSimu, FUN = one.simu.getD,
                   n = my_n,
                   D = 10,
                   type = "isostep",
                   noise = "corrupted",
                   sigma = 1,
                   pcorr = 0.3,
                   mc.cores = cores)

df_D_gauss <- do.call(rbind, lres1_stepD)
df_D_student <- do.call(rbind, lres2_stepD)
df_D_corrupted <- do.call(rbind, lres3_stepD)


######################################################################################

############
############ TABLE 1 ISOLINEAR
############

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


############
############ TABLE 2 ISOSTEP
############

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

############
############ TABLE 3 D hat
############

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

#####################################################################
############
############ Figure 17: Violin plots of the MSE for iso-step simulations
############

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
