
############
############ ALL MEAN
############

meth <- c(df_linear_gauss[10:14,2],"gfpop1","gfpop2","gfpop3","gfpop4")
df_linear_gauss[,2] <- rep(meth, 100)
df_linear_student[,2] <-  rep(meth, 100)
df_linear_corrupted[,2] <-  rep(meth, 100)
df_step_gauss[,2] <-  rep(meth, 100)
df_step_student[,2] <-  rep(meth, 100)
df_step_corrupted[,2] <-  rep(meth, 100)

df_linear_gauss$MSE <- as.numeric(df_linear_gauss$MSE)
df_linear_student$MSE <- as.numeric(df_linear_student$MSE)
df_linear_corrupted$MSE <- as.numeric(df_linear_corrupted$MSE)
df_step_gauss$MSE <- as.numeric(df_step_gauss$MSE)
df_step_student$MSE <- as.numeric(df_step_student$MSE)
df_step_corrupted$MSE <- as.numeric(df_step_corrupted$MSE)

with(df_linear_gauss, tapply(MSE, method, mean))
with(df_linear_student, tapply(MSE, method, mean))
with(df_linear_corrupted, tapply(MSE, method, mean))
with(df_step_gauss, tapply(MSE, method, mean))
with(df_step_student, tapply(MSE, method, mean))
with(df_step_corrupted, tapply(MSE, method, mean))

###

with(df_linear_gauss, tapply(MSE, method, sd))
with(df_linear_student, tapply(MSE, method, sd))
with(df_linear_corrupted, tapply(MSE, method, sd))
with(df_step_gauss, tapply(MSE, method, sd))
with(df_step_student, tapply(MSE, method, sd))
with(df_step_corrupted, tapply(MSE, method, sd))


########## 1 algo mean

df_step_gauss[df_step_gauss$method == "gfpop4",3]

##########



library(ggplot2)

#rename method
meth <- c(df_linear_gauss[10:14,2],"gfpop1","gfpop2","gfpop3","gfpop4")
df_linear_gauss[,2] <- rep(meth, 100)

#create factor
df_linear_gauss$method <- as.factor(df_linear_gauss$method)

#force numeric values
df_linear_gauss$MSE <- as.numeric(df_linear_gauss$MSE)

#plot
p <- ggplot(df_linear_gauss, aes(x = method, y = MSE, color=method)) + geom_violin()
p

#partial plot
me <- meth[2:9]
df <- df_linear_gauss[which(df_linear_gauss$method %in% me),]
p <- ggplot(df, aes(x = method, y = MSE, color=method)) + geom_violin()
p

mean(df$MSE)

#########



library(ggplot2)

#rename method
meth <- c(df_linear_corrupted[10:14,2],"gfpop1","gfpop2","gfpop3","gfpop4")
df_linear_corrupted[,2] <- rep(meth, 100)

#create factor
df_linear_corrupted$method <- as.factor(df_linear_corrupted$method)

#force numeric values
df_linear_corrupted$MSE <- as.numeric(df_linear_corrupted$MSE)

#plot
p <- ggplot(df_linear_corrupted, aes(x = method, y = MSE, color=method)) + geom_violin()
p

#partial plot
me <- meth[c(5,7,9)]
df <- df_linear_corrupted[which(df_linear_corrupted$method %in% me & df_linear_corrupted$MSE < 10),]
p <- ggplot(df, aes(x = method, y = MSE, color=method)) + geom_violin()
p





#rename method
meth <- c(df_step_student[1:5,2],"gfpop1","gfpop2","gfpop3","gfpop4")
df_step_student[,2] <- rep(meth, 100)

#create factor
df_step_student$method <- as.factor(df_step_student$method)

#force numeric values
df_step_student$MSE <- as.numeric(df_step_student$MSE)

#plot
p <- ggplot(df_step_student, aes(x = method, y = MSE, color=method)) + geom_violin()
p

#partial plot
me <- meth[3:9]
df <- df_step_student[which(df_step_student$method %in% me),]
p <- ggplot(df, aes(x = method, y = MSE, color=method)) + geom_violin()
p

#####################################################################
#####################################################################

df_linear_gauss$method <- as.factor(df_linear_gauss$method)
df_linear_gauss$MSE <- as.numeric(df_linear_gauss$MSE)
p <- ggplot(df_linear_gauss, aes(x = method, y = MSE, color=method)) + geom_violin() 
p+scale_y_continuous(trans='log2')

#####################################################################
#####################################################################

library(ggplot2)
library(cowplot)
library(gridExtra)

df_gauss
method <- unique(df_gauss[,2])
met <- df_gauss[,2]
MSE <- df_gauss[,3]

p <-  ggplot(df_gauss, aes(x=method, y=df_gauss[,3], color=method))
    + geom_violin() + geom_boxplot(width=0.05)
    + theme(legend.position = "none")
    + ggtitle("Cas 1")
p <- p + theme(legend.position = "none")
p


#######################

library(ggplot2)
library(cowplot)
library(gridExtra)

d1 <- data.frame(sc=rnorm(1000), met=rep(c("met1", "met2"), each=500))
d2 <- data.frame(sc=rnorm(1000), met=rep(c("met1", "met2"), each=500))
d3 <- data.frame(sc=rnorm(1000), met=rep(c("met1", "met2"), each=500))


p1 <-  ggplot(d1, aes(x=met, y=sc, color=met)) + geom_violin() + geom_boxplot(width=0.1)+ theme(legend.position = "none")+ ggtitle("Cas 1")
p2 <-  ggplot(d2, aes(x=met, y=sc, color=met)) + geom_violin() + geom_boxplot(width=0.1)+ theme(legend.position = "none")+ ggtitle("Cas 2")
p3 <-  ggplot(d3, aes(x=met, y=sc, color=met)) + geom_violin() + geom_boxplot(width=0.1)+ ggtitle("Cas 3")
legend <- cowplot::get_legend(p3)
p3 <- p3 + theme(legend.position = "none")

p <- grid.arrange(p1, p2, p3, legend, ncol = 4, nrow = 1)
