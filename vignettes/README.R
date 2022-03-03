## -----------------------------------------------------------------------------
#devtools::install_github("vrunge/gfpop")
library(gfpop)

## -----------------------------------------------------------------------------
n <- 1000
myData <- dataGenerator(n, c(0.1,0.3,0.5,0.8,1), c(1,2,1,3,1), sigma = 1)

## -----------------------------------------------------------------------------
myGraph <- graph(penalty = 2*log(n), type = "updown")

## -----------------------------------------------------------------------------
gfpop(data = myData, mygraph = myGraph, type = "mean")

## -----------------------------------------------------------------------------
n <- 1000
mydata <- dataGenerator(n, c(0.1, 0.2, 0.3, 0.4, 0.6, 0.8, 1), c(0, 0.5, 1, 1.5, 2, 2.5, 3), sigma = 1)
myGraphIso <- graph(penalty = 2*log(n), type = "isotonic")
gfpop(data =  mydata, mygraph = myGraphIso, type = "mean")

## -----------------------------------------------------------------------------
n <- 1000
mydata <- dataGenerator(n, c(0.1, 0.2, 0.3, 0.4, 0.6, 0.8, 1), c(0, 0.5, 1, 1.5, 2, 2.5, 3), sigma = 1)
beta <- 0
myGraph <- graph(
  Edge(0, 1,"up", beta),
  Edge(1, 2, "up", beta),
  Edge(0, 0, "null"),
  Edge(1, 1, "null"),
  Edge(2, 2, "null"),
  StartEnd(start = 0, end = 2))

gfpop(data =  mydata, mygraph = myGraph, type = "mean")

## -----------------------------------------------------------------------------
n <- 1000
chgtpt <- c(0.1, 0.3, 0.5, 0.8, 1)
myData <- dataGenerator(n, chgtpt, c(0, 1, 0, 1, 0), sigma = 1)
myData <- myData + 5 * rbinom(n, 1, 0.05) - 5 * rbinom(n, 1, 0.05)
beta <- 2 * log(n)
myGraph <- graph(
         Edge("Dw", "Up", type = "up", penalty = beta, gap = 1, K = 3),
         Edge("Up", "Dw", type = "down", penalty = beta, gap = 1, K = 3),
         Edge("Dw", "Dw", type = "null", K = 3),
         Edge("Up", "Up", type = "null", K = 3),
         StartEnd(start = "Dw", end = "Dw"))
gfpop(data =  myData, mygraph = myGraph, type = "mean")

## -----------------------------------------------------------------------------
myGraphStd <- graph(penalty = 2*log(n), type = "std")
gfpop(data =  myData, mygraph = myGraphStd, type = "mean")

## -----------------------------------------------------------------------------
n <- 10000
myData <- dataGenerator(n, c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1), c(0, 1, 0, 2, 1, 2, 0, 1, 0, 1), sigma = 0.5)
beta <- 2*log(n)
myGraph <- graph(
  Edge(0, 0,"abs", penalty = beta, gap = 1),
  Edge(0, 0,"null"))
gfpop(data =  myData, mygraph = myGraph, type = "mean")

## -----------------------------------------------------------------------------
n <- 1000
mydata <- dataGenerator(n, c(0.2, 0.5, 0.8, 1), c(5, 10, 15, 20), sigma = 1, gamma = 0.966)
beta <- 2*log(n)
myGraphDecay <- graph(
  Edge(0, 0, "up", penalty = beta),
  Edge(0, 0, "null", 0, decay = 0.966)
  )
g <- gfpop(data =  mydata, mygraph = myGraphDecay, type = "mean")
g

## ---- expdecay----------------------------------------------------------------
gamma <- 0.966
len <- diff(c(0, g$changepoints))
signal <- NULL
for(i in length(len):1)
  {signal <- c(signal, g$parameters[i]*c(1, cumprod(rep(1/gamma,len[i]-1))))}
signal <- rev(signal)

ylimits <- c(min(mydata), max(mydata))
plot(mydata, type ='p', pch ='+', ylim = ylimits)
par(new = TRUE)
plot(signal, type ='l', col = 4, ylim = ylimits, lwd = 3)

## -----------------------------------------------------------------------------
emptyGraph <- graph()
emptyGraph

## -----------------------------------------------------------------------------
myGraph <- graph(
  Edge("E1", "E1", "null"),
  Edge("E1", "E2", "down", 3.1415, gap = 1.5)
)
myGraph

## -----------------------------------------------------------------------------
beta <- 2 * log(1000)
myGraph <- graph(
  Edge("Dw", "Dw", "null"),
  Edge("Up", "Up", "null"),
  Edge("Dw", "Up", "up", penalty = beta, gap = 1),
  Edge("Dw", "Dw", "down", penalty = beta),
  Edge("Up", "Dw", "down", penalty = beta),
  StartEnd(start = "Dw", end = "Dw"))
myGraph

## -----------------------------------------------------------------------------
myGraphIso <- graph(penalty = 12, type = "isotonic")
myGraphIso

## -----------------------------------------------------------------------------
myGraph <- graph(
  Edge("Up", "Up", "up", penalty = 3.1415),
  Edge("Up", "Up"),
  Node("Up", min = 0, max = 1)
  )
myGraph

