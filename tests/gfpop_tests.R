
#devtools::install_github("vrunge/gfpop")
library(gfpop)

n <- 100
myData <- dataGenerator(n, c(0.1,0.3,0.5,0.8,1), c(1,2,1,3,1), sigma = 1)

myGraph <- graph(penalty = 2*log(n), type = "updown")
myGraph
gfpop(data = myData, mygraph = myGraph, type = "mean")


beta <- 2 * log(1000)
myGraph2 <- graph(
  Node("U", min = 0, max = 10),
  Node("D", min = -10, max = 0),
  Edge("D", "D", "null", decay = 0.9),
  Edge("U", "U", "null"),
  Edge("D", "U", "up", gap = 1, penalty = beta, K = 3, a = 5),
  Edge("D", "D", "down", penalty = beta),
  Edge("U", "D", "down", penalty = beta),
  StartEnd(start = "D", end = c("U","U")))
myGraph2
gfpop(data = myData, mygraph = myGraph2, type = "mean")
