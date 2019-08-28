
#devtools::install_github("vrunge/gfpop")
library(gfpop)

n <- 33
myData <- dataGenerator(n, c(0.1,0.3,0.5,0.8,1), c(1,2,1,3,1), sigma = 1)

myGraph <- graph(penalty = 1, type = "updown")
myGraph

penalty <- 1.45
mygraph <- graph(
  Edge("Dw", "Dw", type = "null", K = 3),
  Edge("Up", "Up", type = "null"),
  Edge("Dw", "Up", type = "up", penalty = 1.4),
  Edge("Up", "Dw", type = "down", penalty = 1.2),
  Node("Dw", min = 0, max = 1),
  StartEnd(start = "Dw", end = "Dw"))

gfpop(data = myData, mygraph = mygraph, type = "mean")




beta <- 2 * log(1000)
myGraph2 <- graph(
  Node("U", min = 0, max = 10),
  Node("D", min = -10, max = 0),
  Edge("D", "D", "null", decay = 0.9),
  Edge("U", "U", "null"),
  Edge("D", "U", "up", gap = 1, penalty = beta, K = 3, a = 5),
  Edge("D", "D", "down", penalty = beta),
  Edge("U", "D", "down", penalty = beta),
  StartEnd(start = "D", end = c("U")))
myGraph2
gfpop(data = myData, mygraph = myGraph2, type = "mean")
