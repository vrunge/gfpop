

### How to get the gfpop PACKAGE
#devtools::install_github("vrunge/gfpop", force = TRUE)
library(gfpop) #PACKAGE AVAILABLE ON CRAN


### ### ### ### ### ### ### ### 
###  4. The gfpop package   ### 
### ### ### ### ### ### ### ### 

## chunk 1 : An example of edge
e1 <- Edge(state1 = "Dw", state2 = "Up", type = "up", penalty = 10, gap = 0.5)
e1

## chunk 2 : An example of graph
graph(
  Edge(state1 = "mu0",state2 = "mu0", penalty = 0, K = 3),
  Edge(state1 = "mu0",state2 = "Coll", penalty = 10, type = "std"),
  Edge(state1 = "Coll",state2 = "Coll", penalty = 0),
  Edge(state1 = "Coll",state2 = "mu0", penalty = 0, type = "std", K = 3),
  StartEnd(start = "mu0", end = c("mu0", "Coll")),
  Node(state = "mu0", min = 0, max = 0)
)

## chunk 3 : Some default graphs.
graph(type = "isotonic", penalty = 12)

## chunk 4 : Gaussian model with an up-down graph
set.seed(1)
n <- 1000
myData <- dataGenerator(n, c(0.1, 0.3, 0.5, 0.8, 1), c(1, 2, 1, 3, 1), sigma = 1)
myGraph <- graph(penalty = 2 * log(n), type = "updown")
gfpop(data = myData, mygraph = myGraph, type = "mean")

## chunk 5 : Gaussian Robust biweight model with an up-down graph
n <- 1000
chgtpt <- c(0.1, 0.3, 0.5, 0.8, 1)
mydata <- dataGenerator(n, chgtpt, c(0, 1, 0, 1, 0), sigma = 1)
mydata <- mydata + 5*(rbinom(n, 1, 0.05)) - 5*(rbinom(n, 1, 0.05))
beta <- 2*log(n)
myGraph <- graph(
  Edge("Dw", "Up", type = "up", penalty = beta, gap = 1, K = 3),
  Edge("Up", "Dw", type = "down", penalty = beta, gap = 1, K = 3),
  Edge("Dw", "Dw", type = "null", K = 3),
  Edge("Up", "Up", type = "null", K = 3),
  StartEnd(start = "Dw", end = "Dw"))
gfpop(data =  mydata, mygraph = myGraph, type = "mean")

## chunk 6 : Poisson model with isotonic up graph
n <- 1000
chgtpt <- c(0.1, 0.3, 0.5, 0.8, 1)
mydata <- dataGenerator(n, chgtpt, c(1, 3, 5, 7, 12), sigma = 1, type = "poisson")
beta <- 2*log(n)
myGraph <- graph(type = "isotonic", gap = 2)
gfpop(data =  mydata, mygraph = myGraph, type = "poisson")

## chunk 7 : Negative binomial model with 3-segment graph
mygraph <- graph(
  Edge("1", "2", type = "std", penalty = 0),
  Edge("2", "3", type = "std", penalty = 0),
  StartEnd(start = "1", end = "3"), 
  all.null.edges = TRUE)
data <- dataGenerator(n = 1000, changepoints = c(0.3,0.7,1),
                      parameters = c(0.2,0.25,0.3), type = "negbin", sigma = 1)
gfpop(data, mygraph, type = "negbin")


## chunk 8 : A plotting function
data <- dataGenerator(1000, c(0.4, 0.8, 1), c(1, 2, 1), "mean", sigma = 3)
g <- gfpop(data, 
           graph(type = "std", penalty = 2*sdDiff(data)^2*(log(1000))), 
           type = "mean")
plot(x = g, data = data)

