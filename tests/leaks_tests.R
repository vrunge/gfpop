#devtools::install_github("vrunge/gfpop")
library(gfpop)


n <- 100
myData <- dataGenerator(n, c(0.1,0.3,0.5,0.8,1), c(1,2,1,3,1), sigma = 1)

myGraph <- graph(penalty = 2*log(n), type = "updown")
myGraph
gfpop(data = myData, mygraph = myGraph, type = "mean")

#CONSOLE TEST CONSOLE TEST CONSOLE TEST CONSOLE TEST CONSOLE TEST CONSOLE TEST CONSOLE TEST CONSOLE TEST CONSOLE TEST
#  R -d valgrind -f Test_Package.R > log.txt 2>&1
#R -d "valgrind --leak-check=full --show-reachable=yes" -f Test_Package.R > log_Test.Rout 2>&1
