<a id="top"></a>

<!--[![Build Status](http://travis-ci.com/vrunge/gfpop.svg?branch=master)](http://travis-ci.com/vrunge/gfpop)
--> 
[![](https://img.shields.io/badge/docs-vignettes-blue.svg)](https://github.com/vrunge/gfpop)

<!-- 
%\VignetteEngine{knitr::rmarkdown} 
%\VignetteIndexEntry{An Introduction to gfpop}
--> 


# gfpop Vignette
### Vincent Runge
#### LaMME, Evry University
### March 2, 2022

> [Quick Start](#qs)

> [Some examples](#se)

> [Graph construction](#gc)

> [Supplementary R functions](#suppl)


## Quick Start

we present a basic use of the main functions of the `gfpop` package. More details about the theory of graph-cosntrained multiple change-point detection can be found [here](https://arxiv.org/abs/2002.03646)

We install the package from Github:

```r
#devtools::install_github("vrunge/gfpop")
library(gfpop)
```

We simulate some univariate gaussian data (`n = 1000` points) with relative change-point positions `0.1, 0.3, 0.5, 0.8, 1` and means `1, 2, 1, 3, 1` with a variance equal to `1`.


```r
n <- 1000
myData <- dataGenerator(n, c(0.1,0.3,0.5,0.8,1), c(1,2,1,3,1), sigma = 1)
```

We define the graph of constraints to use for the dynamic programming algorithm. A simple case is the up-down constraint with a penalty here equal to a classic `2 log(n)`.


```r
myGraph <- graph(penalty = 2*log(n), type = "updown")
```

The gfpop function gives the result of the segmentation using `myData` and `myGraph` as parameters. We choose a gaussian cost.


```r
gfpop(data = myData, mygraph = myGraph, type = "mean")
```

```
## $changepoints
## [1]  100  298  500  800 1000
## 
## $states
## [1] "Dw" "Up" "Dw" "Up" "Dw"
## 
## $forced
## [1] FALSE FALSE FALSE FALSE
## 
## $parameters
## [1] 1.0317856 1.9865414 0.9697470 3.0359938 0.8527903
## 
## $globalCost
## [1] 1010.133
## 
## attr(,"class")
## [1] "gfpop" "mean"
```

The vector `changepoints` gives the last index of each segment. It always ends with the length of the vector `vectData`.

The vector `states` contains the states in which lies each mean. The length of this vector is the same as the length of `changepoint`.

The vector `forced` is a boolean vector. A forced element means that two consecutive means have been forced to satisfy the constraint. For example, the "up" edge with parameter c is forced if m(i+1) - m(i) = c.

The vector `parameters` contains the inferred means/parameters of the successive segments. 
 
The number `globalCost` is equal to the non-penalized cost, that is the value of the fit to the data ignoring the penalties for adding changes.

<a id="se"></a>

## Some examples

### Isotonic regression

The isotonic regression infers a sequence of nondecreasing means. 



```r
n <- 1000
mydata <- dataGenerator(n, c(0.1, 0.2, 0.3, 0.4, 0.6, 0.8, 1), c(0, 0.5, 1, 1.5, 2, 2.5, 3), sigma = 1)
myGraphIso <- graph(penalty = 2*log(n), type = "isotonic")
gfpop(data =  mydata, mygraph = myGraphIso, type = "mean")
```

```
## $changepoints
## [1]  186  361  713 1000
## 
## $states
## [1] "Iso" "Iso" "Iso" "Iso"
## 
## $forced
## [1] FALSE FALSE FALSE
## 
## $parameters
## [1] 0.1250638 1.2483475 2.1161519 2.8636266
## 
## $globalCost
## [1] 992.7441
## 
## attr(,"class")
## [1] "gfpop" "mean"
```

In this example, we use in `gfpop` function a robust biweight gaussian cost with `K = 1` and the `min` parameter in order to infer means greater than `0.5`.

### Fixed number of change-points

This algorithm is called segment neighborhood in the change-point litterature. In this example, we fixed the number of segments at 3 with an isotonic constraint. The graph contains two "up" edges with no cycling.



```r
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
```

```
## $changepoints
## [1]  323  705 1000
## 
## $states
## [1] "0" "1" "2"
## 
## $forced
## [1] FALSE FALSE
## 
## $parameters
## [1] 0.5031993 2.1549691 2.9441510
## 
## $globalCost
## [1] 1102.695
## 
## attr(,"class")
## [1] "gfpop" "mean"
```


### Robust up-down with constrained starting and ending states

In presence of outliers we need a robust loss (biweight). We can also force the starting and ending state and a minimal gap between the means (here equal to `1`)


```r
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
```

```
## $changepoints
## [1]  100  312  500  800 1000
## 
## $states
## [1] "Dw" "Up" "Dw" "Up" "Dw"
## 
## $forced
## [1]  TRUE FALSE FALSE FALSE
## 
## $parameters
## [1]  0.0456763883  1.0456763883 -0.0603308658  1.0383495117 -0.0003792976
## 
## $globalCost
## [1] 1110.176
## 
## attr(,"class")
## [1] "gfpop" "mean"
```



### Robust up-down with constrained starting and ending states

In presence of outliers we need a robust loss (biweight). We can also force the starting and ending state and a minimal gap between the means (here equal to `1`)


```r
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
```

```
## $changepoints
## [1]  113  300  500  796 1000
## 
## $states
## [1] "Dw" "Up" "Dw" "Up" "Dw"
## 
## $forced
## [1] FALSE FALSE  TRUE FALSE
## 
## $parameters
## [1] -0.17825898  1.09464865  0.00177342  1.00177342 -0.20725183
## 
## $globalCost
## [1] 1065.075
## 
## attr(,"class")
## [1] "gfpop" "mean"
```


If we skip all these constraints and use a standard fpop algorithm, the result is the following



```r
myGraphStd <- graph(penalty = 2*log(n), type = "std")
gfpop(data =  myData, mygraph = myGraphStd, type = "mean")
```

```
## $changepoints
##  [1]   29   30   51   56   58   65   66   99  101  143  144  145  197  199  207
## [16]  208  242  245  246  282  283  306  307  350  351  356  357  378  379  392
## [31]  393  426  427  429  430  446  447  494  496  509  510  529  530  556  557
## [46]  570  571  575  577  605  606  616  617  620  621  676  677  718  722  723
## [61]  746  747  769  770  780  781  821  822  829  830  892  893  898  899  908
## [76]  909  912  913  914  921  926  930  931  975  977  998 1000
## 
## $states
##  [1] "Std" "Std" "Std" "Std" "Std" "Std" "Std" "Std" "Std" "Std" "Std" "Std" "Std"
## [14] "Std" "Std" "Std" "Std" "Std" "Std" "Std" "Std" "Std" "Std" "Std" "Std" "Std"
## [27] "Std" "Std" "Std" "Std" "Std" "Std" "Std" "Std" "Std" "Std" "Std" "Std" "Std"
## [40] "Std" "Std" "Std" "Std" "Std" "Std" "Std" "Std" "Std" "Std" "Std" "Std" "Std"
## [53] "Std" "Std" "Std" "Std" "Std" "Std" "Std" "Std" "Std" "Std" "Std" "Std" "Std"
## [66] "Std" "Std" "Std" "Std" "Std" "Std" "Std" "Std" "Std" "Std" "Std" "Std" "Std"
## [79] "Std" "Std" "Std" "Std" "Std" "Std" "Std" "Std" "Std"
## 
## $forced
##  [1] FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE
## [14] FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE
## [27] FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE
## [40] FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE
## [53] FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE
## [66] FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE
## [79] FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE
## 
## $parameters
##  [1]  0.013157272  5.852235983  0.397219746 -2.487036978  4.306750060 -0.459086261
##  [7] -6.020547264 -0.580351710  5.041691711  0.767738682  5.961869168 -4.888515188
## [13]  1.104751718  6.573977997  0.319794566 -4.738167431  0.833793130 -2.899282191
## [19]  6.688410838  1.228639284 -4.269924074  0.404707460 -6.793861272  0.003291056
## [25]  6.506591287  0.629007228 -5.756058094 -0.259927168 -5.736196028  0.362430172
## [31]  5.983851635 -0.083938587 -6.728110389 -0.201933380 -5.596212715  0.515240779
## [37] -5.724921497 -0.233094949 -3.954707046  1.376006649  6.850232664  0.807237359
## [43] -4.772661389  1.572871812  7.044477776  0.419697210  6.576539303  0.668004319
## [49] -4.089779401  1.165036227  6.569345168  0.950959210  6.681891964  1.106775879
## [55]  7.174884098  0.932498620  7.593847560  0.593401906  3.024204730 -4.575301756
## [61]  1.347643390 -4.167946713  0.779419604  6.188344490  1.295509948  6.468646144
## [67]  0.008155765  5.288100272 -0.317542994  5.547977558  0.073240269  5.784884068
## [73]  0.018864107 -5.515696794  0.218333250 -6.987266782  0.513760015 -5.093780701
## [79]  5.643352084 -0.154778961 -4.128204915 -0.544777860  5.729519458 -0.404145572
## [85]  3.440665204  0.088765491  2.860110512
## 
## $globalCost
## [1] 1467.21
## 
## attr(,"class")
## [1] "gfpop" "mean"
```


### abs edge

With a unique `"abs"` edge, we impose a difference between the means of size at least 1.  



```r
n <- 10000
myData <- dataGenerator(n, c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1), c(0, 1, 0, 2, 1, 2, 0, 1, 0, 1), sigma = 0.5)
beta <- 2*log(n)
myGraph <- graph(
  Edge(0, 0,"abs", penalty = beta, gap = 1),
  Edge(0, 0,"null"))
gfpop(data =  myData, mygraph = myGraph, type = "mean")
```

```
## $changepoints
##  [1]   999  2000  3000  4000  4999  6000  7000  8000  8999 10000
## 
## $states
##  [1] "0" "0" "0" "0" "0" "0" "0" "0" "0" "0"
## 
## $forced
## [1]  TRUE  TRUE FALSE FALSE  TRUE FALSE FALSE FALSE  TRUE
## 
## $parameters
##  [1] -0.001621797  0.998378203 -0.001621797  2.018687659  1.000371240  2.000371240
##  [7]  0.007557667  1.039364965 -0.007114006  0.992885994
## 
## $globalCost
## [1] 2417.046
## 
## attr(,"class")
## [1] "gfpop" "mean"
```

Notice that some of the edges are forced, the vector `forced` contains non-zero values.


### Exponential decay

The null edge corresponds to an exponential decay state if its parameter is not equal to 1. 


```r
n <- 1000
mydata <- dataGenerator(n, c(0.2, 0.5, 0.8, 1), c(5, 10, 15, 20), sigma = 1, gamma = 0.966)
beta <- 2*log(n)
myGraphDecay <- graph(
  Edge(0, 0, "up", penalty = beta),
  Edge(0, 0, "null", 0, decay = 0.966)
  )
g <- gfpop(data =  mydata, mygraph = myGraphDecay, type = "mean")
g
```

```
## $changepoints
## [1]  200  322  500  800 1000
## 
## $states
## [1] "0" "0" "0" "0" "0"
## 
## $forced
## [1] FALSE FALSE FALSE FALSE
## 
## $parameters
## [1] 0.0049454559 0.1510457319 0.0028802127 0.0004721947 0.0194129929
## 
## $globalCost
## [1] 976.4308
## 
## attr(,"class")
## [1] "gfpop" "mean"
```


and we plot the result 


```r
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
```

![plot of chunk unnamed-chunk-12](figure/unnamed-chunk-12-1.png)


<a id="gc"></a>

## Graph construction

In the `gfpop` package, graphs are represented by a dataframe with 9 features and build with the R functions `Edge`, `Node`, `StartEnd` and `graph`.



```r
emptyGraph <- graph()
emptyGraph
```

```
## [1] state1    state2    type      parameter penalty   K         a        
## [8] min       max      
## <0 lignes> (ou 'row.names' de longueur nulle)
```


`state1` is the starting node of an edge, `state2` its ending node. `type` is one of the available edge type (`"null"`, `"std"`, `"up"`, `"down"`, `"abs"`). `penalty` is a nonnegative parameter: the additional cost $\beta_i$ to consider when we move within the graph using a edge (or stay on the same node). `parameter` is annother nonnegative parameter, a characteristics of the edge, depending of its type (it is a decay if type is "null" and a gap otherwise). `K` and `a` are robust parameters. `min` and `max` are used to constrain the rang of value for the node parameter.

We add edges into a graph as follows


```r
myGraph <- graph(
  Edge("E1", "E1", "null"),
  Edge("E1", "E2", "down", 3.1415, gap = 1.5)
)
myGraph
```

```
##   state1 state2 type parameter penalty   K a min max
## 1     E1     E1 null       1.0       0 Inf 0  NA  NA
## 2     E1     E2 down       1.5       0 Inf 0  NA  NA
```

we can only add edges to this dataframe using the object `Edge`.

The graph can contain information on the starting and/or ending edge to use with the `StartEnd` function. 


```r
beta <- 2 * log(1000)
myGraph <- graph(
  Edge("Dw", "Dw", "null"),
  Edge("Up", "Up", "null"),
  Edge("Dw", "Up", "up", penalty = beta, gap = 1),
  Edge("Dw", "Dw", "down", penalty = beta),
  Edge("Up", "Dw", "down", penalty = beta),
  StartEnd(start = "Dw", end = "Dw"))
myGraph
```

```
##   state1 state2  type parameter  penalty   K  a min max
## 1     Dw     Dw  null         1  0.00000 Inf  0  NA  NA
## 2     Up     Up  null         1  0.00000 Inf  0  NA  NA
## 3     Dw     Up    up         1 13.81551 Inf  0  NA  NA
## 4     Dw     Dw  down         0 13.81551 Inf  0  NA  NA
## 5     Up     Dw  down         0 13.81551 Inf  0  NA  NA
## 6     Dw   <NA> start        NA       NA  NA NA  NA  NA
## 7     Dw   <NA>   end        NA       NA  NA NA  NA  NA
```


Some graphs are often used: they are defined by default in the `graph` function. To use these graphs, we specify a string `type` equal to `"std"`, `"isotonic"`, `"updown"` or `"relevant"`.
For example,




```r
myGraphIso <- graph(penalty = 12, type = "isotonic")
myGraphIso
```

```
##   state1 state2 type parameter penalty   K a min max
## 1    Iso    Iso null         1       0 Inf 0  NA  NA
## 2    Iso    Iso   up         0      12 Inf 0  NA  NA
```

The function `Node` can be used to restrict the range of value for parameter associated to a node (called also a vertex). For example the following graph is an isotonic graph with inferred parameters between 0 et 1 only.



```r
myGraph <- graph(
  Edge("Up", "Up", "up", penalty = 3.1415),
  Edge("Up", "Up"),
  Node("Up", min = 0, max = 1)
  )
myGraph
```

```
##   state1 state2 type parameter penalty   K  a min max
## 1     Up     Up   up         0  3.1415 Inf  0  NA  NA
## 2     Up     Up null         1  0.0000 Inf  0  NA  NA
## 3     Up     Up node        NA      NA  NA NA   0   1
```

<a id="suppl"></a>

## Supplementary R functions

### Data generator function

the `dataGenerator` function is used to simulate `n` data-points from a distribution of `type` equal to `"mean"`, `"poisson"`, `"exp"`, `"variance"` or `"negbin"`. Standard deviation parameter `sigma` and decay `gamma` are specific to the Gaussian mean model. `size` is linked to the R `rnbinom` function from R stats package.

### Standard deviation estimation

We often need to estimate the standard deviation from the observed data to normalize the data or choose the edge penalties. The `sdDiff` returns such an estimation with the default HALL method [Hall et al., 1990] well suited for time series with change-points.


[Back to Top](#top)

