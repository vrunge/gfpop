<a id="top"></a>

[![Build Status](https://travis-ci.com/vrunge/gfpop.svg?branch=master)](https://travis-ci.com/vrunge/gfpop)
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

we present a basic use of the main functions of the `gfpop` package. More details about optional arguments are given in later sections.

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
## [1]  102  294  499  800 1000
## 
## $states
## [1] "Dw" "Up" "Dw" "Up" "Dw"
## 
## $forced
## [1] FALSE FALSE FALSE FALSE
## 
## $parameters
## [1] 1.0288625 2.1116082 0.9491628 3.0156831 1.0250846
## 
## $globalCost
## [1] 985.455
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
## [1]  142  288  391  773 1000
## 
## $states
## [1] "Iso" "Iso" "Iso" "Iso" "Iso"
## 
## $forced
## [1] FALSE FALSE FALSE FALSE
## 
## $parameters
## [1] 0.1364233 0.7537029 1.4215007 2.2474450 3.0758085
## 
## $globalCost
## [1] 1012.929
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
## [1]  302  603 1000
## 
## $states
## [1] "0" "1" "2"
## 
## $forced
## [1] FALSE FALSE
## 
## $parameters
## [1] 0.5888892 1.8474974 2.7282221
## 
## $globalCost
## [1] 1051.892
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
## [1]   98  315  499  800 1000
## 
## $states
## [1] "Dw" "Up" "Dw" "Up" "Dw"
## 
## $forced
## [1] FALSE  TRUE FALSE  TRUE
## 
## $parameters
## [1] -0.23357606  0.96815229 -0.03184771  1.06690908  0.06690908
## 
## $globalCost
## [1] 1098.621
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
##  [1]    2    3    5    6    7    8   84   85   86   91   93   94   98  106  107
## [16]  113  114  131  134  143  146  155  156  160  161  181  183  237  238  239
## [31]  240  270  271  273  274  282  283  285  287  304  306  312  360  361  365
## [46]  366  375  376  405  406  425  426  430  431  449  450  452  453  472  473
## [61]  516  519  523  524  529  530  564  584  585  631  632  633  634  636  637
## [76]  686  687  795  796  803  804  837  838  848  852  918  922  926  927  929
## [91]  930  942  943  971  972  988  995 1000
## 
## $states
##  [1] "Std" "Std" "Std" "Std" "Std" "Std" "Std" "Std" "Std" "Std" "Std" "Std" "Std"
## [14] "Std" "Std" "Std" "Std" "Std" "Std" "Std" "Std" "Std" "Std" "Std" "Std" "Std"
## [27] "Std" "Std" "Std" "Std" "Std" "Std" "Std" "Std" "Std" "Std" "Std" "Std" "Std"
## [40] "Std" "Std" "Std" "Std" "Std" "Std" "Std" "Std" "Std" "Std" "Std" "Std" "Std"
## [53] "Std" "Std" "Std" "Std" "Std" "Std" "Std" "Std" "Std" "Std" "Std" "Std" "Std"
## [66] "Std" "Std" "Std" "Std" "Std" "Std" "Std" "Std" "Std" "Std" "Std" "Std" "Std"
## [79] "Std" "Std" "Std" "Std" "Std" "Std" "Std" "Std" "Std" "Std" "Std" "Std" "Std"
## [92] "Std" "Std" "Std" "Std" "Std" "Std" "Std"
## 
## $forced
##  [1] FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE
## [14] FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE
## [27] FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE
## [40] FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE
## [53] FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE
## [66] FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE
## [79] FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE
## [92] FALSE FALSE FALSE FALSE FALSE FALSE
## 
## $parameters
##  [1]  0.06178639 -8.27208104  1.71168821 -4.81783933  1.15278300 -6.09072014
##  [7] -0.19140174  5.11567555 -5.78618752 -0.94242960  2.50575369 -6.83988745
## [13] -1.51565324  0.96243414  7.33692635  1.63036178  6.79917175  0.99848840
## [19]  4.37798186  1.30604558  5.76389093  0.20919818  7.37326310  0.76445922
## [25] -4.73666573  0.75754484 -3.50371839  0.99829912 -4.43633515  1.26313233
## [31]  7.94805696  0.83822608 -5.74696340  1.24127254 -5.10155941  1.40974770
## [37]  7.36998502  0.06788884  6.85701888  0.65828648 -5.59560686  2.28746692
## [43] -0.20049915 -5.55693509  0.80184365 -6.67442302  0.86457193  5.57421326
## [49] -0.22889545  5.22753154 -0.86865757  5.24146781 -1.03967567  6.16088158
## [55] -0.30807063  7.45001471 -0.16629008  6.09317372 -0.41061074  7.25682805
## [61]  0.73174815  4.30945492 -0.06389578  6.52082997  0.82133650 -4.57843202
## [67]  1.46415514  0.29382924 -3.86091167  1.11571070 -4.63358568  2.19402906
## [73] -4.24684618  1.89236680  6.88972000  0.92512320  6.65403597  1.12708120
## [79]  6.91043697  0.49331707  6.20559138  0.13319133  5.79442413 -0.09008908
## [85] -4.14091960 -0.02051769 -2.50858493  1.40673283 -6.54638603  0.44283023
## [91] -5.79390525  0.40883417 -5.58817124  0.04452492  6.84495380  0.17285350
## [97] -2.96421299  0.68521942
## 
## $globalCost
## [1] 1773.813
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
##  [1]  1000  2000  3000  4001  5000  6000  7000  8000  9000 10000
## 
## $states
##  [1] "0" "0" "0" "0" "0" "0" "0" "0" "0" "0"
## 
## $forced
## [1] FALSE FALSE FALSE  TRUE FALSE FALSE  TRUE  TRUE FALSE
## 
## $parameters
##  [1] -0.006179849  1.000929491 -0.011641923  1.987545524  0.987545524  1.993071129
##  [7] -0.008145532  0.991854468 -0.008145532  1.008023030
## 
## $globalCost
## [1] 2508.85
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
## [1]  200  500  800 1000
## 
## $states
## [1] "0" "0" "0" "0"
## 
## $forced
## [1] FALSE FALSE FALSE
## 
## $parameters
## [1] 0.0046423682 0.0003159481 0.0004786882 0.0199706094
## 
## $globalCost
## [1] 970.1886
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

![plot of chunk expdecay](figure/expdecay-1.png)


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

