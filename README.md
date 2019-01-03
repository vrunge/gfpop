---
output:
  html_document: default
    keep_md: true
    keep_tex: true
  pdf_document: default
  word_document: default
---

<!-- 
%\VignetteEngine{knitr::rmarkdown} 
%\VignetteIndexEntry{An Introduction to gfpop}
--> 

\DeclareMathOperator*{\argmin}{arg\,min}

<a id="top"></a>

# gfpop Vignette
### Vincent Runge
#### LaMME, Evry University
### January 3, 2019

> [Introduction](#intro)

> [Quick Start](#qs)

> [Some examples](#se)

> [Graph construction](#gc)

<a id="intro"></a>

## Introduction

`gfpop` is an R package for penalized parametric changepoint detection using functional pruning optimal partitioning (fpop) algorithm with additional constraints encoded in a graph (gfpop). Precisely, the successive infered means are constrained to follow edges of a graph (a path), with edges of type "up", "down", "std", "absInf" or "absSup" with the following meaning:

- "up" edge : the next segment has a greater mean (we can also force the size of the gap to be greater than a minimal value)
- "down" edge : the next segment has a lower mean (wan also can force the size of the gap to be greater that a minimal value)
- "std" edge : no contraint, the next segment can have any mean
- "absSup" edge : the absolute value of the difference of two consecutive means is greater than a given parameter
- "absInf" edge : the absolute value of the difference of two consecutive means is lower than a given parameter

A nonnegative internal parameter can thus be associated to an edge (in "up" and "down") or is mandatory ("absSup" and "absInf"). More details on graph construction are given in the last [section](#gc). The graph definition can also include a starting and/or an ending edge and each edge has its own $\beta$ penalty (see further).

The User has the possibility to constraint the infered means to lie between some minimal and maximal values.

Data is modelized by a quadratic cost with possible use of a robust loss, biweight and Huber. In a next version of this package, other parametric costs will be available (L1, Poisson, binomial). 

The package `gfpop` is designed to segment univariate data $y_{1:n} = \{y_1,...,y_n\}$ obeying to a graph structure on segment means. The changepoint vector
$\tau = (\tau_0 < \cdots < \tau_{k+1}) \in \mathbb{N}^{k+2}$ defines the segments $\{\tau_i+1,...,\tau_{i+1}\}$, $i = 0,...,k$ with fixed $\tau_0 = 0$ and  $\tau_{k+1} = n$. We define the set $S_n = \{\hbox{changepoint vector } \tau \in \mathbb{N}^{k+2}\}\,,$ the nonconstrained minimal global cost is given by
\begin{equation}
\label{cost}
Q_n = \min_{\tau \in S_n}\left[ \sum_{i=0}^{k}\lbrace \mathcal{C}(y_{(\tau_i+1):\tau_{i+1}}) + \beta \rbrace \right] \,,
\end{equation}
where $\beta > 0$ is a penalty parameter and $\mathcal{C}(y_{u:v})$ is the minimal cost over the segment $\{u,...,v\}$. The penalty $\beta$ is an additional cost when introducing a new segment. The argmin of this quantity gives us a vector $\tau^*$ containing the last position of each segment (if we do not consider $\tau_0 = 0$).

With any cost, we have
$$\mathcal{C}(y_{(\tau_i+1):\tau_{i+1}}) = \min_{\theta}\mathcal{C}(y_{(\tau_i+1):\tau_{i+1}}, \theta)\,,\quad \hbox{and} \quad m_i = \argmin_{\theta}\mathcal{C}(y_{(\tau_i+1):\tau_{i+1}}, \theta)$$
defining the infered mean of the i+1-th segment $\{\tau_i+1,...,\tau_{i+1}\}$.


The graph $\mathcal{G}$ is defined by its vertices $V = \{1,...,v\} \subset \mathbb{N}$ and edges $E = \{e(i,j)\}$ where $i,j \in V$. Thus $\mathcal{G} = (V,E)$. The successive means $(m_1,...,m_{k+1})$ are constrainted to follow a feasible graph path $P = (e_1,...,e_k)$ with $e_i \in E$. Among all possible paths $P$, the path minimizing the cost is the result $Q_n(\mathcal{G})$ of our algorithm, that is, 
$$Q_n(\mathcal{G}) = \argmin_{P \in \mathcal{G}} (Q_n(P))$$
and for a given path $P = (e_1,...,e_k)$ we definie the path-constrained cost

$$Q_n(P) = \min_{\tau \in S_n}\quad \min_{(\theta_0,\theta_1) \in e_1,..., (\theta_{k-1},\theta_k) \in e_k}\left[ \sum_{i=0}^{k}\lbrace \mathcal{C}(y_{(\tau_i+1):\tau_{i+1}}, \theta_i) + \beta_{e_i} \rbrace \right] $$
with $(\theta_{i-1},\theta_i) \in e_i$ meaning that the two consecutive means have to satisfy the edge constraint. $\beta_{e_i}$ is the penalty associated to edge $e_i$. In many cases, we simply take $\beta_{e_i} = \beta$ for all $i$. The penalty $\beta$ was a constant positive cost we have to pay when adding a new segment, that is when we move on an edge. Thus, we define this penalty within the edge and enable this penalty to be different for each edge.


For each path, the differences $\Delta m_i = m_{i+1} - m_i$ can be or not be constrained depending of the nature of the edge. For example, with an "up" edge with parameter $l > 0$, $\Delta m_i \ge l$,  with an "absInf" edge with parameter $l > 0$, $\Delta |m_i| \le l$, etc... 


<a id="qs"></a>

## Quick Start

we present a basic use of the main functions of the `gfpop` package. More details about optional arguments are given in later sections.

We install the package from Github:

```r
#devtools::install_github("vrunge/gfpop")
library(gfpop)
```

We simulate some univariate gaussian data (`n = 1000` points) with relative changepoint positions `0.1, 0.3, 0.5, 0.8, 1` and means `1, 2, 1, 3, 1` with a variance equal to `1`.

```r
n <- 1000
myData <- dataGenerator(n, c(0.1,0.3,0.5,0.8,1), c(1,2,1,3,1), 1)
```

We define the graph of constraints to use for the dynamic programming algorithm. A simple case is the up-down constraint with a penalty here equal to a classic `2 log(n)`.

```r
myGraph <- graph(penalty = 2*log(n), "updown")
```

The gfpop function gives the result of the segmentation using `myData` and `myGraph` as parameters. We choose a gaussian cost.

```r
gfpop(vectData = myData, mygraph = myGraph, type = "gauss")
```

```
## $changepoints
## [1]  100  302  500  800 1000
## 
## $states
## [1] 0 1 0 1 0
## 
## $forced
## [1] 0 0 0 0
## 
## $means
## [1] 1.0153233 2.0259772 0.9826006 3.0156150 1.0137614
## 
## $cost
## [1] 1135.162
## 
## attr(,"class")
## [1] "gfpop"
```

The vector `changepoints` gives the last index of each segment. It always ends with the length of the vector `vectData`.

The vector `states` contains the states in which lies each mean. The length of this vector is the same as the length of `changepoint`.

The vector `forced` is a boolean vector. A forced element means that two consecutive means have been forced to satisfy the constraint. For example, the "up" edge with parameter $l$ is forced if $m_{i+1} - m_i = l$.

The vector `means` contains the infered means of the successive segments. 
 
The number `cost` is equal to $Q_n(\mathcal{G})$, the overall cost of the segmented data. 


<a id="se"></a>

## Some examples

### Isotonic regression

The isotonic regression infer a sequence of nondecreasing means. 


```r
n <- 1000
mydata <- dataGenerator(n, c(0.1, 0.2, 0.3, 0.4, 0.6, 0.8, 1), c(0, 0.5, 1, 1.5, 2, 2.5, 3), 1)
myGraphIso <- graph(penalty = 2*log(n), type = "isotonic")
gfpop(vectData =  mydata, mygraph = myGraphIso, type = "gauss", K = 1, min = 0.5)
```

```
## $changepoints
## [1]  318  619 1000
## 
## $states
## [1] 0 0 0
## 
## $forced
## [1] 0 0
## 
## $means
## [1] 0.6571523 1.9468952 2.9720038
## 
## $cost
## [1] 569.0306
## 
## attr(,"class")
## [1] "gfpop"
```


In this example, we use in `gfpop` function a robust biweight gaussian cost with `K = 1` and the `min` parameter in order to infer means greater than `0.5`.


### Robust up-down with constrained starting and ending states

We use a robust loss (biweight) in presence of outliers. We can also force the starting and ending state and a minimal gap between the means (here equal to `1`)

```r
n <- 1000
mydata <- dataGenerator(n, c(0.1,0.3,0.5,0.8,1), c(0,1,0,1,0), 1) + 5*(rbinom(n, 1, 0.05)) - 5*(rbinom(n, 1, 0.05))

myGraph <- graph()
beta <- 2*log(n)
myGraph <- addEdge(myGraph, edge(0, 1, "up", beta, 1))
myGraph <- addEdge(myGraph, edge(1, 0, "down", beta, 1))
myGraph <- addStartEnd(myGraph, start = 0, end = 0)

gfpop(vectData =  mydata, mygraph = myGraph, type = "gauss", K = 3.0)
```

```
## $changepoints
## [1]   94  297  503  799 1000
## 
## $states
## [1] 0 1 0 1 0
## 
## $forced
## [1] 0 0 0 0
## 
## $means
## [1] -0.13580915  1.07256017  0.06444424  1.10837831 -0.05891668
## 
## $cost
## [1] 1979.513
## 
## attr(,"class")
## [1] "gfpop"
```

If we skip all these constraints, the result is the following


```r
myGraphStd <- graph(penalty = 2*log(n), type = "std")
gfpop(vectData =  mydata, mygraph = myGraphStd, type = "gauss")
```

```
## $changepoints
##   [1]   43   44   83   84   94   95  115  116  250  252  256  258  271  272
##  [15]  281  282  292  293  336  337  343  344  347  348  367  368  379  384
##  [29]  419  420  446  447  468  469  492  493  505  506  507  514  515  526
##  [43]  527  539  540  548  549  561  563  566  567  574  575  600  601  655
##  [57]  656  658  662  681  682  688  689  691  692  723  724  741  742  746
##  [71]  747  756  757  781  782  788  790  794  796  799  807  808  810  811
##  [85]  851  852  872  873  899  900  907  908  916  925  926  932  933  974
##  [99]  975 1000
## 
## $states
## integer(0)
## 
## $forced
## integer(0)
## 
## $means
##   [1] -0.143162078  5.681582505 -0.374412993  6.481580471 -0.082712576
##   [6] -4.575638081  1.535030098 -4.179953241  1.199456002  5.542332926
##  [11]  1.675896313 -4.332735136  1.243577102 -5.955714934  1.136310266
##  [16]  6.946385659  0.646021365  6.453158465  0.033141666  5.291921146
##  [21] -0.768889500  4.903157973 -0.685278913  5.259606772 -0.086624385
##  [26] -6.009980733 -0.069809896  2.460072361 -0.060433212 -5.650868793
##  [31]  0.168562134 -5.852337176  1.210253708  5.389748602 -0.351723218
##  [36] -6.076462603  0.055433948 -5.445959580  5.950487328  1.582337169
##  [41]  6.829579178  0.257703235  6.015315800  0.699066109  7.801720113
##  [46]  0.836255680 -4.827657838  0.518539225  5.708268923  0.711926153
##  [51]  5.723802278  0.051885756  6.573579749  0.909512294  7.884143356
##  [56]  1.498380627 -4.139747278  4.409443642 -1.385658044  1.334506759
##  [61]  7.045654747  0.856920040 -3.753626095  3.284033583 -3.530561531
##  [66]  1.435488981 -5.022745449  1.407247486 -4.994840329  1.832718264
##  [71] -5.402385919  1.217658039 -5.051499936  1.073685367 -5.673354164
##  [76]  1.797252961 -2.882980592  1.307361754  6.191579264  2.113971379
##  [81] -0.578679232  6.041791196 -0.460563737  5.746843297 -0.308339182
##  [86]  5.612642282 -0.463366644  4.928701749 -0.155486120 -5.626783457
##  [91] -0.002518706  6.098570662 -0.694478618  1.162029297 -5.207695607
##  [96]  0.482855462  5.935026651 -0.122424140  7.330087660  0.378236953
## 
## $cost
## [1] 3246.98
## 
## attr(,"class")
## [1] "gfpop"
```



### absInf and absSup edges

With a unique "absInf" edge, the differences between the means are at most of size 1.  

```r
n <- 10000
mydata <- dataGenerator(n, c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1), c(0, 1, 0, 2, 1, 2, 0, 1, 0, 1), 0.5)
myGraph <- graph()
beta <- 2*log(n)
myGraph <- addEdge(myGraph, edge(0, 0,"absInf", beta, 1))
gfpop(vectData =  mydata, mygraph = myGraph, type = "gauss", K = 3)
```

```
## $changepoints
##  [1]  1000  2000  3000  3003  3999  5000  6000  6001  7001  8000  9000
## [12] 10000
## 
## $states
##  [1] 0 0 0 0 0 0 0 0 0 0 0 0
## 
## $forced
##  [1] 1 0 1 0 0 1 1 0 0 0 1
## 
## $means
##  [1] -0.0157691992  0.9842308008  0.0011168142  1.0011168142  1.9999892007
##  [6]  1.0066564339  2.0066564339  1.0066564339  0.0156792489  0.9941819033
## [11]  0.0001429589  1.0001429589
## 
## $cost
## [1] 2681.673
## 
## attr(,"class")
## [1] "gfpop"
```

With a unique "absSup" edge, the differences between the means are at least of size 1.  

```r
n <- 10000
mydata <- dataGenerator(n, c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1), c(0, 1, 0, 2, 1, 2, 0, 1, 0, 1), 0.5)
myGraph <- graph()
beta <- 2*log(n)
myGraph <- addEdge(myGraph, edge(0, 0,"absSup", beta, 1))
gfpop(vectData =  mydata, mygraph = myGraph, type = "gauss", K = 3)
```

```
## $changepoints
##  [1]   991  2002  3000  4000  5000  6000  7000  8000  9000 10000
## 
## $states
##  [1] 0 0 0 0 0 0 0 0 0 0
## 
## $forced
## [1] 0 0 0 0 0 0 0 0 0
## 
## $means
##  [1] -0.016490648  1.036461482  0.008163164  2.015004335  0.995922493
##  [6]  2.006372821 -0.019650509  1.001105833 -0.022831840  0.997738686
## 
## $cost
## [1] 2643.852
## 
## attr(,"class")
## [1] "gfpop"
```


Notice that some of the edges are forced, the vector `forced` contains non-zero values.

<a id="gc"></a>

## Graph construction

In the `gfpop` package, graphs are represented by a dataframe with 5 features. 


```r
emptyGraph <- graph()
emptyGraph
```

```
## [1] state1    state2    type      penalty   parameter
## <0 rows> (or 0-length row.names)
```

`state1` is the starting vertex of an edge, `state2` its ending vertex. `type` is one of the available edge type ("up", "down", "std", "absInf", "absSup"). `penalty` is a nonnegative parameter: the additional cost $\beta_i$ to consider when we move within the graph using this edge. `parameter` is annother nonnegative parameter, a characteristics of the edge, depending of its type.

We add edges as follows

```r
emptyGraph <- addEdge(emptyGraph, edge(0, 0, "down", 3.1415, 1))
emptyGraph
```

```
##   state1 state2 type penalty parameter
## 1      0      0 down  3.1415         1
```

we can only add edges to this dataframe using the object `edge`. With the example `edge(0, 0, "down", 3.1415, 1)` we give the 5 aforementioned features.

The graph can contain information on the starting and/or ending edge to use. 


```r
myGraph <- graph()
beta <- 4*log(n)
myGraph <- addEdge(myGraph, edge(1, 0,"down", beta, 1))
myGraph <- addEdge(myGraph, edge(0, 0, "down", beta, 0))
myGraph <- addEdge(myGraph, edge(0, 1, "up", beta, 0))
myGraph <- addStartEnd(myGraph, start = 0, end = 0)
myGraph
```

```
##   state1 state2  type  penalty parameter
## 1      1      0  down 36.84136         1
## 2      0      0  down 36.84136         0
## 3      0      1    up 36.84136         0
## 4      0     NA start       NA        NA
## 5      0     NA   end       NA        NA
```


Some graphs are often used: they are defined by default in the `graph` function. To use these graphs, we specify a string `type` equal to "std", "isotonic", "updown" or "infsup".
For example,


```r
myGraphIso <- graph(penalty = 12, type = "isotonic")
myGraphIso
```

```
##   state1 state2 type penalty parameter
## 1      0      0   up      12         0
```



[Back to Top](#top)
