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

Data is modelized by a quadratic cost with possible use of a robust loss, biweight and Huber. In a next version of this package, other parametric costs will be available (L1, Poisson, Binomial). 

The package `gfpop` is designed to segment univariate data $y_{1:n} = \{y_1,...,y_n\}$ obeying to a graph structure on segment means. The changepoint vector
$\tau = (\tau_0 < \cdots < \tau_{k+1}) \in \mathbb{N}^{k+2}$ defines the segments $\{\tau_i+1,...,\tau_{i+1}\}$, $i = 0,...,k$ with fixed $\tau_0 = 0$ and  $\tau_{k+1} = n$. We define the set $S_n = \{\hbox{changepoint vector } \tau \in \mathbb{N}^{k+2}\}\,,$ the nonconstrained minimal global cost is given by

$$Q_n = \min_{\tau \in S_n}\left[ \sum_{i=0}^{k}\lbrace \mathcal{C}(y_{(\tau_i+1):\tau_{i+1}}) + \beta \rbrace \right] \,,$$

where $\beta > 0$ is a penalty parameter and $\mathcal{C}(y_{u:v})$ is the minimal cost over the segment $\{u,...,v\}$. The penalty $\beta$ is an additional cost when introducing a new segment. The argmin of this quantity gives us a vector $\tau^*$ containing the last position of each segment (if we do not consider $\tau_0 = 0$).

With any cost, we have

$$\mathcal{C}(y_{(\tau_i+1):\tau_{i+1}}) = \min_{\theta}\mathcal{C}(y_{(\tau_i+1):\tau_{i+1}}, \theta)\,,\quad \hbox{and} \quad m_i = \mathrm{argmin}_{\theta}\mathcal{C}(y_{(\tau_i+1):\tau_{i+1}}, \theta)$$

defining the infered mean of the i+1-th segment $\{\tau_i+1,...,\tau_{i+1}\}$.

The graph $\mathcal{G}$ is defined by its vertices $V = \{1,...,v\} \subset \mathbb{N}$ and edges $E = \{e(i,j)\}$ where $i,j \in V$. Thus $\mathcal{G} = (V,E)$. The successive means $(m_1,...,m_{k+1})$ are constrainted to follow a feasible graph path $P = (e_1,...,e_k)$ with $e_i \in E$. Among all possible paths $P$, the path minimizing the cost is the result $Q_n(\mathcal{G})$ of our algorithm, that is, 

$$Q_n(\mathcal{G}) = \mathrm{argmin}_{P \in \mathcal{G}} (Q_n(P))$$

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
## changepoints
## [1]  100  294  499  800 1000
## 
## states
## [1] 0 1 0 1 0
## 
## forced
## [1] 0 0 0 0
## 
## means
## [1] 0.9506133 1.9537043 0.9972742 3.0211795 0.9845028
## 
## cost
## [1] 1007.106
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
## changepoints
## [1]  301  608 1000
## 
## states
## [1] 0 0 0
## 
## forced
## [1] 0 0
## 
## means
## [1] 0.500000 1.787239 2.785360
## 
## cost
## [1] 552.7358
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
## changepoints
## [1]   99  299  502  797 1000
## 
## states
## [1] 0 1 0 1 0
## 
## forced
## [1] 1 1 1 0
## 
## means
## [1]  0.04529466  1.04529466  0.04529466  1.04529466 -0.04628134
## 
## cost
## [1] 1785.262
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
## changepoints
##  [1]   42   43   63   64   83   84  103  104  139  140  168  169  194  195
## [15]  196  232  233  240  242  247  248  279  280  294  295  307  308  376
## [29]  377  396  397  417  418  419  420  465  469  482  483  484  485  499
## [43]  501  502  503  531  533  547  548  573  574  586  587  598  599  602
## [57]  603  630  631  633  634  723  725  750  753  767  771  789  790  801
## [71]  802  823  824  843  844  849  850  861  862  866  867  992  993  995
## [85]  997 1000
## 
## states
## integer(0)
## 
## forced
## integer(0)
## 
## means
##  [1]  0.58235050  5.95973617  0.11040567 -5.54332269  0.01597223
##  [6] -5.61655702  0.12396157  6.57021330  0.89009482  6.27698464
## [11]  1.81590845 -4.96846612  0.46173272 -4.01582950  6.07309505
## [16]  1.27510910 -5.23201631  2.00428563 -3.89357755  1.18392035
## [21]  6.13702410  0.32436918  7.31728642  1.21933710  8.18060251
## [26]  0.39171368  5.88250471  0.20508268  6.09000927 -0.14914062
## [31]  6.36188853  0.06691261 -5.88077680  0.82939686 -5.83859468
## [36]  0.24215220 -2.88861258  0.46322050 -5.37509130  5.99463190
## [41] -6.94589881 -0.27494679  6.32749038  0.42509342 -6.09242332
## [46]  1.91362165  5.49980317  0.24100133  6.29686995  1.08355394
## [51]  6.55502695  0.77659393 -4.68638846  1.90402514 -5.28117369
## [56]  0.84892591 -4.78662548  1.27959714 -5.01674930  0.39100256
## [61]  6.42241827  0.88404398 -3.18646047  1.33360243 -3.08100662
## [66]  0.87493475  3.98005805  0.99448371 -5.51246151  0.69887905
## [71] -6.11159623  0.35774709  6.40251185 -0.21096416 -6.13049289
## [76]  0.05862639 -7.23456763  0.74318809 -5.25649072 -0.28539193
## [81] -6.58885517 -0.02340459  5.44173244 -1.54571630  3.37271520
## [86] -4.34779099
## 
## cost
## [1] 2921.03
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
## changepoints
##  [1]  1000  2000  2998  3000  4002  5000  5999  6000  7003  8000  8999
## [12] 10000
## 
## states
##  [1] 0 0 0 0 0 0 0 0 0 0 0 0
## 
## forced
##  [1] 1 1 0 1 0 1 1 1 1 0 0
## 
## means
##  [1]  2.884239e-03  1.002884e+00  2.884239e-03  9.990353e-01  1.999035e+00
##  [6]  9.999281e-01  1.999928e+00  9.999281e-01 -7.188997e-05  9.999281e-01
## [11]  1.059444e-02  9.807759e-01
## 
## cost
## [1] 2736.378
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
## changepoints
##  [1]  1000  2000  3000  4001  5000  6000  7000  8000  9001 10000
## 
## states
##  [1] 0 0 0 0 0 0 0 0 0 0
## 
## forced
## [1] 0 1 0 0 1 0 1 0 0
## 
## means
##  [1] -0.025268821  1.002791848  0.002791848  2.024821363  0.995529680
##  [6]  1.995529680 -0.005607720  0.994392280 -0.022339504  1.020682552
## 
## cost
## [1] 2617.921
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
