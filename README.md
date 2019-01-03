---
output:
  html_document: default
    keep_md: true
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
## [1]  109  298  500  800 1000
## 
## $states
## [1] 0 1 0 1 0
## 
## $forced
## [1] 0 0 0 0
## 
## $means
## [1] 1.0608797 1.9996130 1.0340044 3.0336784 0.9909662
## 
## $cost
## [1] 1019.085
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
## [1]  284  594 1000
## 
## $states
## [1] 0 0 0
## 
## $forced
## [1] 0 0
## 
## $means
## [1] 0.6690761 1.7578968 2.8112646
## 
## $cost
## [1] 543.2958
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
## [1]  103  300  493  800 1000
## 
## $states
## [1] 0 1 0 1 0
## 
## $forced
## [1] 1 0 0 0
## 
## $means
## [1] -0.03258525  0.96741475 -0.13062033  1.07001294 -0.02777633
## 
## $cost
## [1] 1789.352
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
##  [1]    8   11   27   28   42   43   47   50   62   63  124  125  129  130
## [15]  162  163  165  166  167  200  202  205  210  217  219  220  234  235
## [29]  257  258  289  290  345  347  353  354  355  358  359  364  365  386
## [43]  387  409  410  414  415  428  429  470  471  491  535  537  564  565
## [57]  577  578  582  613  614  630  631  712  714  740  741  742  750  751
## [71]  790  791  820  821  873  874  875  891  892  905  906  932  933  934
## [85]  951  952  958  959  966  968  985  986 1000
## 
## $states
## integer(0)
## 
## $forced
## integer(0)
## 
## $means
##  [1] -0.781500279  3.779854671 -0.327578891  5.557008418  0.109480486
##  [6]  6.645229565 -0.305005416 -3.685530003  0.034889832 -5.418562484
## [11]  0.578214349 -5.347844094  1.117297218 -5.355136650  1.372299465
## [16]  7.546492405  1.231653406  6.805972386 -3.828533461  0.647762280
## [21] -5.461701043  2.760017054 -0.584432372  3.286562306 -0.050683534
## [26]  7.321523765  1.382605271 -4.527459248  0.969439366 -6.457723939
## [31]  0.494168219  6.631519140  0.241201921  5.934222902 -0.626102695
## [36] -5.569629080  4.271220503 -0.570818756 -6.429823833  0.512474326
## [41] -6.534349622 -0.070012238 -5.369512378  0.006643091  5.948396198
## [46]  0.099504702 -6.201517125 -0.226840617 -5.773408018 -0.066728528
## [51] -5.343085725 -0.534562447  1.075167629  5.187292290  1.618081377
## [56]  7.141281488  1.380197564 -4.183104832  3.976069595  1.219938431
## [61]  8.326243540  1.357163505 -4.138448781  1.051710677 -3.567050780
## [66]  1.267680060 -4.531871397  6.272091874  0.690479681 -5.276877360
## [71]  1.015572088 -4.653126541  0.130520037  6.072058754 -0.053786559
## [76] -5.032613903  6.529592065 -0.640117097  6.657996125 -0.587149300
## [81] -6.420441808 -0.488898770  5.034746047 -5.178975019  0.022080004
## [86]  5.489831968  0.152277147  6.085726098  0.152972942  4.981150804
## [91] -0.114816419 -6.189479391  0.333670143
## 
## $cost
## [1] 2960.465
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
##  [1]  1000  2000  2999  3000  4000  5000  6000  6001  7000  8001  8999
## [12] 10000
## 
## $states
##  [1] 0 0 0 0 0 0 0 0 0 0 0 0
## 
## $forced
##  [1] 1 0 0 1 0 1 1 1 0 0 0
## 
## $means
##  [1]  0.002876827  1.002876827  0.009253125  0.982244641  1.982244641
##  [6]  0.997217686  1.997217686  0.997217686 -0.002782314  0.996845171
## [11] -0.002054608  0.984199556
## 
## $cost
## [1] 2680.238
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
##  [1]  1001  2002  3000  4000  5000  6000  7002  8000  9000 10000
## 
## $states
##  [1] 0 0 0 0 0 0 0 0 0 0
## 
## $forced
## [1] 0 0 0 0 0 0 0 0 0
## 
## $means
##  [1]  0.00864413  1.03313497  0.01021457  2.01073808  1.00198514
##  [6]  2.01830073  0.01387377  1.01754926 -0.02043876  0.98441879
## 
## $cost
## [1] 2615.079
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
