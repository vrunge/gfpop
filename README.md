---
output:
  html_document: default
    variant: markdown_github
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

A nonnegative internal parameter can thus be associated to an edge (in "up" and "down") or is mandatory ("absSup" and "absInf"). More details on graph construction are given in the last [section](#gc). The graph definition can also include a starting and/or an ending edge and each edge has its own <a href="https://www.codecogs.com/eqnedit.php?latex=\beta" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\beta" title="\beta" /></a> penalty (see further).

The User has the possibility to constraint the infered means to lie between some minimal and maximal values.

Data is modelized by a quadratic cost with possible use of a robust loss, biweight and Huber. In a next version of this package, other parametric costs will be available (L1, Poisson, binomial). 

The package `gfpop` is designed to segment univariate data <a href="https://www.codecogs.com/eqnedit.php?latex=y_{1:n}&space;=&space;\{y_1,...,y_n\}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?y_{1:n}&space;=&space;\{y_1,...,y_n\}" title="y_{1:n} = \{y_1,...,y_n\}" /></a> obeying to a graph structure on segment means. The changepoint vector
<a href="https://www.codecogs.com/eqnedit.php?latex=\tau&space;=&space;(\tau_0&space;<&space;\cdots&space;<&space;\tau_{k&plus;1})&space;\in&space;\mathbb{N}^{k&plus;2}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\tau&space;=&space;(\tau_0&space;<&space;\cdots&space;<&space;\tau_{k&plus;1})&space;\in&space;\mathbb{N}^{k&plus;2}" title="\tau = (\tau_0 < \cdots < \tau_{k+1}) \in \mathbb{N}^{k+2}" /></a> defines the segments <a href="https://www.codecogs.com/eqnedit.php?latex=\{\tau_i&plus;1,...,\tau_{i&plus;1}\}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\{\tau_i&plus;1,...,\tau_{i&plus;1}\}" title="\{\tau_i+1,...,\tau_{i+1}\}" /></a>, <a href="https://www.codecogs.com/eqnedit.php?latex=i&space;=&space;0,...,k" target="_blank"><img src="https://latex.codecogs.com/gif.latex?i&space;=&space;0,...,k" title="i = 0,...,k" /></a> with fixed <a href="https://www.codecogs.com/eqnedit.php?latex=\tau_0&space;=&space;0" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\tau_0&space;=&space;0" title="\tau_0 = 0" /></a> and  <a href="https://www.codecogs.com/eqnedit.php?latex=\tau_{k&plus;1}&space;=&space;n" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\tau_{k&plus;1}&space;=&space;n" title="\tau_{k+1} = n" /></a>. We define the set <a href="https://www.codecogs.com/eqnedit.php?latex=S_n&space;=&space;\{\hbox{changepoint&space;vector&space;}&space;\tau&space;\in&space;\mathbb{N}^{k&plus;2}\}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?S_n&space;=&space;\{\hbox{changepoint&space;vector&space;}&space;\tau&space;\in&space;\mathbb{N}^{k&plus;2}\}" title="S_n = \{\hbox{changepoint vector } \tau \in \mathbb{N}^{k+2}\}" /></a>, the nonconstrained minimal global cost is given by

<a href="https://www.codecogs.com/eqnedit.php?latex=\begin{equation}&space;\label{cost}&space;Q_n&space;=&space;\min_{\tau&space;\in&space;S_n}\left[&space;\sum_{i=0}^{k}\lbrace&space;\mathcal{C}(y_{(\tau_i&plus;1):\tau_{i&plus;1}})&space;&plus;&space;\beta&space;\rbrace&space;\right]&space;\,,&space;\end{equation}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\begin{equation}&space;\label{cost}&space;Q_n&space;=&space;\min_{\tau&space;\in&space;S_n}\left[&space;\sum_{i=0}^{k}\lbrace&space;\mathcal{C}(y_{(\tau_i&plus;1):\tau_{i&plus;1}})&space;&plus;&space;\beta&space;\rbrace&space;\right]&space;\,,&space;\end{equation}" title="\begin{equation} \label{cost} Q_n = \min_{\tau \in S_n}\left[ \sum_{i=0}^{k}\lbrace \mathcal{C}(y_{(\tau_i+1):\tau_{i+1}}) + \beta \rbrace \right] \,, \end{equation}" /></a>

where <a href="https://www.codecogs.com/eqnedit.php?latex=\beta&space;>&space;0" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\beta&space;>&space;0" title="\beta > 0" /></a> is a penalty parameter and <a href="https://www.codecogs.com/eqnedit.php?latex=\mathcal{C}(y_{u:v})" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\mathcal{C}(y_{u:v})" title="\mathcal{C}(y_{u:v})" /></a> is the minimal cost over the segment <a href="https://www.codecogs.com/eqnedit.php?latex=\{u,...,v\}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\{u,...,v\}" title="\{u,...,v\}" /></a>. The penalty <a href="https://www.codecogs.com/eqnedit.php?latex=\beta" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\beta" title="\beta" /></a> is an additional cost when introducing a new segment. The argmin of this quantity gives us a vector $\tau^*$ containing the last position of each segment (if we do not consider <a href="https://www.codecogs.com/eqnedit.php?latex=\tau_0&space;=&space;0" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\tau_0&space;=&space;0" title="\tau_0 = 0" /></a>).

With any cost, we have

<a href="https://www.codecogs.com/eqnedit.php?latex=$$\mathcal{C}(y_{(\tau_i&plus;1):\tau_{i&plus;1}})&space;=&space;\min_{\theta}\mathcal{C}(y_{(\tau_i&plus;1):\tau_{i&plus;1}},&space;\theta)\,,\quad&space;\hbox{and}&space;\quad&space;m_i&space;=&space;\argmin_{\theta}\mathcal{C}(y_{(\tau_i&plus;1):\tau_{i&plus;1}},&space;\theta)$$" target="_blank"><img src="https://latex.codecogs.com/gif.latex?$$\mathcal{C}(y_{(\tau_i&plus;1):\tau_{i&plus;1}})&space;=&space;\min_{\theta}\mathcal{C}(y_{(\tau_i&plus;1):\tau_{i&plus;1}},&space;\theta)\,,\quad&space;\hbox{and}&space;\quad&space;m_i&space;=&space;\argmin_{\theta}\mathcal{C}(y_{(\tau_i&plus;1):\tau_{i&plus;1}},&space;\theta)$$" title="$$\mathcal{C}(y_{(\tau_i+1):\tau_{i+1}}) = \min_{\theta}\mathcal{C}(y_{(\tau_i+1):\tau_{i+1}}, \theta)\,,\quad \hbox{and} \quad m_i = \argmin_{\theta}\mathcal{C}(y_{(\tau_i+1):\tau_{i+1}}, \theta)$$" /></a>

defining the infered mean of the i+1-th segment <a href="https://www.codecogs.com/eqnedit.php?latex=\{\tau_i&plus;1,...,\tau_{i&plus;1}\}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\{\tau_i&plus;1,...,\tau_{i&plus;1}\}" title="\{\tau_i+1,...,\tau_{i+1}\}" /></a>.


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
## [1]  103  302  500  800 1000
## 
## $states
## [1] 0 1 0 1 0
## 
## $forced
## [1] 0 0 0 0
## 
## $means
## [1] 0.9551194 1.9972674 0.9983557 2.9706495 0.9207470
## 
## $cost
## [1] 972.2915
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
## [1]  331  794 1000
## 
## $states
## [1] 0 0 0
## 
## $forced
## [1] 0 0
## 
## $means
## [1] 0.5319267 2.2162633 3.3649432
## 
## $cost
## [1] 536.1175
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
## [1]  113  298  500  802 1000
## 
## $states
## [1] 0 1 0 1 0
## 
## $forced
## [1] 0 0 0 0
## 
## $means
## [1] -0.04169195  1.05685805 -0.07933690  1.05906530 -0.03644652
## 
## $cost
## [1] 1775.219
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
##  [1]   34   35   54   55  106  109  118  119  120  122  124  162  163  175
## [15]  176  178  179  202  203  215  216  279  280  293  295  311  312  397
## [29]  398  399  402  403  409  410  450  451  455  459  502  546  547  561
## [43]  562  581  583  584  598  599  654  657  680  681  711  712  757  758
## [57]  781  782  796  797  801  810  811  826  827  835  836  844  845  848
## [71]  849  892  893  900  901  910  927  928  938  939 1000
## 
## $states
## integer(0)
## 
## $forced
## integer(0)
## 
## $means
##  [1] -0.202899780 -5.252245894  0.442987108 -5.296136103 -0.239180059
##  [6]  3.680115465  0.644057191 -3.361670580  6.758106042  1.067800445
## [11] -4.852239914  0.947034781 -4.783622790  0.960542298 -4.046887975
## [16]  3.999659228 -4.708194707  1.473122831 -4.701830519  0.956939255
## [21] -5.318879701  1.418343422 -4.489795914  0.475824088  4.400712679
## [26]  0.484278547  6.577740091  0.003575062 -6.878841241  5.266309508
## [31] -0.501683304  5.419626537 -0.313664886  5.729254104 -0.217138645
## [36] -5.641325693  0.083074174 -3.454955674 -0.109649754  1.220078872
## [41] -5.090350094  0.506225505  6.278012254  0.780286836  5.159355290
## [46] -3.832445085  0.757025111  6.703107112  0.866826885  5.182727157
## [51]  0.886123706 -5.432512155  0.728464424  6.471864649  0.800756865
## [56]  7.493835860  1.421325245 -5.635293734  1.305615612 -4.721374998
## [61]  2.335801230 -0.033805335 -5.362515494  0.039620346  5.718571795
## [66]  0.848413833  6.463912778 -0.119862294 -5.470270168  1.303223748
## [71] -5.250252573 -0.704228466 -7.557327609  0.091767817 -4.036631928
## [76]  1.590200212 -0.477567740  6.015416093 -0.316810279  5.936119091
## [81] -0.158227561
## 
## $cost
## [1] 2863.719
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
##  [1]   999  2000  3000  3002  4000  5000  6000  6001  7001  8000  9002
## [12] 10000
## 
## $states
##  [1] 0 0 0 0 0 0 0 0 0 0 0 0
## 
## $forced
##  [1] 0 1 1 0 0 0 1 1 1 0 1
## 
## $means
##  [1]  0.009113768  1.008539723  0.008539723  1.008539723  2.001906635
##  [6]  1.021822820  1.992751113  0.992751113 -0.007248887  0.992751113
## [11]  0.012334452  1.012334452
## 
## $cost
## [1] 2663.17
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
##  [1]  1001  2000  3001  4000  5000  6000  6999  8000  9000 10000
## 
## $states
##  [1] 0 0 0 0 0 0 0 0 0 0
## 
## $forced
## [1] 1 0 0 1 0 0 0 0 1
## 
## $means
##  [1]  0.015553907  1.015553907 -0.021128362  2.000306219  1.000306219
##  [6]  2.001602952 -0.005200685  1.009990903  0.003067433  1.003067433
## 
## $cost
## [1] 2600.7
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
