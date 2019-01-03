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


Some graphs are often used: they are defined by default in the `graph` function. To use these graphs, we specify a string `type` equal to "std", "isotonic", "updown"<!-- 
%\VignetteEngine{knitr::rmarkdown} 
%\VignetteIndexEntry{An Introduction to gfpop}
--> 


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
$$\DeclareMathOperator*{\argmin}{arg\,min}\mathcal{C}(y_{(\tau_i+1):\tau_{i+1}}) = \min_{\theta}\mathcal{C}(y_{(\tau_i+1):\tau_{i+1}}, \theta)\,,\quad \hbox{and} \quad m_i = \argmin_{\theta}\mathcal{C}(y_{(\tau_i+1):\tau_{i+1}}, \theta)$$
defining the infered mean of the i+1-th segment $\{\tau_i+1,...,\tau_{i+1}\}$.


The graph $\mathcal{G}$ is defined by its vertices $V = \{1,...,v\} \subset \mathbb{N}$ and edges $E = \{e(i,j)\}$ where $i,j \in V$. Thus $\mathcal{G} = (V,E)$. The successive means $(m_1,...,m_{k+1})$ are constrainted to follow a feasible graph path $P = (e_1,...,e_k)$ with $e_i \in E$. Among all possible paths $P$, the path minimizing the cost is the result $Q_n(\mathcal{G})$ of our algorithm, that is, 
$$\DeclareMathOperator*{\argmin}{arg\,min} Q_n(\mathcal{G}) = \argmin_{P \in \mathcal{G}} (Q_n(P))$$
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
## [1]  100  305  500  800 1000
## 
## states
## [1] 0 1 0 1 0
## 
## forced
## [1] 0 0 0 0
## 
## means
## [1] 0.9125526 1.9913644 1.0253301 3.0987031 0.9915935
## 
## cost
## [1] 1065.235
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
## [1]  378  630 1000
## 
## states
## [1] 0 0 0
## 
## forced
## [1] 0 0
## 
## means
## [1] 0.6298355 1.9696933 2.7447847
## 
## cost
## [1] 555.2229
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
## [1]  100  294  475  789 1000
## 
## states
## [1] 0 1 0 1 0
## 
## forced
## [1] 0 0 1 0
## 
## means
## [1] -0.10024695  1.03253136  0.01495186  1.01495186 -0.07439726
## 
## cost
## [1] 1888.009
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
##   [1]    9   26   27   38   39   40   55   56   57   81   87   90  100  136
##  [15]  137  198  199  201  210  212  219  220  226  229  233  234  237  238
##  [29]  247  248  263  271  272  298  299  302  304  313  314  319  322  334
##  [43]  345  347  348  352  354  355  403  404  422  424  453  454  491  492
##  [57]  493  535  536  537  538  546  547  558  559  562  563  565  566  589
##  [71]  590  599  600  601  609  610  611  623  624  629  631  652  653  681
##  [85]  685  698  699  724  725  730  734  794  854  855  865  866  901  905
##  [99]  938  939  983  984  985  988  989 1000
## 
## states
## integer(0)
## 
## forced
## integer(0)
## 
## means
##   [1]  1.58127963 -0.30133475 -6.82550836 -0.21944968  5.21216926
##   [6] -6.42331003 -0.70544639  4.94902262 -5.10064308 -0.55018502
##  [11]  1.45518811 -3.86747653 -0.43428973  1.21946215 -5.30510805
##  [16]  0.89898997 -5.60786647  5.77188700  1.10727775 -4.27899076
##  [21]  2.19469166  7.55529453  0.77832510  4.21559869 -0.07408956
##  [26] -5.20603852  1.05230828  7.15877316  0.70305774  8.09334797
##  [31]  1.24489767 -0.77860566  5.54677215  0.39792742  8.68705675
##  [36]  0.25501563 -5.28525562 -0.25106447  5.38090390  0.33301263
##  [41] -3.50214411  1.09714736 -0.58701084  5.66014614  1.09133086
##  [46] -3.91541980  0.72988269 -5.81503653  0.35867575 -5.26036022
##  [51]  0.04248718  4.78996865  0.37985875 -5.21352732  0.62016378
##  [56]  4.94663057 -5.10629748  0.40712449  7.40232496 -0.84334965
##  [61]  7.10601418  0.39077966  7.13247644  0.32073889  6.27849558
##  [66]  0.91911694 -4.79170207  4.00520924 -5.33687692  1.03925152
##  [71] -5.08942463  1.07468589  6.55249185 -4.62420242  0.51192863
##  [76] -4.69758940  7.68601264  1.04977328 -4.52711056  0.80988475
##  [81] -4.40029388  1.21954483 -4.61507748  0.30203481  5.29673735
##  [86]  0.75923284  7.46490515  1.10540875 -6.51244384  0.23263435
##  [91]  3.45960421  0.81278703 -0.09557443  5.42307019 -0.04376894
##  [96]  6.19570093 -0.16519583  3.28349501  0.03270829 -6.23504763
## [101] -0.19421287  5.78646170 -5.46293128  0.36699394  5.16269801
## [106] -0.39636117
## 
## cost
## [1] 3122.935
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
##  [1]   999  1999  3000  3001  4000  4998  6000  6001  7000  8001  9000
## [12] 10000
## 
## states
##  [1] 0 0 0 0 0 0 0 0 0 0 0 0
## 
## forced
##  [1] 1 0 1 0 0 0 1 0 0 0 0
## 
## means
##  [1] 0.001623115 1.001623115 0.002556100 1.002556100 1.995161386
##  [6] 1.011486554 1.992875983 0.992875983 0.013908776 0.989742839
## [11] 0.027096438 0.995170324
## 
## cost
## [1] 2687.888
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
##  [1]  1001  1998  3000  4000  5000  6000  7000  8000  9000 10000
## 
## states
##  [1] 0 0 0 0 0 0 0 0 0 0
## 
## forced
## [1] 0 1 0 1 1 0 1 0 0
## 
## means
##  [1]  0.016744546  1.017446553  0.017446553  2.006344283  1.006344283
##  [6]  2.006344283  0.001067912  1.001067912 -0.016846564  0.991117753
## 
## cost
## [1] 2610.746
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



[Back to Top](#top) or "infsup".
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
