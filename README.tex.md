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
## [1]   97  306  499  800 1000
## 
## $states
## [1] 0 1 0 1 0
## 
## $forced
## [1] 0 0 0 0
## 
## $means
## [1] 0.9867554 1.9348412 0.9280818 2.9946186 0.8912949
## 
## $cost
## [1] 1000.553
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
## [1]  204  392 1000
## 
## $states
## [1] 0 0 0
## 
## $forced
## [1] 0 0
## 
## $means
## [1] 0.500000 1.296115 2.615312
## 
## $cost
## [1] 555.1794
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
## [1]  101  304  499  798 1000
## 
## $states
## [1] 0 1 0 1 0
## 
## $forced
## [1] 0 0 0 0
## 
## $means
## [1] -0.08860356  0.93617052 -0.10382370  1.09366956 -0.03921620
## 
## $cost
## [1] 1923.031
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
##   [1]    5    6    8    9   16   17   29   30   40   41   76   77   83   84
##  [15]  113  114  119  120  121  127  128  131  132  196  197  221  222  239
##  [29]  241  242  260  261  301  302  313  314  382  384  388  390  403  404
##  [43]  405  423  424  437  439  453  458  476  477  483  484  499  523  524
##  [57]  545  546  576  579  604  605  608  609  626  627  644  645  647  648
##  [71]  657  683  684  690  691  698  699  722  723  744  745  761  762  779
##  [85]  780  799  801  808  809  817  826  827  834  835  870  872  873  902
##  [99]  903  917  919  978  982  989  990 1000
## 
## $states
## integer(0)
## 
## $forced
## integer(0)
## 
## $means
##   [1] -0.28074641  6.82163355 -0.02711986  5.24340232 -1.14618822
##   [6] -5.68412774  0.06673554 -6.93351640 -0.61191049  5.36710790
##  [11]  0.14795004 -7.47655289  0.08336177 -6.29007759  0.34872751
##  [16] -4.24054619  1.66986423 -4.84178594  5.11555627  0.23077468
##  [21] -5.07028141  0.82342571  6.18710935  0.69692495  6.03873406
##  [26]  1.07030093 -4.63278137  2.20740632  6.32637182 -4.24762995
##  [31]  1.00045186 -4.90410319  0.92084769 -5.27659192  0.28542054
##  [36]  7.49756419  0.19404173  6.20001996  0.16441652  5.00587253
##  [41]  0.04647774  5.66070337 -4.90327974 -0.64111457  4.79226082
##  [46] -0.44654382 -5.68755843 -0.09746062  3.23854807  0.10212617
##  [51] -5.11868963  0.77501230  6.83027570  0.02136138  1.49954387
##  [56] -4.09622745  0.97858252 -6.68102402  1.14898242 -3.27202747
##  [61]  1.29992305 -4.64670865  0.78487316 -5.22377573  0.54916860
##  [66]  6.32217606  0.97437415  8.18196294  0.15868638 -5.45003150
##  [71]  2.73760485  1.16000703 -4.40038455  0.36041985 -4.58571322
##  [76]  1.61194173 -4.03600581  0.91019685  6.40097941  1.06336928
##  [81]  6.27742014  1.75976583 -5.66185483  1.03752881  7.48112163
##  [86]  1.25636363 -4.53500606 -0.04021723  5.79170834 -1.11699875
##  [91]  0.82887325 -6.28863666  1.06834175 -4.90507919 -0.11259664
##  [96] -4.34284718  6.32173004  0.44475074 -5.87861430 -0.02652423
## [101] -5.32901210  0.11967497  3.37706156 -0.32303418  6.59982300
## [106] -0.56357031
## 
## $cost
## [1] 3160.045
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
##  [1]  1000  2000  3000  3002  4000  5000  6000  6002  7000  8000  9000
## [12] 10000
## 
## $states
##  [1] 0 0 0 0 0 0 0 0 0 0 0 0
## 
## $forced
##  [1] 0 1 1 1 1 1 1 0 1 0 1
## 
## $means
##  [1] 0.010915715 1.003567530 0.003567530 1.003567530 2.003567530
##  [6] 1.003567530 2.003567530 1.003567530 0.004911061 1.004911061
## [11] 0.008377556 1.008377556
## 
## $cost
## [1] 2675.34
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
##  [1]  1000  2000  3000  3999  5004  6000  7000  8001  9001 10000
## 
## $states
##  [1] 0 0 0 0 0 0 0 0 0 0
## 
## $forced
## [1] 0 0 0 0 1 0 1 0 1
## 
## $means
##  [1]  0.006853368  1.011663891 -0.009694916  1.999076254  0.992060124
##  [6]  1.992060124 -0.001894111  0.998105889 -0.019339104  0.980660896
## 
## $cost
## [1] 2575.895
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
## [1] 1.084572 2.080527 0.980150 2.962705 1.018889
## 
## $cost
## [1] 1098.127
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
## [1]  289  716 1000
## 
## $states
## [1] 0 0 0
## 
## $forced
## [1] 0 0
## 
## $means
## [1] 0.500000 1.834487 2.858796
## 
## $cost
## [1] 544.2706
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
## [1]   98  299  500  805 1000
## 
## $states
## [1] 0 1 0 1 0
## 
## $forced
## [1] 1 0 0 0
## 
## $means
## [1]  0.066528175  1.066528175  0.007549165  1.105678226 -0.147898255
## 
## $cost
## [1] 1852.263
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
##  [1]    5    6   47   48   62   63   71   72  133  134  161  162  202  203
## [15]  287  288  293  294  300  301  316  317  380  381  432  433  434  453
## [29]  454  460  477  478  482  484  496  497  512  513  534  535  540  541
## [43]  558  560  565  567  574  575  584  586  587  591  594  651  652  657
## [57]  658  699  700  712  715  735  736  752  753  805  822  823  844  847
## [71]  889  890  938  939  997 1000
## 
## $states
## integer(0)
## 
## $forced
## integer(0)
## 
## $means
##  [1] -0.04305469  5.99010419  0.40702925 -6.15197617 -0.21114237
##  [6] -6.38803314  0.54828862 -4.87329323  0.57496404  8.26783147
## [11]  1.34134053  6.90361975  1.23880699  6.98354981  0.98574893
## [16] -4.57235711  1.26705550  6.54456140  0.73943324 -5.05229307
## [21]  0.18048973  5.89234199  0.33024900 -4.63610381 -0.25014586
## [26] -5.09347160  5.49740514 -0.34263172  7.02174418 -1.21948266
## [31]  1.03014059 -5.67148132 -0.09411212 -3.98280146  0.66520986
## [36]  6.36356708  0.78772253  6.83688116  0.65383333 -5.04171026
## [41]  1.76174390  7.74427225  1.05104948  5.21476303  0.78786096
## [46]  6.49764203  1.93576841 -3.88251028  1.27699368  4.66061463
## [51] -5.69168680  1.55459977 -3.19365248  0.99260920 -4.22733517
## [56]  0.46377523  6.30730774  0.88590144  7.19570146  1.39342784
## [61] -2.89735479  0.88629482 -6.68716586  1.19265333  6.92207636
## [66]  1.35460814  0.07853272 -5.56576294 -0.08622355  2.98788318
## [71] -0.23751024  5.34640520 -0.02403903  7.12546471 -0.26685535
## [76]  2.79638073
## 
## $cost
## [1] 2900.493
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
##  [1]  1006  1996  3000  3001  4001  5000  5996  6000  7003  8000  9000
## [12] 10000
## 
## $states
##  [1] 0 0 0 0 0 0 0 0 0 0 0 0
## 
## $forced
##  [1] 1 0 1 0 0 1 0 1 0 0 1
## 
## $means
##  [1] -0.0113752805  0.9886247195  0.0064298046  1.0064298046  1.9873619847
##  [6]  1.0087838773  2.0087838773  1.0234045829  0.0234045829  0.9694855950
## [11]  0.0007930579  1.0007930579
## 
## $cost
## [1] 2714.772
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
##  [1]  1000  1999  3000  4000  5000  6000  7000  8001  9000 10000
## 
## $states
##  [1] 0 0 0 0 0 0 0 0 0 0
## 
## $forced
## [1] 0 1 0 1 0 0 1 1 0
## 
## $means
##  [1] -0.0190730826  1.0034706113  0.0034706113  1.9906950312  0.9906950312
##  [6]  1.9940465640  0.0008343378  1.0008343378  0.0008343378  1.0080318609
## 
## $cost
## [1] 2630.11
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
