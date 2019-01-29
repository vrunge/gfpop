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
### January 29, 2019

> [Introduction](#intro)

> [Quick Start](#qs)

> [Some examples](#se)

> [Graph construction](#gc)

<a id="intro"></a>

## Introduction

`gfpop` is an `R` package developed to perform parametric changepoint detection in univariate time series constrained to a graph structure. The constraints are imposed to the sequence of infered means of consecutive segments and are related to the direction and/or the magnitude of the mean changes. Changepoint detection is performed using the functional pruning optimal partitioning method (fpop) based on an exact dynamic programming algorithm.  The user can specify some other constraints on the graph (starting and ending vertices) and constraint the range of means. 

For each data point, the algorithm updates a function and we go througth an edge. The edges of the graph can be of type "null", up", "down", "std", "absInf" or "absSup" with the following meaning:

- "null" edge : there is no changepoint introduced. We stay on the same segment
- "up" edge : the next segment has a greater mean (we can also force the size of the gap to be greater than a minimal value)
- "down" edge : the next segment has a lower mean (wan also can force the size of the gap to be greater that a minimal value)
- "std" edge : no contraint, the next segment can have any mean
- "absSup" edge : the absolute value of the difference of two consecutive means is greater than a given parameter
- "absInf" edge : the absolute value of the difference of two consecutive means is lower than a given parameter

A nonnegative internal parameter can thus be associated to an edge (in "up" and "down") or is mandatory ("absSup" and "absInf"). The "null" edge refers to an absence of changepoint. If this edge is not present, we are able to constraint segment lengths. Our package includes the segment neighborhood algorithm for which the number of changepoints is fixed. More details on graph construction are given in the last [section](#gc).

Data is modelized by a quadratic cost with possible use of a robust loss, biweight and Huber. In a next version of this package, other parametric costs will be available (L1, Poisson). 

The package `gfpop` is designed to segment univariate data $y_{1:n} = \{y_1,...,y_n\}$ obeying to a graph structure on segment means. The changepoint vector $\overline{\tau} = (\tau_0 < \cdots < \tau_{k+1}) \in \mathbb{N}^{k+2}$ defines the $k+1$ segments $\{\tau_i+1,...,\tau_{i+1}\}$, $i = 0,...,k$ with fixed bounds $\tau_0 = 0$ and  $\tau_{k+1} = n$. We use the set $S_n = \{\hbox{changepoint vector } \overline{\tau} \in \mathbb{N}^{k+2}\}$ to define the nonconstrained minimal global cost given by

$$Q_n = \min_{\overline{\tau} \in S_n}\left[ \sum_{i=0}^{k}\lbrace \mathcal{C}(y_{(\tau_i+1):\tau_{i+1}}) + \beta \rbrace \right]\,,$$

where $\beta > 0$ is a penalty parameter and $\mathcal{C}(y_{u:v})$ is the minimal cost over the segment $\{u,...,v\}$. The penalty $\beta$ is understood as an additional cost when introducing a new segment. The argminimum of this quantity gives a vector $\tau^*$ containing the last position of each segment (if we do not consider $\tau_0 = 0$). The quantity $Q_n$ is the solution of the nonconstrained optimal partitioning method. 

In our setting, the cost $\mathcal{C}$ is the result of the minimization of a cost function with additive property:

$$ \left\{
\begin{array}{l}
  \mathcal{C}(y_{(\tau_i+1):\tau_{i+1}}) = \min_{\theta_i}\mathcal{C}(y_{(\tau_i+1):\tau_{i+1}}, \theta_i) = \min_{\theta_i} \sum_{j = \tau_i+1}^{\tau_{i+1}}\gamma(y_j, \theta_i)\,,\\
  m_i = \argmin_{\theta_i}\mathcal{C}(y_{(\tau_i+1):\tau_{i+1}}, \theta_i)\,,
\end{array}
\right.$$

with the argminimum defining the infered mean $m_i$ of the i+1-th segment $\{\tau_i+1,...,\tau_{i+1}\}$ with $i \in \{0,...,k\}$. Additivity of the cost (the $\gamma$ decomposition) is guaranteed as we will use costs deriving from a likelihood. 

More generally, we can associate a current mean $\mu_i$ to each data point $y_i$ and we write a cost on a segment $\{u,...,v\}$ as a result of a constaint minimization:

$$\mathcal{C}(y_{u:v}) = \min_{\overset{(\mu_u,...,\mu_v)\in \mathbb{R}^{{v-u+1}}}{\mu_u = \cdots = \mu_v}}\sum_{j = u}^{v}\gamma(y_j, \mu_j)\,,$$

so that we get another description of the objective function :

$$Q_n = \min_{(\mu_1,...,\mu_n)\in \mathbb{R}^{n}\,,\,\mu_{n+1} = +\infty}\quad\sum_{i=1}^n\gamma(y_i,\mu_i) + \beta I(\mu_{i} \ne \mu_{i+1})\,,$$

with $I$ the indicator function. This approach can be generalized to more complex constraints on consecutive means. We define a transition graph $\mathsf{G} = (\mathcal{S}, \mathcal{T}(1),...,\mathcal{T}(n-1))$ as a set of states $\mathcal{S} =  \{0,...,S\} \subset \mathbb{N}$ and transition sets $\mathcal{T}(t)$ constraining the adding of data $y_{t+1}$:

$$\mathcal{T}(t) \subset \mathcal{S} \times \mathcal{S} \times \mathcal{B}\,,$$

with

$$\mathcal{B} = \{\hbox{constraint function}\,\,  g : \mathbb{R}\times \mathbb{R} \mapsto \mathbb{R} \,,\, \hbox{penalty}\,\beta \ge 0 \}\,.$$

An element $e(t) \in \mathcal{T}(t)$ is described by the elements $e(t) = (s_t,s_{t+1},g_{e(t)},\beta^{e(t)})$. The transition function $g_{e(t)}$ defines the kind of constraint we use in the transition from state $s_t$ to state $s_{t+1}$ and penalized by $\beta^{e(t)}$. The following able summarizes all the possible constraints encoded into the `gfpop` package.


| constraints | $g_{e(t)} : \mathbb{R}\times \mathbb{R} \mapsto \mathbb{R}$, $c \in \mathbb{R}^+$ |
|---------:|-----:|
| no changepoint | $g(\mu_{t},\mu_{t+1}) = (\mu_{t} - \mu_{t+1})^2 \le 0$ |
| up | $g(\mu_{t},\mu_{t+1}) = \mu_{t} - \mu_{t+1}  + c \le 0$ |
| down | $g(\mu_{t},\mu_{t+1}) = -\mu_{t} + \mu_{t+1} + c \le 0$ |
| absSup | $g(\mu_{t},\mu_{t+1}) = -\ell_1(\mu_{t} - \mu_{t+1}) + c \le 0$ |
| absInf | $g(\mu_{t},\mu_{t+1}) = \ell_1(\mu_{t} - \mu_{t+1}) - c \le 0$ |
| no constraint | $g(\mu_{t},\mu_{t+1}) = 0$ | 

We can now write with $\beta^{e(1)} = 0$,

$$ \begin{aligned}
        Q_n(\mathsf{G}) = &\min_{(\mu_1,...,\mu_n)\,\in\, \mathbb{R}^{n}} \sum_{i=1}^n\gamma(y_i,\mu_i) + \beta^{e(i)}\,,\\
        &\hbox{subject to} \,\, e(i) = (s_{i},s_{i+1},g_{e(i)},\beta^{e(i)}) \,\in\,  \mathcal{T}(i)\,,\, i = 1,...,n-1\,,\\
        &\quad\quad\quad\quad\quad(s_1,...,s_n) \in \mathcal{S}^n\,.\\
      \end{aligned}$$

In the `gfpop` package, the transition set $\mathcal{T}(t)$ is not time dependent. In this case, we define the graph $\mathcal{G}= (V,E)$  by the label of its vertices $V = \{0,...,S\} \subset \mathbb{N}$ and its set of oriented edges $E \subset V \times V \times \mathcal{B}$ with $e = (i,j,b) \in E$ defining an edge going from vertex $i$ to vertex $j$ with constraint option $b$. We get

$$Q_n(\mathcal{G}) = \min_{(e(1),...,e(n-1)) \,\in \, P} Q_n(e(1),...,e(n-1))\,,$$

with 

$$\begin{aligned}
       Q_n(e(1),...,e(n-1)) = &\min_{(\mu_1,...,\mu_n)\,\in\, \mathbb{R}^{n}} \sum_{i=1}^n\gamma(y_i,\mu_i) + \beta^{e(i)}\,,\\
        &\hbox{subject to} \,\, g_{e(1)}(\mu_1,\mu_2)\le 0,...,g_{e(n-1)}(\mu_{n-1},\mu_n)\le 0\,,\\
      \end{aligned}$$

and $P \subset \mathcal{T}^{n-1}$ being the set of valid paths of the graph. We say that $(e(1),...,e(n-1))$ is a valid path if $e(1) = (s_1,s_2,b_2) \in E$, $e(2) = (s_2,s_3,b_3) \in E$,..., $e(n-1) = (s_{n-1},s_n,b_n) \in E$. 

Notice that in most applications, the set of edges always contains edges of type $(i,i,g(\mu,\nu)=(\mu-\nu)^2,0)$ for all $i \in V$ as it corresponds to an absence of constraint on segment lengths. 

The main algorithm of the package, the `gfpop` function, returns the changepoint vector $\tau^*$ defined as $\{i \in \{1,...,n\}\,,\, m_{i}\ne m_{i+1}\}$, with $(m_1,...,m_n)$ the argminimum of $(\mu_1,...,\mu_n)$ and $m_{n+1} = +\infty$.

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
## [1]  110  304  498  799 1000
## 
## states
## [1] 0 1 0 1 0
## 
## forced
## [1] 0 0 0 0
## 
## means
## [1] 1.1078968 1.9624378 0.9154199 2.9230037 0.9634287
## 
## cost
## [1] 1008.885
## 
## attr(,"class")
## [1] "gfpop"
```

The vector `changepoints` gives the last index of each segment. It always ends with the length of the vector `vectData`.

The vector `states` contains the states in which lies each mean. The length of this vector is the same as the length of `changepoint`.

The vector `forced` is a boolean vector. A forced element means that two consecutive means have been forced to satisfy the constraint. For example, the "up" edge with parameter $c$ is forced if $m_{i+1} - m_i = c$.

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
## [1]  342  778 1000
## 
## states
## [1] 0 0 0
## 
## forced
## [1] 0 0
## 
## means
## [1] 0.500000 2.118050 3.089515
## 
## cost
## [1] 559.5639
## 
## attr(,"class")
## [1] "gfpop"
```

In this example, we use in `gfpop` function a robust biweight gaussian cost with `K = 1` and the `min` parameter in order to infer means greater than `0.5`.

### Fixed number of changepoints

This algorithm is called segment neighborhood in the changepoint litterature. In this example, we fixed the number of segments at $3$ with an isotonic constraint.

```r
n <- 1000
mydata <- dataGenerator(n, c(0.1, 0.2, 0.3, 0.4, 0.6, 0.8, 1), c(0, 0.5, 1, 1.5, 2, 2.5, 3), 1)
myGraph <- graph()
beta <- 0
myGraph <- addEdge(myGraph, edge(0, 1,"up", beta, 0))
myGraph <- addEdge(myGraph, edge(1, 2, "up", beta, 0))
myGraph <- addEdge(myGraph, edge(0, 0, "null", 0, 0))
myGraph <- addEdge(myGraph, edge(1, 1, "null", 0, 0))
myGraph <- addEdge(myGraph, edge(2, 2, "null", 0, 0))
myGraph <- addStartEnd(myGraph, start = 0, end = 2)

gfpop(vectData =  mydata, mygraph = myGraph, type = "gauss")
```

```
## changepoints
## [1]  215  599 1000
## 
## states
## [1] 0 1 2
## 
## forced
## [1] 0 0
## 
## means
## [1] 0.2627238 1.6544963 2.7962799
## 
## cost
## [1] 1175.894
## 
## attr(,"class")
## [1] "gfpop"
```

### Robust up-down with constrained starting and ending states

We use a robust loss (biweight) in presence of outliers. We can also force the starting and ending state and a minimal gap between the means (here equal to `1`)

```r
n <- 1000
mydata <- dataGenerator(n, c(0.1,0.3,0.5,0.8,1), c(0,1,0,1,0), 1) + 5*(rbinom(n, 1, 0.05)) - 5*(rbinom(n, 1, 0.05))

myGraph <- graph()
beta <- 2*log(n)
myGraph <- addEdge(myGraph, edge(0, 1, "up", beta, 1))
myGraph <- addEdge(myGraph, edge(1, 0, "down", beta, 1))
myGraph <- addEdge(myGraph, edge(0, 0, "null", 0, 0))
myGraph <- addEdge(myGraph, edge(1, 1, "null", 0, 0))
myGraph <- addStartEnd(myGraph, start = 0, end = 0)

gfpop(vectData =  mydata, mygraph = myGraph, type = "gauss", K = 3.0)
```

```
## changepoints
## [1]  108  303  501  801 1000
## 
## states
## [1] 0 1 0 1 0
## 
## forced
## [1] 0 0 0 0
## 
## means
## [1] -0.21300948  1.05732964  0.02170257  1.06765299 -0.03098637
## 
## cost
## [1] 1861.727
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
##  [1]   35   36   50   51   53   54   90   91   92  114  115  143  144  175
## [15]  176  205  206  210  211  217  218  219  220  239  240  267  268  281
## [29]  282  311  312  327  328  362  363  386  391  413  415  425  426  484
## [43]  485  508  509  521  522  524  526  549  551  587  588  617  618  629
## [57]  630  640  641  642  643  659  660  669  671  672  680  681  735  736
## [71]  809  810  827  828  875  876  877  928  929  955  956  961  962 1000
## 
## states
## integer(0)
## 
## forced
## integer(0)
## 
## means
##  [1] -0.29509651  6.14978739 -0.29257691 -6.18092662  0.09186331
##  [6] -6.91223263 -0.28468862 -4.10017658  5.61583421  0.33027204
## [11]  6.36411357  1.24884240 -5.56451864  0.88529301 -5.51867616
## [16]  1.61122876  7.68849704  1.54174283 -5.24650640  2.12902748
## [21] -5.77915378  1.32897784  7.24821231  0.64751660  6.48033277
## [26]  1.08337874 -4.07794736  1.44412114  6.14650103  0.81186246
## [31] -6.62753427 -0.35393244  5.08917536 -0.36195058 -5.74330002
## [36] -0.35324118  2.64265935 -0.11267510 -4.09348505  1.30770067
## [41] -4.83030005  0.16617456  5.68427733 -0.11478326  5.97633452
## [46]  0.77189078  6.73481539 -3.37350596  4.42770993  0.34826774
## [51] -3.93944021  1.43911168 -3.40505111  0.59344743  6.59593720
## [56]  0.83816338 -6.22437273  1.06074636 -4.56425142  2.08474558
## [61]  7.40572544  1.05797066  7.42808738  2.28196750 -1.07778118
## [66]  8.03498653  1.46615722  7.22957074  1.11017494 -5.54145065
## [71]  1.07972737 -5.15729922  0.55174328 -5.68680955 -0.02376106
## [76]  4.27489514 -4.94437411 -0.29078347  6.26653210  0.18239916
## [81] -5.41287533 -0.16586203  5.40762504  0.08107869
## 
## cost
## [1] 2950.915
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
myGraph <- addEdge(myGraph, edge(0, 0,"null", 0, 0))
gfpop(vectData =  mydata, mygraph = myGraph, type = "gauss", K = 3)
```

```
## changepoints
##  [1]  1001  2003  2998  3000  4000  4999  6000  6002  7000  8001  9000
## [12] 10000
## 
## states
##  [1] 0 0 0 0 0 0 0 0 0 0 0 0
## 
## forced
##  [1] 0 0 0 1 1 1 1 0 1 0 0
## 
## means
##  [1] -0.022133631  0.957062160  0.012464823  0.987525332  1.987525332
##  [6]  0.987525332  1.987525332  0.987525332 -0.003673743  0.996326257
## [11]  0.007073716  0.976821519
## 
## cost
## [1] 2678.014
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
myGraph <- addEdge(myGraph, edge(0, 0,"null", 0, 0))
gfpop(vectData =  mydata, mygraph = myGraph, type = "gauss", K = 3)
```

```
## changepoints
##  [1]  1001  2000  3001  3996  4999  6000  6999  8000  9002 10000
## 
## states
##  [1] 0 0 0 0 0 0 0 0 0 0
## 
## forced
## [1] 0 1 0 0 1 0 0 0 1
## 
## means
##  [1] -0.001545194  1.004876055  0.004876055  1.997319887  0.989647733
##  [6]  1.989647733 -0.010992911  1.022179317 -0.013045297  0.986954703
## 
## cost
## [1] 2633.055
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

`state1` is the starting vertex of an edge, `state2` its ending vertex. `type` is one of the available edge type ("null", up", "down", "std", "absInf", "absSup"). `penalty` is a nonnegative parameter: the additional cost $\beta$ to consider when we move within the graph using this edge. `parameter` is annother nonnegative parameter, a characteristics of the edge, depending of its type.

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
beta <- 2 * log(n)
myGraph <- addEdge(myGraph, edge(0, 0, "null", 0, 0))
myGraph <- addEdge(myGraph, edge(1, 1, "null", 0, 0))
myGraph <- addEdge(myGraph, edge(0, 1, "up", beta, 1))
myGraph <- addEdge(myGraph, edge(0, 0, "down", beta, 0))
myGraph <- addEdge(myGraph, edge(1, 0, "down", beta, 0))
myGraph <- addStartEnd(myGraph, start = 0, end = 0)
myGraph
```

```
##   state1 state2  type  penalty parameter
## 1      0      0  null  0.00000         0
## 2      1      1  null  0.00000         0
## 3      0      1    up 18.42068         1
## 4      0      0  down 18.42068         0
## 5      1      0  down 18.42068         0
## 6      0     NA start       NA        NA
## 7      0     NA   end       NA        NA
```

Some graphs are often used: they are defined by default in the `graph` function. To use these graphs, we specify a string `type` equal to "std", "isotonic", "updown" or "infsup".
For example,

```r
myGraphIso <- graph(penalty = 12, type = "isotonic")
myGraphIso
```

```
##   state1 state2 type penalty parameter
## 1      0      0 null       0         0
## 2      0      0   up      12         0
```

[Back to Top](#top)
