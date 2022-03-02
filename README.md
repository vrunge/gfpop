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

> [The Graph-constrained Change-point problem](#intro)

> [Quick Start](#qs)

> [Some examples](#se)

> [Graph construction](#gc)

> [Supplementary R functions](#suppl)

> [More on the gfpop function and its C++ structure](#gfpop)


<a id="intro"></a>

## The Graph-constrained Change-point problem

### General description

`gfpop` is an `R` package developed to complete parametric change-point detection in univariate time series constrained to a graph structure. Constraints are imposed to the sequence of inferred parameters (most of the time a mean parameter) of consecutive segments and related to the direction and/or the magnitude of the mean changes. 

Change-point detection is performed using the functional pruning optimal partitioning method (fpop) based on an exact dynamic programming algorithm. [Cf paper. On optimal multiple changepoint algorithms for large data](https://link.springer.com/article/10.1007/s11222-016-9636-3). This algorithmic strategy can be seen as an extented Viterbi algorithm for a Hidden Markov Model with continuous parameter states. A presentation of the generalized fpop (gfpop) algorithm is available in the paper [A log-linear time algorithm for constrained changepoint detection.](https://arxiv.org/pdf/1703.03352.pdf)

In the main `gfpop` function the user chooses a global parametric model for the change-point problem to solve (change in "mean", "variance", "exp", "poisson" or "negbin" distribution in `type` variable) and a `graph` structure. To run the function the user also gives some `data` to segment potentially associated with a `weights` vector of same length.

```r
gfpop <- function(data, mygraph, type = "mean", weights = NULL)
```

The algorithm of our package is designed to consider a large variety of constraints. These constraints are modelled by a graph. At each time $t$ a number of states is accessible, these states are the nodes of the graph. Possible transitions between states at time $t$ and $t+1$ are represented by the edges of the graph. Each edge is associated to $4$ main elements: a **constraint** (for example $\mu_t \le \mu_{t+1}$), some **parameters** associated to the constraint, a **penalty** (possibly null) and **robust parameters** to use Huber or biweight losses only on the considered edge.

In addition to the edge definition, the nodes can be constrained to a range of values (means) in the inference. We also can specify starting and ending nodes. In our implementation the transitions are not time-dependent except for starting and ending steps.

The edges of the graph can be of type "null", "std",  "up", "down" or "abs" with the following meaning:

- "null" edge : there is no change-point introduced. We stay on the same segment. "null" corresponds to the constraint $I_{\mu_t = \alpha\mu_{t+1}}$. The value does not change (or we have an exponential decay if $0 < \alpha < 1$);
- "std" edge : no contraint, the next segment can have any mean;
- "up" edge : the next segment has a greater mean (we can also force the size of the gap to be greater than a minimal value). "up" corresponds to the constraint $I_{\mu_t \leq \mu_{t+1}}$;
- "down" edge : the next segment has a lower mean (wan also can force the size of the gap to be greater that a minimal value). "down" corresponds to the constraint $I_{\mu_t \geq \mu_{t+1}}$;
- "abs" edge : the absolute value of the difference of two consecutive means is greater than a given parameter.


A nonnegative internal parameter can thus be associated to an edge (in "up", "down" and "abs""). The "null" edge refers to an absence of change-point. This edge can be used between different states to constraint segment lengths. Thus our package includes the segment neighborhood algorithm for which the number of change-points is fixed. More details on graph construction are given in the section [Graph construction](#gc).

### The non-constrained changepoint detection problem

The change-point vector $\overline{\tau} = (\tau_0 < \cdots < \tau_{k+1}) \in \mathbb{N}^{k+2}$ defines the $k+1$ segments $\{\tau_i+1,...,\tau_{i+1}\}$, $i = 0,...,k$ with fixed bounds $\tau_0 = 0$ and  $\tau_{k+1} = n$. We use the set $S_n = \{\hbox{change-point vector } \overline{\tau} \in \mathbb{N}^{k+2}\}$ to define the non-constrained minimal global cost given by

<img src="https://render.githubusercontent.com/render/math?math=Q_n = \min_{\overline{\tau} \in S_n}\left[ \sum_{i=0}^{k}\lbrace \mathcal{C}(y_{(\tau_i+1):\tau_{i+1}}) + \beta \rbrace \right]\,,">"

where $\beta > 0$ is a penalty parameter and $\mathcal{C}(y_{u:v})$ is the minimal cost over the segment $\{u,...,v\}$. The penalty $\beta$ is understood as an additional cost when introducing a new segment. The argminimum of this quantity gives a vector $\tau^*$ containing the last position of each segment (if we do not consider $\tau_0 = 0$). The quantity $Q_n$ is the solution of the non-constrained optimal partitioning method. 

In our setting, the cost $\mathcal{C}$ is the result of the minimization of a cost function with additive property:

$$ \left\{
\begin{array}{l}
  \mathcal{C}(y_{(\tau_i+1):\tau_{i+1}}) = \min_{\mu_i}\mathcal{C}(y_{(\tau_i+1):\tau_{i+1}}, \mu_i) = \min_{\mu_i} \sum_{j = \tau_i+1}^{\tau_{i+1}}\gamma(y_j, \mu_i)\,,\\
  m_i = \mathrm{argmin}_{\mu_i}\mathcal{C}(y_{(\tau_i+1):\tau_{i+1}}, \mu_i)\,,
\end{array}
\right.$$

with the argminimum defining the infered mean $m_i$ of the i+1-th segment $\{\tau_i+1,...,\tau_{i+1}\}$ with $i \in \{0,...,k\}$. Additivity of the cost (the $\gamma$ decomposition) is guaranteed as we will use costs deriving from a likelihood. $\gamma$ has in our package $3$ possible basis decompositions:

\begin{description}
\item[Gauss decomposition.] $f_1 : \mu \mapsto 1$, $f_2 : \mu \mapsto \mu$ and $f_3 : \mu \mapsto \mu^2$. This decomposition allows to consider $\ell_2$, biweight and Huber loss functions; 
\item[Poisson decomposition.] $f_1 : \mu \mapsto 1$, $f_2 : \mu \mapsto \mu$ and $f_3 : \mu \mapsto \log(\mu)$. This decomposition allows to consider the Poisson, Exponential, as well as the fixed mean and variable variance Normal likelihoods;

\item[Binomial decomposition.] $f_1 : \mu \mapsto 1$, $f_2 : \mu \mapsto \log(\mu)$ and $f_3 : \mu \mapsto \log(1-\mu)$. This decomposition allows to consider the Binomial and Negative binomial likelihoods.
\end{description}

We can associate a current mean $\mu_i$ to each data point $y_i$ and we write a cost on a segment $\{u,...,v\}$ as a result of a constrained minimization:

$$\mathcal{C}(y_{u:v}) = \min_{\overset{(\mu_u,...,\mu_v)\in \mathbb{R}^{{v-u+1}}}{\mu_u = \cdots = \mu_v}}\sum_{j = u}^{v}\gamma(y_j, \mu_j)\,,$$

so that we get another description of the objective function :

$$Q_n = \min_{\overset{(\mu_1,...,\mu_n)\in \mathbb{R}^{n}}{\mu_{n+1} = +\infty}}\quad\sum_{i=1}^n\gamma(y_i,\mu_i) + \beta I(\mu_{i} \ne \mu_{i+1})\,,$$

From this expression, it will be possible to impose some constraints on the vector of means $\mu$.

### Graph-constrained problem

$I \in \{0,1\}$ is the indicator function. This previous approach can be generalized to complex constraints on consecutive means using a graph structure. We define the transition graph $\mathcal{G}_n$ as a directed acyclic graph with the following properties:

1. Nodes are indexed by time and state. $v = (t,s)$ for node $v$ with time $t$ and state $s$. The states are elements of the set $\mathcal{S} =  \{0,...,S\} \subset \mathbb{N}$;

2. All the nodes have time $t$ in $\{1,...,n\}$ except for the unique starting node $v_0= (0,\emptyset)$ and the unique ending node $v_{n+1} = (n+1,\emptyset)$, where $\emptyset$ denotes an undefinite state;

3. Edges are transitions between "time-consecutive" vertices of type $v = (t,s)$ and $v' = (t+1,s')$ which gives us the edge $e = (t,s,s')$ for $t \in \{0,...,n\}$;

4. Each edge $e$ is associated to a function $I_e : \mathbb{R}\times \mathbb{R} \mapsto \{0,1\}$ constraining consecutive means $\mu_t$ and $\mu_{t+1}$ by the indicator function $(\mu_t, \mu_{t+1}) \mapsto I_e(\mu_t, \mu_{t+1})$ (for example $I_e(\mu_t, \mu_{t+1}) = I(\mu_t \le \mu_{t+1})$) and a penalty $\beta_e \ge 0$.

The next table summarizes all the possible constraints encoded into the `gfpop` package.

| constraints | $I_e : \mathbb{R}\times \mathbb{R} \mapsto \mathbb{R}$, $c \in \mathbb{R}^+$ |
|---------:|-----:|
| null | $I_e(\mu_{t},\mu_{t+1}) =  I_{\mu_t = \mu_{t+1}}$ |
| std | $I_e(\mu_{t},\mu_{t+1}) =  I_{\mu_t \ne \mu_{t+1}}$ | 
| up | $I_e(\mu_{t},\mu_{t+1}) = I_{\mu_{t}  + c \le \mu_{t+1}}$ |
| down | $I_e(\mu_{t},\mu_{t+1}) = I_{\mu_{t+1} + c \le \mu_{t}}$ |
| abs | $I_e(\mu_{t},\mu_{t+1}) = I_{c \le \ell_1(\mu_{t} - \mu_{t+1})}$|


We define a path $p \in \mathcal{G}_n$ of the graph as a collection of $n+2$ vertices $(v_0,...,v_{n+1})$ with $v_0 = (0,\emptyset)$ and $v_{n+1} = (n+1,\emptyset)$ and $v_t = (t,s_t)$ for $t \in \{1,...,n\}$ and $s_t \in \mathcal{S}$. Morever, the path is made of $n+1$ edges denoted $e_0,...,e_n$ with $\beta_{e_n} = 0$. A vector $\mu \in \mathbb{R}^n$ verifies the path $p$ if for all $t \in \{1,...,n-1\}$, we have $I_{e_t}(\mu_t,\mu_{t+1}) = 1$ (valid constraint). We write $p(\mu)$ to say that the vector $\mu$ verifies the path $p$. The formulation of our graph-constrained problem is then the following:

$$
Q_n(\mathcal{G}_n)  = \underset{\underset{\mu \in \mathbb{R}^n\,,\,p(\mu)}{p \in \mathcal{G}_n} }{\min} 
\left\{\sum_{t=1}^n (\gamma(y_t, \mu_t) + \beta_{e_t}) \right\}\,.
$$


In the `gfpop` package, the edges $(t,s,s')$ for $t\in\{1,...,n\}$ are not time-dependent. In this case, we redefine a graph $\mathcal{G}= (V,E)$  by the label of its vertices $V = \{0,...,S\} \subset \mathbb{N}$ and its set of oriented edges $E \subset V \times V$ with $e = (s,s') \in E$ defining an edge going from node $s$ to node $s'$ with indicator function $I_e$ and penalty $\beta_e$. 

In most applications, the set of edges always contains edges of type $(s,s)$ for all $s \in V$ with indicator function $I_e(\mu,\nu)= I(\mu \ne \nu)$ as it corresponds to an absence of constraint on segment lengths. 

The main algorithm of the package, the `gfpop` function, returns the change-point vector $\tau^*$ defined as $\{i \in \{1,...,n\}\,,\, m_{i}\ne m_{i+1}\}$, with $(m_1,...,m_n)$ the argminimum of $(\mu_1,...,\mu_n)$ in $Q_n(\mathcal{G}_n)$ and $m_{n+1} = +\infty$. we also give the possibility to restrict the set of valid paths by imposing a starting and/or an ending ,node and contraint the range of infered means, replacing $\mu \in \mathbb{R}^n$ by $\mu \in [A,B]$ in the definition of $Q_n(\mathcal{G}_n)$.

<a id="qs"></a>

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
## changepoints
## [1]  100  299  500  801 1000
## 
## states
## [1] "Dw" "Up" "Dw" "Up" "Dw"
## 
## forced
## [1] FALSE FALSE FALSE FALSE
## 
## parameters
## [1] 0.9767950 1.8985302 0.9593281 3.0247773 1.0370582
## 
## globalCost
## [1] 958.077
## 
## attr(,"class")
## [1] "gfpop" "mean"
```

The vector `changepoints` gives the last index of each segment. It always ends with the length of the vector `vectData`.

The vector `states` contains the states in which lies each mean. The length of this vector is the same as the length of `changepoint`.

The vector `forced` is a boolean vector. A forced element means that two consecutive means have been forced to satisfy the constraint. For example, the "up" edge with parameter $c$ is forced if $m_{i+1} - m_i = c$.

The vector `parameters` contains the inferred means/parameters of the successive segments. 
 
The number `globalCost` is equal to $Q_n(\mathcal{G})$, the overall cost of the segmented data. 


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
## changepoints
## [1]   58  235  400  762 1000
## 
## states
## [1] "Iso" "Iso" "Iso" "Iso" "Iso"
## 
## forced
## [1] FALSE FALSE FALSE FALSE
## 
## parameters
## [1] -0.1822423  0.3884548  1.2961028  2.2785127  2.8569659
## 
## globalCost
## [1] 996.0674
## 
## attr(,"class")
## [1] "gfpop" "mean"
```

In this example, we use in `gfpop` function a robust biweight gaussian cost with `K = 1` and the `min` parameter in order to infer means greater than `0.5`.

### Fixed number of change-points

This algorithm is called segment neighborhood in the change-point litterature. In this example, we fixed the number of segments at $3$ with an isotonic constraint. The graph contains two "up" edges with no cycling.



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
## changepoints
## [1]  301  609 1000
## 
## states
## [1] "0" "1" "2"
## 
## forced
## [1] FALSE FALSE
## 
## parameters
## [1] 0.5214622 1.7950420 2.7039511
## 
## globalCost
## [1] 1087.774
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
## changepoints
## [1]  100  297  493  798 1000
## 
## states
## [1] "Dw" "Up" "Dw" "Up" "Dw"
## 
## forced
## [1] FALSE FALSE FALSE  TRUE
## 
## parameters
## [1] -0.24757635  1.18202359 -0.02892791  1.04627242  0.04627242
## 
## globalCost
## [1] 1076.384
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
## changepoints
##  [1]   24   26   73   74   83   84  101  102  115  119  129  147  148  171  172
## [16]  216  219  220  238  239  288  289  297  305  306  307  334  335  340  341
## [31]  393  394  396  397  398  419  420  441  442  496  497  541  543  544  585
## [46]  586  603  604  628  629  663  664  696  697  706  707  718  719  747  748
## [61]  756  761  772  773  792  793  814  815  844  845  848  849  867  869  874
## [76]  875  888  889  950  952  955 1000
## 
## states
##  [1] "Std" "Std" "Std" "Std" "Std" "Std" "Std" "Std" "Std" "Std" "Std" "Std" "Std"
## [14] "Std" "Std" "Std" "Std" "Std" "Std" "Std" "Std" "Std" "Std" "Std" "Std" "Std"
## [27] "Std" "Std" "Std" "Std" "Std" "Std" "Std" "Std" "Std" "Std" "Std" "Std" "Std"
## [40] "Std" "Std" "Std" "Std" "Std" "Std" "Std" "Std" "Std" "Std" "Std" "Std" "Std"
## [53] "Std" "Std" "Std" "Std" "Std" "Std" "Std" "Std" "Std" "Std" "Std" "Std" "Std"
## [66] "Std" "Std" "Std" "Std" "Std" "Std" "Std" "Std" "Std" "Std" "Std" "Std" "Std"
## [79] "Std" "Std" "Std" "Std"
## 
## forced
##  [1] FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE
## [14] FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE
## [27] FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE
## [40] FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE
## [53] FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE
## [66] FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE
## [79] FALSE FALSE FALSE
## 
## parameters
##  [1]  0.09531845 -3.85703556 -0.13019394  5.16267205  0.44447672  5.20321468
##  [7] -0.39426869 -5.08788196  1.17547781 -1.61382166  2.51844573  0.66099190
## [13] -6.13502552  1.16677323 -5.25144245  1.19790595 -1.10033854  7.66759284
## [19]  1.47574267 -3.98368417  1.05260128 -5.38083210  1.73201851 -0.26765007
## [25] -5.61166029  6.06497166 -0.05299352  4.45206610 -1.49900686 -6.24316215
## [31]  0.15086007  5.61934257  0.52500592 -4.51976413  5.24207543 -0.54093753
## [37] -5.42663094  0.03931173 -6.71419815  0.13592288  5.56618535  0.66407336
## [43]  4.20858670 -3.97925015  1.24780049  6.59153592  0.77567253  7.12261221
## [49]  1.46381318  6.31997462  0.79516336  6.14695385  0.55942101 -4.40451015
## [55]  1.21243384 -5.01827488  1.94092687  6.15208912  0.71233302  6.82105479
## [61]  1.18033057  4.03182114  0.35281200  7.50141174  1.54135300 -4.28366480
## [67]  0.33847303 -5.35427462 -0.22562237  5.23822258  0.56197995 -5.59772017
## [73]  0.36645125 -4.80571764  1.81510038  7.04674283  0.60538203  6.12114295
## [79] -0.17145432  5.19794411 -3.11723378  0.37162072
## 
## globalCost
## [1] 1824.284
## 
## attr(,"class")
## [1] "gfpop" "mean"
```



### abs edge

With a unique "abs" edge, we impose a difference between the means of size at least 1.  


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
## changepoints
##  [1]  1000  2000  3000  3999  5000  6000  6999  8000  8998 10000
## 
## states
##  [1] "0" "0" "0" "0" "0" "0" "0" "0" "0" "0"
## 
## forced
## [1] FALSE  TRUE FALSE FALSE FALSE FALSE FALSE FALSE  TRUE
## 
## parameters
##  [1] -0.008507931  0.993648981 -0.006351019  1.998147387  0.978719714  2.017369226
##  [7]  0.005955304  1.006008367 -0.005573578  0.994426422
## 
## globalCost
## [1] 2492.677
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
## changepoints
## [1]  200  500  800 1000
## 
## states
## [1] "0" "0" "0" "0"
## 
## forced
## [1] FALSE FALSE FALSE
## 
## parameters
## [1] 0.0049237269 0.0003130684 0.0004639405 0.0201732198
## 
## globalCost
## [1] 981.2278
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

![plot of chunk unnamed-chunk-15](figure/unnamed-chunk-15-1.png)

![plot](figure/expDecay.png)


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

`state1` is the starting node of an edge, `state2` its ending node. `type` is one of the available edge type ("null", "std", "up", "down", "abs"). `penalty` is a nonnegative parameter: the additional cost $\beta_i$ to consider when we move within the graph using a edge (or stay on the same node). `parameter` is annother nonnegative parameter, a characteristics of the edge, depending of its type (it is a decay if type is "null" and a gap otherwise). `K` and `a` are robust parameters. `min` and `max` are used to constrain the rang of value for the node parameter.

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

Some graphs are often used: they are defined by default in the `graph` function. To use these graphs, we specify a string `type` equal to "std", "isotonic", "updown" or "relevant".
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

the `dataGenerator` function is used to simulate $n$ data-points from a distribution of `type` equal to “mean”, “poisson”, “exp”, “variance” or “negbin”. Standard deviation parameter `sigma` and decay `gamma` are specific to the Gaussian mean model. `size` is linked to the R “rnbinom” function from R stats package.

### Standard deviation estimation

We often need to estimate the standard deviation from the observed data to normalize the data or choose the edge penalties. The `sdDiff` returns such an estimation with the default HALL method [Hall et al., 1990] well suited for time series with change-points.

<a id="gfpop"></a>

## More on the gfpop function and its C++ structure

To be done.

[Back to Top](#top)
