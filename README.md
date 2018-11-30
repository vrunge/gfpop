---
output:
  html_document: default
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
### November 30, 2018

> [Introduction](#intro)

> [Quick Start](#qs)

> [Some examples](#se)

<a id="intro"></a>

<img src="https://latex.codecogs.com/svg.latex?\Large&space;x=\frac{-b\pm\sqrt{b^2-4ac}}{2a}"/>

## Introduction

`gfpop` is an R package for penalized parametric changepoint detection using functional pruning dynamic programming. The successive means can be constrained using a graph structure with edge of type "up", "down", "std", "absInf" or "absSup". To each edge we can use an additional nonnegative parameter allowing us to force a minimal gap between two successive means. The user can also constraint the infered means to lie between some minimal and maximal values. Data is modelized by a quadratic cost with possible use of a robust loss, biweight and Huber. In a next version of this package, other parametric losses will be available (L1, Poisson, binomial). 

The package gfpop is designed to segment univariate data $y_{1:n} = \{y_1,...,y_n\}$ obeying to a graph structure on segment means. The changepoint vector
$\tau = (\tau_1 < \cdots < \tau_k) \in \mathbb{N}^k$ defines the segments $\{\tau_i,...,\tau_{i+1}-1\}$, $i = 0,...,k$ using notations $\tau_0 = 1$ and  $\tau_{k+1} = n+1$. With the set $S_n = \{\tau=(\tau_0,...,\tau_{k+1}) \in \mathbb{N}^{k+2} \,|\, 1 = \tau_0 < \tau_1 < \cdots < \tau_k < \tau_{k+1} = n+1\}\,,$ the nonconstrained minimal global cost is given by
\begin{equation}
\label{cost}
Q_n = \min_{\tau \in S_n}\left[ \sum_{i=0}^{k}\lbrace \mathcal{C}(y_{\tau_i:\tau_{i+1}-1}) + \beta \rbrace \right] \,,
\end{equation}
where $\beta > 0$ is a penalty parameter and $\mathcal{C}(y_{u:v})$ is the minimal cost over the segement $\{u,...,v\}$. The penalty $\beta$ is an additional cost where introducing a new segment. \\


The graph $\mathcal{G}$ is defined by its vertices $V = \{1,...,v\} \subset \mathbb{N}^*$ and edges $E = \{e(i,j)\}$ where $i,j \in V$. $\mathcal{G} = (V,E)$. The successive means in $Q_n$ can be constrainted to follow the graph structure.

<a id="qs"></a>

## Quick Start

we present a basic use of the main functions of the package. More details about optional arguments are given in later sections.

We install the package from Github:
```{r}
#devtools::install_github("vrunge/gfpop")
library(gfpop)
```

We simulate some univariate gaussian data ($n = 1000$ points) with changepoints $0.1, 0.3, 0.5, 0.8, 1$ and means $1, 2, 1, 3, 1$ with a variance equal to $1$.
```{r}
n <- 1000
myData <- dataGenerator(n, c(0.1,0.3,0.5,0.8,1), c(1,2,1,3,1), 1)
```

We then define the graph of constraints to use for the dynamic programming algorithm. A simple case is the up-down constraint with a penalty here equal to $2 \log(n)$.
```{r}
n <- 1000
myGraph <- graph(penalty = 2*log(n), "updown")
```

The gfpop function gives the result of the segmentation using $myData$ and $myGraph$ as parameters.
```{r}
gfpop(vectData = myData, mygraph = myGraph, type = "gauss")
```

The vector `changepoints` gives the first index of each segment. It always begins with the element $1$.

The vector `states` contains the states in which lies each mean. The length of this vector is the same as the length of changepoint.

The vector `forced` is a boolean vector. A forced element means that two consecutive means have been constrained by the parameter of the edge to be at least greater that the value of the parameter. This is the case for edges "up", "down" and "absSup". With "absInf" we constraint the absolute value of the difference of means to be less that the parameter.  

The vector `means` contains the means of the successive segments. 
 
The number `cost` is equal to $Q_n$, the overall cost of the segmented data. 


<a id="se"></a>

## Some examples

### Robust loss and constrained starting and ending states

We use a robust loss (biweight) in presence of outliers. We can also force the start and end state and a minimal gap between the means (here equal to $1$)
```{r}
n <- 1000
mydata <- dataGenerator(n, c(0.1,0.3,0.5,0.8,1), c(0,1,0,1,0), 1) + 5*(rbinom(n, 1, 0.05)) - 5*(rbinom(n, 1, 0.05))

myGraph <- graph()
beta <- 2*log(n)
myGraph <- addEdge(myGraph, edge(0, 1, "up", beta, 1))
myGraph <- addEdge(myGraph, edge(1, 0, "down", beta, 1))
myGraph <- addStartEnd(myGraph, start = 0, end = 0)
gfpop(vectData =  mydata, mygraph = myGraph, type = "gauss", K = 3.0)
```

If we skip all these constraints, the result is the following

```{r}
myGraphStd <- graph(penalty = 2*log(n), "std")
gfpop(vectData =  mydata, mygraph = myGraphStd, type = "gauss")
```



### absInf and absSup edges

```{r}
n <- 10000
myGraph <- graph()
beta <- 2*log(n)
myGraph <- addEdge(myGraph, edge(0, 0,"absInf", beta, 1))
mydata <- dataGenerator(n, c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1), c(0, 1, 0, 2, 1, 2, 0, 1, 0, 1), 0.5)
gfpop(vectData =  mydata, mygraph = myGraph, type = "gauss", K = 3)
```



```{r}
n <- 10000
myGraph <- graph()
beta <- 2*log(n)
myGraph <- addEdge(myGraph, edge(0, 0,"absSup", beta, 1))
mydata <- dataGenerator(n, c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1), c(0, 1, 0, 2, 1, 2, 0, 1, 0, 1), 0.5)
gfpop(vectData =  mydata, mygraph = myGraph, type = "gauss", K = 3)
```








[Back to Top](#top)
