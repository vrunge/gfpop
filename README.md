<a id="top"></a>

# gfpop Vignette
### Vincent Runge
#### LaMME, Evry University
### November 30, 2018

> [Introduction](#intro)

> [Quick Start](#qs)

> [Some examples](#se)

<a id="intro"></a>

## Introduction

`gfpop` is an R package for penalized parametric changepoint detection using functional pruning dynamic programming. The successive means can be constrained using a graph structure with edge of type "up", "down", "std", "absInf" or "absSup". To each edge we can associate an additional nonnegative parameter allowing us to force a minimal gap length between two successive means. The user can also constraint the infered means to lie between some minimal and maximal values. Data is modelized by a quadratic cost with possible use of a robust loss, biweight and Huber. In a next version of this package, other parametric losses will be available (L1, Poisson, binomial). 

The package gfpop is designed to segment univariate data <a href="https://www.codecogs.com/eqnedit.php?latex=\inline&space;y_{1:n}&space;=&space;\{y_1,...,y_n\}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\inline&space;y_{1:n}&space;=&space;\{y_1,...,y_n\}"/></a> obeying to a graph structure on segment means. The changepoint vector <a href="https://www.codecogs.com/eqnedit.php?latex=\inline&space;\tau&space;=&space;(\tau_1&space;<&space;\cdots&space;<&space;\tau_k)&space;\in&space;\mathbb{N}^k" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\inline&space;\tau&space;=&space;(\tau_1&space;<&space;\cdots&space;<&space;\tau_k)&space;\in&space;\mathbb{N}^k" /></a> defines the segments <a href="https://www.codecogs.com/eqnedit.php?latex=\inline&space;\{\tau_i,...,\tau_{i&plus;1}-1\}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\inline&space;\{\tau_i,...,\tau_{i&plus;1}-1\}" /></a>, <a href="https://www.codecogs.com/eqnedit.php?latex=\inline&space;i&space;=&space;0,...,k" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\inline&space;i&space;=&space;0,...,k"/></a> using notations <a href="https://www.codecogs.com/eqnedit.php?latex=\inline&space;\tau_0&space;=&space;1" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\inline&space;\tau_0&space;=&space;1" /></a> and <a href="https://www.codecogs.com/eqnedit.php?latex=\inline&space;\tau_{k&plus;1}&space;=&space;n&plus;1" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\inline&space;\tau_{k&plus;1}&space;=&space;n&plus;1"/></a>. With the set 

<a href="https://www.codecogs.com/eqnedit.php?latex=S_n&space;=&space;\{\tau=(\tau_0,...,\tau_{k&plus;1})&space;\in&space;\mathbb{N}^{k&plus;2}&space;\,|\,&space;1&space;=&space;\tau_0&space;<&space;\tau_1&space;<&space;\cdots&space;<&space;\tau_k&space;<&space;\tau_{k&plus;1}&space;=&space;n&plus;1\}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?S_n&space;=&space;\{\tau=(\tau_0,...,\tau_{k&plus;1})&space;\in&space;\mathbb{N}^{k&plus;2}&space;\,|\,&space;1&space;=&space;\tau_0&space;<&space;\tau_1&space;<&space;\cdots&space;<&space;\tau_k&space;<&space;\tau_{k&plus;1}&space;=&space;n&plus;1\}"/></a>

the nonconstrained minimal global cost is given by

<a href="https://www.codecogs.com/eqnedit.php?latex=Q_n&space;=&space;\min_{\tau&space;\in&space;S_n}\left[&space;\sum_{i=0}^{k}\lbrace&space;\mathcal{C}(y_{\tau_i:\tau_{i&plus;1}-1})&space;&plus;&space;\beta&space;\rbrace&space;\right]" target="_blank"><img src="https://latex.codecogs.com/gif.latex?Q_n&space;=&space;\min_{\tau&space;\in&space;S_n}\left[&space;\sum_{i=0}^{k}\lbrace&space;\mathcal{C}(y_{\tau_i:\tau_{i&plus;1}-1})&space;&plus;&space;\beta&space;\rbrace&space;\right]"/></a>

where <a href="https://www.codecogs.com/eqnedit.php?latex=\beta&space;>&space;0" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\beta&space;>&space;0" title="\beta > 0" /></a> is a penalty parameter and <a href="https://www.codecogs.com/eqnedit.php?latex=\inline&space;\mathcal{C}(y_{u:v})" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\inline&space;\mathcal{C}(y_{u:v})"/></a> is the minimal cost over the segement <a href="https://www.codecogs.com/eqnedit.php?latex=\inline&space;\{u,...,v\}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\inline&space;\{u,...,v\}"/></a>. The penalty <a href="https://www.codecogs.com/eqnedit.php?latex=\inline&space;\beta" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\inline&space;\beta" /></a> is an additional cost where introducing a new segment.


The graph <a href="https://www.codecogs.com/eqnedit.php?latex=\inline&space;\mathcal{G}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\inline&space;\mathcal{G}"/></a> is defined by its vertices <a href="https://www.codecogs.com/eqnedit.php?latex=\inline&space;V&space;=&space;\{1,...,v\}&space;\subset&space;\mathbb{N}^*" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\inline&space;V&space;=&space;\{1,...,v\}&space;\subset&space;\mathbb{N}^*"/></a> and edges <a href="https://www.codecogs.com/eqnedit.php?latex=\inline&space;E&space;=&space;\{e(i,j)\}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\inline&space;E&space;=&space;\{e(i,j)\}"/></a> where <a href="https://www.codecogs.com/eqnedit.php?latex=\inline&space;i,j&space;\in&space;V" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\inline&space;i,j&space;\in&space;V"/></a>. <a href="https://www.codecogs.com/eqnedit.php?latex=\inline&space;\mathcal{G}&space;=&space;(V,E)" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\inline&space;\mathcal{G}&space;=&space;(V,E)"/></a>. The successive means in <a href="https://www.codecogs.com/eqnedit.php?latex=\inline&space;Q_n" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\inline&space;Q_n" /></a> can be constrainted to follow the graph structure.

<a id="qs"></a>

## Quick Start

we present a basic use of the main functions of the package. More details about optional arguments are given in later sections.

We install the package from Github:
```{r}
#devtools::install_github("vrunge/gfpop")
library(gfpop)
```

We simulate some univariate gaussian data (n = 1000 points) with changepoints 0.1, 0.3, 0.5, 0.8, 1 and means 1, 2, 1, 3, 1 with a variance equal to 1.
```{r}
n <- 1000
myData <- dataGenerator(n, c(0.1,0.3,0.5,0.8,1), c(1,2,1,3,1), 1)
```

We then define the graph of constraints to use for the dynamic programming algorithm. A simple case is the up-down constraint with a penalty here equal to <a href="https://www.codecogs.com/eqnedit.php?latex=\inline&space;2&space;\log(n)" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\inline&space;2&space;\log(n)" title="2 \log(n)" /></a>.
```{r}
n <- 1000
myGraph <- graph(penalty = 2*log(n), "updown")
```

The gfpop function gives the result of the segmentation using myData and myGraph as parameters.
```{r}
gfpop(vectData = myData, mygraph = myGraph, type = "gauss")
```

The vector `changepoints` gives the first index of each segment. It always begins with the element 1.

The vector `states` contains the states in which lies each mean. The length of this vector is the same as the length of changepoints.

The vector `forced` is a boolean vector. A forced element means that two consecutive means have been constrained by the parameter of the edge to be at least greater that the value of the parameter. This is the case for edges "up", "down" and "absSup". With "absInf" we constraint the absolute value of the difference of means to be less that the parameter.  

The vector `means` contains the means of the successive segments. 
 
The number `cost` is equal to <a href="https://www.codecogs.com/eqnedit.php?latex=\inline&space;Q_n" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\inline&space;Q_n" title="Q_n" /></a>, the overall cost of the segmented data. 


<a id="se"></a>

## Some examples

### Robust loss and constrained starting and ending states

We use a robust loss (biweight) in presence of outliers. We can also force the start and end state and a minimal gap between the means (here equal to 1)
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
