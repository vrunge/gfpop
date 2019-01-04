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

A nonnegative internal parameter can thus be associated to an edge (in "up" and "down") or is mandatory ("absSup" and "absInf"). More details on graph construction are given in the last [section](#gc). The graph definition can also include a starting and/or an ending edge and each edge has its own <img src="/tex/8217ed3c32a785f0b5aad4055f432ad8.svg?invert_in_darkmode&sanitize=true" align=middle width=10.16555099999999pt height=22.831056599999986pt/> penalty (see further).

The User has the possibility to constraint the infered means to lie between some minimal and maximal values.

Data is modelized by a quadratic cost with possible use of a robust loss, biweight and Huber. In a next version of this package, other parametric costs will be available (L1, Poisson, Binomial). 

The package `gfpop` is designed to segment univariate data <img src="/tex/b704a970e46e8418f3ef56718438b122.svg?invert_in_darkmode&sanitize=true" align=middle width=126.38913869999998pt height=24.65753399999998pt/> obeying to a graph structure on segment means. The changepoint vector
<img src="/tex/72b0e057a8a58ba5a30cbae5e8731343.svg?invert_in_darkmode&sanitize=true" align=middle width=209.11494284999998pt height=27.91243950000002pt/> defines the segments <img src="/tex/97a268c17395aace06ce389334ba5322.svg?invert_in_darkmode&sanitize=true" align=middle width=115.02097364999997pt height=24.65753399999998pt/>, <img src="/tex/59680c09e0e2d723d0bcf2005047b028.svg?invert_in_darkmode&sanitize=true" align=middle width=73.18587374999998pt height=22.831056599999986pt/> with fixed <img src="/tex/370a29c873e3d269a6111aa219085d0b.svg?invert_in_darkmode&sanitize=true" align=middle width=44.697406049999984pt height=21.18721440000001pt/> and  <img src="/tex/afdb85da7c3e7b7c3d226050994dbf5f.svg?invert_in_darkmode&sanitize=true" align=middle width=63.70246739999999pt height=14.15524440000002pt/>. We define the set <img src="/tex/8ec36d117f5c48fb6fe6fb0e54d9df29.svg?invert_in_darkmode&sanitize=true" align=middle width=271.57361054999996pt height=27.91243950000002pt/> the nonconstrained minimal global cost is given by

<p align="center"><img src="/tex/2724e825a3f45c095a9569ff9e19081b.svg?invert_in_darkmode&sanitize=true" align=middle width=277.1488236pt height=49.315569599999996pt/></p>

where <img src="/tex/99751e94989c68f9be0f6aa442bc80d5.svg?invert_in_darkmode&sanitize=true" align=middle width=40.302373649999986pt height=22.831056599999986pt/> is a penalty parameter and <img src="/tex/a44ff4154fc3bc708e9e752a14051324.svg?invert_in_darkmode&sanitize=true" align=middle width=49.762892849999986pt height=24.65753399999998pt/> is the minimal cost over the segment <img src="/tex/add1478513cabbadcd5004323f01b74c.svg?invert_in_darkmode&sanitize=true" align=middle width=62.71697189999998pt height=24.65753399999998pt/>. The penalty <img src="/tex/8217ed3c32a785f0b5aad4055f432ad8.svg?invert_in_darkmode&sanitize=true" align=middle width=10.16555099999999pt height=22.831056599999986pt/> is an additional cost when introducing a new segment. The argmin of this quantity gives us a vector <img src="/tex/fa1931c6895840344a67aa4df9cf3f59.svg?invert_in_darkmode&sanitize=true" align=middle width=15.782028899999991pt height=22.63846199999998pt/> containing the last position of each segment (if we do not consider <img src="/tex/370a29c873e3d269a6111aa219085d0b.svg?invert_in_darkmode&sanitize=true" align=middle width=44.697406049999984pt height=21.18721440000001pt/>).

With any cost, we have

<p align="center"><img src="/tex/ca6a9827bab0d298f28a4c7a533691bc.svg?invert_in_darkmode&sanitize=true" align=middle width=547.45790385pt height=23.0593242pt/></p>

defining the infered mean of the i+1-th segment <img src="/tex/97a268c17395aace06ce389334ba5322.svg?invert_in_darkmode&sanitize=true" align=middle width=115.02097364999997pt height=24.65753399999998pt/>.

The graph <img src="/tex/68a463cbf8842017bbbab8ca879333c7.svg?invert_in_darkmode&sanitize=true" align=middle width=10.75346414999999pt height=22.465723500000017pt/> is defined by its vertices <img src="/tex/17723a15d7438a3eb1184f1afb7c5eb5.svg?invert_in_darkmode&sanitize=true" align=middle width=130.47537855pt height=24.65753399999998pt/> and edges <img src="/tex/ca2470eb23ef165907bddf445541668d.svg?invert_in_darkmode&sanitize=true" align=middle width=92.55732584999998pt height=24.65753399999998pt/> where <img src="/tex/73e18d026de1f20e9f779ef5f1e0eeb3.svg?invert_in_darkmode&sanitize=true" align=middle width=54.01270214999999pt height=22.465723500000017pt/>. Thus <img src="/tex/39af0bee918c59081f4ecf8c2a4cc8b8.svg?invert_in_darkmode&sanitize=true" align=middle width=76.34688434999998pt height=24.65753399999998pt/>. The successive means <img src="/tex/ba97fea5b98d1b95ef2bfdfbba39044b.svg?invert_in_darkmode&sanitize=true" align=middle width=102.06838785pt height=24.65753399999998pt/> are constrainted to follow a feasible graph path <img src="/tex/f88ef5828577ac5b460c213ad9f91505.svg?invert_in_darkmode&sanitize=true" align=middle width=106.62095399999998pt height=24.65753399999998pt/> with <img src="/tex/f656003fdf20f1bf1dc44afea451ff83.svg?invert_in_darkmode&sanitize=true" align=middle width=46.30026884999999pt height=22.465723500000017pt/>. Among all possible paths <img src="/tex/df5a289587a2f0247a5b97c1e8ac58ca.svg?invert_in_darkmode&sanitize=true" align=middle width=12.83677559999999pt height=22.465723500000017pt/>, the path minimizing the cost is the result <img src="/tex/e96673f14f7d2280cfd060ac442da0b9.svg?invert_in_darkmode&sanitize=true" align=middle width=45.482259899999995pt height=24.65753399999998pt/> of our algorithm, that is, 

<p align="center"><img src="/tex/9653b6cc1b5bb51607ddfb3a1f13589b.svg?invert_in_darkmode&sanitize=true" align=middle width=206.416914pt height=17.46572685pt/></p>

and for a given path <img src="/tex/f88ef5828577ac5b460c213ad9f91505.svg?invert_in_darkmode&sanitize=true" align=middle width=106.62095399999998pt height=24.65753399999998pt/> we definie the path-constrained cost

<p align="center"><img src="/tex/f5104b69464610c67cf95e8ccccc404d.svg?invert_in_darkmode&sanitize=true" align=middle width=502.2140508pt height=49.315569599999996pt/></p>

with <img src="/tex/a7b70cf2530bd437854a7672a5d27827.svg?invert_in_darkmode&sanitize=true" align=middle width=95.69340989999999pt height=24.65753399999998pt/> meaning that the two consecutive means have to satisfy the edge constraint. <img src="/tex/14eb94e11e4e7df75df5246a34c67891.svg?invert_in_darkmode&sanitize=true" align=middle width=19.92040709999999pt height=22.831056599999986pt/> is the penalty associated to edge <img src="/tex/b95c2b0aab2482e5bebd25332a4bbde0.svg?invert_in_darkmode&sanitize=true" align=middle width=12.30503669999999pt height=14.15524440000002pt/>. In many cases, we simply take <img src="/tex/c13ff1e7a1e99be22dd58e6f7ea398e7.svg?invert_in_darkmode&sanitize=true" align=middle width=53.647382249999986pt height=22.831056599999986pt/> for all <img src="/tex/77a3b857d53fb44e33b53e4c8b68351a.svg?invert_in_darkmode&sanitize=true" align=middle width=5.663225699999989pt height=21.68300969999999pt/>. The penalty <img src="/tex/8217ed3c32a785f0b5aad4055f432ad8.svg?invert_in_darkmode&sanitize=true" align=middle width=10.16555099999999pt height=22.831056599999986pt/> was a constant positive cost we have to pay when adding a new segment, that is when we move on an edge. Thus, we define this penalty within the edge and enable this penalty to be different for each edge.

For each path, the differences <img src="/tex/d3a11478349f31da44cf7954e1f4f367.svg?invert_in_darkmode&sanitize=true" align=middle width=131.24719905pt height=22.465723500000017pt/> can be or not be constrained depending of the nature of the edge. For example, with an "up" edge with parameter <img src="/tex/b55131511b15b0ca60eee9efe6f939fc.svg?invert_in_darkmode&sanitize=true" align=middle width=35.36518424999999pt height=22.831056599999986pt/>, <img src="/tex/7221766392c47dfeaa41ac8d8f29b35c.svg?invert_in_darkmode&sanitize=true" align=middle width=60.75055139999999pt height=22.831056599999986pt/>,  with an "absInf" edge with parameter <img src="/tex/b55131511b15b0ca60eee9efe6f939fc.svg?invert_in_darkmode&sanitize=true" align=middle width=35.36518424999999pt height=22.831056599999986pt/>, <img src="/tex/d17a3fd6963758b83a8744679af4ddc0.svg?invert_in_darkmode&sanitize=true" align=middle width=69.88299779999998pt height=24.65753399999998pt/>, etc... 

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

The vector `forced` is a boolean vector. A forced element means that two consecutive means have been forced to satisfy the constraint. For example, the "up" edge with parameter <img src="/tex/2f2322dff5bde89c37bcae4116fe20a8.svg?invert_in_darkmode&sanitize=true" align=middle width=5.2283516999999895pt height=22.831056599999986pt/> is forced if <img src="/tex/40da0c26af3ff48ec5e1122cbabcaadf.svg?invert_in_darkmode&sanitize=true" align=middle width=103.69287719999998pt height=22.831056599999986pt/>.

The vector `means` contains the infered means of the successive segments. 
 
The number `cost` is equal to <img src="/tex/e96673f14f7d2280cfd060ac442da0b9.svg?invert_in_darkmode&sanitize=true" align=middle width=45.482259899999995pt height=24.65753399999998pt/>, the overall cost of the segmented data. 


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

`state1` is the starting vertex of an edge, `state2` its ending vertex. `type` is one of the available edge type ("up", "down", "std", "absInf", "absSup"). `penalty` is a nonnegative parameter: the additional cost <img src="/tex/3d13090ef3ed1448f3c4dc166d06ab4d.svg?invert_in_darkmode&sanitize=true" align=middle width=13.948864049999989pt height=22.831056599999986pt/> to consider when we move within the graph using this edge. `parameter` is annother nonnegative parameter, a characteristics of the edge, depending of its type.

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
