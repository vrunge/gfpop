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

The package `gfpop` is designed to segment univariate data <img src="/tex/b704a970e46e8418f3ef56718438b122.svg?invert_in_darkmode&sanitize=true" align=middle width=126.38913869999998pt height=24.65753399999998pt/> obeying to a graph structure on segment means. The changepoint vector <img src="/tex/ed535ca37308856c5f0c6d85b4d4c676.svg?invert_in_darkmode&sanitize=true" align=middle width=209.11492305pt height=27.91243950000002pt/> defines the <img src="/tex/33359de825e43daa97171e27f6558ae9.svg?invert_in_darkmode&sanitize=true" align=middle width=37.38576269999999pt height=22.831056599999986pt/> segments <img src="/tex/97a268c17395aace06ce389334ba5322.svg?invert_in_darkmode&sanitize=true" align=middle width=115.02097364999997pt height=24.65753399999998pt/>, <img src="/tex/59680c09e0e2d723d0bcf2005047b028.svg?invert_in_darkmode&sanitize=true" align=middle width=73.18587374999998pt height=22.831056599999986pt/> with fixed bounds <img src="/tex/370a29c873e3d269a6111aa219085d0b.svg?invert_in_darkmode&sanitize=true" align=middle width=44.697406049999984pt height=21.18721440000001pt/> and  <img src="/tex/afdb85da7c3e7b7c3d226050994dbf5f.svg?invert_in_darkmode&sanitize=true" align=middle width=63.70246739999999pt height=14.15524440000002pt/>. We use the set <img src="/tex/36dd5900c84c0ddd4a48f4858bbd6e8f.svg?invert_in_darkmode&sanitize=true" align=middle width=264.26770754999995pt height=27.91243950000002pt/> to define the nonconstrained minimal global cost given by

<p align="center"><img src="/tex/55e63ede9d3bc1ac8c64fa74af425a24.svg?invert_in_darkmode&sanitize=true" align=middle width=277.1488236pt height=49.315569599999996pt/></p>

where <img src="/tex/99751e94989c68f9be0f6aa442bc80d5.svg?invert_in_darkmode&sanitize=true" align=middle width=40.302373649999986pt height=22.831056599999986pt/> is a penalty parameter and <img src="/tex/a44ff4154fc3bc708e9e752a14051324.svg?invert_in_darkmode&sanitize=true" align=middle width=49.762892849999986pt height=24.65753399999998pt/> is the minimal cost over the segment <img src="/tex/add1478513cabbadcd5004323f01b74c.svg?invert_in_darkmode&sanitize=true" align=middle width=62.71697189999998pt height=24.65753399999998pt/>. The penalty <img src="/tex/8217ed3c32a785f0b5aad4055f432ad8.svg?invert_in_darkmode&sanitize=true" align=middle width=10.16555099999999pt height=22.831056599999986pt/> is understood as an additional cost when introducing a new segment. The argminimum of this quantity gives a vector <img src="/tex/fa1931c6895840344a67aa4df9cf3f59.svg?invert_in_darkmode&sanitize=true" align=middle width=15.782028899999991pt height=22.63846199999998pt/> containing the last position of each segment (if we do not consider <img src="/tex/370a29c873e3d269a6111aa219085d0b.svg?invert_in_darkmode&sanitize=true" align=middle width=44.697406049999984pt height=21.18721440000001pt/>). The quantity <img src="/tex/0ce745c733fc96f40b3901e008412ce9.svg?invert_in_darkmode&sanitize=true" align=middle width=21.121448699999988pt height=22.465723500000017pt/> is the solution of the nonconstrained optimal partitioning method. 

In our setting, the cost <img src="/tex/db5f7b3e9934fbc5a2859d88c4ba84a3.svg?invert_in_darkmode&sanitize=true" align=middle width=9.614228249999991pt height=22.465723500000017pt/> is the result of the minimization of a cost function with additive property:

<p align="center"><img src="/tex/62b773158b0e824a56f2d40e326b99fe.svg?invert_in_darkmode&sanitize=true" align=middle width=496.7556429pt height=40.23961095pt/></p>

with the argminimum defining the infered mean <img src="/tex/47b592a798cd56ccf668b67abad36a61.svg?invert_in_darkmode&sanitize=true" align=middle width=19.083998999999988pt height=14.15524440000002pt/> of the i+1-th segment <img src="/tex/97a268c17395aace06ce389334ba5322.svg?invert_in_darkmode&sanitize=true" align=middle width=115.02097364999997pt height=24.65753399999998pt/> with <img src="/tex/0794cee3e13634952ef9036a12e8fb9e.svg?invert_in_darkmode&sanitize=true" align=middle width=87.79779359999999pt height=24.65753399999998pt/>. Additivity of the cost (the <img src="/tex/11c596de17c342edeed29f489aa4b274.svg?invert_in_darkmode&sanitize=true" align=middle width=9.423880949999988pt height=14.15524440000002pt/> decomposition) is guaranteed as we will use costs deriving from a likelihood. 

More generally, we can associate a current mean <img src="/tex/ce9c41bf6906ffd46ac330f09cacc47f.svg?invert_in_darkmode&sanitize=true" align=middle width=14.555823149999991pt height=14.15524440000002pt/> to each data point <img src="/tex/2b442e3e088d1b744730822d18e7aa21.svg?invert_in_darkmode&sanitize=true" align=middle width=12.710331149999991pt height=14.15524440000002pt/> and we write a cost on a segment <img src="/tex/add1478513cabbadcd5004323f01b74c.svg?invert_in_darkmode&sanitize=true" align=middle width=62.71697189999998pt height=24.65753399999998pt/> as a result of a constaint minimization:

<p align="center"><img src="/tex/4060746cfc34359c50eb0701ada25324.svg?invert_in_darkmode&sanitize=true" align=middle width=276.95500799999996pt height=48.7101186pt/></p>

so that we get another description of the objective function :

<p align="center"><img src="/tex/c1ba660018d6595156881adfbf0c7bc3.svg?invert_in_darkmode&sanitize=true" align=middle width=433.3941711pt height=44.89738935pt/></p>

with <img src="/tex/21fd4e8eecd6bdf1a4d3d6bd1fb8d733.svg?invert_in_darkmode&sanitize=true" align=middle width=8.515988249999989pt height=22.465723500000017pt/> the indicator function. This approach can be generalized to more complex constraints on consecutive means. We define a transition graph <img src="/tex/cc2fa70b4141dfd7a7fb8fe015ae1fb9.svg?invert_in_darkmode&sanitize=true" align=middle width=190.69513649999996pt height=24.65753399999998pt/> as a set of states <img src="/tex/4b924de9e932698c85d5af3791bea6bf.svg?invert_in_darkmode&sanitize=true" align=middle width=130.89004995pt height=24.65753399999998pt/> and transition sets <img src="/tex/3a0f2719f2ef8b08fde68fce48b78eef.svg?invert_in_darkmode&sanitize=true" align=middle width=31.852679099999992pt height=24.65753399999998pt/> constraining the adding of data <img src="/tex/da853d6c3fbed5f19188402bbce562e7.svg?invert_in_darkmode&sanitize=true" align=middle width=29.669143349999988pt height=14.15524440000002pt/>:

<p align="center"><img src="/tex/6560dfb0cb78500d2196c3fa76c1e272.svg?invert_in_darkmode&sanitize=true" align=middle width=134.92969875pt height=16.438356pt/></p>

with

<p align="center"><img src="/tex/c4d902f213fe885cc28bdbef2933e7a9.svg?invert_in_darkmode&sanitize=true" align=middle width=411.30755325pt height=16.438356pt/></p>

An element <img src="/tex/5808ad42fe58946d47eccfa1dc12fc3e.svg?invert_in_darkmode&sanitize=true" align=middle width=78.31948574999998pt height=24.65753399999998pt/> is described by the elements <img src="/tex/a0076764b98171f70edddddf925d32ab.svg?invert_in_darkmode&sanitize=true" align=middle width=189.22996124999997pt height=29.190975000000005pt/>. The transition function <img src="/tex/e16ea4600d976009912c7053126daad5.svg?invert_in_darkmode&sanitize=true" align=middle width=29.317325399999987pt height=14.15524440000002pt/> defines the kind of constraint we use in the transition from state <img src="/tex/1f1c28e0a1b1708c6889fb006c886784.svg?invert_in_darkmode&sanitize=true" align=middle width=12.67127234999999pt height=14.15524440000002pt/> to state <img src="/tex/b02dd36b8a10566f2a0ad9cbb2e74858.svg?invert_in_darkmode&sanitize=true" align=middle width=29.31519194999999pt height=14.15524440000002pt/> and penalized by <img src="/tex/8ffefdc6ce671e81194aa2df66aed0eb.svg?invert_in_darkmode&sanitize=true" align=middle width=31.64227604999999pt height=29.190975000000005pt/>. Table \ref{constraintstable} summarizes all the possible constraints encoded into the \proglang{gfpop} package.\\


| constraints | <img src="/tex/fb61bc5592a383b9b7bb157ba07377cc.svg?invert_in_darkmode&sanitize=true" align=middle width=125.1159426pt height=22.648391699999998pt/>, <img src="/tex/69e4724bb8055a140759e9aac47561fb.svg?invert_in_darkmode&sanitize=true" align=middle width=49.168495199999995pt height=26.17730939999998pt/> |
|---------:|-----:|
| no changepoint | <img src="/tex/3d1df3fcec1b16d211b0015ccd0aed45.svg?invert_in_darkmode&sanitize=true" align=middle width=216.88549905000002pt height=26.76175259999998pt/> |
| up | <img src="/tex/ac10c005566bdd6740023c3d10bb5f6b.svg?invert_in_darkmode&sanitize=true" align=middle width=223.93060139999994pt height=24.65753399999998pt/> |
| down | <img src="/tex/946c1e0b7c83cf57740dd328af89efd4.svg?invert_in_darkmode&sanitize=true" align=middle width=236.7160356pt height=24.65753399999998pt/> |
| absSup | <img src="/tex/0fd835f301d104615a7c6152a06d794c.svg?invert_in_darkmode&sanitize=true" align=middle width=263.7252948pt height=24.65753399999998pt/> |
| absInf | <img src="/tex/812e0170eb25ce3aa49ff182bce235ff.svg?invert_in_darkmode&sanitize=true" align=middle width=250.93986225pt height=24.65753399999998pt/> |
| no constraint | <img src="/tex/641b7841d9bcb93149b51c149843c1ee.svg?invert_in_darkmode&sanitize=true" align=middle width=106.68764864999997pt height=24.65753399999998pt/> | 

We can now write with <img src="/tex/4e6a93c1f807fa86474f40b37b73651e.svg?invert_in_darkmode&sanitize=true" align=middle width=64.18775384999998pt height=29.190975000000005pt/>,

<p align="center"><img src="/tex/5ac0ab2e47b3596655878afae57c4ab0.svg?invert_in_darkmode&sanitize=true" align=middle width=507.52615484999995pt height=97.30109565pt/></p>

In the `gfpop` package, the transition set <img src="/tex/3a0f2719f2ef8b08fde68fce48b78eef.svg?invert_in_darkmode&sanitize=true" align=middle width=31.852679099999992pt height=24.65753399999998pt/> is not time dependent. In this case, we define the graph <img src="/tex/9f1efdaeed741c2a46969cdb880cabbf.svg?invert_in_darkmode&sanitize=true" align=middle width=76.34688434999998pt height=24.65753399999998pt/>  by the label of its vertices <img src="/tex/768af24250bc547f14842daac9e5f22b.svg?invert_in_darkmode&sanitize=true" align=middle width=132.94492035pt height=24.65753399999998pt/> and its set of oriented edges <img src="/tex/deaf67742213f7a8f1c73a4583687abc.svg?invert_in_darkmode&sanitize=true" align=middle width=112.96305899999999pt height=22.465723500000017pt/> with <img src="/tex/475ea0a9e4a5111deca8a9175421df98.svg?invert_in_darkmode&sanitize=true" align=middle width=109.65749849999997pt height=24.65753399999998pt/> defining an edge going from vertex <img src="/tex/77a3b857d53fb44e33b53e4c8b68351a.svg?invert_in_darkmode&sanitize=true" align=middle width=5.663225699999989pt height=21.68300969999999pt/> to vertex <img src="/tex/36b5afebdba34564d884d347484ac0c7.svg?invert_in_darkmode&sanitize=true" align=middle width=7.710416999999989pt height=21.68300969999999pt/> with constraint option <img src="/tex/4bdc8d9bcfb35e1c9bfb51fc69687dfc.svg?invert_in_darkmode&sanitize=true" align=middle width=7.054796099999991pt height=22.831056599999986pt/>. We get

<p align="center"><img src="/tex/6dd5ac9cb97a5886486b59693ae2c347.svg?invert_in_darkmode&sanitize=true" align=middle width=345.56907494999996pt height=26.5753257pt/></p>

with 

<p align="center"><img src="/tex/3136235feeec674b990d24f57b5db978.svg?invert_in_darkmode&sanitize=true" align=middle width=551.9422012499999pt height=72.76950285pt/></p>

and <img src="/tex/42bd324589782e0f763bed8fcdc9ac3d.svg?invert_in_darkmode&sanitize=true" align=middle width=72.83814284999998pt height=26.76175259999998pt/> being the set of valid paths of the graph. We say that <img src="/tex/4f652592f0dea416b0d672295cbd1cd7.svg?invert_in_darkmode&sanitize=true" align=middle width=128.37149984999996pt height=24.65753399999998pt/> is a valid path if <img src="/tex/4a7a7d5388dddffa76a2bd3b4931c100.svg?invert_in_darkmode&sanitize=true" align=middle width=155.73607995pt height=24.65753399999998pt/>, <img src="/tex/a2b7fa301fe08a63c504d1e1a99e9207.svg?invert_in_darkmode&sanitize=true" align=middle width=155.73607995pt height=24.65753399999998pt/>,..., <img src="/tex/b6b5c37b83df12ca7ef1445aa324979a.svg?invert_in_darkmode&sanitize=true" align=middle width=207.24112035pt height=24.65753399999998pt/>. 

Notice that in most applications, the set of edges always contains edges of type <img src="/tex/3a21f5ebc13c1bc653e588fd35038b39.svg?invert_in_darkmode&sanitize=true" align=middle width=183.08234835pt height=26.76175259999998pt/> for all <img src="/tex/9b2a3c5dc55668a26467cbc9b734a6e6.svg?invert_in_darkmode&sanitize=true" align=middle width=38.996401949999985pt height=22.465723500000017pt/> as it corresponds to an absence of constraint on segment lengths. 

The main algorithm of the package, the `gfpop` function, returns the changepoint vector <img src="/tex/fa1931c6895840344a67aa4df9cf3f59.svg?invert_in_darkmode&sanitize=true" align=middle width=15.782028899999991pt height=22.63846199999998pt/> defined as <img src="/tex/0b35e94a3c1848e17627a2ad4d8fe4d8.svg?invert_in_darkmode&sanitize=true" align=middle width=196.18626389999997pt height=24.65753399999998pt/>, with <img src="/tex/48162c2e6d8e19213e506579c4cccc51.svg?invert_in_darkmode&sanitize=true" align=middle width=86.28446969999999pt height=24.65753399999998pt/> the argminimum of <img src="/tex/05839a0807c0b4eb51ed66d3f64b7ac0.svg?invert_in_darkmode&sanitize=true" align=middle width=77.228118pt height=24.65753399999998pt/> and <img src="/tex/ac4f2e6484bfba9965ba05d82ffb5cce.svg?invert_in_darkmode&sanitize=true" align=middle width=91.16641709999999pt height=19.1781018pt/>.

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

The vector `forced` is a boolean vector. A forced element means that two consecutive means have been forced to satisfy the constraint. For example, the "up" edge with parameter <img src="/tex/3e18a4a28fdee1744e5e3f79d13b9ff6.svg?invert_in_darkmode&sanitize=true" align=middle width=7.11380504999999pt height=14.15524440000002pt/> is forced if <img src="/tex/541dd8a10e3c8d8e898e9bf1d9b2dd56.svg?invert_in_darkmode&sanitize=true" align=middle width=105.57833054999998pt height=19.1781018pt/>.

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

This algorithm is called segment neighborhood in the changepoint litterature. In this example, we fixed the number of segments at <img src="/tex/5dc642f297e291cfdde8982599601d7e.svg?invert_in_darkmode&sanitize=true" align=middle width=8.219209349999991pt height=21.18721440000001pt/> with an isotonic constraint.

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

`state1` is the starting vertex of an edge, `state2` its ending vertex. `type` is one of the available edge type ("null", up", "down", "std", "absInf", "absSup"). `penalty` is a nonnegative parameter: the additional cost <img src="/tex/8217ed3c32a785f0b5aad4055f432ad8.svg?invert_in_darkmode&sanitize=true" align=middle width=10.16555099999999pt height=22.831056599999986pt/> to consider when we move within the graph using this edge. `parameter` is annother nonnegative parameter, a characteristics of the edge, depending of its type.

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
