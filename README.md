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
### February 20, 2019

> [Introduction](#intro)

> [Quick Start](#qs)

> [Some examples](#se)

> [Graph construction](#gc)

<a id="intro"></a>

## Introduction

`gfpop` is an `R` package developed to perform parametric changepoint detection in univariate time series constrained to a graph structure. The constraints are imposed to the sequence of infered means of consecutive segments and are related to the direction and/or the magnitude of the mean changes. Changepoint detection is performed using the functional pruning optimal partitioning method (fpop) based on an exact dynamic programming algorithm.  The user can specify some other constraints on the graph (starting and ending vertices) and constraint the range of means. 

For each data point, the algorithm updates a function (a functional cost) and we go through an edge. The edges of the graph can be of type "null", up", "down", "std", "absInf" or "absSup" with the following meaning:

- "null" edge : there is no changepoint introduced. We stay on the same segment
- "up" edge : the next segment has a greater mean (we can also force the size of the gap to be greater than a minimal value)
- "down" edge : the next segment has a lower mean (wan also can force the size of the gap to be greater that a minimal value)
- "std" edge : no contraint, the next segment can have any mean
- "absSup" edge : the absolute value of the difference of two consecutive means is greater than a given parameter
- "absInf" edge : the absolute value of the difference of two consecutive means is lower than a given parameter

A nonnegative internal parameter can thus be associated to an edge (in "up" and "down") or is mandatory ("absSup" and "absInf"). The "null" edge refers to an absence of changepoint. This edge can be used between different states to constraint segment lengths. Our package includes the segment neighborhood algorithm for which the number of changepoints is fixed. More details on graph construction are given in the last [section](#gc).

Data is modelized by a quadratic cost with possible use of a robust loss, biweight and Huber. In a next version of this package, other parametric costs will be available (L1, Poisson). 

The package `gfpop` is designed to segment univariate data <img src="/tex/b704a970e46e8418f3ef56718438b122.svg?invert_in_darkmode&sanitize=true" align=middle width=126.38913869999998pt height=24.65753399999998pt/> obeying to a graph structure on segment means. The changepoint vector <img src="/tex/ed535ca37308856c5f0c6d85b4d4c676.svg?invert_in_darkmode&sanitize=true" align=middle width=209.11492305pt height=27.91243950000002pt/> defines the <img src="/tex/33359de825e43daa97171e27f6558ae9.svg?invert_in_darkmode&sanitize=true" align=middle width=37.38576269999999pt height=22.831056599999986pt/> segments <img src="/tex/97a268c17395aace06ce389334ba5322.svg?invert_in_darkmode&sanitize=true" align=middle width=115.02097364999997pt height=24.65753399999998pt/>, <img src="/tex/59680c09e0e2d723d0bcf2005047b028.svg?invert_in_darkmode&sanitize=true" align=middle width=73.18587374999998pt height=22.831056599999986pt/> with fixed bounds <img src="/tex/370a29c873e3d269a6111aa219085d0b.svg?invert_in_darkmode&sanitize=true" align=middle width=44.697406049999984pt height=21.18721440000001pt/> and  <img src="/tex/afdb85da7c3e7b7c3d226050994dbf5f.svg?invert_in_darkmode&sanitize=true" align=middle width=63.70246739999999pt height=14.15524440000002pt/>. We use the set <img src="/tex/36dd5900c84c0ddd4a48f4858bbd6e8f.svg?invert_in_darkmode&sanitize=true" align=middle width=264.26770754999995pt height=27.91243950000002pt/> to define the nonconstrained minimal global cost given by

<p align="center"><img src="/tex/55e63ede9d3bc1ac8c64fa74af425a24.svg?invert_in_darkmode&sanitize=true" align=middle width=277.1488236pt height=49.315569599999996pt/></p>

where <img src="/tex/99751e94989c68f9be0f6aa442bc80d5.svg?invert_in_darkmode&sanitize=true" align=middle width=40.302373649999986pt height=22.831056599999986pt/> is a penalty parameter and <img src="/tex/a44ff4154fc3bc708e9e752a14051324.svg?invert_in_darkmode&sanitize=true" align=middle width=49.762892849999986pt height=24.65753399999998pt/> is the minimal cost over the segment <img src="/tex/add1478513cabbadcd5004323f01b74c.svg?invert_in_darkmode&sanitize=true" align=middle width=62.71697189999998pt height=24.65753399999998pt/>. The penalty <img src="/tex/8217ed3c32a785f0b5aad4055f432ad8.svg?invert_in_darkmode&sanitize=true" align=middle width=10.16555099999999pt height=22.831056599999986pt/> is understood as an additional cost when introducing a new segment. The argminimum of this quantity gives a vector <img src="/tex/fa1931c6895840344a67aa4df9cf3f59.svg?invert_in_darkmode&sanitize=true" align=middle width=15.782028899999991pt height=22.63846199999998pt/> containing the last position of each segment (if we do not consider <img src="/tex/370a29c873e3d269a6111aa219085d0b.svg?invert_in_darkmode&sanitize=true" align=middle width=44.697406049999984pt height=21.18721440000001pt/>). The quantity <img src="/tex/0ce745c733fc96f40b3901e008412ce9.svg?invert_in_darkmode&sanitize=true" align=middle width=21.121448699999988pt height=22.465723500000017pt/> is the solution of the nonconstrained optimal partitioning method. 

In our setting, the cost <img src="/tex/db5f7b3e9934fbc5a2859d88c4ba84a3.svg?invert_in_darkmode&sanitize=true" align=middle width=9.614228249999991pt height=22.465723500000017pt/> is the result of the minimization of a cost function with additive property:

<p align="center"><img src="/tex/aca209e810ec973a24c9e1fb4d2828a0.svg?invert_in_darkmode&sanitize=true" align=middle width=496.7556429pt height=40.23961095pt/></p>

with the argminimum defining the infered mean <img src="/tex/47b592a798cd56ccf668b67abad36a61.svg?invert_in_darkmode&sanitize=true" align=middle width=19.083998999999988pt height=14.15524440000002pt/> of the i+1-th segment <img src="/tex/97a268c17395aace06ce389334ba5322.svg?invert_in_darkmode&sanitize=true" align=middle width=115.02097364999997pt height=24.65753399999998pt/> with <img src="/tex/0794cee3e13634952ef9036a12e8fb9e.svg?invert_in_darkmode&sanitize=true" align=middle width=87.79779359999999pt height=24.65753399999998pt/>. Additivity of the cost (the <img src="/tex/11c596de17c342edeed29f489aa4b274.svg?invert_in_darkmode&sanitize=true" align=middle width=9.423880949999988pt height=14.15524440000002pt/> decomposition) is guaranteed as we will use costs deriving from a likelihood. 

More generally, we can associate a current mean <img src="/tex/ce9c41bf6906ffd46ac330f09cacc47f.svg?invert_in_darkmode&sanitize=true" align=middle width=14.555823149999991pt height=14.15524440000002pt/> to each data point <img src="/tex/2b442e3e088d1b744730822d18e7aa21.svg?invert_in_darkmode&sanitize=true" align=middle width=12.710331149999991pt height=14.15524440000002pt/> and we write a cost on a segment <img src="/tex/add1478513cabbadcd5004323f01b74c.svg?invert_in_darkmode&sanitize=true" align=middle width=62.71697189999998pt height=24.65753399999998pt/> as a result of a constrained minimization:

<p align="center"><img src="/tex/4060746cfc34359c50eb0701ada25324.svg?invert_in_darkmode&sanitize=true" align=middle width=276.95500799999996pt height=48.7101186pt/></p>

so that we get another description of the objective function :

<p align="center"><img src="/tex/c1ba660018d6595156881adfbf0c7bc3.svg?invert_in_darkmode&sanitize=true" align=middle width=433.3941711pt height=44.89738935pt/></p>

with <img src="/tex/37ca9e919ee5a54f28ed9f70e6ed1277.svg?invert_in_darkmode&sanitize=true" align=middle width=68.78982494999998pt height=24.65753399999998pt/> the indicator function. This approach can be generalized to more complex constraints on consecutive means using a graph structure. We define the transition graph <img src="/tex/0553770832b41fddeb39cbfa76cf036a.svg?invert_in_darkmode&sanitize=true" align=middle width=17.90463839999999pt height=22.465723500000017pt/> as a directed acyclic graph with the following properties:

1. Vertices are indexed by time and state. <img src="/tex/6c8ced74f5b3bcbfbd287c91837d0b64.svg?invert_in_darkmode&sanitize=true" align=middle width=64.20836564999999pt height=24.65753399999998pt/> for vertex <img src="/tex/6c4adbc36120d62b98deef2a20d5d303.svg?invert_in_darkmode&sanitize=true" align=middle width=8.55786029999999pt height=14.15524440000002pt/> with time <img src="/tex/4f4f4e395762a3af4575de74c019ebb5.svg?invert_in_darkmode&sanitize=true" align=middle width=5.936097749999991pt height=20.221802699999984pt/> and state <img src="/tex/6f9bad7347b91ceebebd3ad7e6f6f2d1.svg?invert_in_darkmode&sanitize=true" align=middle width=7.7054801999999905pt height=14.15524440000002pt/>. The states are elements of the set <img src="/tex/4b924de9e932698c85d5af3791bea6bf.svg?invert_in_darkmode&sanitize=true" align=middle width=130.89004995pt height=24.65753399999998pt/>;

2. All the vertices have time <img src="/tex/4f4f4e395762a3af4575de74c019ebb5.svg?invert_in_darkmode&sanitize=true" align=middle width=5.936097749999991pt height=20.221802699999984pt/> in <img src="/tex/868738293e2285473ac5c2b2c7d7f7a8.svg?invert_in_darkmode&sanitize=true" align=middle width=62.83494359999998pt height=24.65753399999998pt/> except for the unique starting vertex <img src="/tex/ce927a5e91802f349da081165382e133.svg?invert_in_darkmode&sanitize=true" align=middle width=73.78989254999999pt height=24.65753399999998pt/> and the unique ending vertex <img src="/tex/f508dd141f08e86bda0a55c9e2815f37.svg?invert_in_darkmode&sanitize=true" align=middle width=121.96534844999998pt height=24.65753399999998pt/>, where <img src="/tex/53fe7f87db69e0ed1312d865111c131f.svg?invert_in_darkmode&sanitize=true" align=middle width=8.219209349999991pt height=24.65753399999998pt/> denotes an undefinite state;

3. Edges are transitions between "time-consecutive" vertices of type <img src="/tex/6c8ced74f5b3bcbfbd287c91837d0b64.svg?invert_in_darkmode&sanitize=true" align=middle width=64.20836564999999pt height=24.65753399999998pt/> and <img src="/tex/997bd6ae1076fb8e024546d2d8852578.svg?invert_in_darkmode&sanitize=true" align=middle width=101.74251614999997pt height=24.7161288pt/> which gives us the edge <img src="/tex/7d3f780dafd2bd5f98d6febd6cc41f02.svg?invert_in_darkmode&sanitize=true" align=middle width=82.92789945pt height=24.7161288pt/> for <img src="/tex/c8d3e8f08b1d18c8f9bee7f42c6d6797.svg?invert_in_darkmode&sanitize=true" align=middle width=88.86217889999998pt height=24.65753399999998pt/>;

4. Each edge <img src="/tex/8cd34385ed61aca950a6b06d09fb50ac.svg?invert_in_darkmode&sanitize=true" align=middle width=7.654137149999991pt height=14.15524440000002pt/> is associated to a function <img src="/tex/71d7882ef649a7537f2e3c662f98092a.svg?invert_in_darkmode&sanitize=true" align=middle width=137.5721589pt height=24.65753399999998pt/> constraining consecutive means <img src="/tex/87eefe082e181864d1321025c2705ecd.svg?invert_in_darkmode&sanitize=true" align=middle width=14.870715749999988pt height=14.15524440000002pt/> and <img src="/tex/1790d105076a15b61f11a9d6e6296025.svg?invert_in_darkmode&sanitize=true" align=middle width=31.51463699999999pt height=14.15524440000002pt/> by the indicator function <img src="/tex/c261ab74762f29919bc44aa0c519c8e4.svg?invert_in_darkmode&sanitize=true" align=middle width=176.09637045pt height=24.65753399999998pt/> (for example <img src="/tex/4919996dd46d961cb23ba7804cbf9d71.svg?invert_in_darkmode&sanitize=true" align=middle width=195.57111419999998pt height=24.65753399999998pt/>) and a penalty <img src="/tex/2111bbb895243afe1ffa4d502cfa57ed.svg?invert_in_darkmode&sanitize=true" align=middle width=46.493620799999995pt height=22.831056599999986pt/>.

The next table summarizes all the possible constraints encoded into the `gfpop` package.


| constraints | <img src="/tex/cc7c3637db1e725c477f85fd4d9903cd.svg?invert_in_darkmode&sanitize=true" align=middle width=109.26161729999998pt height=22.648391699999998pt/>, <img src="/tex/69e4724bb8055a140759e9aac47561fb.svg?invert_in_darkmode&sanitize=true" align=middle width=49.168495199999995pt height=26.17730939999998pt/> |
|---------:|-----:|
| no changepoint | <img src="/tex/fd7c00e2a1d565658787e8000b2f94ba.svg?invert_in_darkmode&sanitize=true" align=middle width=195.57111419999998pt height=24.65753399999998pt/> |
| up | <img src="/tex/8e4b447e027d0fad46fd5e86d22ceece.svg?invert_in_darkmode&sanitize=true" align=middle width=222.7761096pt height=24.65753399999998pt/> |
| down | <img src="/tex/5d084c33cc699e5de28eb57d5eeef4d0.svg?invert_in_darkmode&sanitize=true" align=middle width=222.7761096pt height=24.65753399999998pt/> |
| absSup | <img src="/tex/c305f8d88f6b5b252ea9903d92b74329.svg?invert_in_darkmode&sanitize=true" align=middle width=249.78536879999996pt height=24.65753399999998pt/> |
| absInf | <img src="/tex/36e1f046428d5a9dd03e10dfec598742.svg?invert_in_darkmode&sanitize=true" align=middle width=249.78536879999996pt height=24.65753399999998pt/> |
| no constraint | <img src="/tex/aa6acfdc302d20cce40de97d9a0d6bdc.svg?invert_in_darkmode&sanitize=true" align=middle width=195.57111419999998pt height=24.65753399999998pt/> | 


We define a path <img src="/tex/f26713629da615784507a3edf48fb24b.svg?invert_in_darkmode&sanitize=true" align=middle width=46.26634319999999pt height=22.465723500000017pt/> of the graph as a collection of <img src="/tex/014ae02bf2b677c4c225a7dfd00b8420.svg?invert_in_darkmode&sanitize=true" align=middle width=38.17727759999999pt height=21.18721440000001pt/> vertices <img src="/tex/54d840342c0d4e1f4770f2572852818c.svg?invert_in_darkmode&sanitize=true" align=middle width=89.99831279999998pt height=24.65753399999998pt/> with <img src="/tex/93729215cb2ecc9aea94ffd3ff1ed5bf.svg?invert_in_darkmode&sanitize=true" align=middle width=73.78989254999999pt height=24.65753399999998pt/> and <img src="/tex/f508dd141f08e86bda0a55c9e2815f37.svg?invert_in_darkmode&sanitize=true" align=middle width=121.96534844999998pt height=24.65753399999998pt/> and <img src="/tex/73eb8592e5cb867e19391571ae12bfb0.svg?invert_in_darkmode&sanitize=true" align=middle width=75.19396665pt height=24.65753399999998pt/> for <img src="/tex/985dbad85d61b07e704840368824ee09.svg?invert_in_darkmode&sanitize=true" align=middle width=88.86217889999998pt height=24.65753399999998pt/> and <img src="/tex/308736d40f48c7c0c2e28c99b956051e.svg?invert_in_darkmode&sanitize=true" align=middle width=44.77148444999999pt height=22.465723500000017pt/>. Morever, the path is made of <img src="/tex/3f18d8f60c110e865571bba5ba67dcc6.svg?invert_in_darkmode&sanitize=true" align=middle width=38.17727759999999pt height=21.18721440000001pt/> edges denoted <img src="/tex/892f838396f392c3a6c9da53214ef424.svg?invert_in_darkmode&sanitize=true" align=middle width=59.119196399999986pt height=14.15524440000002pt/> with <img src="/tex/b8e1a11c627257800b07edffe7fbe8a9.svg?invert_in_darkmode&sanitize=true" align=middle width=54.554976299999986pt height=22.831056599999986pt/>. A vector <img src="/tex/ea48643af203a0e7ad01bcafddbc8f82.svg?invert_in_darkmode&sanitize=true" align=middle width=49.99426409999999pt height=22.648391699999998pt/> verifies the path <img src="/tex/2ec6e630f199f589a2402fdf3e0289d5.svg?invert_in_darkmode&sanitize=true" align=middle width=8.270567249999992pt height=14.15524440000002pt/> if for all <img src="/tex/875f0d0d37d49d24b63ba045d9e7491a.svg?invert_in_darkmode&sanitize=true" align=middle width=117.1725786pt height=24.65753399999998pt/>, we have <img src="/tex/a9f7bf6538d048b2d8fae862a45abeb5.svg?invert_in_darkmode&sanitize=true" align=middle width=117.99692024999999pt height=24.65753399999998pt/> (valid constraint). We write <img src="/tex/5a2bcf863517582f776c5598b8a1c6b2.svg?invert_in_darkmode&sanitize=true" align=middle width=30.960925049999993pt height=24.65753399999998pt/> to say that the vector <img src="/tex/07617f9d8fe48b4a7b3f523d6730eef0.svg?invert_in_darkmode&sanitize=true" align=middle width=9.90492359999999pt height=14.15524440000002pt/> verifies the path <img src="/tex/2ec6e630f199f589a2402fdf3e0289d5.svg?invert_in_darkmode&sanitize=true" align=middle width=8.270567249999992pt height=14.15524440000002pt/>. The formulation of our graph-constrained problem is then the following:

<p align="center"><img src="/tex/25a089a47a987094b647297062ee8766.svg?invert_in_darkmode&sanitize=true" align=middle width=310.64316525pt height=51.74449499999999pt/></p>


In the `gfpop` package, the edges <img src="/tex/46ed17799525427705742ef69fb71fd9.svg?invert_in_darkmode&sanitize=true" align=middle width=53.356130849999985pt height=24.7161288pt/> for <img src="/tex/2da0ae8ede87918239efbfd469567202.svg?invert_in_darkmode&sanitize=true" align=middle width=88.86217889999998pt height=24.65753399999998pt/> are not time-dependent. In this case, we redefine a graph <img src="/tex/9f1efdaeed741c2a46969cdb880cabbf.svg?invert_in_darkmode&sanitize=true" align=middle width=76.34688434999998pt height=24.65753399999998pt/>  by the label of its vertices <img src="/tex/768af24250bc547f14842daac9e5f22b.svg?invert_in_darkmode&sanitize=true" align=middle width=132.94492035pt height=24.65753399999998pt/> and its set of oriented edges <img src="/tex/526e9ac1fb94812f38332e979abb135a.svg?invert_in_darkmode&sanitize=true" align=middle width=81.57506609999999pt height=22.465723500000017pt/> with <img src="/tex/88879fa4b41455559f35230fb5d4dcf1.svg?invert_in_darkmode&sanitize=true" align=middle width=102.85925264999997pt height=24.7161288pt/> defining an edge going from vertex <img src="/tex/6f9bad7347b91ceebebd3ad7e6f6f2d1.svg?invert_in_darkmode&sanitize=true" align=middle width=7.7054801999999905pt height=14.15524440000002pt/> to vertex <img src="/tex/675c2f5707a1fa7050c12adc1872ba32.svg?invert_in_darkmode&sanitize=true" align=middle width=11.49544109999999pt height=24.7161288pt/> with indicator function <img src="/tex/629cbe1de238666b8fbb23255f22ea6b.svg?invert_in_darkmode&sanitize=true" align=middle width=13.46296874999999pt height=22.465723500000017pt/> and penalty <img src="/tex/8a51f4aa3be26df848de9e1e89f6a4a6.svg?invert_in_darkmode&sanitize=true" align=middle width=15.534883649999992pt height=22.831056599999986pt/>. 

In most applications, the set of edges always contains edges of type <img src="/tex/f308eb98d3204881e73de0e34ff40cfa.svg?invert_in_darkmode&sanitize=true" align=middle width=35.502276149999986pt height=24.65753399999998pt/> for all <img src="/tex/e6bfc1a9751a1dda5621c949c34e793c.svg?invert_in_darkmode&sanitize=true" align=middle width=41.03865479999999pt height=22.465723500000017pt/> with indicator function <img src="/tex/34c893ff78debf4a9f1763de90d3a70a.svg?invert_in_darkmode&sanitize=true" align=middle width=137.6560581pt height=24.65753399999998pt/> as it corresponds to an absence of constraint on segment lengths. 

The main algorithm of the package, the `gfpop` function, returns the changepoint vector <img src="/tex/fa1931c6895840344a67aa4df9cf3f59.svg?invert_in_darkmode&sanitize=true" align=middle width=15.782028899999991pt height=22.63846199999998pt/> defined as <img src="/tex/0b35e94a3c1848e17627a2ad4d8fe4d8.svg?invert_in_darkmode&sanitize=true" align=middle width=196.18626389999997pt height=24.65753399999998pt/>, with <img src="/tex/48162c2e6d8e19213e506579c4cccc51.svg?invert_in_darkmode&sanitize=true" align=middle width=86.28446969999999pt height=24.65753399999998pt/> the argminimum of <img src="/tex/05839a0807c0b4eb51ed66d3f64b7ac0.svg?invert_in_darkmode&sanitize=true" align=middle width=77.228118pt height=24.65753399999998pt/> in <img src="/tex/63dee4f83c8372c83bbf2389f1c3a964.svg?invert_in_darkmode&sanitize=true" align=middle width=53.455350299999985pt height=24.65753399999998pt/> and <img src="/tex/ac4f2e6484bfba9965ba05d82ffb5cce.svg?invert_in_darkmode&sanitize=true" align=middle width=91.16641709999999pt height=19.1781018pt/>. we also give the possibility to restrict the set of valid paths by imposing a starting and/or an ending vertex and contraint the range of infered means, replacing <img src="/tex/ea48643af203a0e7ad01bcafddbc8f82.svg?invert_in_darkmode&sanitize=true" align=middle width=49.99426409999999pt height=22.648391699999998pt/> by <img src="/tex/502c144cca7929a6af3ab94ff27edc73.svg?invert_in_darkmode&sanitize=true" align=middle width=72.0565956pt height=24.65753399999998pt/> in the definition of <img src="/tex/63dee4f83c8372c83bbf2389f1c3a964.svg?invert_in_darkmode&sanitize=true" align=middle width=53.455350299999985pt height=24.65753399999998pt/>.

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
myGraph <- graph(penalty = 2*log(n), type = "updown")
```

The gfpop function gives the result of the segmentation using `myData` and `myGraph` as parameters. We choose a gaussian cost.

```r
gfpop(vectData = myData, mygraph = myGraph, type = "gauss")
```

```
## changepoints
## [1]  103  301  504  800 1000
## 
## states
## [1] 0 1 0 1 0
## 
## forced
## [1] 0 0 0 0
## 
## means
## [1] 1.0510341 2.0324955 0.9894555 2.9390439 1.1213021
## 
## cost
## [1] 1003.595
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

The isotonic regression infers a sequence of nondecreasing means. 


```r
n <- 1000
mydata <- dataGenerator(n, c(0.1, 0.2, 0.3, 0.4, 0.6, 0.8, 1), c(0, 0.5, 1, 1.5, 2, 2.5, 3), 1)
myGraphIso <- graph(penalty = 2*log(n), type = "isotonic")
gfpop(vectData =  mydata, mygraph = myGraphIso, type = "gauss", K = 1, min = 0.5)
```

```
## changepoints
## [1]  286  622 1000
## 
## states
## [1] 0 0 0
## 
## forced
## [1] 0 0
## 
## means
## [1] 0.500000 1.949243 2.839941
## 
## cost
## [1] 546.5847
## 
## attr(,"class")
## [1] "gfpop"
```


In this example, we use in `gfpop` function a robust biweight gaussian cost with `K = 1` and the `min` parameter in order to infer means greater than `0.5`.

### Fixed number of changepoint

This algorithm is called segment neighborhood in the changepoint litterature. In this example, we fixed the number of segments at <img src="/tex/5dc642f297e291cfdde8982599601d7e.svg?invert_in_darkmode&sanitize=true" align=middle width=8.219209349999991pt height=21.18721440000001pt/> with an isotonic constraint. The graph contains two "up" edges with no cycling.


```r
n <- 1000
mydata <- dataGenerator(n, c(0.1, 0.2, 0.3, 0.4, 0.6, 0.8, 1), c(0, 0.5, 1, 1.5, 2, 2.5, 3), 1)
beta <- 0
myGraph <- graph(
  edge(0, 1,"up", beta, 0),
  edge(1, 2, "up", beta, 0),
  edge(0, 0, "null", 0, 0),
  edge(1, 1, "null", 0, 0),
  edge(2, 2, "null", 0, 0),
  StartEnd(start = 0, end = 2))

gfpop(vectData =  mydata, mygraph = myGraph, type = "gauss")
```

```
## changepoints
## [1]  319  572 1000
## 
## states
## [1] 0 1 2
## 
## forced
## [1] 0 0
## 
## means
## [1] 0.5484494 1.7347737 2.7097749
## 
## cost
## [1] 1040.883
## 
## attr(,"class")
## [1] "gfpop"
```


### Robust up-down with constrained starting and ending states

In presence of outliers we need a robust loss (biweight). We can also force the starting and ending state and a minimal gap between the means (here equal to `1`)

```r
n <- 1000
mydata <- dataGenerator(n, c(0.1,0.3,0.5,0.8,1), c(0,1,0,1,0), 1) + 5*(rbinom(n, 1, 0.05)) - 5*(rbinom(n, 1, 0.05))
beta <- 2*log(n)
myGraph <- graph(
  edge(0, 1, "up", beta, 1),
  edge(1, 0, "down", beta, 1),
  edge(0, 0, "null", 0, 0),
  edge(1, 1, "null", 0, 0),
  StartEnd(start = 0, end = 0))
gfpop(vectData =  mydata, mygraph = myGraph, type = "gauss", K = 3.0)
```

```
## changepoints
## [1]  100  303  509  802 1000
## 
## states
## [1] 0 1 0 1 0
## 
## forced
## [1] 1 1 0 0
## 
## means
## [1] -0.02179054  0.97820946 -0.02179054  1.04794919 -0.08340491
## 
## cost
## [1] 1735.151
## 
## attr(,"class")
## [1] "gfpop"
```

If we skip all these constraints and use a standard fpop algorithm, the result is the following


```r
myGraphStd <- graph(penalty = 2*log(n), type = "std")
gfpop(vectData =  mydata, mygraph = myGraphStd, type = "gauss")
```

```
## changepoints
##  [1]   15   16   97   99  176  177  191  192  193  206  207  231  232  239
## [15]  240  243  244  253  254  291  292  303  313  314  337  341  353  354
## [29]  382  383  400  401  422  423  436  438  501  502  553  559  560  567
## [43]  576  577  607  608  633  634  665  666  685  686  713  714  761  762
## [57]  774  775  802  833  834  836  855  857  910  914  974  975 1000
## 
## states
## integer(0)
## 
## forced
## integer(0)
## 
## means
##  [1]  0.35516059  7.31539702  0.14105579 -2.94273330  0.94784112
##  [6]  7.47034210  0.03027980  6.44775993 -3.56469578  1.12363695
## [11]  6.79551386  0.58513289  5.88494768  0.08772338 -4.59107925
## [16]  2.25543420  6.59442166  0.95170787 -5.05427726  0.90129719
## [21]  6.44712402  1.15188801 -0.65090828 -5.35417484  0.15556922
## [26] -3.08864520  0.30447161 -6.30623891 -0.24484582  6.33510727
## [31] -0.22640480  5.22219708  0.14319821  5.33118850 -0.80120905
## [36] -4.51670060  0.06083162  5.68602427  0.79619509 -1.66803663
## [41]  5.39348658  0.71371503  3.22587516 -3.32261711  1.14748334
## [46] -4.30720760  1.12006466  6.56211080  0.94195841 -4.32791321
## [51]  0.28925888  6.78461467  1.81052567 -4.27711226  1.41438988
## [56]  6.55629102  0.46431458  6.86045040  1.60554148 -0.71720586
## [61] -4.79349869  2.95168993 -0.20751388  3.77722757 -0.15332081
## [66]  3.22710032 -0.31051196 -6.14093927  0.49413508
## 
## cost
## [1] 2699.994
## 
## attr(,"class")
## [1] "gfpop"
```


### absInf and absSup edges

With a unique "absInf" edge, we impose a difference between the means of size at most 1.  

```r
n <- 10000
mydata <- dataGenerator(n, c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1), c(0, 1, 0, 2, 1, 2, 0, 1, 0, 1), 0.5)
beta <- 2*log(n)
myGraph <- graph(
  edge(0, 0,"absInf", beta, 1),
  edge(0, 0,"null", 0, 0))
gfpop(vectData =  mydata, mygraph = myGraph, type = "gauss", K = 3)
```

```
## changepoints
##  [1]   999  1999  3000  3001  4003  5000  5999  6000  7000  8000  9001
## [12] 10000
## 
## states
##  [1] 0 0 0 0 0 0 0 0 0 0 0 0
## 
## forced
##  [1] 1 0 1 0 0 0 1 1 1 0 1
## 
## means
##  [1] -0.0155661567  0.9844338433 -0.0067182663  0.9932817337  1.9829459867
##  [6]  1.0239684275  1.9966187678  0.9966187678 -0.0033812322  0.9966187678
## [11]  0.0007690086  1.0007690086
## 
## cost
## [1] 2651.255
## 
## attr(,"class")
## [1] "gfpop"
```

With a unique "absSup" edge, we impose a difference between the means of size at least 1.  

```r
n <- 10000
mydata <- dataGenerator(n, c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1), c(0, 1, 0, 2, 1, 2, 0, 1, 0, 1), 0.5)
myGraph <- graph()
beta <- 2*log(n)
myGraph <- graph(
  edge(0, 0,"absSup", beta, 1),
  edge(0, 0,"null", 0, 0))
gfpop(vectData =  mydata, mygraph = myGraph, type = "gauss", K = 3)
```

```
## changepoints
##  [1]  1000  1997  3000  4000  5000  6000  7000  7999  9000 10000
## 
## states
##  [1] 0 0 0 0 0 0 0 0 0 0
## 
## forced
## [1] 0 1 0 0 1 0 0 0 1
## 
## means
##  [1] -0.0007571859  1.0132097237  0.0132097237  1.9970387390  0.9901447187
##  [6]  1.9901447187 -0.0014977899  1.0032445639  0.0002328453  1.0002328453
## 
## cost
## [1] 2619.784
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

`state1` is the starting vertex of an edge, `state2` its ending vertex. `type` is one of the available edge type ("null", up", "down", "std", "absInf", "absSup"). `penalty` is a nonnegative parameter: the additional cost <img src="/tex/3d13090ef3ed1448f3c4dc166d06ab4d.svg?invert_in_darkmode&sanitize=true" align=middle width=13.948864049999989pt height=22.831056599999986pt/> to consider when we move within the graph using this edge. `parameter` is annother nonnegative parameter, a characteristics of the edge, depending of its type.

We add edges as follows

```r
myGraph <- graph(
  edge(0, 0, "down", 3.1415, 1),
  edge(0, 0))
myGraph
```

```
##   state1 state2 type penalty parameter
## 1      0      0 down  3.1415         1
## 2      0      0 null  0.0000         0
```

we can only add edges to this dataframe using the object `edge`. With the example `edge(0, 0, "down", 3.1415, 1)` we give the 5 aforementioned features.

The graph can contain information on the starting and/or ending edge to use with the `StartEnd` function. 


```r
beta <- 2 * log(n)
myGraph <- graph(
  edge(0, 0, "null", 0, 0),
  edge(1, 1, "null", 0, 0),
  edge(0, 1, "up", beta, 1),
  edge(0, 0, "down", beta, 0),
  edge(1, 0, "down", beta, 0),
  StartEnd(start = 0, end = 0))
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

