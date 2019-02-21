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
### February 21, 2019

> [Introduction](#intro)

> [Quick Start](#qs)

> [Some examples](#se)

> [Graph construction](#gc)

<a id="intro"></a>

## Introduction

`gfpop` is an `R` package developed to perform parametric changepoint detection in univariate time series constrained to a graph structure. The constraints are imposed to the sequence of infered means of consecutive segments and are related to the direction and/or the magnitude of the mean changes. Changepoint detection is performed using the functional pruning optimal partitioning method (fpop) based on an exact dynamic programming algorithm.  The user can specify some other constraints on the graph (starting and ending vertices) and constraint the range of means. 

For each data point, the algorithm updates a function (a functional cost) and we go through an edge. The edges of the graph can be of type "null", up", "down", "std", "absInf" or "absSup" with the following meaning:

- "null" edge : there is no changepoint introduced. We stay on the same segment (if its parameter < 1 there is an exponetional decay with gamma = parameter)
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
## [1]  101  305  501  800 1000
## 
## states
## [1] 0 1 0 1 0
## 
## forced
## [1] 0 0 0 0
## 
## means
## [1] 0.9412859 1.9042115 1.0064542 2.9324436 1.0024089
## 
## cost
## [1] 1080.428
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
gfpop(vectData =  mydata, mygraph = myGraphIso, type = "gauss", K = 1, min = 0)
```

```
## changepoints
## [1]  293  594 1000
## 
## states
## [1] 0 0 0
## 
## forced
## [1] 0 0
## 
## means
## [1] 0.5339437 1.8341719 2.7824167
## 
## cost
## [1] 561.8481
## 
## attr(,"class")
## [1] "gfpop"
```


In this example, we use in `gfpop` function a robust biweight gaussian cost with `K = 1` and the `min` parameter in order to infer means greater than `0.5`.

### Fixed number of changepoints

This algorithm is called segment neighborhood in the changepoint litterature. In this example, we fixed the number of segments at <img src="/tex/5dc642f297e291cfdde8982599601d7e.svg?invert_in_darkmode&sanitize=true" align=middle width=8.219209349999991pt height=21.18721440000001pt/> with an isotonic constraint. The graph contains two "up" edges with no cycling.


```r
n <- 1000
mydata <- dataGenerator(n, c(0.1, 0.2, 0.3, 0.4, 0.6, 0.8, 1), c(0, 0.5, 1, 1.5, 2, 2.5, 3), 1)
beta <- 0
myGraph <- graph(
  edge(0, 1,"up", beta, 0),
  edge(1, 2, "up", beta, 0),
  edge(0, 0, "null"),
  edge(1, 1, "null"),
  edge(2, 2, "null"),
  StartEnd(start = 0, end = 2))

gfpop(vectData =  mydata, mygraph = myGraph, type = "gauss")
```

```
## changepoints
## [1]  194  413 1000
## 
## states
## [1] 0 1 2
## 
## forced
## [1] 0 0
## 
## means
## [1] 0.1924616 1.3769578 2.5875370
## 
## cost
## [1] 1034.578
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
  edge(0, 0, "null"),
  edge(1, 1, "null"),
  StartEnd(start = 0, end = 0))
gfpop(vectData =  mydata, mygraph = myGraph, type = "gauss", K = 3.0)
```

```
## changepoints
## [1]   96  295  489  802 1000
## 
## states
## [1] 0 1 0 1 0
## 
## forced
## [1] 0 0 0 0
## 
## means
## [1] -0.08351631  1.04045198 -0.15651782  0.93319389 -0.11748538
## 
## cost
## [1] 1740.751
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
##  [1]    5    6   24   27   29   30   41   42   96  144  145  171  173  176
## [15]  178  205  206  236  238  262  263  289  291  309  312  320  321  370
## [29]  371  411  416  443  445  446  451  452  487  488  508  509  520  521
## [43]  529  530  531  541  542  584  585  586  628  629  687  689  690  722
## [57]  723  729  731  800  807  808  834  835  856  857  858  894  895  935
## [71]  936  939  940  941  943  968  969 1000
## 
## states
## integer(0)
## 
## forced
## integer(0)
## 
## means
##  [1] -0.06013445  7.11421933 -0.59869423  3.62388914 -0.41207528
##  [6]  6.48511471 -0.20776566 -5.92732191 -0.03661876  1.34082032
## [11]  6.86803702  0.93634094  4.70619384 -0.08778119  5.57030161
## [16]  1.30777679 -4.88591074  0.99432297  5.67969496  0.75588296
## [21]  6.57865137  1.03417121 -3.51513212  0.23617510 -3.82238332
## [26]  0.33942023 -6.32922241 -0.03165873  5.97623876 -0.33528805
## [31]  2.30657546 -0.45914429 -4.14007789  4.63979848 -0.23235047
## [36]  5.75254104  0.04376510  5.66081320  0.99441630 -4.65957926
## [41]  1.04924942 -5.86571587  1.09354768  7.20938416 -4.30372472
## [46]  1.78287878 -4.38531646  1.15980381  6.50141769 -4.64412405
## [51]  0.91219510 -5.59853085  0.88897307 -3.86200315  4.65701787
## [56]  0.45249591  6.07602562  1.18201496  5.76491831  0.84653415
## [61] -1.21518778 -6.50467305 -0.22839960 -6.40144962 -0.10881665
## [66]  7.04677615 -4.71194613 -0.18989655  6.11902886  0.12685642
## [71]  6.17752706 -0.05709973 -4.70792489  5.29613535 -3.05890514
## [76]  0.01279152  6.29002388 -0.15696088
## 
## cost
## [1] 2652.346
## 
## attr(,"class")
## [1] "gfpop"
```


### absSup edges

With a unique "absSup" edge, we impose a difference between the means of size at least 1.  

```r
n <- 10000
mydata <- dataGenerator(n, c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1), c(0, 1, 0, 2, 1, 2, 0, 1, 0, 1), 0.5)
beta <- 2*log(n)
myGraph <- graph(
  edge(0, 0,"absSup", beta, 1),
  edge(0, 0,"null"))
gfpop(vectData =  mydata, mygraph = myGraph, type = "gauss", K = 3)
```

```
## changepoints
##  [1]  1000  2000  3000  4000  5000  6000  7005  7999  9000 10000
## 
## states
##  [1] 0 0 0 0 0 0 0 0 0 0
## 
## forced
## [1] 0 0 0 0 1 0 1 0 1
## 
## means
##  [1]  0.0110928807  1.0327701200 -0.0024735656  2.0088922942  0.9866282085
##  [6]  1.9866282085 -0.0001655706  0.9998344294 -0.0040009817  0.9959990183
## 
## cost
## [1] 2732.468
## 
## attr(,"class")
## [1] "gfpop"
```

Notice that some of the edges are forced, the vector `forced` contains non-zero values.


### Exponential decay

The null edge corresponds to an exponential decay state if its parameter is not equal to 1. 


```r
n <- 1000
mydata <- dataGenerator(n, c(0.2, 0.5, 0.8, 1), c(5, 10, 15, 20), 1, gamma = 0.966)
beta <- 2*log(n)
myGraphDecay <- graph(
  edge(0, 0, "up", beta, 0),
  edge(0, 0, "null", 0, 0.966)
  )
g <- gfpop(vectData =  mydata, mygraph = myGraphDecay, type = "gauss", min = 0)
g
```

```
## changepoints
## [1]  200  500  800 1000
## 
## states
## [1] 0 0 0 0
## 
## forced
## [1] 0 0 0
## 
## means
## [1] 0.0048666904 0.0003157411 0.0004712629 0.0198766992
## 
## cost
## [1] 1078.436
## 
## attr(,"class")
## [1] "gfpop"
```

and we plot the result 

```r
gamma <- 0.966
len <- diff(c(0, g<img src="/tex/46024bfe5198e5e1ff626011648c9edf.svg?invert_in_darkmode&sanitize=true" align=middle width=558.8873268pt height=24.65753399999998pt/>means[i]*c(1, cumprod(rep(1/gamma,len[i]-1))))}
signal <- rev(signal)

ylimits <- c(min(mydata), max(mydata))
plot(mydata, type ='p', pch ='+', ylim = ylimits)
par(new = TRUE)
plot(signal, type ='l', col = 4, ylim = ylimits, lwd = 3)
```

![plot of chunk unnamed-chunk-11](figure/unnamed-chunk-11-1.png)


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
## 2      0      0 null  0.0000         1
```

we can only add edges to this dataframe using the object `edge`. With the example `edge(0, 0, "down", 3.1415, 1)` we give the 5 aforementioned features.

The graph can contain information on the starting and/or ending edge to use with the `StartEnd` function. 


```r
beta <- 2 * log(n)
myGraph <- graph(
  edge(0, 0, "null"),
  edge(1, 1, "null"),
  edge(0, 1, "up", beta, 1),
  edge(0, 0, "down", beta, 0),
  edge(1, 0, "down", beta, 0),
  StartEnd(start = 0, end = 0))
myGraph
```

```
##   state1 state2  type  penalty parameter
## 1      0      0  null  0.00000         1
## 2      1      1  null  0.00000         1
## 3      0      1    up 13.81551         1
## 4      0      0  down 13.81551         0
## 5      1      0  down 13.81551         0
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
## 1      0      0 null       0         1
## 2      0      0   up      12         0
```


[Back to Top](#top)
