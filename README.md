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
### February 21, 2022

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

The algorithm of our package is designed to consider a large variety of constraints. These constraints are modelled by a graph. At each time <img src="/tex/4f4f4e395762a3af4575de74c019ebb5.svg?invert_in_darkmode&sanitize=true" align=middle width=5.936097749999991pt height=20.221802699999984pt/> a number of states is accessible, these states are the nodes of the graph. Possible transitions between states at time <img src="/tex/4f4f4e395762a3af4575de74c019ebb5.svg?invert_in_darkmode&sanitize=true" align=middle width=5.936097749999991pt height=20.221802699999984pt/> and <img src="/tex/628783099380408a32610228991619a8.svg?invert_in_darkmode&sanitize=true" align=middle width=34.24649744999999pt height=21.18721440000001pt/> are represented by the edges of the graph. Each edge is associated to <img src="/tex/ecf4fe2774fd9244b4fd56f7e76dc882.svg?invert_in_darkmode&sanitize=true" align=middle width=8.219209349999991pt height=21.18721440000001pt/> main elements: a **constraint** (for example <img src="/tex/d0b08429a8092571e5683ca898c7e0a0.svg?invert_in_darkmode&sanitize=true" align=middle width=69.12487889999998pt height=20.908638300000003pt/>), some **parameters** associated to the constraint, a **penalty** (possibly null) and **robust parameters** to use Huber or biweight losses only on the considered edge.

In addition to the edge definition, the nodes can be constrained to a range of values (means) in the inference. We also can specify starting and ending nodes. In our implementation the transitions are not time-dependent except for starting and ending steps.

The edges of the graph can be of type "null", "std",  "up", "down" or "abs" with the following meaning:

- "null" edge : there is no change-point introduced. We stay on the same segment. "null" corresponds to the constraint <img src="/tex/7a973d7d4f51344eadbf9e25b2c18d88.svg?invert_in_darkmode&sanitize=true" align=middle width=65.97801374999999pt height=22.465723500000017pt/>. The value does not change (or we have an exponential decay if <img src="/tex/3622f66f0df4e0212870fbe5f45e395b.svg?invert_in_darkmode&sanitize=true" align=middle width=70.8501783pt height=21.18721440000001pt/>);
- "std" edge : no contraint, the next segment can have any mean;
- "up" edge : the next segment has a greater mean (we can also force the size of the gap to be greater than a minimal value). "up" corresponds to the constraint <img src="/tex/3fb11313b3f76ebb5dc499c8caa916f3.svg?invert_in_darkmode&sanitize=true" align=middle width=57.61474949999999pt height=22.465723500000017pt/>;
- "down" edge : the next segment has a lower mean (wan also can force the size of the gap to be greater that a minimal value). "down" corresponds to the constraint <img src="/tex/10f9abe364a7546a343e9ca6d717ed34.svg?invert_in_darkmode&sanitize=true" align=middle width=57.61474949999999pt height=22.465723500000017pt/>;
- "abs" edge : the absolute value of the difference of two consecutive means is greater than a given parameter.


A nonnegative internal parameter can thus be associated to an edge (in "up", "down" and "abs""). The "null" edge refers to an absence of change-point. This edge can be used between different states to constraint segment lengths. Thus our package includes the segment neighborhood algorithm for which the number of change-points is fixed. More details on graph construction are given in the section [Graph construction](#gc).

### The non-constrained changepoint detection problem

The change-point vector <img src="/tex/ed535ca37308856c5f0c6d85b4d4c676.svg?invert_in_darkmode&sanitize=true" align=middle width=209.11492305pt height=27.91243950000002pt/> defines the <img src="/tex/33359de825e43daa97171e27f6558ae9.svg?invert_in_darkmode&sanitize=true" align=middle width=37.38576269999999pt height=22.831056599999986pt/> segments <img src="/tex/97a268c17395aace06ce389334ba5322.svg?invert_in_darkmode&sanitize=true" align=middle width=115.02097364999997pt height=24.65753399999998pt/>, <img src="/tex/59680c09e0e2d723d0bcf2005047b028.svg?invert_in_darkmode&sanitize=true" align=middle width=73.18587374999998pt height=22.831056599999986pt/> with fixed bounds <img src="/tex/370a29c873e3d269a6111aa219085d0b.svg?invert_in_darkmode&sanitize=true" align=middle width=44.697406049999984pt height=21.18721440000001pt/> and  <img src="/tex/afdb85da7c3e7b7c3d226050994dbf5f.svg?invert_in_darkmode&sanitize=true" align=middle width=63.70246739999999pt height=14.15524440000002pt/>. We use the set <img src="/tex/de3cb25372b9c18c99bc9f21a0efab11.svg?invert_in_darkmode&sanitize=true" align=middle width=269.7471512999999pt height=27.91243950000002pt/> to define the non-constrained minimal global cost given by

<p align="center"><img src="/tex/55e63ede9d3bc1ac8c64fa74af425a24.svg?invert_in_darkmode&sanitize=true" align=middle width=277.1488236pt height=49.315569599999996pt/></p>

where <img src="/tex/99751e94989c68f9be0f6aa442bc80d5.svg?invert_in_darkmode&sanitize=true" align=middle width=40.302373649999986pt height=22.831056599999986pt/> is a penalty parameter and <img src="/tex/a44ff4154fc3bc708e9e752a14051324.svg?invert_in_darkmode&sanitize=true" align=middle width=49.762892849999986pt height=24.65753399999998pt/> is the minimal cost over the segment <img src="/tex/add1478513cabbadcd5004323f01b74c.svg?invert_in_darkmode&sanitize=true" align=middle width=62.71697189999998pt height=24.65753399999998pt/>. The penalty <img src="/tex/8217ed3c32a785f0b5aad4055f432ad8.svg?invert_in_darkmode&sanitize=true" align=middle width=10.16555099999999pt height=22.831056599999986pt/> is understood as an additional cost when introducing a new segment. The argminimum of this quantity gives a vector <img src="/tex/fa1931c6895840344a67aa4df9cf3f59.svg?invert_in_darkmode&sanitize=true" align=middle width=15.782028899999991pt height=22.63846199999998pt/> containing the last position of each segment (if we do not consider <img src="/tex/370a29c873e3d269a6111aa219085d0b.svg?invert_in_darkmode&sanitize=true" align=middle width=44.697406049999984pt height=21.18721440000001pt/>). The quantity <img src="/tex/0ce745c733fc96f40b3901e008412ce9.svg?invert_in_darkmode&sanitize=true" align=middle width=21.121448699999988pt height=22.465723500000017pt/> is the solution of the non-constrained optimal partitioning method. 

In our setting, the cost <img src="/tex/db5f7b3e9934fbc5a2859d88c4ba84a3.svg?invert_in_darkmode&sanitize=true" align=middle width=9.614228249999991pt height=22.465723500000017pt/> is the result of the minimization of a cost function with additive property:

<p align="center"><img src="/tex/a2c69706ac0380a606c182ae3c52eea7.svg?invert_in_darkmode&sanitize=true" align=middle width=504.5946576pt height=40.51346475pt/></p>

with the argminimum defining the infered mean <img src="/tex/47b592a798cd56ccf668b67abad36a61.svg?invert_in_darkmode&sanitize=true" align=middle width=19.083998999999988pt height=14.15524440000002pt/> of the i+1-th segment <img src="/tex/97a268c17395aace06ce389334ba5322.svg?invert_in_darkmode&sanitize=true" align=middle width=115.02097364999997pt height=24.65753399999998pt/> with <img src="/tex/0794cee3e13634952ef9036a12e8fb9e.svg?invert_in_darkmode&sanitize=true" align=middle width=87.79779359999999pt height=24.65753399999998pt/>. Additivity of the cost (the <img src="/tex/11c596de17c342edeed29f489aa4b274.svg?invert_in_darkmode&sanitize=true" align=middle width=9.423880949999988pt height=14.15524440000002pt/> decomposition) is guaranteed as we will use costs deriving from a likelihood. <img src="/tex/11c596de17c342edeed29f489aa4b274.svg?invert_in_darkmode&sanitize=true" align=middle width=9.423880949999988pt height=14.15524440000002pt/> has in our package <img src="/tex/5dc642f297e291cfdde8982599601d7e.svg?invert_in_darkmode&sanitize=true" align=middle width=8.219209349999991pt height=21.18721440000001pt/> possible basis decompositions:

<p align="center"><img src="/tex/16ee85ae1690e8b380c396dd91218e25.svg?invert_in_darkmode&sanitize=true" align=middle width=700.2745497pt height=161.23476104999997pt/></p>

We can associate a current mean <img src="/tex/ce9c41bf6906ffd46ac330f09cacc47f.svg?invert_in_darkmode&sanitize=true" align=middle width=14.555823149999991pt height=14.15524440000002pt/> to each data point <img src="/tex/2b442e3e088d1b744730822d18e7aa21.svg?invert_in_darkmode&sanitize=true" align=middle width=12.710331149999991pt height=14.15524440000002pt/> and we write a cost on a segment <img src="/tex/add1478513cabbadcd5004323f01b74c.svg?invert_in_darkmode&sanitize=true" align=middle width=62.71697189999998pt height=24.65753399999998pt/> as a result of a constrained minimization:

<p align="center"><img src="/tex/4060746cfc34359c50eb0701ada25324.svg?invert_in_darkmode&sanitize=true" align=middle width=276.95500799999996pt height=48.7101186pt/></p>

so that we get another description of the objective function :

<p align="center"><img src="/tex/d069c19a2f5fdecd518c4a0f5c8f9dd7.svg?invert_in_darkmode&sanitize=true" align=middle width=354.07496355pt height=49.931608649999994pt/></p>

From this expression, it will be possible to impose some constraints on the vector of means <img src="/tex/07617f9d8fe48b4a7b3f523d6730eef0.svg?invert_in_darkmode&sanitize=true" align=middle width=9.90492359999999pt height=14.15524440000002pt/>.

### Graph-constrained problem

<img src="/tex/37ca9e919ee5a54f28ed9f70e6ed1277.svg?invert_in_darkmode&sanitize=true" align=middle width=68.78982494999998pt height=24.65753399999998pt/> is the indicator function. This previous approach can be generalized to complex constraints on consecutive means using a graph structure. We define the transition graph <img src="/tex/0553770832b41fddeb39cbfa76cf036a.svg?invert_in_darkmode&sanitize=true" align=middle width=17.90463839999999pt height=22.465723500000017pt/> as a directed acyclic graph with the following properties:

1. Nodes are indexed by time and state. <img src="/tex/6c8ced74f5b3bcbfbd287c91837d0b64.svg?invert_in_darkmode&sanitize=true" align=middle width=64.20836564999999pt height=24.65753399999998pt/> for node <img src="/tex/6c4adbc36120d62b98deef2a20d5d303.svg?invert_in_darkmode&sanitize=true" align=middle width=8.55786029999999pt height=14.15524440000002pt/> with time <img src="/tex/4f4f4e395762a3af4575de74c019ebb5.svg?invert_in_darkmode&sanitize=true" align=middle width=5.936097749999991pt height=20.221802699999984pt/> and state <img src="/tex/6f9bad7347b91ceebebd3ad7e6f6f2d1.svg?invert_in_darkmode&sanitize=true" align=middle width=7.7054801999999905pt height=14.15524440000002pt/>. The states are elements of the set <img src="/tex/4b924de9e932698c85d5af3791bea6bf.svg?invert_in_darkmode&sanitize=true" align=middle width=130.89004995pt height=24.65753399999998pt/>;

2. All the nodes have time <img src="/tex/4f4f4e395762a3af4575de74c019ebb5.svg?invert_in_darkmode&sanitize=true" align=middle width=5.936097749999991pt height=20.221802699999984pt/> in <img src="/tex/868738293e2285473ac5c2b2c7d7f7a8.svg?invert_in_darkmode&sanitize=true" align=middle width=62.83494359999998pt height=24.65753399999998pt/> except for the unique starting node <img src="/tex/ce927a5e91802f349da081165382e133.svg?invert_in_darkmode&sanitize=true" align=middle width=73.78989254999999pt height=24.65753399999998pt/> and the unique ending node <img src="/tex/f508dd141f08e86bda0a55c9e2815f37.svg?invert_in_darkmode&sanitize=true" align=middle width=121.96534844999998pt height=24.65753399999998pt/>, where <img src="/tex/53fe7f87db69e0ed1312d865111c131f.svg?invert_in_darkmode&sanitize=true" align=middle width=8.219209349999991pt height=24.65753399999998pt/> denotes an undefinite state;

3. Edges are transitions between "time-consecutive" vertices of type <img src="/tex/6c8ced74f5b3bcbfbd287c91837d0b64.svg?invert_in_darkmode&sanitize=true" align=middle width=64.20836564999999pt height=24.65753399999998pt/> and <img src="/tex/997bd6ae1076fb8e024546d2d8852578.svg?invert_in_darkmode&sanitize=true" align=middle width=101.74251614999997pt height=24.7161288pt/> which gives us the edge <img src="/tex/7d3f780dafd2bd5f98d6febd6cc41f02.svg?invert_in_darkmode&sanitize=true" align=middle width=82.92789945pt height=24.7161288pt/> for <img src="/tex/c8d3e8f08b1d18c8f9bee7f42c6d6797.svg?invert_in_darkmode&sanitize=true" align=middle width=88.86217889999998pt height=24.65753399999998pt/>;

4. Each edge <img src="/tex/8cd34385ed61aca950a6b06d09fb50ac.svg?invert_in_darkmode&sanitize=true" align=middle width=7.654137149999991pt height=14.15524440000002pt/> is associated to a function <img src="/tex/71d7882ef649a7537f2e3c662f98092a.svg?invert_in_darkmode&sanitize=true" align=middle width=137.5721589pt height=24.65753399999998pt/> constraining consecutive means <img src="/tex/87eefe082e181864d1321025c2705ecd.svg?invert_in_darkmode&sanitize=true" align=middle width=14.870715749999988pt height=14.15524440000002pt/> and <img src="/tex/1790d105076a15b61f11a9d6e6296025.svg?invert_in_darkmode&sanitize=true" align=middle width=31.51463699999999pt height=14.15524440000002pt/> by the indicator function <img src="/tex/c261ab74762f29919bc44aa0c519c8e4.svg?invert_in_darkmode&sanitize=true" align=middle width=176.09637045pt height=24.65753399999998pt/> (for example <img src="/tex/4919996dd46d961cb23ba7804cbf9d71.svg?invert_in_darkmode&sanitize=true" align=middle width=195.57111419999998pt height=24.65753399999998pt/>) and a penalty <img src="/tex/2111bbb895243afe1ffa4d502cfa57ed.svg?invert_in_darkmode&sanitize=true" align=middle width=46.493620799999995pt height=22.831056599999986pt/>.

The next table summarizes all the possible constraints encoded into the `gfpop` package.

| constraints | <img src="/tex/cc7c3637db1e725c477f85fd4d9903cd.svg?invert_in_darkmode&sanitize=true" align=middle width=109.26161729999998pt height=22.648391699999998pt/>, <img src="/tex/69e4724bb8055a140759e9aac47561fb.svg?invert_in_darkmode&sanitize=true" align=middle width=49.168495199999995pt height=26.17730939999998pt/> |
|---------:|-----:|
| null | <img src="/tex/99f6e4a4989dea6b5e336f8f326ca8b6.svg?invert_in_darkmode&sanitize=true" align=middle width=161.75505059999998pt height=24.65753399999998pt/> |
| std | <img src="/tex/f21b890a3381a0da31c3f229444dba9f.svg?invert_in_darkmode&sanitize=true" align=middle width=161.75505059999998pt height=24.65753399999998pt/> | 
| up | <img src="/tex/6de9c4956a1bef11d38879f6f1b94eae.svg?invert_in_darkmode&sanitize=true" align=middle width=177.90372269999997pt height=24.65753399999998pt/> |
| down | <img src="/tex/da6a35ac039e3397f4da159c56918967.svg?invert_in_darkmode&sanitize=true" align=middle width=177.90371609999997pt height=24.65753399999998pt/> |
| abs | <img src="/tex/84f0c45ef917380c122efe039401e111.svg?invert_in_darkmode&sanitize=true" align=middle width=201.100251pt height=24.65753399999998pt/>|


We define a path <img src="/tex/f26713629da615784507a3edf48fb24b.svg?invert_in_darkmode&sanitize=true" align=middle width=46.26634319999999pt height=22.465723500000017pt/> of the graph as a collection of <img src="/tex/014ae02bf2b677c4c225a7dfd00b8420.svg?invert_in_darkmode&sanitize=true" align=middle width=38.17727759999999pt height=21.18721440000001pt/> vertices <img src="/tex/54d840342c0d4e1f4770f2572852818c.svg?invert_in_darkmode&sanitize=true" align=middle width=89.99831279999998pt height=24.65753399999998pt/> with <img src="/tex/93729215cb2ecc9aea94ffd3ff1ed5bf.svg?invert_in_darkmode&sanitize=true" align=middle width=73.78989254999999pt height=24.65753399999998pt/> and <img src="/tex/f508dd141f08e86bda0a55c9e2815f37.svg?invert_in_darkmode&sanitize=true" align=middle width=121.96534844999998pt height=24.65753399999998pt/> and <img src="/tex/73eb8592e5cb867e19391571ae12bfb0.svg?invert_in_darkmode&sanitize=true" align=middle width=75.19396665pt height=24.65753399999998pt/> for <img src="/tex/985dbad85d61b07e704840368824ee09.svg?invert_in_darkmode&sanitize=true" align=middle width=88.86217889999998pt height=24.65753399999998pt/> and <img src="/tex/308736d40f48c7c0c2e28c99b956051e.svg?invert_in_darkmode&sanitize=true" align=middle width=44.77148444999999pt height=22.465723500000017pt/>. Morever, the path is made of <img src="/tex/3f18d8f60c110e865571bba5ba67dcc6.svg?invert_in_darkmode&sanitize=true" align=middle width=38.17727759999999pt height=21.18721440000001pt/> edges denoted <img src="/tex/892f838396f392c3a6c9da53214ef424.svg?invert_in_darkmode&sanitize=true" align=middle width=59.119196399999986pt height=14.15524440000002pt/> with <img src="/tex/b8e1a11c627257800b07edffe7fbe8a9.svg?invert_in_darkmode&sanitize=true" align=middle width=54.554976299999986pt height=22.831056599999986pt/>. A vector <img src="/tex/ea48643af203a0e7ad01bcafddbc8f82.svg?invert_in_darkmode&sanitize=true" align=middle width=49.99426409999999pt height=22.648391699999998pt/> verifies the path <img src="/tex/2ec6e630f199f589a2402fdf3e0289d5.svg?invert_in_darkmode&sanitize=true" align=middle width=8.270567249999992pt height=14.15524440000002pt/> if for all <img src="/tex/875f0d0d37d49d24b63ba045d9e7491a.svg?invert_in_darkmode&sanitize=true" align=middle width=117.1725786pt height=24.65753399999998pt/>, we have <img src="/tex/a9f7bf6538d048b2d8fae862a45abeb5.svg?invert_in_darkmode&sanitize=true" align=middle width=117.99692024999999pt height=24.65753399999998pt/> (valid constraint). We write <img src="/tex/5a2bcf863517582f776c5598b8a1c6b2.svg?invert_in_darkmode&sanitize=true" align=middle width=30.960925049999993pt height=24.65753399999998pt/> to say that the vector <img src="/tex/07617f9d8fe48b4a7b3f523d6730eef0.svg?invert_in_darkmode&sanitize=true" align=middle width=9.90492359999999pt height=14.15524440000002pt/> verifies the path <img src="/tex/2ec6e630f199f589a2402fdf3e0289d5.svg?invert_in_darkmode&sanitize=true" align=middle width=8.270567249999992pt height=14.15524440000002pt/>. The formulation of our graph-constrained problem is then the following:

<p align="center"><img src="/tex/25a089a47a987094b647297062ee8766.svg?invert_in_darkmode&sanitize=true" align=middle width=310.64316525pt height=51.74449499999999pt/></p>


In the `gfpop` package, the edges <img src="/tex/46ed17799525427705742ef69fb71fd9.svg?invert_in_darkmode&sanitize=true" align=middle width=53.356130849999985pt height=24.7161288pt/> for <img src="/tex/2da0ae8ede87918239efbfd469567202.svg?invert_in_darkmode&sanitize=true" align=middle width=88.86217889999998pt height=24.65753399999998pt/> are not time-dependent. In this case, we redefine a graph <img src="/tex/9f1efdaeed741c2a46969cdb880cabbf.svg?invert_in_darkmode&sanitize=true" align=middle width=76.34688434999998pt height=24.65753399999998pt/>  by the label of its vertices <img src="/tex/768af24250bc547f14842daac9e5f22b.svg?invert_in_darkmode&sanitize=true" align=middle width=132.94492035pt height=24.65753399999998pt/> and its set of oriented edges <img src="/tex/526e9ac1fb94812f38332e979abb135a.svg?invert_in_darkmode&sanitize=true" align=middle width=81.57506609999999pt height=22.465723500000017pt/> with <img src="/tex/88879fa4b41455559f35230fb5d4dcf1.svg?invert_in_darkmode&sanitize=true" align=middle width=102.85925264999997pt height=24.7161288pt/> defining an edge going from node <img src="/tex/6f9bad7347b91ceebebd3ad7e6f6f2d1.svg?invert_in_darkmode&sanitize=true" align=middle width=7.7054801999999905pt height=14.15524440000002pt/> to node <img src="/tex/675c2f5707a1fa7050c12adc1872ba32.svg?invert_in_darkmode&sanitize=true" align=middle width=11.49544109999999pt height=24.7161288pt/> with indicator function <img src="/tex/629cbe1de238666b8fbb23255f22ea6b.svg?invert_in_darkmode&sanitize=true" align=middle width=13.46296874999999pt height=22.465723500000017pt/> and penalty <img src="/tex/8a51f4aa3be26df848de9e1e89f6a4a6.svg?invert_in_darkmode&sanitize=true" align=middle width=15.534883649999992pt height=22.831056599999986pt/>. 

In most applications, the set of edges always contains edges of type <img src="/tex/f308eb98d3204881e73de0e34ff40cfa.svg?invert_in_darkmode&sanitize=true" align=middle width=35.502276149999986pt height=24.65753399999998pt/> for all <img src="/tex/e6bfc1a9751a1dda5621c949c34e793c.svg?invert_in_darkmode&sanitize=true" align=middle width=41.03865479999999pt height=22.465723500000017pt/> with indicator function <img src="/tex/34c893ff78debf4a9f1763de90d3a70a.svg?invert_in_darkmode&sanitize=true" align=middle width=137.6560581pt height=24.65753399999998pt/> as it corresponds to an absence of constraint on segment lengths. 

The main algorithm of the package, the `gfpop` function, returns the change-point vector <img src="/tex/fa1931c6895840344a67aa4df9cf3f59.svg?invert_in_darkmode&sanitize=true" align=middle width=15.782028899999991pt height=22.63846199999998pt/> defined as <img src="/tex/0b35e94a3c1848e17627a2ad4d8fe4d8.svg?invert_in_darkmode&sanitize=true" align=middle width=196.18626389999997pt height=24.65753399999998pt/>, with <img src="/tex/48162c2e6d8e19213e506579c4cccc51.svg?invert_in_darkmode&sanitize=true" align=middle width=86.28446969999999pt height=24.65753399999998pt/> the argminimum of <img src="/tex/05839a0807c0b4eb51ed66d3f64b7ac0.svg?invert_in_darkmode&sanitize=true" align=middle width=77.228118pt height=24.65753399999998pt/> in <img src="/tex/63dee4f83c8372c83bbf2389f1c3a964.svg?invert_in_darkmode&sanitize=true" align=middle width=53.455350299999985pt height=24.65753399999998pt/> and <img src="/tex/ac4f2e6484bfba9965ba05d82ffb5cce.svg?invert_in_darkmode&sanitize=true" align=middle width=91.16641709999999pt height=19.1781018pt/>. we also give the possibility to restrict the set of valid paths by imposing a starting and/or an ending ,node and contraint the range of infered means, replacing <img src="/tex/ea48643af203a0e7ad01bcafddbc8f82.svg?invert_in_darkmode&sanitize=true" align=middle width=49.99426409999999pt height=22.648391699999998pt/> by <img src="/tex/502c144cca7929a6af3ab94ff27edc73.svg?invert_in_darkmode&sanitize=true" align=middle width=72.0565956pt height=24.65753399999998pt/> in the definition of <img src="/tex/63dee4f83c8372c83bbf2389f1c3a964.svg?invert_in_darkmode&sanitize=true" align=middle width=53.455350299999985pt height=24.65753399999998pt/>.

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
## [1]  100  300  500  799 1000
## 
## states
## [1] "Dw" "Up" "Dw" "Up" "Dw"
## 
## forced
## [1] 0 0 0 0
## 
## parameters
## [1] 1.031384 2.130883 1.047181 2.961262 1.005459
## 
## globalCost
## [1] 984.3485
## 
## attr(,"class")
## [1] "gfpop"
```

The vector `changepoints` gives the last index of each segment. It always ends with the length of the vector `vectData`.

The vector `states` contains the states in which lies each mean. The length of this vector is the same as the length of `changepoint`.

The vector `forced` is a boolean vector. A forced element means that two consecutive means have been forced to satisfy the constraint. For example, the "up" edge with parameter <img src="/tex/3e18a4a28fdee1744e5e3f79d13b9ff6.svg?invert_in_darkmode&sanitize=true" align=middle width=7.11380504999999pt height=14.15524440000002pt/> is forced if <img src="/tex/541dd8a10e3c8d8e898e9bf1d9b2dd56.svg?invert_in_darkmode&sanitize=true" align=middle width=105.57833054999998pt height=19.1781018pt/>.

The vector `parameters` contains the inferred means/parameters of the successive segments. 
 
The number `globalCost` is equal to <img src="/tex/e96673f14f7d2280cfd060ac442da0b9.svg?invert_in_darkmode&sanitize=true" align=middle width=45.482259899999995pt height=24.65753399999998pt/>, the overall cost of the segmented data. 


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
## [1]  211  383  713 1000
## 
## states
## [1] "Iso" "Iso" "Iso" "Iso"
## 
## forced
## [1] 0 0 0
## 
## parameters
## [1] 0.2898816 1.2191021 2.1403436 2.8936415
## 
## globalCost
## [1] 1025.616
## 
## attr(,"class")
## [1] "gfpop"
```


In this example, we use in `gfpop` function a robust biweight gaussian cost with `K = 1` and the `min` parameter in order to infer means greater than `0.5`.

### Fixed number of change-points

This algorithm is called segment neighborhood in the change-point litterature. In this example, we fixed the number of segments at <img src="/tex/5dc642f297e291cfdde8982599601d7e.svg?invert_in_darkmode&sanitize=true" align=middle width=8.219209349999991pt height=21.18721440000001pt/> with an isotonic constraint. The graph contains two "up" edges with no cycling.


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
## [1]  298  600 1000
## 
## states
## [1] "0" "1" "2"
## 
## forced
## [1] 0 0
## 
## parameters
## [1] 0.5195712 1.7562499 2.6745127
## 
## globalCost
## [1] 1061.791
## 
## attr(,"class")
## [1] "gfpop"
```


### Robust up-down with constrained starting and ending states

In presence of outliers we need a robust loss (biweight). We can also force the starting and ending state and a minimal gap between the means (here equal to `1`)

```r
n <- 1000
mydata <- dataGenerator(n, c(0.1,0.3,0.5,0.8,1), c(0,1,0,1,0), sigma = 1) + 5*(rbinom(n, 1, 0.05)) - 5*(rbinom(n, 1, 0.05))
beta <- 2*log(n)
myGraph <- graph(
  Edge(0, 1, "up", penalty = beta, gap = 1, K = 3),
  Edge(1, 0, "down", penalty = beta, gap = 1, K = 3),
  Edge(0, 0, "null"),
  Edge(1, 1, "null"),
  StartEnd(start = 0, end = 0))
gfpop(data =  mydata, mygraph = myGraph, type = "mean")
```

```
## changepoints
## [1]  138  305  493  818 1000
## 
## states
## [1] "0" "1" "0" "1" "0"
## 
## forced
## [1] 1 1 0 1
## 
## parameters
## [1] 0.01850678 1.01850678 0.01850678 1.02272647 0.02272647
## 
## globalCost
## [1] 1079.697
## 
## attr(,"class")
## [1] "gfpop"
```

If we skip all these constraints and use a standard fpop algorithm, the result is the following


```r
myGraphStd <- graph(penalty = 2*log(n), type = "std")
gfpop(data =  mydata, mygraph = myGraphStd, type = "mean")
```

```
## changepoints
##  [1]   17   18   31   32   33   48   51   52   65   66   88   89  111  112  129  130  149  150  174  176  194  195  205  206  220  221  227  229
## [29]  284  285  340  341  358  359  407  408  424  425  428  442  449  450  473  474  479  480  486  487  510  512  549  550  551  565  566  570
## [57]  571  574  590  591  595  596  613  614  627  629  639  640  701  702  778  779  796  797  818  847  848  873  875  878  882  928  929  933
## [85]  934  953  955  957  958  962  963 1000

## states
##  [1] "Std" "Std" "Std" "Std" "Std" "Std" "Std" "Std" "Std" "Std" "Std" "Std" "Std" "Std" "Std" "Std" "Std" "Std" "Std" "Std" "Std" "Std" "Std"
## [24] "Std" "Std" "Std" "Std" "Std" "Std" "Std" "Std" "Std" "Std" "Std" "Std" "Std" "Std" "Std" "Std" "Std" "Std" "Std" "Std" "Std" "Std" "Std"
## [47] "Std" "Std" "Std" "Std" "Std" "Std" "Std" "Std" "Std" "Std" "Std" "Std" "Std" "Std" "Std" "Std" "Std" "Std" "Std" "Std" "Std" "Std" "Std"
## [70] "Std" "Std" "Std" "Std" "Std" "Std" "Std" "Std" "Std" "Std" "Std" "Std" "Std" "Std" "Std" "Std" "Std" "Std" "Std" "Std" "Std" "Std" "Std"
## 
## forced
## [1] 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ## 0 0
## [72] 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
## 
## parameters
##  [1]  0.37883666 -5.60737627  0.14265729 -4.92600986  5.76837573  0.46081049 -2.07025038  6.85184085 -0.08186373  6.07962475  0.16536163
## [12] -6.72115878  0.11779333 -6.03851292  0.38897505  6.36318506  1.22328631 -4.84795400  1.37996121  4.93788071  0.28409636  6.30089144
## [23]  0.79739606 -5.77855389  1.43330020  6.66018030  0.79843808 -4.43936467  1.17189869  7.11750913  0.30445508  5.63498352  0.07145084
## [34] -6.13223532 -0.10524605  6.53335020  0.29761434  4.42059093 -3.52301587 -0.41869828  1.60922633 -6.59086009  0.07051877  4.87531990
## [45] -2.14802188  7.64909748 -0.03687840 -7.05437136  0.62557676  5.83147112  1.26743138  6.64006018 -4.30260774  0.94530626 -4.62237889
## [56]  1.24394216  8.62245556 -1.51272772  1.14038368  6.62490433  0.97797273  6.91505783  0.45741709 -5.00510212  0.94368809 -4.22010595
## [67]  0.82243663  6.88636862  1.06272055 -4.38864051  0.97078589 -4.92447920  1.46891367 -3.97447863  1.27138347 -0.25178268  6.50770084
## [78]  0.66843446  4.86486110 -0.36798139 -3.52991185  0.28266676 -5.37296951  1.26764318  6.32132912 -0.29253887  5.10483206 -0.24247171
## [89]  7.11097890 -0.62541613 -5.62987116  0.27378843
## 
## globalCost
## [1] 1747.779
## 
## attr(,"class")
## [1] "gfpop"
```


### abs edge

With a unique "abs" edge, we impose a difference between the means of size at least 1.  

```r
n <- 10000
mydata <- dataGenerator(n, c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1), c(0, 1, 0, 2, 1, 2, 0, 1, 0, 1), sigma = 0.5)
beta <- 2*log(n)
myGraph <- graph(
  Edge(0, 0,"abs", penalty = beta, gap = 1),
  Edge(0, 0,"null"))
gfpop(data =  mydata, mygraph = myGraph, type = "mean")
```

```
## changepoints
##  [1]  1000  1999  3000  4000  5000  6000  7000  8002  9000 10000
## 
## states
## [1] "0" "0" "0" "0" "0" "0" "0" "0" "0" "0"
## 
## forced
## [1] 0 1 0 0 1 0 1 1 0
## 
## parameters
## [1] 0.008077166 1.023073852 0.023073852 1.994958839 0.990916305 1.990916305 0.001495140 1.001495140 0.001495140 1.024786144
## 
## globalCost
## [1] 2494.115
## 
## attr(,"class")
## [1] "gfpop"

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
## [1] 0 0 0
## 
## parameters
## [1] 0.0052340381 0.0003157659 0.0004575503 0.0198921241
## 
## globalCost
## [1] 935.4619
## 
## attr(,"class")
## [1] "gfpop"
```

and we plot the result 

```r
gamma <- 0.966
len <- diff(c(0, g<img src="/tex/46024bfe5198e5e1ff626011648c9edf.svg?invert_in_darkmode&sanitize=true" align=middle width=558.8873268pt height=24.65753399999998pt/>parameters[i]*c(1, cumprod(rep(1/gamma,len[i]-1))))}
signal <- rev(signal)

ylimits <- c(min(mydata), max(mydata))
plot(mydata, type ='p', pch ='+', ylim = ylimits)
par(new = TRUE)
plot(signal, type ='l', col = 4, ylim = ylimits, lwd = 3)
```

![plot](figure/expDecay.png)


<a id="gc"></a>

## Graph construction

In the `gfpop` package, graphs are represented by a dataframe with 9 features and build with the R functions `Edge`, `Node`, `StartEnd` and `graph`.


```r
emptyGraph <- graph()
emptyGraph
```

```
## [1] state1    state2    type      penalty   parameter K         a         min       max
## <0 rows> (or 0-length row.names)
```

`state1` is the starting node of an edge, `state2` its ending node. `type` is one of the available edge type ("null", "std", "up", "down", "abs"). `penalty` is a nonnegative parameter: the additional cost <img src="/tex/3d13090ef3ed1448f3c4dc166d06ab4d.svg?invert_in_darkmode&sanitize=true" align=middle width=13.948864049999989pt height=22.831056599999986pt/> to consider when we move within the graph using a edge (or stay on the same node). `parameter` is annother nonnegative parameter, a characteristics of the edge, depending of its type (it is a decay if type is "null" and a gap otherwise). `K` and `a` are robust parameters. `min` and `max` are used to constrain the rang of value for the node parameter.

We add edges into a graph as follows

```r
myGraph <- graph(
  Edge("E1", "E1", "null"),
  Edge("E1", "E2", "down", 3.1415, gap = 1.5)
)
myGraph
```

```
##  state1 state2 type parameter penalty   K   a min max
## 1     E1     E1 null       1.0       0 Inf Inf  NA  NA
## 2     E1     E2 down       1.5       0 Inf Inf  NA  NA
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
## state1 state2  type parameter  penalty   K   a min max
## 1     Dw     Dw  null         1  0.00000 Inf Inf  NA  NA
## 2     Up     Up  null         1  0.00000 Inf Inf  NA  NA
## 3     Dw     Up    up         1 13.81551 Inf Inf  NA  NA
## 4     Dw     Dw  down         0 13.81551 Inf Inf  NA  NA
## 5     Up     Dw  down         0 13.81551 Inf Inf  NA  NA
## 6     Dw   <NA> start        NA       NA  NA  NA  NA  NA
## 7     Dw   <NA>   end        NA       NA  NA  NA  NA  NA
```

Some graphs are often used: they are defined by default in the `graph` function. To use these graphs, we specify a string `type` equal to "std", "isotonic", "updown" or "relevant".
For example,


```r
myGraphIso <- graph(penalty = 12, type = "isotonic")
myGraphIso
```

```
##   state1 state2 type parameter penalty   K   a min max
## 1    Iso    Iso null         1       0 Inf Inf  NA  NA
## 2    Iso    Iso   up         0      12 Inf Inf  NA  NA
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
##   state1 state2 type parameter penalty   K   a min max
## 1     Up     Up   up         0  3.1415 Inf Inf  NA  NA
## 2     Up     Up null         1  0.0000 Inf Inf  NA  NA
## 3     Up     Up node        NA      NA  NA  NA   0   1
```

<a id="suppl"></a>

## Supplementary R functions

### Data generator function

the `dataGenerator` function is used to simulate <img src="/tex/55a049b8f161ae7cfeb0197d75aff967.svg?invert_in_darkmode&sanitize=true" align=middle width=9.86687624999999pt height=14.15524440000002pt/> data-points from a distribution of `type` equal to “mean”, “poisson”, “exp”, “variance” or “negbin”. Standard deviation parameter `sigma` and decay `gamma` are specific to the Gaussian mean model. `size` is linked to the R “rnbinom” function from R stats package.

### Standard deviation estimation

We often need to estimate the standard deviation from the observed data to normalize the data or choose the edge penalties. The `sdDiff` returns such an estimation with the default HALL method [Hall et al., 1990] well suited for time series with change-points.

<a id="gfpop"></a>

## More on the gfpop function and its C++ structure

To be done.

[Back to Top](#top)
