<a id="top"></a>

[![Build
Status](https://travis-ci.com/vrunge/gfpop.svg?branch=master)](https://travis-ci.com/vrunge/gfpop)
[![](https://img.shields.io/badge/docs-vignettes-blue.svg)](https://github.com/vrunge/gfpop)

<!-- 
%\VignetteEngine{knitr::rmarkdown} 
%\VignetteIndexEntry{An Introduction to gfpop}
-->

# gfpop Vignette

### Vincent Runge

#### LaMME, Evry University

### March 2, 2022

> [Quick Start](#qs)

> [Some examples](#se)

> [Graph construction](#gc)

> [Supplementary R functions](#suppl)

> [More on the gfpop function and its C++ structure](#gfpop)

## Quick Start

we present a basic use of the main functions of the `gfpop` package.
More details about optional arguments are given in later sections.

We install the package from Github:

    #devtools::install_github("vrunge/gfpop")
    library(gfpop)

Details about the optimization problem solved by gfpop can be found [here](https://arxiv.org/abs/2002.03646)

We simulate some univariate gaussian data (`n = 1000` points) with
relative change-point positions `0.1, 0.3, 0.5, 0.8, 1` and means
`1, 2, 1, 3, 1` with a variance equal to `1`.

    n <- 1000
    myData <- dataGenerator(n, c(0.1,0.3,0.5,0.8,1), c(1,2,1,3,1), sigma = 1)

We define the graph of constraints to use for the dynamic programming
algorithm. A simple case is the up-down constraint with a penalty here
equal to a classic `2 log(n)`.

    myGraph <- graph(penalty = 2*log(n), type = "updown")

The gfpop function gives the result of the segmentation using `myData`
and `myGraph` as parameters. We choose a gaussian cost.

    gfpop(data = myData, mygraph = myGraph, type = "mean")

    ## $changepoints
    ## [1]  100  299  500  800 1000
    ## 
    ## $states
    ## [1] "Dw" "Up" "Dw" "Up" "Dw"
    ## 
    ## $forced
    ## [1] FALSE FALSE FALSE FALSE
    ## 
    ## $parameters
    ## [1] 0.9781741 2.0585601 1.0311449 2.9681224 0.9012283
    ## 
    ## $globalCost
    ## [1] 1048.816
    ## 
    ## attr(,"class")
    ## [1] "gfpop" "mean"

The vector `changepoints` gives the last index of each segment. It
always ends with the length of the vector `vectData`.

The vector `states` contains the states in which lies each mean. The
length of this vector is the same as the length of `changepoint`.

The vector `forced` is a boolean vector. A forced element means that two
consecutive means have been forced to satisfy the constraint. For
example, the “up” edge with parameter *c* is forced if
*m*<sub>*i* + 1</sub> − *m*<sub>*i*</sub> = *c*.

The vector `parameters` contains the inferred means/parameters of the
successive segments.

The number `globalCost` is equal to *Q*<sub>*n*</sub>(𝒢), the overall
cost of the segmented data.

<a id="se"></a>

## Some examples

### Isotonic regression

The isotonic regression infers a sequence of nondecreasing means.

    n <- 1000
    mydata <- dataGenerator(n, c(0.1, 0.2, 0.3, 0.4, 0.6, 0.8, 1), c(0, 0.5, 1, 1.5, 2, 2.5, 3), sigma = 1)
    myGraphIso <- graph(penalty = 2*log(n), type = "isotonic")
    gfpop(data =  mydata, mygraph = myGraphIso, type = "mean")

    ## $changepoints
    ## [1]   91  215  400  616 1000
    ## 
    ## $states
    ## [1] "Iso" "Iso" "Iso" "Iso" "Iso"
    ## 
    ## $forced
    ## [1] FALSE FALSE FALSE FALSE
    ## 
    ## $parameters
    ## [1] -0.2566487  0.5429243  1.2257508  2.0048839  2.7082123
    ## 
    ## $globalCost
    ## [1] 869.5067
    ## 
    ## attr(,"class")
    ## [1] "gfpop" "mean"

In this example, we use in `gfpop` function a robust biweight gaussian
cost with `K = 1` and the `min` parameter in order to infer means
greater than `0.5`.

### Fixed number of change-points

This algorithm is called segment neighborhood in the change-point
litterature. In this example, we fixed the number of segments at 3 with
an isotonic constraint. The graph contains two “up” edges with no
cycling.

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

    ## $changepoints
    ## [1]  301  601 1000
    ## 
    ## $states
    ## [1] "0" "1" "2"
    ## 
    ## $forced
    ## [1] FALSE FALSE
    ## 
    ## $parameters
    ## [1] 0.4525906 1.9078955 2.7251198
    ## 
    ## $globalCost
    ## [1] 1071.464
    ## 
    ## attr(,"class")
    ## [1] "gfpop" "mean"

### Robust up-down with constrained starting and ending states

In presence of outliers we need a robust loss (biweight). We can also
force the starting and ending state and a minimal gap between the means
(here equal to `1`)

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

    ## $changepoints
    ## [1]  104  300  507  802 1000
    ## 
    ## $states
    ## [1] "Dw" "Up" "Dw" "Up" "Dw"
    ## 
    ## $forced
    ## [1] FALSE FALSE  TRUE FALSE
    ## 
    ## $parameters
    ## [1] -0.06055366  1.02881347 -0.01852851  0.98147149 -0.03055950
    ## 
    ## $globalCost
    ## [1] 1093.429
    ## 
    ## attr(,"class")
    ## [1] "gfpop" "mean"
