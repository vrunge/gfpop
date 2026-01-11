\name{neuroSpike}
\alias{neuroSpike}
\docType{data}
\title{
  Neuro spike train data
}
\description{
  The goal in these data s to detect the spikes = abrupt up
  changes in calcium concentration, which correspond to neurons firing
  in laboratory animals.
}
\usage{data("neuroSpike")}
\format{
  A data frame with 20000 observations on the following 3 variables.
  \describe{
    \item{\code{seconds}}{time of measurement}
    \item{\code{calcium}}{noisy data measured}
    \item{\code{AR1}}{autoregressive model mean, for comparison}
  }
}
\source{
  Spikefinder challenge, dataset number 7, cell number 13, subset
  1:20000, 100 fps, for AR1: 0.98 exp decay parameter, 3.6 penalty parameter.
}
\keyword{datasets}
