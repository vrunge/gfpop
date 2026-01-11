\name{profile614chr2}
\alias{profile614chr2}
\docType{data}
\title{
 Labeled DNA copy number profile
}
\description{
  This DNA copy number data set has changepoint labels that are impossible to all
  predict correctly using the maximum likelihood gaussian changepoint model.
}
\usage{data("profile614chr2")}
\format{
Named list of two elements: probes is a data.table with 153663 noisy
logratio observations in which we would like to detect changepoints;
labels is a data.table with 8 rows that define regions with 1 or 0 breakpoints.
}
\source{
  http://members.cbio.mines-paristech.fr/~thocking/neuroblastoma/test-data-01/614-regions.csv.gz
  http://members.cbio.mines-paristech.fr/~thocking/neuroblastoma/test-data-01/614-probes.csv.gz
}
\references{
  Labels created with SegAnnDB software, Hocking et al, Bioinformatics 2014. https://www.ncbi.nlm.nih.gov/pubmed/24493034
}
\keyword{datasets}
