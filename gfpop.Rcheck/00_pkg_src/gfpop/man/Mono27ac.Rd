\name{Mono27ac}
\alias{Mono27ac}
\docType{data}
\title{
  A small ChIP-seq data set in which peaks can be found using PeakSegFPOP
}
\description{
  The data come from an H3K27ac ChIP-seq experiment which was aligned to
  the human reference genome (hg19), aligned read counts were used to
  produce the coverage data; looking at these data in a genome browser
  was used to produce the labels.
}
\usage{data("Mono27ac")}
\format{
A list of 2 data.tables: coverage has 4 columns (chrom, chromStart,
chromEnd, count=number of aligned reads at each position on
chrom:chromStart-chromEnd); labels has 4 columns (chrom, chromStart,
chromEnd, annotation=label at chrom:chromStart-chromEnd).
}
\source{
  https://github.com/tdhock/feature-learning-benchmark, prob.dir=
  H3K27ac-H3K4me3_TDHAM_BP/samples/Mono1_H3K27ac/S001YW_NCMLS/problems/chr11:60000-580000
}
\keyword{datasets}
