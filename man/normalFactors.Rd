\name{normalFactors}
\alias{normalFactors}
\title{Estimating size factors from the reads count table}
\usage{normalFactors(object)}
\description{Function for esitmating size factors}
\details{Given a matrix of count data, this function esitmates the size
factors by selected method.
It provides four different methods for normalizing according to user-defined size factors,
total reads, up quantile (75\%) or DESeq (See \code{\link{estimateSizeFactorsForMatrix}}).}
\value{an ABSDataSet object with the estimates size factors, one element per column. Use the \code{\link{sizeFactors}}
to show it.}
\arguments{\item{object}{an ABSDataSet object with element of 'counts' and 'normMethod', see the constructor functions
\code{\link{ABSDataSet}}.}}
\examples{data(simuN5)
obj <- ABSDataSet(counts=simuN5$counts, groups=factor(simuN5$groups))
obj <- normalFactors(obj)
sizeFactors(obj)
}
