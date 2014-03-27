\name{estimateSizeFactorsForMatrix}
\alias{estimateSizeFactorsForMatrix}
\title{Low-level function to estimate size factors with robust regression.}
\usage{estimateSizeFactorsForMatrix(counts, locfunc=median)}
\description{This function is borrowed from DESeq.}
\details{Given a matrix or data frame of count data, this function estimates the size
factors as follows: Each column is divided by the geometric means of the
rows. The median (or, if requested, another location estimator) of these
ratios (skipping the genes with a geometric mean of zero) is used as the size
factor for this column. Typically, you will not call this function directly.}
\value{a vector with the estimates size factors, one element per column}
\author{Simon Anders}
\references{Simon Anders, Wolfgang Huber: Differential expression analysis for sequence count data. Genome Biology 11 (2010) R106, \url{http://dx.doi.org/10.1186/gb-2010-11-10-r106}}
\arguments{\item{counts}{a matrix or data frame of counts, i.e., non-negative integer
values}
\item{locfunc}{a function to compute a location for a sample. By default, the
median is used.}}
\examples{data(simuN5)
dat=simuN5
estimateSizeFactorsForMatrix(dat$counts)}
