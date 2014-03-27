\name{calPara}
\alias{calPara}
\title{Calculate parameters for generalized Poisson test (GPTest)}
\usage{calPara(object)}
\description{Calculate a set of parameters from normalized counts table before \code{\link{GPTest}}}
\value{An ABSDataSet object with absolute differences, basemean, mean of each group, variance, 
log2 of foldchange and mean of each group, named as 'absD', 'baseMean', 'Amean', 'Bmean', 
'Variance' and 'foldChange', respectively. Use the \code{\link{results}} to get access it}
\note{This function should run after \code{\link{normalFactors}} or providing size factors.}
\arguments{\item{object}{an \code{\link{ABSDataSet}} object.}}
\examples{data(simuN5)
obj <- ABSDataSet(counts=simuN5$counts, groups=factor(simuN5$groups))
obj <- normalFactors(obj)
obj <- calPara(obj)
plotDifftoBase(obj)}
