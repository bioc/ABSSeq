\name{GPTest}
\alias{GPTest}
\title{Testing the differential expression}
\usage{GPTest(object, pPrior=0.3, dcut=0)}
\description{Fitting generalized Poisson model and calculating p-value for each gene}
\details{This function firstly fits the generalized Possion model used absolute differences between two
groups, then calculates the pvalue for each gene and finally adjusts the pvalues via BH method.}
\value{an \code{\link{ABSDataSet}} object with additional elements: pvalue and adjusted p-value,
denoted by pvalue and adj-pvalue, respectively. Use the \code{\link{results}} to get access it.}
\note{This function should run after \code{\link{calPara}}}
\arguments{\item{object}{an \code{\link{ABSDataSet}} object.}
\item{pPrior}{the parameter for estimating rates of generalized Possion model for each gene, 
default is 0.30. The predefined parameter for estimating GP rate according to expression level of each gene. It is estimated by reducing the 
type I error rate in a acceptable level, i.e. 10\% for perPrior=0.30 or 5\% for perPrior=0.40. For details, please see the references,}
\item{dcut}{the parameter for fitting GP model, default is 0.}}
\examples{data(simuN5)
obj <- ABSDataSet(counts=simuN5$counts, groups=factor(simuN5$groups))
obj <- normalFactors(obj)
obj <- calPara(obj)
obj <- GPTest(obj)
head(results(obj))}
