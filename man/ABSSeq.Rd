\name{ABSSeq}
\alias{ABSSeq}
\title{Differential expression analysis based on the absolute difference.}
\usage{ABSSeq(object, perPrior=0.3, dcut=0, quiet=FALSE)}
\description{This function performs a default analysis by calling, in order, the functions:
\code{\link{normalFactors}},
\code{\link{calPara}},
\code{\link{GPTest}}.}
\details{The differential expression analysis uses a generalized generalized Poisson model of the form:
\deqn{ GP(D=k|\lambda,\theta)=(\lambda/(1+\lambda * \theta)^k * (1+\theta * k)^(k-1)/k! * exp(-\lambda * (1+\theta * k) / (1+\lambda * \theta))}
for \eqn{ \lambda > 0, \theta >= 0, k=0, 1, 2, \ldots}.}
\value{a table with return a table with elements:
\item{Amean and Bmean}{mean of reads count for group A and B} 
\item{foldChange}{log 2 of fold-change, B vs. A} 
\item{pvalue}{pvalue from GP model} 
\item{adj.pvalue}{adjuested p-value used BH method}
}
\author{Wentao Yang}
\references{Wentao Yang, Philip Rosenstiel & Hinrich Schulenburg: ABSSeq: a new RNA-Seq analysis method based on absolute expression differences and generalized Poisson model}
\arguments{\item{object}{an \code{\link{ABSDataSet}} object, contains the reads count matrix, groups and normalization method.}
\item{perPrior}{default is 0.30, see \code{\link{GPTest}} for description.}
\item{dcut}{default is 0, see \code{\link{GPTest}} for description.}
\item{quiet}{default is FALSE, whether to print messages at each step}}
\examples{data(simuN5)
obj <- ABSDataSet(counts=simuN5$counts, groups=factor(simuN5$groups))
res <- ABSSeq(obj)
head(res)}
