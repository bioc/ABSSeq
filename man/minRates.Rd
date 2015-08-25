\docType{methods}
\name{minRates}
\title{Accessors for the 'minRates' slot of a ABSDataSet object.}
\usage{\S4method{minRates}{ABSDataSet}(object)

\S4method{minRates}{ABSDataSet,numeric}(object)<-value}
\description{Accessor functions for the 'minRates' slot of a ABSDataSet
object.}
\details{The 'minRates' is the lower bound of rate for baseline of counts difference esitimation.}
\alias{minRates}
\alias{minRates<-}
\alias{minRates,ABSDataSet-method}
\alias{minRates<-,ABSDataSet,numeric-method}
\seealso{\code{\link{callParameter}},\code{\link{ABSDataSet}}}
\arguments{\item{object}{an \code{ABSDataSet} object.}
\item{value}{a positive numeric object}}
\examples{data(simuN5)
obj <- ABSDataSet(counts=simuN5$counts, groups=factor(simuN5$groups))
minRates(obj)
minRates(obj) <- 0.3
minRates(obj)}