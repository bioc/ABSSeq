\docType{methods}
\name{maxRates}
\title{Accessors for the 'maxRates' slot of a ABSDataSet object.}
\usage{\S4method{maxRates}{ABSDataSet}(object)

\S4method{maxRates}{ABSDataSet,numeric}(object)<-value}
\description{Accessor functions for the 'maxRates' slot of a ABSDataSet
object.}
\details{The 'maxRates' is the upper bound of rate for baseline of counts difference esitimation.}
\alias{maxRates}
\alias{maxRates<-}
\alias{maxRates,ABSDataSet-method}
\alias{maxRates<-,ABSDataSet,numeric-method}
\seealso{\code{\link{callParameter}},\code{\link{ABSDataSet}}}
\arguments{\item{object}{an \code{ABSDataSet} object.}
\item{value}{a positive numeric object}}
\examples{data(simuN5)
obj <- ABSDataSet(counts=simuN5$counts, groups=factor(simuN5$groups))
maxRates(obj)
maxRates(obj) <- 0.4
maxRates(obj)}