\docType{methods}
\name{minimalDispersion}
\title{Accessors for the 'minDispersion' slot of a ABSDataSet object.}
\usage{\S4method{minimalDispersion}{ABSDataSet}(object)

\S4method{minimalDispersion}{ABSDataSet,numeric}(object)<-value}
\description{Accessor functions for the 'minDispersion' slot of a ABSDataSet
object.}
\details{The 'minimalDispersion' is the penalty of dispersion estimation.
User can set the penalty of dispersion by this function}
\alias{minimalDispersion}
\alias{minimalDispersion<-}
\alias{minimalDispersion,ABSDataSet-method}
\alias{minimalDispersion<-,ABSDataSet,numeric-method}
\seealso{\code{\link{callParameter}},\code{\link{ABSDataSet}}}
\arguments{\item{object}{an \code{ABSDataSet} object.}
\item{value}{a positive numeric object}}
\examples{data(simuN5)
obj <- ABSDataSet(counts=simuN5$counts, groups=factor(simuN5$groups))
minimalDispersion(obj)
minimalDispersion(obj) <- 0.2
minimalDispersion(obj)}