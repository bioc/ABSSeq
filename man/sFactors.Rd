\docType{methods}
\name{sFactors}
\title{Accessors for the 'sizeFactor' slot of a ABSDataSet object.}
\usage{\S4method{sFactors}{ABSDataSet}(object)

\S4method{sFactors}{ABSDataSet,numeric}(object)<-value}
\description{Accessor functions for the 'sizeFactor' slot of a ABSDataSet
object.}
\details{The sizeFactors vector assigns to each sample a value, used to normalize the
counts in each sample according to selected normMethod.}
\alias{sFactors}
\alias{sFactors<-}
\alias{sFactors,ABSDataSet-method}
\alias{sFactors<-,ABSDataSet,numeric-method}
\seealso{\code{\link{normalFactors}}}
\arguments{\item{object}{an \code{ABSDataSet} object.}
\item{value}{a numeric object, one for each sample}}
\examples{data(simuN5)
obj <- ABSDataSet(counts=simuN5$counts, groups=factor(simuN5$groups))
obj <- normalFactors(obj)
sFactors(obj)
sFactors(obj) <- runif(10,1,2)
sFactors(obj)}
