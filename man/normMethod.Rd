\docType{methods}
\name{normMethod}
\title{Accessors for the 'normMethod' slot of a ABSDataSet object.}
\usage{\S4method{normMethod}{ABSDataSet}(object)

\S4method{normMethod}{ABSDataSet,character}(object)<-value}
\description{Accessor functions for the 'normMethod' information in a ABSDataSet
object.}
\details{The 'normMethod' is the method for calculating the size factors.
Currently, Four methods: 'user', 'total', 'quartile' and 'DESeq' are 
available.}
\alias{normMethod}
\alias{normMethod<-}
\alias{normMethod,ABSDataSet-method}
\alias{normMethod<-,ABSDataSet,character-method}
\arguments{\item{object}{an \code{ABSDataSet} object.}
\item{value}{a character object, should  be one of 'user', 'total', 'quartile' and 'geometric'.
See \code{\link{normalFactors}}}}
\examples{data(simuN5)
obj <- ABSDataSet(counts=simuN5$counts, groups=factor(simuN5$groups))
normMethod(obj)
normMethod(obj) <- "geometric"
normMethod(obj)}
