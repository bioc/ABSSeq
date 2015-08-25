\docType{methods}
\name{LevelstoNormFC}
\title{Accessors for the 'LevelstoNormFC' slot of a ABSDataSet object.}
\usage{\S4method{LevelstoNormFC}{ABSDataSet}(object)

\S4method{LevelstoNormFC}{ABSDataSet,numeric}(object)<-value}
\description{Accessor functions for the 'LevelstoNormFC' slot of a ABSDataSet
object.}
\details{The 'LevelstoNormFC' is maximal level of average standard deviation in fold-change normalization 
according to expression level.}
\alias{LevelstoNormFC}
\alias{LevelstoNormFC<-}
\alias{LevelstoNormFC,ABSDataSet-method}
\alias{LevelstoNormFC<-,ABSDataSet,numeric-method}
\seealso{\code{\link{ABSDataSet}}, \code{\link{callParameter}}}
\arguments{\item{object}{an \code{ABSDataSet} object.}
\item{value}{a positive numeric object}}
\examples{data(simuN5)
obj <- ABSDataSet(counts=simuN5$counts, groups=factor(simuN5$groups))
LevelstoNormFC(obj)
LevelstoNormFC(obj) <- 200
LevelstoNormFC(obj)}