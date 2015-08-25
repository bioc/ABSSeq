\docType{methods}
\name{groups}
\title{Accessors for the 'groups' slot of a ABSDataSet object.}
\usage{\S4method{groups}{ABSDataSet}(object)

\S4method{groups}{ABSDataSet,factor}(object)<-value}
\description{Accessor functions for the 'groups' information in a ABSDataSet
object.}
\details{The 'groups' is a factor object, contains the experiment design for differential expression analysis.
Its length should be equal with the sample size.}
\alias{groups}
\alias{groups<-}
\alias{groups,ABSDataSet-method}
\alias{groups<-,ABSDataSet,factor-method}
\arguments{\item{object}{an \code{ABSDataSet} object.}
\item{value}{a \code{factor} object, includes two groups, equal with the number of samples}}
\examples{data(simuN5)
obj <- ABSDataSet(counts=simuN5$counts, groups=factor(simuN5$groups))
groups(obj)
groups(obj) <- factor(rep(c("A","B"),c(5,5)))
groups(obj)
}
