\docType{methods}
\name{paired}
\title{Accessors for the 'paired' slot of a ABSDataSet object.}
\usage{\S4method{paired}{ABSDataSet}(object)

\S4method{paired}{ABSDataSet,logical}(object)<-value}
\description{Accessors for the 'paired' slot of a ABSDataSet object, return a logical value}
\details{The 'paired' is the switch for differential expression detection among paired samples, with a boolean value: TRUE or FALSE (default).
When "paired" is TRUE, the replicates in each group should be equal.}
\alias{paired}
\alias{paired<-}
\alias{paired,ABSDataSet-method}
\alias{paired<-,ABSDataSet,logical-method}
\arguments{\item{object}{a \code{ABSDataSet} object.}
\item{value}{value a boolean object, should  be either TRUE or FALSE.}}
\examples{data(simuN5)
obj <- ABSDataSet(counts=simuN5$counts, groups=factor(simuN5$groups))
paired(obj)
paired(obj) <- TRUE
paired(obj)}