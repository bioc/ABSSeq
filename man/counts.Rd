\docType{methods}
\name{counts}
\title{Accessors for the 'counts' slot of a ABSDataSet object.}
\usage{\S4method{counts}{ABSDataSet}(object,norm=FALSE)

\S4method{counts}{ABSDataSet,matrix}(object)<-value}
\description{Accessors for the 'counts' slot of a ABSDataSet object, return a matrix}
\details{The counts slot holds the count data as a matrix of non-negative integer
count values, rows and columns for genes and samples, respectively.}
\alias{counts}
\alias{counts<-}
\alias{counts,ABSDataSet-method}
\alias{counts<-,ABSDataSet,matrix-method}
\seealso{\code{\link{sizeFactors}}, \code{\link{normalFactors}}}
\arguments{\item{object}{a \code{ABSDataSet} object.}
\item{norm}{logical indicating whether or not to normalize the counts before returning}
\item{value}{an numeric matrix}}
\examples{data(simuN5)
obj <- ABSDataSet(counts=simuN5$counts, groups=factor(simuN5$groups))
head(counts(obj))
counts(obj) <- matrix(1:50,nrow=5,ncol=10)
head(counts(obj))}
