\docType{methods}
\name{excounts}
\title{Accessors for the 'excounts' slot of a ABSDataSet object.}
\usage{
\S4method{excounts}{ABSDataSet,matrix}(object)<-value}
\description{Accessors for the 'excounts' slot of a ABSDataSet object, return a matrix}
\details{The excounts slot holds the nomarlized (trimmed or not) count data as a matrix of non-negative integer
count values, rows and columns for genes and samples, respectively.}
\alias{excounts}
\alias{excounts<-}
\alias{excounts,ABSDataSet-method}
\alias{excounts<-,ABSDataSet,matrix-method}
\seealso{\code{\link{ABSDataSet}}, \code{\link{ReplaceOutliersByMAD}}}
\arguments{\item{object}{a \code{ABSDataSet} object.}
\item{value}{an numeric matrix}}
\examples{data(simuN5)
obj <- ABSDataSet(counts=simuN5$counts, groups=factor(simuN5$groups))
obj <- normalFactors(obj)
obj <- ReplaceOutliersByMAD(obj)
head(excounts(obj))}