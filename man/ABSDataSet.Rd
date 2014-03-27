\name{ABSDataSet}
\title{ABSDataSet object and constructors}
\usage{ABSDataSet(counts, groups, normMethod=c("User", "total","quartile", "DESeq"), sfactors=0)}
\description{ABSDataSet object and constructors}
\details{The function contructs an ABSDataSet object with counts table and groups.
It also checks the structure of counts and groups.The ABSDataSet is a class, used to store the input
  values, intermediate calculations and results of an analysis of differential expression. It also
   contains information for the running time of an analysis.}
\value{An ABSDataSet object.}
\alias{ABSDataSet}
\alias{ABSDataSet-class}
\alias{SumInfo-class}
\alias{[[<-,SumInfo,character,missing-method}
\docType{class}
\arguments{\item{counts}{a matrix or table with at least two columns and one row,}
\item{groups}{a factor with two groups, whose length should be equal  with sample size}
\item{normMethod}{method for estimating the size factors, should be one of 'User', 'total', 'quartile' and 'DESeq'. See \code{normalFactor} for description.}
\item{sfactors}{size factors for 'User' method, self-defined size factors by user.}}
\examples{counts <- matrix(1:4,ncol=2)
groups <- factor(c("a","b"))
obj <- ABSDataSet(counts, groups)}
