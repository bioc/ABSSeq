\name{ABSDataSet}
\title{ABSDataSet object and constructors}
\usage{ABSDataSet(counts, groups, normMethod = c("user", "total", "quartile",
  "geometric"), sizeFactor = 0, minDispersion = NULL, minRates = 0.1,
  maxRates = 0.3, LevelstoNormFC = 100)}
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

\item{normMethod}{method for estimating the size factors, should be one of 'user', 'total', 'quartile' and 'geometric'. See \code{\link{normalFactors}} for description.}

\item{sizeFactor}{ize factors for 'user' method, self-defined size factors by user.}

\item{minDispersion}{a positive double for user-defined penalty of dispersion estimation}

\item{minRates}{low bounder rate of baseline estimation for counts difference, default is 0.1}

\item{maxRates}{up bounder rate of baseline estimation for counts difference, default is 0.3. Setting minRates equal with maxRates will result in a testing on user-define rate,}

\item{LevelstoNormFC}{maximal level of average standard deviation in fold-change normalization according to expression level, default is 100.}
}
\examples{counts <- matrix(1:4,ncol=2)
groups <- factor(c("a","b"))
obj <- ABSDataSet(counts, groups)}