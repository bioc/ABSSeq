\docType{methods}
\name{results}
\title{Accessor functions for the result from a ABSDataSet}
\usage{\S4method{results}{ABSDataSet}(object,cnames)}
\description{Accessor functions for the result from a ABSDataSet by given names}
\details{This function returns the result of ABSSeq as a table or a vector depended on
the given names.
}
\alias{results}
\alias{results,ABSDataSet-method}
\value{a table according to 'canmes'.}
\seealso{\code{\link{ABSSeq}}}
\arguments{\item{object}{an ABSDataSet}
\item{cnames}{a vecotr of names for output, which are among:

'baseMean': average expression for each gene including all samples,

'Amean' and 'Bmean': average expression for each group:A and B,

'absD': absolute expression differences between two groups,

'foldChange': log2 of fold-change for each gene based on B vs. A,

'Variance': variance for each gene

'pvalue':pvalue from \code{\link{GPTest}},

'adj.pvalue':ajusted pvalues by BH method. See \code{\link{p.adjust}}.}}
\examples{data(simuN5)
obj <- ABSDataSet(counts=simuN5$counts, groups=factor(simuN5$groups))
obj <- normalFactors(obj)
obj <- calPara(obj)
obj <- GPTest(obj)
head(results(obj))
head(results(obj,c("pvalue","adj.pvalue")))}
