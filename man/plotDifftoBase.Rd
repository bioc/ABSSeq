\name{plotDifftoBase}
\alias{plotDifftoBase}
\title{Plot absolute differencs}
\usage{plotDifftoBase(object, cols=c("black", "red"), fcut=1.5, pch=16,
    xlab="Expression level", ylab="ABS diffs", ...)}
\description{Plot absolute differencs against expression levels}
\details{Plot absolute differencs against expression levels
and mark the gene with a color at a given cutoff of fold-change}
\arguments{\item{object}{an ABSDataSet object}
\item{fcut}{the cutoff of fold-change}
\item{cols}{the colors to mark the genes which greater or less than fcut}
\item{xlab}{xlab}
\item{ylab}{ylab}
\item{pch}{pch}
\item{...}{further arguments to \code{plot}}}
\examples{data(simuN5)
obj <- ABSDataSet(counts=simuN5$counts, groups=factor(simuN5$groups))
obj <- normalFactors(obj)
obj <- calPara(obj)
plotDifftoBase(obj)}
