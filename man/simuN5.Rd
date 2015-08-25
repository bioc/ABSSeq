\name{simuN5}
\alias{simuN5}
\docType{data}
\title{
Simulated study with random outliers
}
\description{
Simulated study with random outliers, include five samples for two groups.
It contains counts table, groups and defined differential expression genes.
}
\usage{data(simuN5)}
\format{
  The format is:
List of 3

 $ counts: integer, reads count matrix

 $ groups: two groups

 $ DEs   : differential expression genes
}
\details{
Multiple each gene with a value from 5-10 by chance at pvalue of 0.05.
}
\source{
http://bcf.isb-sib.ch/data/compcodeR/
}
\references{
Soneson C, Delorenzi M: A comparison of methods for differential expression analysis of RNA-seq data. BMC Bioinformatics 2013, 14(1):91.
}
\examples{
data(simuN5)
}
\keyword{datasets}
