#' Plot absolute differencs against expression levels
#'
#' Plot absolute differencs against expression levels
#' and mark the gene with a color at a given cutoff of fold-change
#'
#' @title Plot absolute log2 fold-change against base mean of expression
#' @param object a ABSDataSet
#' @param foldname indicates kind of fold-change in plotting, default is 'foldChange', see \code{results}
#' @param adj.pcut cutoff for differential expressed genes, marked by different color, default is 0.05
#' @param cols the colors to mark the non-DE and DE genes, defualt is black and red, respectively
#' @param xlab xlab, default is 'log2 of Expression level'
#' @param ylab ylab, default is 'log2 fold-change'
#' @param pch pch, default is 16
#' @param ..., further arguments to \code{plot}
#'
#' @examples
#' 
#' data(simuN5)
#' obj <- ABSDataSet(counts=simuN5$counts, groups=factor(simuN5$groups))
#' obj <- ABSSeq(obj)
#' plotDifftoBase(obj)
#' 
#' @export
plotDifftoBase = function(object,foldname="foldChange", adj.pcut=0.05, cols = c("black","red"),pch=16, xlab = "log2 of Expression level",ylab = "log2 fold-change", ...)
{ 
  if(!is(object,"ABSDataSet"))
  {
    stop("input is not an ABSDataSet object!")
  }
  if(length(cols)!=2)
  {
    stop("Please provide two colors!")
  }
  if(is.null(object[["Amean"]]) || is.null(object[["Bmean"]]) || is.null(object[[foldname]]) || is.null(object[[foldname]]) )
  {
    stop("Please run ABSSeq firstly!")
  }
  ccols <- rep(cols[1],length(object[["Amean"]]))
  ccols[object[["adj.pvalue"]]<adj.pcut] <- cols[2]
  plot((object[["Amean"]]+object[["Bmean"]])/2,object[[foldname]],col=ccols,pch=pch, xlab=xlab, ylab=ylab,...)
}