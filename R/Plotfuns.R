#' Plot absolute differencs against expression levels
#'
#' Plot absolute differencs against expression levels
#' and mark the gene with a color at a given cutoff of fold-change
#'
#' @title Plot absolute differencs
#' @param object a ABSDataSet
#' @param fcut the cutoff of fold-change
#' @param cols the colors to mark the genes which greater or less than fcut
#' @param xlab xlab
#' @param ylab ylab
#' @param pch pch
#' @param ... further arguments to \code{plot}
#'
#' @examples
#' 
#' data(simuN5)
#' obj <- ABSDataSet(counts=simuN5$counts, groups=factor(simuN5$groups))
#' obj <- normalized(obj)
#' plotDifftoBase(obj)
#' 
#' @export
plotDifftoBase = function( object, cols = c("black","red"), fcut = 1.5, pch=16, xlab = "Expression level",ylab = "ABS diffs", ...)
{ 
  if(!is(object,"ABSDataSet"))
  {
    stop("input is not an ABSDataSet object!")
  }
  if(fcut<=0)
  {
   stop("fcut is not positive!")
  }
  dx=object[["baseMean"]]
  dy=object[["absD"]]
  dfold=object[["foldChange"]]
  if(length(dy)!=length(dx))
  {
   stop("the lengths of dmean and baseMean are not equal!")
  }
  dy=dy[dx>0]
  dfold=dfold[dx>0]
  dx=dx[dx>0]
  iso.reg <- isoreg(dx, dy)

  plot(dx, dy, xlab=xlab, ylab=ylab, pch=16, col=cols[(dfold>log2(fcut))+1], ...)
  
  lines(iso.reg,
        col= 'gray19',
        pch=10,
        do.points=TRUE,
        lwd=.5, cex=.5
       )
}