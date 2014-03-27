setClass("SumInfo",representation(begins = "POSIXct",ends= "POSIXct"),contains = "environment")
setMethod("[[<-", c("SumInfo","character","missing"),
   function(x, i,j,..., value) { 
      cnames=c("baseMean","Amean","Bmean","absD","foldChange","Variance","pvalue","adj.pvalue")
      if(!i %in% cnames)
      {
       stop( "Can't have an element in this name!" )
      }
      if(!is.numeric(value))
      {
       stop( "Can't have an non-numeric value!" )
      }
      ev <- as(x, "environment")
      ev[[i]] <- value  
       x@ends <- Sys.time() # the update time
       x})

setClass("ABSDataSet", representation(counts = "matrix", groups = "factor",normMethod="character",sizeFactor="numeric"),contains="SumInfo")

setValidity( "ABSDataSet", function( object ) {
  if( any( is.na(counts(object))) )
    return( "the count data contains NA values" )
  if( any( is.infinite(counts(object)) ) )
    return( "the count data contains infinite values" )
  if( any( !is.numeric(counts(object)) ) )
    return( "the count data contains non-numeric values" )
  if( any(round(counts(object))!= counts(object)) )
    return( "the count data is not in integer mode" )
  if( any( counts(object) < 0 ) )
    return( "the count data contains negative values" )

  ngroup=unique(object@groups)
  if(length(ngroup)!=2)
  {
   return("the number of group is not equal with 2!")
  }
  if(ncol(object@counts)!=length(object@groups))
  {
   return("the col number of counts table is not equal with length of groups!")
  }
  if(length(object@normMethod) !=1 || !object@normMethod %in% c("User","total","quartile","DESeq"))
  {
     return("Please choose one of the normalization methods as below: 'User', 'total', 'quartile' and 'DESeq'!")
  } 
  if(object@normMethod =="User" && (any(is.na(object@sizeFactor)) || length(object@sizeFactor)!= length(object@groups) || any(is.infinite(object@sizeFactor)) || any(object@sizeFactor<0) ))
  {
     return("Please provide right size factors for each sample if you choose 'User' as normalization method!")
  } 
  TRUE
} )

#' ABSDataSet object and constructors
#'
#' The function contructs an ABSDataSet object with counts table and groups.
#' It also checks the structure of counts and groups. The ABSDataSet is a class, used to store the input
#' values, intermediate calculations and results of an
#' analysis of differential expression.It also contains information for the running time of an analysis.
#'
#' @title ABSDataSet object
#' @param counts a matrix or table with at least two columns and one row,
#' @param groups a factor with two groups, whose length should be equal  with sample size
#' @param normMethod method for estimating the size factors, should be one of 'User', 'total', 'quartile' and 'DESeq'. See \code{normalFactor} for description.
#' @param sfactors size factors for 'User' method, self-defined size factors by user.
#'
#' @return A ABSDataSet object.
#' 
#' @aliases ABSDataSet ABSDataSet-class
#'
#' @docType class
#' @examples
#'
#' counts <- matrix(1:4,ncol=2)
#' groups <- factor(c("a","b"))
#' obj <- ABSDataSet(counts, groups)
#' @export
ABSDataSet <- function(counts, groups, normMethod=c("User","total","quartile","DESeq"),sfactors=0) {
  if (is.null(dim(counts))) {
      stop("'counts' is not like a matrix or a table!")
    }
  if(length(normMethod)!=1)
  {
     normMethod="quartile"
  } 
  if(!normMethod %in% c("User","total","quartile","DESeq"))
  {
     stop("Please use one of the normalization methods as below: 'User', 'total', 'quartile' and 'DESeq'!")
  } 
  if(normMethod=="User")
  {
    obj=new("ABSDataSet",counts=as.matrix(counts),groups=as.factor(groups),normMethod=normMethod,sizeFactor=sfactors,begins=Sys.time(),ends=Sys.time())
  }else
  {
    obj=new("ABSDataSet",counts=as.matrix(counts),groups=as.factor(groups),normMethod=normMethod,begins=Sys.time(),ends=Sys.time())
  }
  return(obj)
}