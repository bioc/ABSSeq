setClass("SumInfo",representation(begins = "POSIXct",ends= "POSIXct"),contains = "environment")
setMethod("[[<-", c("SumInfo","character","missing"),
   function(x, i,j,..., value) { 
      cnames <- c("baseMean","Amean","Bmean","absD","foldChange","rawFC","lowFC","Variance","priors","pvalue","adj.pvalue","m1","mults","trimmed","preab","mts")
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

setClass("ABSDataSet", representation(counts = "matrix",excounts = "matrix" ,groups = "factor",normMethod="character",sizeFactor="numeric",paired="logical",minDispersion="numeric"
                                      ,minRates="numeric",maxRates="numeric",LevelstoNormFC="numeric"),contains="SumInfo")

setValidity( "ABSDataSet", function( object ) {
  if( any( is.na(counts(object))) || any( is.na(excounts(object))))
    return( "the count data contains NA values" )
  if( any( is.infinite(counts(object)) ) || any( is.infinite(excounts(object)) ) )
    return( "the count data contains infinite values" )
  if( any( !is.numeric(counts(object)) ) || any( !is.numeric(excounts(object)) ) )
    return( "the count data contains non-numeric values" )
#  if( any(round(counts(object))!= counts(object)) )
#    return( "the count data is not in integer mode" )
  if( any( counts(object) < 0 ) || any( excounts(object) < 0 ) )
    return( "the count data contains negative values" )

  ngroup=unique(object@groups)
  if(length(ngroup)!=2)
  {
   return("the number of group is not equal with 2!")
  }
  if(paired(object))
  {
    if(sum(object@groups==ngroup[1])!=sum(object@groups==ngroup[2]))
    {
      return("For paired comparison, the replicates in each group should be equal!")
    }
  }
  if(object@maxRates>=1.0 || object@maxRates<=0.0)
  {
    return("'maxRates' is bigger than 1 or less than 0!") 
  }
  if(object@minRates>=1 || object@minRates<=0 || object@minRates>object@maxRates)
  {
    return("'minRates' is bigger than 1 or less than 0 or bigger than 'maxRates'!") 
  }
  if(object@LevelstoNormFC<=0)
  {
    return("'LevelstoNormFC' is less than 0!") 
  }
  if(length(object@minDispersion)>0 && object@minDispersion < 0)
  {
    return("'minDispersion' is less than 0!")
  }
  if(ncol(object@counts)!=length(object@groups))
  {
   return("the col number of counts table is not equal with length of groups!")
  }
  if(length(object@normMethod) !=1 || !object@normMethod %in% c("user","total","quartile","geometric"))
  {
     return("Please choose one of the normalization methods as below: 'user', 'total', 'quartile' and 'geometric'!")
  } 
  if(object@normMethod =="user" && (any(is.na(object@sizeFactor)) || length(object@sizeFactor)!= length(object@groups) || any(is.infinite(object@sizeFactor)) || any(object@sizeFactor<0) ))
  {
   
     return("Please provide right size factors for each sample if you choose 'user' as normalization method!")
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
#' @param normMethod method for estimating the size factors, should be one of 'user', 'total', 'quartile' and 'geometric'. See \code{\link{normalFactors}} for description.
#' @param sizeFactor size factors for 'user' method, self-defined size factors by user.
#' @param paired switch for differential expression detection in paired samples.
#' @param minDispersion a positive double for user-defined penalty of dispersion estimation
#' @param minRates low bounder rate of baseline estimation for counts difference, default is 0.1
#' @param maxRates up bounder rate of baseline estimation for counts difference, default is 0.3. Setting minRates equal with maxRates will result in a testing on user-define rate, 
#' @param LevelstoNormFC maximal level of average standard deviation in fold-change normalization according to expression level, default is 100.
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
#' obj <- ABSDataSet(counts, groups, paired=TRUE)
#' @export
ABSDataSet <- function(counts, groups, normMethod=c("user","total","quartile","geometric"),sizeFactor=0,paired=FALSE,minDispersion=NULL,minRates=0.1,maxRates=0.3,LevelstoNormFC=100) {
  if (is.null(dim(counts))) {
      stop("'counts' is not like a matrix or a table!")
    }
  if(length(normMethod)!=1)
  {
     normMethod <- "quartile"
  }
  if(!normMethod %in% c("user","total","quartile","geometric"))
  {
     stop("Please use one of the normalization methods as below: 'user', 'total', 'quartile' and 'geometric'!")
  } 
  if(normMethod=="user")
  {
    obj=new("ABSDataSet",counts=as.matrix(counts),groups=as.factor(groups),
            minRates=minRates,maxRates=maxRates,paired=paired,LevelstoNormFC=LevelstoNormFC,normMethod=normMethod,sizeFactor=sizeFactor,begins=Sys.time(),ends=Sys.time())
  }else
  {
    obj=new("ABSDataSet",counts=as.matrix(counts),groups=as.factor(groups),
            minRates=minRates,maxRates=maxRates,paired=paired,LevelstoNormFC=LevelstoNormFC,normMethod=normMethod,begins=Sys.time(),ends=Sys.time())
  }
  if(!is.null(minDispersion))
  {
    obj@minDispersion <- minDispersion
  }
    
  return(obj)
}