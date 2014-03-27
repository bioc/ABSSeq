#' Accessors for the 'counts' slot of a ABSDataSet object, return a matrix 
#' 
#' The counts slot holds the count data as a matrix of non-negative integer
#' count values, rows and columns for genes and samples, respectively. 
#'
#' @usage
#' \S4method{counts}{ABSDataSet}(object,norm=FALSE)
#'
#' \S4method{counts}{ABSDataSet,matrix}(object)<-value
#'
#' @docType methods
#' @name counts
#' @title Accessors for the 'counts' slot of a ABSDataSet object.
#' @aliases counts ABSDataSet-method counts<- ABSDataSet matrix-method
#' @param object a \code{ABSDataSet} object.
#' @param normalized logical indicating whether or not to normalize the counts before returning
#' @param value an numeric matrix
#' @seealso \code{\link{sizeFactors}}, \code{\link{normalFactors}}
#'
#' @examples
#' 
#' data(simuN5)
#' obj <- ABSDataSet(counts=simuN5$counts, groups=factor(simuN5$groups))
#' head(counts(obj))
counts.ABSDataSet <- function(object, norm=FALSE) {
            if (!norm) {
              return(object@counts)
            } else {
               message(paste("Normalizing used ",object@normMethod,"!",sep=""))
               object=normalFactors(object)
               return( t( t( object@counts ) / sizeFactors(object) ) )
            }
}
#' @name counts      
#' @rdname counts                                                                
#' @export
setMethod("counts", signature(object="ABSDataSet"), counts.ABSDataSet)
#' @name counts
#' @rdname counts
#' @exportMethod "counts<-"
setReplaceMethod("counts", signature(object="ABSDataSet", value="matrix"),
  function( object, value ) {
   object@counts <- value
   validObject(object)
   object
}) 




#' Accessor functions for the 'sizeFactor' slot of a ABSDataSet
#' object.
#' 
#' The sizeFactors vector assigns to each sample a value, used to normalize the
#' counts in each sample according to selected normMethod.
#' 
#' @usage
#' \S4method{sizeFactors}{ABSDataSet}(object)
#'
#' \S4method{sizeFactors}{ABSDataSet,numeric}(object)<-value
#'
#' @docType methods
#' @title Accessors for the 'sizeFactor' slot of a ABSDataSet object.
#' @name sizeFactors
#' @aliases sizeFactors ABSDataSet-method sizeFactors<- ABSDataSet numeric-method
#' @param object a \code{ABSDataSet} object.
#' @param value a numeric object, one for each sample
#' @seealso \code{\link{normalFactors}}
#' @examples
#' 
#' data(simuN5)
#' obj <- ABSDataSet(counts=simuN5$counts, groups=factor(simuN5$groups))
#' normalFactors(obj)
#' sizeFactors(obj)
sizeFactors.ABSDataSet <- function(object) {
siz=object@sizeFactor
 if(is.null(siz))
 {
   message("Run the normalized funtion firstly to get sizefactor!")
   return(NULL)
 }
 return(object@sizeFactor)
}
#' @name sizeFactors
#' @rdname sizeFactors
#' @export
setMethod("sizeFactors", signature(object="ABSDataSet"),sizeFactors.ABSDataSet)

#' @name sizeFactors
#' @rdname sizeFactors
#' @exportMethod "sizeFactors<-"
setReplaceMethod("sizeFactors", signature(object="ABSDataSet", value="numeric"),
                 function( object, value ) {
                   if (any(value <= 0)) {
                     stop("size factors must be positive")
                   }
                   object@normMethod="User"
                   object@sizeFactor=value
                   validObject(object)
                   object
                 })
                 
#' Accessor functions for the 'groups' information in a ABSDataSet
#' object.
#' 
#' The 'groups' is a factor object, contains the experiment design for differential expression analysis.
#' Its length should be equal with the sample size.
#' 
#' @usage
#' \S4method{groups}{ABSDataSet}(object)
#'
#' \S4method{groups}{ABSDataSet,factor}(object)<-value
#'
#' @docType methods
#' @title Accessors for the 'groups' slot of a ABSDataSet object.
#' @name groups
#' @aliases groups ABSDataSet-method groups<- ABSDataSet factor-method
#' @param object a \code{ABSDataSet} object.
#' @param value a \code{factor} object, includes two groups, equal with the number of samples
#' @examples
#' 
#' data(simuN5)
#' obj <- ABSDataSet(counts=simuN5$counts, groups=factor(simuN5$groups))
#' groups(obj)
groups.ABSDataSet <- function(object) {
gps=object@groups
 if(length(gps)==0)
 {
   message("There is no groups information!")
   return(NULL)
 }
 return(gps)
}
#' @name groups
#' @rdname groups
#' @export
setMethod("groups", signature(object="ABSDataSet"),groups.ABSDataSet)

#' @name groups
#' @rdname groups
#' @exportMethod "groups<-"
setReplaceMethod("groups", signature(object="ABSDataSet", value="factor"),
                 function( object, value ) {
                   object@groups=value
                   validObject(object)
                   object
                 }) 
 
#' Accessor functions for the 'normMethod' information in a ABSDataSet
#' object.
#' 
#' The 'normMethod' is the method for calculating the size factors.
#' Currently, Four methods: 'User', 'total', 'quartile' and 'DESeq' are 
#' available.
#' 
#' @usage
#' \S4method{normMethod}{ABSDataSet}(object)
#'
#' \S4method{normMethod}{ABSDataSet,character}(object)<-value
#'
#' @docType methods
#' @title Accessors for the 'normMethod' slot of a ABSDataSet object.
#' @name normMethod
#' @aliases normMethod ABSDataSet-method normMethod<- ABSDataSet character-method
#' @param object a \code{ABSDataSet} object.
#' @param value a character object, should  be one of 'User', 'total', 'quartile' and 'DESeq'.
#' See \code{\link{ABSDataSet}}
#' @examples
#' 
#' data(simuN5)
#' obj <- ABSDataSet(counts=simuN5$counts, groups=factor(simuN5$groups))
#' normMethod(obj)
normMethod.ABSDataSet <- function(object) {
nm=object@normMethod
 if(length(nm)==0)
 {
   message("There is no 'normMethod' information!")
   return(NULL)
 }
 return(nm)
}
#' @name normMethod
#' @rdname normMethod
#' @export
setMethod("normMethod", signature(object="ABSDataSet"),normMethod.ABSDataSet)

#' @name normMethod
#' @rdname normMethod
#' @exportMethod "normMethod<-"
setReplaceMethod("normMethod", signature(object="ABSDataSet", value="character"),
                 function( object, value ) {
                   object@normMethod=value
                   validObject(object)
                   object
                 }) 
 
 
#' Accessor functions for the result from a ABSDataSet by given names
#' 
#' This function returns the result of ABSSeq as a table or a vector depended on
#' the given names, see \code{\link{ABSSeq}}
#'
#' @usage
#' \S4method{results}{ABSDataSet}(object,cname)
#'
#' @docType methods
#' @name results
#' @title Accessor functions for the result from a ABSDataSet
#' @rdname results
#' @aliases results results,ABSDataSet-method
#' @param object a ABSDataSet
#' @param cname a vecotr of names for output, which are among:
#' 'sizeFactor','baseMean','Amean','Bmean','absD','foldChange','Variance','pvalue','adj.pvalue'
#' \code{\link[genefilter]{shorth}} function from the genefilter package may give better results.
#' @return a table according to canmes.
#' @seealso \code{\link{ABSSeq}}
#' 
#' @examples
#' 
#' data(simuN5)
#' obj <- ABSDataSet(counts=simuN5$counts, groups=factor(simuN5$groups))
#' obj <- normalized(obj)
#' obj <- estimatePara(obj)
#' obj <- GPTest(obj)
#' head(results(obj))

results.ABSDataSet <- function(object, cnames=c("baseMean","Amean","Bmean","absD","foldChange","Variance","pvalue","adj.pvalue") ) {
  if(sum(! cnames %in% c("baseMean","Amean","Bmean","absD","foldChange","Variance","pvalue","adj.pvalue")) > 0)
  {
    stop("input 'cnames' contains names which not in 'baseMean','Amean','Bmean','absD','foldChange','Variance','pvalue','adj.pvalue'!")
  }
  res=c()
  for(each in 1:length(cnames))
  {
    buff=object[[cnames[each]]]
    if(is.null(buff))
    {
     if(cnames[each] %in% c("pvalue","adj.pvalue"))
     {
      stop("Please run GPTest firstly to get pvalues and adjusted pvalues!")
     }
      stop("Please run estimatePara firstly to get general factors!")
    }
    if(each ==1)
    {
      res=buff
    }else
    {
      res=cbind(res,buff)
    }
  }
  if(length(cnames)<2)
  {
    names(res)=rownames(counts(object))
  }else
  {
    rownames(res)=rownames(counts(object))
    colnames(res)=cnames
  }
  return(res)
}
  
#' @rdname results
#' @export
setMethod("results", c("ABSDataSet"),results.ABSDataSet)


