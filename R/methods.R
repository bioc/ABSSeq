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
#' @aliases counts ABSDataSet-method logical-method
#' @aliases counts ABSDataSet-method counts<- ABSDataSet matrix-method
#' @param object a \code{ABSDataSet} object.
#' @param norm logical indicating whether or not to normalize the counts before returning
#' @param value an numeric matrix
#' @seealso \code{\link{sFactors}}, \code{\link{normalFactors}}
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
               return( t( t( object@counts ) / sFactors(object) ) )
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


#' Accessors for the 'excounts' slot of a ABSDataSet object, return a matrix 
#' 
#' The counts slot holds the normalized (trimmed or not) as a matrix of non-negative integer
#' count values, rows and columns for genes and samples, respectively. 
#'
#' @usage
#' \S4method{excounts}{ABSDataSet}(object)
#'
#' \S4method{excounts}{ABSDataSet,matrix}(object)<-value
#'
#' @docType methods
#' @name excounts
#' @title Accessors for the 'excounts' slot of a ABSDataSet object.
#' @aliases excounts ABSDataSet-method
#' @aliases excounts ABSDataSet-method excounts<- ABSDataSet matrix-method
#' @param object a \code{ABSDataSet} object.
#' @param value an numeric matrix
#' @seealso \code{\link{ReplaceOutliersByMAD}}
#'
#' @examples
#' 
#' data(simuN5)
#' obj <- ABSDataSet(counts=simuN5$counts, groups=factor(simuN5$groups))
#' obj <- normalFactors(obj)
#' obj <- ReplaceOutliersByMAD(obj)
#' head(excounts(obj))
excounts.ABSDataSet <- function(object) {
    return(object@excounts)
}
#' @name counts      
#' @rdname counts                                                                
#' @export
setMethod("excounts", signature(object="ABSDataSet"), excounts.ABSDataSet)
#' @name counts
#' @rdname counts
#' @exportMethod "counts<-"
setReplaceMethod("excounts", signature(object="ABSDataSet", value="matrix"),
                 function( object, value ) {
                   object@excounts <- value
                   validObject(object)
                   object
                 }) 

#' Accessor functions for the 'sFactor' slot of a ABSDataSet
#' object.
#' 
#' The sFactor vector assigns to each sample a value, used to normalize the
#' counts in each sample according to selected normMethod.
#' 
#' @usage
#' \S4method{sFactors}{ABSDataSet}(object)
#'
#' \S4method{sFactors}{ABSDataSet,numeric}(object)<-value
#'
#' @docType methods
#' @title Accessors for the 'sizeFactor' slot of a ABSDataSet object.
#' @name sFactors
#' @aliases sFactors ABSDataSet-method sFactors ABSDataSet numeric-method
#' @aliases sFactors ABSDataSet-method sFactors<- ABSDataSet numeric-method
#' @param object a \code{ABSDataSet} object.
#' @param value a numeric object, one for each sample
#' @seealso \code{\link{normalFactors}}
#' @examples
#' 
#' data(simuN5)
#' obj <- ABSDataSet(counts=simuN5$counts, groups=factor(simuN5$groups))
#' pbj <- normalFactors(obj)
#' sFactors(obj)
sizeFactors.ABSDataSet <- function(object) {
siz <- object@sizeFactor
 if(length(siz)==0)
 {
   message("Run the normalized funtion firstly to get sizefactor!")
   return(NULL)
 }
 return(object@sizeFactor)
}
#' @name sFactors
#' @rdname sFactors
#' @export
setMethod("sFactors", signature(object="ABSDataSet"),sizeFactors.ABSDataSet)

#' @name sFactors
#' @rdname sFactors
#' @exportMethod "sFactors<-"
setReplaceMethod("sFactors", signature(object="ABSDataSet", value="numeric"),
                 function( object, value ) {
                   if (any(value <= 0)) {
                     stop("size factors must be positive")
                   }
                   object@normMethod <- "user"
                   object@sizeFactor <- value
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
gps <- object@groups
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
                   object@groups <- value
                   validObject(object)
                   object
                 }) 
 
#' Accessor functions for the 'normMethod' information in a ABSDataSet
#' object.
#' 
#' The 'normMethod' is the method for calculating the size factors.
#' Currently, Four methods: 'user', 'total', 'quartile' and 'geometric' are 
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
#' @param value a character object, should  be one of 'user', 'total', 'quartile' and 'geometric'.
#' See \code{\link{ABSDataSet}}
#' @examples
#' 
#' data(simuN5)
#' obj <- ABSDataSet(counts=simuN5$counts, groups=factor(simuN5$groups))
#' normMethod(obj)
#' normMethod(obj) <- "geometric"
#' normMethod(obj)
normMethod.ABSDataSet <- function(object) {
nm <- object@normMethod
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
                   object@normMethod <- value
                   validObject(object)
                   object
                 }) 
 
#' Accessor functions for the 'minimalDispersion' information in a ABSDataSet
#' object.
#' 
#' The 'minimalDispersion' is the penalty of dispersion estimation. See \code{\link{callParameter}} and \code{\link{ABSDataSet}}.
#' User can set the penalty of dispersion by this function
#' 
#' @usage
#' \S4method{minimalDispersion}{ABSDataSet}(object)
#'
#' \S4method{minimalDispersion}{ABSDataSet,numeric}(object)<-value
#'
#' @docType methods
#' @title Accessors for the 'minimalDispersion' slot of a ABSDataSet object.
#' @name minimalDispersion
#' @aliases minimalDispersion ABSDataSet-method minimalDispersion<- ABSDataSet numeric-method
#' @param object a \code{ABSDataSet} object.
#' @param value a double object, should  be positive.
#' See \code{\link{ABSDataSet}}
#' @examples
#' 
#' data(simuN5)
#' obj <- ABSDataSet(counts=simuN5$counts, groups=factor(simuN5$groups), minDispersion=0.1)
#' minimalDispersion(obj)
#' minimalDispersion(obj) <- 0.2
#' minimalDispersion(obj)
minimalDispersion.ABSDataSet <- function(object) {
  nm <- object@minDispersion
  if(length(nm)==0)
  {
    message("There is no 'minimalDispersion' information!")
    return(NULL)
  }
  return(nm)
}
#' @name minimalDispersion
#' @rdname minimalDispersion
#' @export
setMethod("minimalDispersion", signature(object="ABSDataSet"),minimalDispersion.ABSDataSet)

#' @name minimalDispersion
#' @rdname minimalDispersion
#' @exportMethod "minimalDispersion<-"
setReplaceMethod("minimalDispersion", signature(object="ABSDataSet", value="numeric"),
                 function( object, value ) {
                   object@minDispersion=value
                   validObject(object)
                   object
                 }) 
#' Accessor functions for the 'minRates' information in a ABSDataSet
#' object.
#' 
#' The 'minRates' is the lower bound of rate for baseline of counts difference esitimation. See \code{\link{callParameter}} and \code{\link{ABSDataSet}}.
#' 
#' @usage
#' \S4method{minRates}{ABSDataSet}(object)
#'
#' \S4method{minRates}{ABSDataSet,numeric}(object)<-value
#'
#' @docType methods
#' @title Accessors for the 'minRates' slot of a ABSDataSet object.
#' @name minRates
#' @aliases minRates ABSDataSet-method minRates<- ABSDataSet numeric-method
#' @param object a \code{ABSDataSet} object.
#' @param value a numeric object, should  be greater than 0.
#' See \code{\link{ABSDataSet}}
#' @examples
#' 
#' data(simuN5)
#' obj <- ABSDataSet(counts=simuN5$counts, groups=factor(simuN5$groups))
#' minRates(obj)
#' minRates(obj) <- 0.2
#' minRates(obj)
minRates.ABSDataSet <- function(object) {
  nm <- object@minRates
  if(length(nm)==0)
  {
    message("There is no 'minRates' information!")
    return(NULL)
  }
  return(nm)
}
#' @name minRates
#' @rdname minRates
#' @export
setMethod("minRates", signature(object="ABSDataSet"),minRates.ABSDataSet)

#' @name minRates
#' @rdname minRates
#' @exportMethod "minRates<-"
setReplaceMethod("minRates", signature(object="ABSDataSet", value="numeric"),
                 function( object, value ) {
                   object@minRates=value
                   validObject(object)
                   object
                 }) 
#' Accessor functions for the 'maxRates' information in a ABSDataSet
#' object.
#' 
#' The 'maxRates' is the upper bound of rate for baseline of counts difference esitimation. See \code{\link{callParameter}} and \code{\link{ABSDataSet}}.
#' 
#' @usage
#' \S4method{maxRates}{ABSDataSet}(object)
#'
#' \S4method{maxRates}{ABSDataSet,numeric}(object)<-value
#'
#' @docType methods
#' @title Accessors for the 'maxRates' slot of a ABSDataSet object.
#' @name maxRates
#' @aliases maxRates ABSDataSet-method maxRates<- ABSDataSet numeric-method
#' @param object a \code{ABSDataSet} object.
#' @param value a numeric object, should  be greater than 0.
#' See \code{\link{ABSDataSet}}
#' @examples
#' 
#' data(simuN5)
#' obj <- ABSDataSet(counts=simuN5$counts, groups=factor(simuN5$groups))
#' maxRates(obj)
#' maxRates(obj) <- 0.4
#' maxRates(obj)
maxRates.ABSDataSet <- function(object) {
  nm <- object@maxRates
  if(length(nm)==0)
  {
    message("There is no 'maxRates' information!")
    return(NULL)
  }
  return(nm)
}
#' @name maxRates
#' @rdname maxRates
#' @export
setMethod("maxRates", signature(object="ABSDataSet"),maxRates.ABSDataSet)

#' @name maxRates
#' @rdname maxRates
#' @exportMethod "maxRates<-"
setReplaceMethod("maxRates", signature(object="ABSDataSet", value="numeric"),
                 function( object, value ) {
                   object@maxRates=value
                   validObject(object)
                   object
                 }) 
#' Accessor functions for the 'LevelstoNormFC' information in a ABSDataSet
#' object.
#' 
#' The 'LevelstoNormFC' is maximal level of average standard deviation in fold-change normalization according to expression level. See \code{\link{callParameter}} and \code{\link{ABSDataSet}}.
#' 
#' @usage
#' \S4method{LevelstoNormFC}{ABSDataSet}(object)
#'
#' \S4method{LevelstoNormFC}{ABSDataSet,numeric}(object)<-value
#'
#' @docType methods
#' @title Accessors for the 'LevelstoNormFC' slot of a ABSDataSet object.
#' @name LevelstoNormFC
#' @aliases LevelstoNormFC ABSDataSet-method LevelstoNormFC<- ABSDataSet numeric-method
#' @param object a \code{LevelstoNormFC} object.
#' @param value a numeric object, should  be greater than 0.
#' See \code{\link{ABSDataSet}}
#' @examples
#' 
#' data(simuN5)
#' obj <- ABSDataSet(counts=simuN5$counts, groups=factor(simuN5$groups))
#' LevelstoNormFC(obj)
#' LevelstoNormFC(obj) <- 200
#' LevelstoNormFC(obj)
LevelstoNormFC.ABSDataSet <- function(object) {
  nm <- object@LevelstoNormFC
  if(length(nm)==0)
  {
    message("There is no 'LevelstoNormFC' information!")
    return(NULL)
  }
  return(nm)
}
#' @name LevelstoNormFC
#' @rdname LevelstoNormFC
#' @export
setMethod("LevelstoNormFC", signature(object="ABSDataSet"),LevelstoNormFC.ABSDataSet)

#' @name LevelstoNormFC
#' @rdname LevelstoNormFC
#' @exportMethod "LevelstoNormFC<-"
setReplaceMethod("LevelstoNormFC", signature(object="ABSDataSet", value="numeric"),
                 function( object, value ) {
                   object@LevelstoNormFC=value
                   validObject(object)
                   object
                 }) 
#' Accessor functions for the result from a ABSDataSet by given names
#' 
#' This function returns the result of ABSSeq as a table or a vector depended on
#' the given names, see \code{\link{ABSSeq}}
#'
#' @usage
#' \S4method{results}{ABSDataSet}(object,cnames)
#'
#' @docType methods
#' @name results
#' @title Accessor functions for the result from a ABSDataSet
#' @rdname results
#' @aliases results results,ABSDataSet-method
#' @param object a ABSDataSet
#' @param cnames a vecotr of names for output, which are among:
#' 'Amean','Bmean', log2 of mean counts for group A and B,
#' "baseMean', estimated mean for absolute counts difference (absD), used for mu in \code{\link{pnbinom}}
#' 'absD', absolute counts difference in total
#' 'Variance', pooled Variance for two groups
#' 'rawFC','lowFC', 'foldChange', log2 fold-change of original (Bmean-Amean), corrected by expression level and corrected by both expression level and gene-specific dispersion
#' 'pvalue','adj.pvalue', pvalue and adjusted pvalue
#' 'trimmed', number of trimmed outliers
#' @return a table according to canmes.
#' @seealso \code{\link{ABSSeq}}
#' 
#' @examples
#' 
#' data(simuN5)
#' obj <- ABSDataSet(counts=simuN5$counts, groups=factor(simuN5$groups))
#' obj <- normalFactors(obj)
#' obj <- callParameter(obj)
#' obj <- callDEs(obj)
#' head(results(obj))

results.ABSDataSet <- function(object, cnames=c("Amean","Bmean","baseMean","absD","Variance","rawFC","lowFC","foldChange","pvalue","adj.pvalue","trimmed") ) {
  if(sum(! cnames %in% c("baseMean","Amean","Bmean","absD","foldChange","rawFC","lowFC","Variance","pvalue","adj.pvalue","trimmed")) > 0)
  {
    stop("input 'cnames' contains names which not in 'baseMean','Amean','rawFC','lowFC','Bmean','absD','foldChange','Variance','pvalue','trimmed','adj.pvalue'!")
  }
  res <- c()
  for(each in 1:length(cnames))
  {
    buff=object[[cnames[each]]]
    if(is.null(buff))
    {
     if(cnames[each] %in% c("pvalue","adj.pvalue"))
     {
      stop("Please run callDEs firstly to get pvalues and adjusted pvalues!")
     }
      stop("Please run callParameter firstly to get general factors!")
    }
    if(each ==1)
    {
      res <- buff
    }else
    {
      res <- cbind(res,buff)
    }
  }
  if(length(cnames)<2)
  {
    names(res) <- rownames(counts(object))
  }else
  {
    rownames(res) <- rownames(counts(object))
    colnames(res) <- cnames
  }
  return(res)
}
  
#' @rdname results
#' @export
setMethod("results", c("ABSDataSet"),results.ABSDataSet)


