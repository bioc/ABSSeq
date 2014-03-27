#' @rdname counts
#' @export
setGeneric("counts", function(object, ...) standardGeneric("counts"))

#' @rdname counts<-
#' @export
setGeneric("counts<-", function(object, value) standardGeneric("counts<-"))

#' @rdname sizeFactors
#' @export
setGeneric("sizeFactors", function(object) standardGeneric("sizeFactors"))

#' @rdname sizeFactors<-
#' @export
setGeneric("sizeFactors<-", function(object, value) standardGeneric("sizeFactors<-"))

#' @rdname results
#' @export
setGeneric("results", function(object,...) standardGeneric("results"))

#' @rdname groups
#' @export
setGeneric("groups", function(object) standardGeneric("groups"))

#' @rdname groups<-
#' @export
setGeneric("groups<-", function(object,value) standardGeneric("groups<-"))


#' @rdname normMethod
#' @export
setGeneric("normMethod", function(object) standardGeneric("normMethod"))

#' @rdname normMethod<-
#' @export
setGeneric("normMethod<-", function(object,value) standardGeneric("normMethod<-"))
