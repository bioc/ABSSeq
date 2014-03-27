#' This function performs a default analysis by calling, in order, the functions:
#' \code{\link{normalized}},
#' \code{\link{estimatePara}},
#' \code{\link{GPTest}}.
#'
#' The differential expression analysis uses a generalized generalized Poisson model of the form:
#' \deqn{ GP(D=k|\lambda,\theta)=(\lambda/(1+\lambda * \theta)^k * (1+\theta * k)^(k-1)/k! * exp(-\lambda * (1+\theta * k) / (1+\lambda * \theta))}
#' for \eqn{ \lambda > 0, \theta >= 0, k=0, 1, 2, \ldots}.
#' 
#' @title Differential expression analysis based on the absolute difference.
#' @param object an \code{\link{ABSDataSet}} object, contains the reads count matrix, groups and normalization method.
#' @param perPrior, default is 0.30. See \code{\link{GPTest}} for details.
#' @param dcut default is 0, see \code{\link{GPTest}} for description.
#' @param quiet default is FALSE, whether to print messages at each step
#' @return a table with return a table with elements:
#' Amean and Bmean, mean of reads count for group A and B, 
#' foldChange, log 2 of fold-change, B vs. A, 
#' pvalue, pvalue from GP model, 
#' adj.pvalue, adjuested p-value used BH method.
#' @author Wentao Yang
#' @references Wentao Yang, Philip Rosenstiel & Hinrich Schulenburg: ABSSeq: a new RNA-Seq analysis method based on absolute expression differences and generalized Poisson model
#' @examples
#' 
#' data(simuN5)
#' obj <- ABSDataSet(counts=simuN5$counts, groups=factor(simuN5$groups))
#' res <- ABSSeq(obj)
#' head(res)
#' @export
ABSSeq <- function(object,perPrior=0.30,dcut=0,quiet=FALSE) {
  if (!quiet) message("eistimating size factors....")
  object <- normalFactors(object)
  if (!quiet) message("calculating parameters....")
  object <- calPara(object)
  if (!quiet) message("fitting generalized GP model....")
  object <- GPTest(object, pPrior=perPrior)
  result=results(object,c("Amean","Bmean","foldChange","pvalue","adj.pvalue"))
  return(result)
}
#' This function is borrowed from DESeq.
#' 
#' Given a matrix or data frame of count data, this function estimates the size
#' factors as follows: Each column is divided by the geometric means of the
#' rows. The median (or, if requested, another location estimator) of these
#' ratios (skipping the genes with a geometric mean of zero) is used as the size
#' factor for this column. Typically, you will not call this function directly.
#'
#' @title Low-level function to estimate size factors with robust regression.
#' @param counts a matrix or data frame of counts, i.e., non-negative integer
#' values
#' @param locfunc a function to compute a location for a sample. By default, the
#' median is used.
#' @return a vector with the estimates size factors, one element per column
#' @author Simon Anders
#' @references Simon Anders, Wolfgang Huber: Differential expression analysis for sequence count data. Genome Biology 11 (2010) R106, \url{http://dx.doi.org/10.1186/gb-2010-11-10-r106}
#' @examples
#' 
#' data(simuN5)
#' dat <- simuN5
#' estimateSizeFactorsForMatrix(dat$counts)
#' 
#' @export
estimateSizeFactorsForMatrix <- function( counts, locfunc = median )
{
   loggeomeans <- rowMeans( log( counts ) )
   if (all(is.infinite(loggeomeans))) {
     stop("every gene contains at least one zero, cannot compute log geometric means")
   }
   apply( counts, 2, function(cnts)
      exp( locfunc( ( log(cnts) - loggeomeans )[ is.finite(loggeomeans) ] ) ) )
}

#' Function for esitmating size factors
#' 
#' Given a matrix of count data, this function esitmates the size
#' factors by selected method.
#' It aslo provides four different methods for normalizing according to user-defined size factors,
#' total reads, up quantile (75%) or DESeq (See \code{\link{estimateSizeFactorsForMatrix}}).
#'
#' @title Estimating size factors from the reads count table
#' @param object a ABSSeq object with element of 'counts' and 'normMethod', see the constructor functions
#' \code{\link{ABSDataSet}}. 
#' @return a ABSDataSet object with the estimates size factors, one element per column. Use the \code{\link{sizeFactors}}
#' to show it.  
#' @examples
#' 
#' data(simuN5)
#' obj <- ABSDataSet(counts=simuN5$counts, groups=factor(simuN5$groups))
#' obj <- normalFactors(obj)
#' sizeFactors(obj)
#'
#' @export
normalFactors <- function(object){
  if(!is(object,"ABSDataSet"))
  {
    stop("input is not an ABSDataSet object!")
  }
  validObject(object)
  method=object@normMethod
  if (method == "total") {
    sizefactor=colSums(object@counts)
  
  }
  if (method == "quartile") {
    rowQuar=function(drow){
     drow=sort(drow[drow>0])
     quantile(drow)[4]
    }
    sizefactor=apply(object@counts,2,rowQuar)
  }
  if (method == "DESeq") {
    sizefactor=estimateSizeFactorsForMatrix(object@counts)
  }
  if (method == "User") {
    sizefactor=object@sizeFactor
  }
  sizemax=max(sizefactor)
  if(sizemax<=0)
    {
      stop("The counts table is wrong, with zero row or negative number!")
    }
   sizefactor=sizefactor/mean(sizefactor)
   sizefactor[sizefactor<=0]=1.0
   #sizefactor=1/sizefactor
   #object$counts=t(t(object$counts)*sizefactor)
   #object$sizefactor=sizefactor
   object@sizeFactor=sizefactor
   return(object)
}
#' Calculate mean and variance for a vector
#' not available for user

getGWei=function(dat,rmean,msd,outl)
{
    nmean=mean(dat)
    nvar=var(dat)
    nlen=length(dat)
    if(nlen<2 || nmean==0 || nvar==0 ) return(c(nmean,0))
    nwei=c()
    if(any(outl>3.5))
    {
       if(nlen==2 && (max(dat)<min(rmean) || min(dat)>max(rmean)))
       {
         nmean=mean(c(rmean,dat)) 
       }
       nwei=dnorm((dat-nmean)/(nmean+0.001),0,msd,log=TRUE)       
    }else 
    {
      nwei=dnorm(dat,nmean,sqrt(nvar),log=TRUE)
    }
    nwei=exp(nwei-max(nwei))
    nwei=nwei/sum(nwei)
    if(sum(nwei>0)==0)
    {
      return(c(nmean,nvar))
    }
    nn=sum(nwei>0)
    nmean=sum(dat*nwei)
    nvar=sum(nwei*(dat-nmean)^2)
    nvar=nvar*nn/(nn-1)
    return(c(nmean,nvar)) 
}
#' Calculate mean and variance for a gene
#' not available for user
 getWei=function(dat,msd,mdg,tind1,tind2)
  {
    outl=0.6745*abs(dat-median(dat))/mad(dat)
    if(mad(dat)==0) outl=rep(0,length(outl))
    return(c(getGWei(dat[tind1],dat[tind2],msd,(outl[tind1])),getGWei(dat[tind2],dat[tind1],mdg,(outl[tind2])))) 
  }
#' Calculate a set of parameters from normalized counts table before \code{\link{GPTest}}
#'
#' @param object a \code{\link{ABSDataSet}} object.
#'
#' @return A ABSDataSet object with absolute differences, basemean, mean of each group, variance, 
#' log2 of foldchange and mean of each group, named as 'absD', 'baseMean', 'Amean', 'Bmean', 
#'  'Variance' and 'foldChange', respectively. Use the \code{\link{results}} to get access it
#'
#' @title Calculate parameters for generalized Poisson test (GPTest) 
#' @note This function should run after \code{\link{normalFactors}} or providing size factors.
#' @examples
#'
#' data(simuN5)
#' obj <- ABSDataSet(counts=simuN5$counts, groups=factor(simuN5$groups))
#' obj <- normalFactors(obj)
#' obj <- calPara(obj)
#' plotDifftoBase(obj)
#'
#' @export
calPara <- function(object) {
  if(!is(object,"ABSDataSet"))
  {
    stop("input is not an ABSDataSet object!")
  }
   if(is.null(sizeFactors(object)))
  {
    stop("Please run normalized function before 'calPara'!")
  }
  if(length(sizeFactors(object)) !=length(groups(object)) || sum(sizeFactors(object)<=0)>0) 
  {
   stop("Size factors are wrong, not equal with sample size or with none positive number!")
  }
  
  ncounts=counts(object,TRUE)
  igroups=groups(object)
  ngr1=igroups[1]
  ngr2=igroups[igroups!=igroups[1]][1]
  n1=sum(igroups==ngr1)
  n2=sum(igroups==ngr2)
  
  mdv=0
  if(n1>2)
  {
   mds1= apply(ncounts[,igroups==ngr1],1,mean)
   mds=ncounts[,igroups==ngr1]-mds1
   mdv=mds[mds1>0,]/mds1[mds1>0] 
  }
  if(n2>2)
  {
    mds1= apply(ncounts[,igroups==ngr2],1,mean)
    mds=ncounts[,igroups==ngr2]-mds1
    if(n1>2) 
    {
      mdv=c(mdv,mds[mds1>0,]/mds1[mds1>0])
    }
  }
  if(n1>2 || n2>2) mdv=sd(mdv)
  mdg=0
  if(n1<=2 || n2<=2)
  {
     mds1= apply(ncounts,1,mean)
     mds=ncounts[,igroups==ngr1]-mds1
     mdg=mds[mds1>0,]/mds1[mds1>0]
     mdg=sd(mdg) 
  }
  msd1=mdv
  msd2=mdv
  if(n1<=2) msd1=mdg
  if(n2<=2) msd2=mdg
  wei <- apply(ncounts,1,getWei,msd1,msd2,igroups==ngr1,igroups==ngr2)
  #wei2 <- apply(ncounts,1,getWei,mdssd,igroups==ngr2)
 
  wei=t(wei)
  
  object[["baseMean"]]=(wei[,1]*n1+wei[,3]*n2)/(n1+n2)
  object[["absD"]]=(abs(wei[,1]-wei[,3]))        #ceiling
  object[["Variance"]]=(wei[,2]*(n1-1)+wei[,4]*(n2-1))/(n1+n2-2)
  object[["baseMean"]][is.na(object[["baseMean"]])]=0
  object[["Variance"]][is.na(object[["Variance"]])]=0
  object[["baseMean"]][object[["baseMean"]]<0.5]=0.5
  object[["foldChange"]]=log2((wei[,3]+0.0000001)/(wei[,1]+0.0000001))
  object[["Amean"]]=wei[,1]
  object[["Bmean"]]=wei[,3]
  return(object)
}



#' Fit dispersions for generalized Possion model
#'
#' This function estimates the rate (lambda) and dispersion parameter (theta) for genralized
#' Poisson model used Metropolis method.
#'
#' @param dMean a vector of absolute differences between groups.
#' @param dcut the cutoff of dMean for fitting, default is 0
#' 
#' return a list with elements: rate, theta, vector of log looklihood based on rate and theta
#' @examples
#'
#' data(simuN5)
#' obj <- ABSDataSet(counts=simuN5$counts, groups=factor(simuN5$groups))
#' obj <- normalFactors(obj)
#' obj <- calPara(obj)
#' obj <- fitGP(results(obj,"absD"))
#' names(obj)
#'
#' @note this function is not available for users
fitGP <- function (dMean, dcut=0) {
  dMean=dMean[dMean>dcut]
  .Call(C_fitGP,dMean=dMean)
}
#' Calculating p-value for each gene based on fitted generalized Poisson model
#'
#' @param dmean a vector of absolute differences between groups.
#' @param logm the vector of log looklihood from fitGP
#' \code{\link{fitGP}}.
#' @param gams a vector of rates for each gene.
#' @param alphas the vector of dispersion parameters for each gene
#' 
#' @return a vector of pvalue for each gene
#' @note this function is not available for users
getPvalue <- function (dmean,logm,gams,alphas) {
  .Call(C_getPvalue, dmean=ceiling(dmean),logm=logm,gams=gams,alphas=alphas)
}
#' Fitting generalized Poisson model and calculating p-value for each gene
#'
#' This function firstly fits the generalized Possion model used absolute differences between two
#' groups, then calculates the pvalue for each gene and finally adjusts the pvalues via BH method.
#'
#' @title Testing the differential expression
#' @param object an \code{\link{ABSDataSet}} object.
#' @param pPrior the parameter for estimating rates of generalized Possion model for each gene, 
#' default is 0.30. The predefined parameter for estimating GP rate according to expression level of each gene. It is estimated by reducing the 
#' type I error rate in a acceptable level, i.e. 10% (perPrior=0.30) or 5% (perPrior=0.40). For details, please see the references 
#' @param dcut the parameter for fitting GP model, default is 0.
#' 
#' @return an \code{\link{ABSDataSet}} object with additional elements: pvalue and adjusted p-value,
#' denoted by pvalue and adj-pvalue, respectively. Use the 'results' method to get access it. 
#' @note this function should run after \code{\link{calPara}}
#' @examples
#'
#' data(simuN5)
#' obj <- ABSDataSet(counts=simuN5$counts, groups=factor(simuN5$groups))
#' obj <- normalFactors(obj)
#' obj <- calPara(obj)
#' obj <- GPTest(obj)
#' head(results(obj))
#' @export
GPTest <- function(object,pPrior=0.30,dcut=0) {
  if(!is(object,"ABSDataSet"))
  {
    stop("input is not an ABSDataSet object!")
  }
  if(is.null(object[["absD"]]) | is.null(object[["baseMean"]]) | is.null(object[["Variance"]]))
  {
   stop("Please run 'calPara' function before 'GPTest'!")
  }
  absD=object[["absD"]][object[["baseMean"]]>0]
  logm=fitGP(object[["absD"]],dcut);
  lambdas=object[["baseMean"]]*pPrior;
  thetas=(sqrt(object[["Variance"]]/object[["baseMean"]])-1)/object[["baseMean"]]
  thetas[object[["baseMean"]]<=0]=0
  thetas[is.na(thetas)]=0
  thetas[thetas<0]=0
  object[["pvalue"]]=getPvalue(object[["absD"]],logm[[3]],lambdas,thetas)
  object[["adj.pvalue"]]=p.adjust(object[["pvalue"]],method="BH");
  return(object);
}




