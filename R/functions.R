#' This function performs a default analysis by calling, in order, the functions:
#' \code{\link{normalFactors}},
#' \code{\link{callParameter}},
#' \code{\link{callDEs}}.
#'
#' The differential expression analysis models the total counts difference by a Negative binomal distribution \deqn{NB(\mu,r)}:
#' 
#' @title Differential expression analysis based on the total counts difference.
#' @param object an \code{\link{ABSDataSet}} object, contains the reads count matrix, groups and normalization method.
#' @param replaceOutliers default is TRUE, switch for outlier replacement.
#' @param adjmethod defualt is 'BH', method for p-value adjusted, see \code{\link{p.adjust.methods}} for details
#' @param useaFold defualt is FALSE, switch for DE detection through fold-change, see \code{\link{callDEs}} for details
#' @param quiet default is FALSE, whether to print messages at each step
#' @param ... parameters passed to \code{\link{ReplaceOutliersByMAD}} and \code{\link{genAFold}} from \code{\link{callParameter}}
#' @return an ABSDataSet object with additional elements, which can be retrieved by \code{\link{results}}:
#' Amean and Bmean, mean of log2 normalized reads count for group A and B, 
#' foldChange, shrinked (expression level and gene-specific) log2 of fold-change, B - A, 
#' rawFC, raw log2 of fold-change, B-A (without shrinkage),
#' lowFC, expression level corrected log2 fold-change,
#' pvalue, pvalue from NB distribution model, 
#' adj.pvalue, adjuested p-value used p.adjust method.
#' @author Wentao Yang
#' @references Wentao Yang, Philip Rosenstiel & Hinrich Schulenburg: ABSSeq: a new RNA-Seq analysis method based on modelling absolute expression differences
#' @examples
#' 
#' data(simuN5)
#' obj <- ABSDataSet(counts=simuN5$counts, groups=factor(simuN5$groups))
#' obj <- ABSSeq(obj)
#' res <- results(obj,c("Amean","Bmean","foldChange","pvalue","adj.pvalue"))
#' head(res)
#' @export
ABSSeq <- function(object,adjmethod="BH", replaceOutliers=TRUE, useaFold=FALSE, quiet=FALSE,...) {
  if (!quiet) message("eistimating size factors....")
  object <- normalFactors(object)
  if (!quiet) message("calculating parameters and fitting....")
  if(length(groups(object))==2)
  {
    message("No replicates! switch to 'callParameterwithoutReplicates'!")
    object <- callParameterwithoutReplicates(object)
    useaFold <- FALSE
  }
  else
    object <- callParameter(object,replaceOutliers,...)
  if (!quiet) message("Calling p-value and adjusted it....")
  object <- callDEs(object, adjmethod,useaFold)
  return(object)
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
      exp( locfunc( ( log(cnts) - loggeomeans )[ is.finite(loggeomeans)&cnts>0] ) ) )
}
#' Function of qtotal for esitmating size factors 
#' 
#' Given a matrix of count data, this function esitmates the size
#' factors by qtotal method, which is based on assessing DE (CV) and ranking. The CV is estimated via sliding window.
#'
#' @title Estimating size factors from the reads count table via ranking 
#' @param ma a count matrix
#' @param qper quantile for assessing dispersion of data, default is 0.95, which serves to avoid outliers, should in (0,1]
#' @param qst start of quantile for estimating cv ratio, should be in [0,1], default is 0.1
#' @param qend end of quantile for estimating cv ratio, should be in [qbound,1-qbound], default is .95
#' @param qstep step of quantile for estimating cv ratio (sliding window), should be in (0,1], default is 0.01
#' @param qbound window size for estimating cv and shifted size factor, default is 0.05, a smaller window size is suitable if number of genes is large
#' @param mcut cutoff of mean from sliding window to avoid abnormal cv, should >=0, default is 5
#' @param qcl scale for outlier detection, should >=0, default is 0.2
#' @return a vector with the estimates size factors, one element per column
#' @examples
#' 
#' data(simuN5)
#' counts <- simuN5$counts
#' qtotalNormalized(counts)
#'
#' @export

qtotalNormalized=function(ma,qper=0.95,qst=0.1,qend=.95,qstep=0.01,qbound=0.05,mcut=5,qcl=0.2)
{
  if(qper>1||qper<=0||qst>1||qst>qend||qend>1||qend<0||qbound<0||qbound>1||mcut<0) 
    stop("quartile and boundary should be in (0,1]")
  
  rowQuar=function(x,y,qper=0.95,qst=0.1,qend=.95,qstep=0.01,qbound=0.05,mcut=5,qcl=0.2) {
    if(length(x)!=length(y)) stop("Please input samples with equal gene number!")
    #in case outliers
    indx <- x<=quantile(x,qper)&y<=quantile(y,qper)
    x <-x[indx]
    y <-y[indx]
    alens <- length(x)
    scl <- sum(x)/sum(y)
    if(!is.numeric(scl)||is.infinite(scl)||scl==0) return(1)
    qstep <- max(qstep, 1/alens)
    if(qbound < 10/alens)
    {
      #message("qbound for normalization is too small! reset as 10 divided by gene number!")
      qbound <- min(1,10/alens)
    }
    qend <- min(qend,1-qbound)
    qst <- min(qst,qend)
    yt <-x
    ys <- y
    
    aas <- seq(qst,qend,qstep)
    ra <- c()
    axs <- quantile(x,aas)
    axe <- quantile(x,aas+qbound)
    ays <- quantile(y,aas)
    aye <- quantile(y,aas+qbound)
    for(i in 1:length(aas))
    {
      indx <- x>=axs[i]&x<=axe[i]
      indy <- y>=ays[i]&y<=aye[i]
      mx <- x[indx]
      my <- y[indy]
      ma <-mean(mx)
      mb <-mean(my)
      ##use CV instead of sd of log data to avoid influence of zero counts, slightly different
      sa <- sd(mx)
      sb <- sd(my)
      if(ma>mcut&&mb>mcut) ra <- c(ra, sa/ma/(sb/mb))
    }
    if(length(ra)<1) ra <-1 
    indx <- is.numeric(ra)||is.finite(ra)||ra>0||is.na(ra)
    if(sum(ind)>0) ra <- ra[indx]
    if(length(ra)<1) ra <-1 
    ##remove outliers
    if(length(ra)<0)
    {
      rdif <- abs(diff(log(ra)))
      qcl <- max(qcl,mean(rdif)*2)
      rind <- which(rdif>qcl)
      if(length(rind)>0)
      {
        ra[rind] <- 1
        ra[rind+1] <- 1
        rind <- rind-1
        rind <- rind[rind>0]
        if(length(rind)>0) ra[rind] <-1
      }
    }
    rvm <- max(ra)
    rvi <- min(ra)
    indx <- x>=quantile(x,1-qbound)
    mr <- sum(x[indx])/sum(y[indx])
    indx <- y>=quantile(y,1-qbound)
    sr <- sum(x[indx])/sum(y[indx])
    m1 <- mr/(mr/sr)^(1/(1+rvm))
    m2 <- mr/(mr/sr)^(1/(1+rvi))
    
    if(abs(log(m2/scl))>abs(log(m1/scl))) m1 <- m2
    m1
  }
  sif <- colSums(ma)
  ind <- 1:ncol(ma)
  cind <- ind[which.max(sif)]
  sfac <- rep(1,ncol(ma))
  ay <- ma[,cind]
  for(i in 1:ncol(ma))
  {
    if(i==cind) next
    indx <- ma[,i]>0|ay>0
    sfac[i] <- rowQuar(ma[indx,i],ay[indx],qper,qst,qend,qstep,qbound,mcut)
  }
  return(sfac)
}

#' Function for esitmating size factors
#' 
#' Given a matrix of count data, this function esitmates the size
#' factors by selected method.
#' It aslo provides four different methods for normalizing according to user-defined size factors,
#' total reads, up quantile (75%), qtotal (\code{\link{qtotalNormalized}}), TMM (edgeR) or geometric from DESeq (See \code{\link{estimateSizeFactorsForMatrix}}).
#'
#' @title Estimating size factors from the reads count table
#' @param object a ABSSeq object with element of 'counts' and 'normMethod', see the constructor functions
#' \code{\link{ABSDataSet}}. 
#' @return a ABSDataSet object with the estimates size factors, one element per column. Use the \code{\link{sFactors}}
#' to show it.  
#' @examples
#' 
#' data(simuN5)
#' obj <- ABSDataSet(counts=simuN5$counts, groups=factor(simuN5$groups))
#' obj <- normalFactors(obj)
#' sFactors(obj)
#'
#' @export
normalFactors <- function(object){
  if(!is(object,"ABSDataSet"))
  {
    stop("input is not an ABSDataSet object!")
  }
  validObject(object)
  sizefactor <-rep(1,ncol(object@counts))
  method=object@normMethod
  if (method == "total") {
    sizefactor<- colSums(object@counts)
  
  }
  if (method == "quartile") {
    rowQuar=function(z) {
      x <- z[z > 0]
      sum(x[x <= quantile(x, 0.75, na.rm = TRUE)], na.rm = TRUE)
    }
    cbuf <- object@counts
    cbuf <- cbuf[colSums(cbuf)>0,]
    sizefactor <- apply(cbuf,2,rowQuar)
  }
  if (method == "qtotal") {
    sizefactor <- qtotalNormalized(object@counts)
  }
  if (method == "TMM") {
    require(edgeR)
    nf <- calcNormFactors(object@counts)
    sizefactor=colSums(object@counts) * nf
  }
  if (method == "geometric") {
    sizefactor <- estimateSizeFactorsForMatrix(object@counts)
  }
  if (method == "user") {
    sizefactor<- object@sizeFactor
  }
  #message(sizefactor)
  sizemax=max(sizefactor)
  if(sizemax<=0)
    {
      stop("The counts table is wrong, with zero row or negative number!")
    }
   sizefactor=sizefactor/mean(sizefactor)
  #in case a zero sizefactor
   sizefactor[sizefactor<=0]<- 1.0
   object@sizeFactor<- sizefactor
   return(object)
}

#' Function for replacing outliers, this function is not available for user
#' 
#' @param rowd, a vector of data for trimming
#' @param madpriors, priors for moderating MAD, estimated by local regression and set to the highest standard deviation
#' @param group1 and group2, group index for conditions
#' @param baselevel1 and baselevel2, baselevel for trimming (avoid over-trimming)
#' @param len1 and len2, length of groups
#' @param spriors the weight size for madpriors
#' @param cutoff cutoff of MAD for trimming, default is 2
#' @param caseon switch for specifically dealing with sample size of 2, default is TRUE
#' Given a vector of count data, this function replaces the outliers with trimmed mean
#'
#' @export
replaceByrow=function(rowd,madpriors,group1,group2,baselevel1,baselevel2,len1,len2,spriors=2,cutoff=2.0,caseon=TRUE)
{
  b1 <- baselevel1
  b2 <- baselevel2

  indexs=0
  am <- max(median(rowd[group1]),b1)
  bm <- max(median(rowd[group2]),b2)
  indA <- any(rowd[group1]>bm)
  indB <- any(rowd[group2]>am)
  
  if(caseon&&len1<=2)
  {
    am <- max(min(rowd[group1]),b1)
    indA <- all(rowd[group1]>=max(rowd[group2]))
  }
  if(caseon && len2<=2)
  {
    bm <- max(min(rowd[group2]),b2)
    indB <- all(rowd[group2]>=max(rowd[group1]))
  }
  if(indA)
  {
    ##posterior of MAD, moderated by Poisson Dispersion
    postmad <- sqrt((mad(rowd[group1])^2*len1+madpriors^2)/(len1+spriors)+1/am)
    if(len1<=2 && caseon) postmad <- sqrt(madpriors^2+1/am)
    ab <- rowd[group1]-am-postmad*cutoff
    abuf <- rowd[group1]
    
    ##compared to median of other group, avoid over-trimmed
    indexs=indexs+sum(ab>0)
    abuf[ab>0] <- max(bm,am+postmad)
    rowd[group1]=abuf
  }
  if(indB)
  {
    postmad <- sqrt((mad(rowd[group2])^2*len1+madpriors^2)/(len2+spriors)+1/bm)
    if(len2<=2 && caseon) postmad <- sqrt(madpriors^2+1/bm)
    ab <- rowd[group2]-bm-postmad*cutoff
    indexs=indexs+sum(ab>0)
    abuf <- rowd[group2]
    abuf[ab>0] <- max(am,bm+postmad)
    rowd[group2]=abuf
  }
  return(c(rowd,indexs))
}
#' Function for replacing the outliers by MAD
#' 
#' Given a matrix of count data, this function replacing the outliers by MAD. Noticely, this function also provides part of parameters for DEs calling.
#' It is called by \code{\link{callParameter}}
#'
#' @title Replacing outliers by moderated MAD
#' @param object a ABSSeq object with element of 'counts' and 'normMethod', see the constructor functions
#' \code{\link{ABSDataSet}}.
#' @param replaceOutlier switch for replacing, default is TRUE.
#' @param cutoff cutoff of moderating MAD for outliers, default is 2
#' @param baseMean parameter for limiting the trimming at low expression level by baseMean/(sample size), default is 100.
#' @param limitMad the minimal prior for moderating MAD, default is set to 0.707, which is usually the highest standard deviation at expression level of 1
#' @param spriors prior weight size for prior MAD, default is 2
#' @param Caseon switch for dealing with outlier trimming at sample size of 2
#' @param ... reserved parameters
#' @return a ABSDataSet object with normalized counts after trimming (replaceOutlier=TRUE) or not (replaceOutlier=FALSE). Use the \code{\link{excounts}}
#' to show it. Use \code{\link{results}} with name 'trimmed' to view the trimming status.
#' @examples
#' 
#' data(simuN5)
#' obj <- ABSDataSet(counts=simuN5$counts, groups=factor(simuN5$groups))
#' obj <- normalFactors(obj)
#' obj <- ReplaceOutliersByMAD(obj)
#' head(excounts(obj))
#' head(results(obj,c("trimmed")))
#'
#' @export
ReplaceOutliersByMAD=function(object, replaceOutlier=TRUE,cutoff=2.0,baseMean=100,limitMad=0.707,spriors=2,Caseon=TRUE,...)
{
  if(!is(object,"ABSDataSet"))
  {
    stop("input is not an ABSDataSet object!")
  }
  if(is.null(sFactors(object)))
  {
    stop("Please run normalized function before 'ReplaceOutliersByMAD'!")
  }
  if(length(sFactors(object)) !=length(groups(object)) || sum(sFactors(object)<=0)>0) 
  {
    stop("Size factors are wrong, not equal with sample size or with none positive number!")
  }
  igroups <- groups(object)
  ngr1 <- igroups[1]
  ngr2 <- igroups[igroups!=igroups[1]][1]
  gr1 <- igroups==ngr1
  gr2 <- igroups==ngr2
  n1 <- sum(igroups==ngr1)
  n2 <- sum(igroups==ngr2)
  ncounts <-log2(counts(object,TRUE)+1)

  
  ##estimate highest standard deviation by local fit
  design <- model.matrix(~0+igroups)
  colnames(design)  <- levels(igroups)
  fit <- lmFit(ncounts,design)
  tvar <-(fit$sigma)*sqrt(fit$stdev.unscaled[,ngr1]^2+fit$stdev.unscaled[,ngr2]^2)
  sx <- (fit$Amean)
  sy <- tvar
  allzero <- fit$Amean <=1
  if(any(allzero)) {
    sx <- sx[!allzero]
    sy <- sy[!allzero]
  }
  
  ab <- locfit(sy[sy>0]~lp(sx[sy>0],deg=1,scale=F,nn=.5),family="gamma")
  mults <- min(2,max(predict(ab,1)*sqrt(n1+n2-2)/sqrt(2)/1,1))
  m1 <- predict(ab,1)
  object[["m1"]] <- m1
  object[["mults"]] <- mults
  object[["Amean"]] <- fit$coefficients[,ngr1]
  object[["Bmean"]] <- fit$coefficients[,ngr2]
  object[["rawFC"]] <- object[["Bmean"]]-object[["Amean"]]
  baselevel1 <- log2(baseMean/max(n1,1))
  baselevel2 <- log2(baseMean/max(n2,1))
  trimmed <- rep(0,nrow(ncounts))
  ##replace outliers
  if(replaceOutlier)
  {
    ebuf=t(apply(ncounts,1,replaceByrow,max(m1*sqrt(n1+n2-2)/sqrt(2),limitMad),gr1,gr2,baselevel1,baselevel2,n1,n2,spriors,cutoff,Caseon))
    trimmed=ebuf[,n1+n2+1]
    ebuf=ebuf[,1:(n1+n2)]
    excounts(object) <- 2^(ebuf)-1
  }else
  {
    excounts(object) <- 2^ncounts-1
  }
  object[["trimmed"]] <- trimmed
  return(object)
}

#' Calculate parameters for each gene (the moderating basemean and dispersions), without replicates
#' 
#' buliding a pseudo group to esitimate parameter by mean difference. shifted and calculate a set of parameters from normalized counts table before \code{\link{callDEs}}
#'
#' @param object a \code{\link{ABSDataSet}} object.
#'
#'
#' @return A ABSDataSet object with absolute differences, basemean, mean of each group, variance, 
#' log2 of foldchange, named as 'absD', 'baseMean', 'Amean', 'Bmean', 
#'  'Variance' and 'foldChange', respectively. Use the \code{\link{results}} to get access it
#'
#' @title Calculate parameters for differential expression test base on absolute counts differences without replicates
#' @note This function should run after \code{\link{normalFactors}} or providing size factors. This function firstly constructs an expression level depended fold-change cutoffs
#' and then separate the data into two groups. The group with fold-change less than cutoffs is used to training the dispersion. However, the cutoff might be too small when applied
#' on data set without or with less DEs. To avoid it, we set a prior value (0.5) to it.
#' @examples
#'
#' data(simuN5)
#' obj <- ABSDataSet(counts=(simuN5$counts)[,c(1,2)], groups=factor(c(1,2)))
#' obj <- normalFactors(obj)
#' obj <- callParameterwithoutReplicates(obj)
#' obj <- callDEs(obj)
#' head(results(obj))
#'
#' @export
callParameterwithoutReplicates <- function(object) {
  if(!is(object,"ABSDataSet"))
  {
    stop("input is not an ABSDataSet object!")
  }
  if(is.null(sFactors(object)))
  {
    stop("Please run normalized function before 'callParameter'!")
  }
  if(length(sFactors(object)) !=length(groups(object)) || sum(sFactors(object)<=0)>0) 
  {
    stop("Size factors are wrong, not equal with sample size or with none positive number!")
  }
  igroups <- groups(object)
  ngr1 <- igroups[1]
  ngr2 <- igroups[igroups!=igroups[1]][1]
  gr1 <- igroups==ngr1
  gr2 <- igroups==ngr2
  
  ncounts <- log2(counts(object,TRUE)+1)
  object[["Amean"]] <- ncounts[,gr1]
  object[["Bmean"]] <- ncounts[,gr2]
  object[["rawFC"]] <- object[["Bmean"]]-object[["Amean"]]
  object[["lowFC"]] <- object[["Bmean"]]-object[["Amean"]]
  trimmed <- rep(0,nrow(ncounts))
  object[["trimmed"]] <- trimmed
  lvar <- apply(ncounts,1,sd)
  lmean <- apply(ncounts,1,mean) 
  
  sx <- lmean
  sy <- lvar/sqrt(2)
  allzero <- sx <=0 
  if(any(allzero)) {
    sx <- sx[!allzero]
    sy <- sy[!allzero]
  }
  ab <- locfit(sy[sy>0]~lp(sx[sy>0],deg=1,scale=F,nn=.5),family="gamma")
  mults <- min(2,max(predict(ab,1)/1,1))
  m1 <- max(0.707,predict(ab,1))
  object[["m1"]] <- m1
  object[["mults"]] <- mults
  inds <- abs(object[["rawFC"]]) < pmax(0.5,2*pmax(m1/sqrt(pmax(lmean,1)),predict(ab, lmean)))
  
  
  
  excounts(object) <- 2^ncounts-1
  
  ex=excounts(object)

  nncounts <- ex
  
  AAmean <- nncounts[,gr1]
  BBmean <- nncounts[,gr2]

  Mmax <- pmax(AAmean,BBmean)
  
  object[["absD"]] <- round(abs(AAmean-BBmean)+0.1)

  mdis <- mean(object[["absD"]])
  LevelstoNormFC(object) <- pmin(mdis,LevelstoNormFC(object))

  mvar <- apply(nncounts[inds,],1,var)
  mmean <- apply(nncounts[inds,],1,mean)
  sx <- mmean
  sy <- (mvar-mmean)/mmean^2/2
  sy[is.na(sy)] <- 0
  if(all(sy<=0))
  {
    stop("All dispersions is <= 0! The program stops!")
  }
  ## estimate penalty of dispersion
  if(length(object@minDispersion)==0)
  {
    mcut <- 10
    scut <- (sum(sy[sx>mcut]<0)/length(sy[sx>mcut]))
    
    mults1 <- quantile(sy[sy>0.005&sx>0],sqrt(scut))
    if(scut<0.05 || scut>0.9) mults1 <- scut*0.2
    minimalDispersion(object) <- mults1
  }
  mults1 <- minimalDispersion(object)
  ## in case abnormal of scut
  if(is.na(mults1)) mults1 <- 0
  allzero <- sx <=1
  if(any(allzero)) {
    sx <- sx[!allzero]
    sy <- sy[!allzero]
  }
  
  ab1 <- locfit(sqrt(sy[sy>0])~lp(log2(sx[sy>0]),deg=1,scale=F,nn=.5),family="gamma")
  
  
  a <- Mmax
  a[a==0] <- 1
  preabs <- predict(ab1,log2(a))
  baseadd <- sqrt((a+a^2*preabs^2)*LevelstoNormFC(object)/a)

  ###shift counts and recalculate the mean and variances
  ncounts <- (nncounts+baseadd)
  
  Mmax <- Mmax+baseadd
  Mmax[Mmax==0] <- 1
  
  predis <- predict(ab1,log2(Mmax))^2
  mults1 <- max(0,mults1-min(predis)) 
  minimalDispersion(object) <- mults1
  object[["priors"]] <- 1/(pmax(predis+mults1+(m1/sqrt(8))^2,0))
  
  
  design <- model.matrix(~0+igroups)
  colnames(design)  <- levels(igroups)
  
  #in case zero counts
  #ncounts[ncounts==0] <- 1
  
  object[["foldChange"]] <- log2(ncounts[,gr2])-log2(ncounts[,gr1])
  ncounts[ncounts==0] <- 1
  tvar <- mean((apply(log2(ncounts),1,sd)/sqrt(2))[inds])
  Mmax <- pmax(Mmax,1)
  rats <- max(minRates(object),min(m1/mults/sqrt(8)/2+tvar/2,maxRates(object)))
  object[["mts"]] <- rats
  object[["baseMean"]] <- (Mmax)*rats
  object[["Variance"]] <- rep(0,length(object[["baseMean"]]))
  object[["preab"]] <- Mmax
  return(object)
}
#' Calculate baseline and penalty for uncertainty
#' 
#' 
#' Called by \code{\link{genAFold}}
#'
#' @param svar pooled variance of read count.
#' @param smean pooled mean of read count.
#' @param amean mean of group A.
#' @param bmean mean of group B.
#' @param sunc pooled uncertainty.
#' @param size sample size
#' @param preval control of variance scale in case over-scaled, default is 0.05.
#' @param qforkappa quantile for estimating kappa(>=qforkappa), default is 0 (no trimming of data).
#' @param ... reserved paramaters
#' @return A list object with baseline of uncertainty & cv and penalty of sample size
#' 
#' @title Calculate parameters for differential expression test base on absolute counts differences
#' @note Not available for user.
#' @examples
#'
preAFold <- function(svar,smean,amean,bmean,sunc,size,preval=0.05,qforkappa=0,...) {
  #estimating kappa factor (assessing variances stability)
  ind <-  smean>=1
  kam <- smean[ind]
  kvar <- svar[ind]
  ka <- sqrt(kvar)/(kam+sunc[ind])
  qka<- kam>=quantile(kam,qforkappa)
  ka <- max(1,mean(ka[qka])^2/var(ka[qka]))
  #no replicates
  if(!is.numeric(ka)) ka <-1
  #stable CVs, no need to borrow information
  if(is.infinite(ka)) 
  {
    message("The CVs of all genes are identical!")
    ka <- 1.0e9
  }
  #kappa for Poisson dis.,if closed to Poisson dis, then use sample size as kappa
  qper <- sum(kvar<=2*kam)/sum(ind)*(size-1)
  # exclude small variances for fitting
  ind <- smean>=1&svar>smean
  AAvar <- svar[ind]
  AAmean <- smean[ind]
  sx <-log(AAmean)
  sunct <-sunc[ind]
  sy <- sqrt(AAvar)/(AAmean+sunct)
  
  amin <- 1
  mmean <- pmax(amean,bmean)
  cpra <- 0
  absd <- pmax(abs(amean-bmean),amin)
  thres <- max(preval,1/sqrt(max(mmean)))
  tmax <- preval
  tmin <- preval
  if(sum(ind)>0.1*length(ind))
  {
    preb <- locfit(sy~lp(sx,deg=1,scale=F,nn=.5),family="gamma")
    cpra <- predict(preb,log(pmax(mmean,amin)))
    pre0 <- predict(preb,0)
    tmin <- min(cpra)
    tmax <- max(cpra)
    if(tmax>pre0) tmin <- tmax
    else tmax <- pre0
    ##as cv in log2
    scl <-max(1,tmin/tmax/log(2))
    tmax <- tmax*scl
    tmin <- tmin*scl
    #actual level of DE is sqrt(absd)
    cpra <- predict(preb,0.5*log(pmax(absd,amin)))
  }else
  {
    message("All variances are less than mean! Use Poisson as default!")
    cpra <- sqrt(absd)/(absd+sqrt(absd)+1)
  }
  amin <- 1/size
  apra <- (cpra/preval)^2/max(qper,ka)*(sqrt(pmax(amean,amin))+sqrt(pmax(bmean,amin)))
  thres <- min(sqrt(max(size-1,0))*preval,max(tmax/2,tmin))

  return(list(apra,thres,min(1,max(ka-size+1,0))))
}
#' Calculate aFold for each gene and general sd
#' 
#' 
#' shifted and calculate a set of parameters from normalized counts table before \code{\link{callDEs}}
#'
#' @param nncounts matrix for read count.
#' @param cond factor for conditions. If provide only one condition, fold-change estimation will be suppressed. 
#' @param preval pre-defined scale control for variance normalization, default is 0.05, a large value generally increases the fold-changes (decreases penalty of variances) under low expression.
#' @param qforkappa quantile for estimating kappa(>=qforkappa), default is 0 (without trimming of data). Please set up a value in [0,1) if you want to trim the low expressed data.
#' @param priorgenesd prior value  for general SD of fold change, if provided, the estimation of general SD will be replaced by this value.
#' @return A list with log2 foldchange, general SD for calculating pvalue and variance stablized counts
#' 
#' @title Calculate parameters for differential expression test base on absolute counts differences
#' @note This function should run after \code{\link{normalFactors}}.
#' @examples
#'
#' data(simuN5)
#' obj <- ABSDataSet(counts=simuN5$counts, groups=factor(simuN5$groups))
#' mtx <- counts(obj,TRUE)
#' aFold <- genAFold(mtx,factor(simuN5$groups))
#' hist(aFold[[1]])
#'
#' @export
genAFold <- function(nncounts,cond,preval=0.05,qforkappa=0,priorgenesd) {
  if(is.null(ncol(nncounts)))
  {
    stop("Please input a matrix as nncounts!")
  }
  if(!missing(priorgenesd)&&(!is.numeric(priorgenesd)||priorgenesd<=0)) 
    stop("Please positive value for priorgenesd!")
  if(ncol(nncounts)!=length(cond))
  {
    stop("The columns of nncounts and length of cond differ!")
  }
  if(!is.numeric(qforkappa)||qforkappa<0||qforkappa>1)
    stop("qforkappa should be in [0,1)!")
  if(!is.numeric(preval)||preval<0)
    stop("preval should be positive!")
  igroups <- cond
  ngr1 <- igroups[1]
  ngr2 <- NULL
  if(!all(igroups==ngr1)) ngr2 <- igroups[igroups!=igroups[1]][1]
  gr1 <- igroups==ngr1
  gr2 <- FALSE
  if(!all(igroups==ngr1)) gr2 <- igroups==ngr2
  n1 <- sum(gr1)
  n2 <- sum(gr2)
  if(n1<2 &&n2<2)
  {
    message("No replicates! Use Poisson as default!")
  }
  AAmean <- 0
  AAvar <- 0
  ind <- rowSums(nncounts)>0
  tmean <- c()
  tvar <- 0
  if(n1==1) AAmean <- nncounts[,gr1]
  else if(n1>1)
  {
    AAmean <- apply(nncounts[,gr1],1,mean)
    AAvar <- apply(nncounts[,gr1],1,var)
    if(all(AAvar==0)) message("Counts in sample identical to each other in group:", ngr1)
  }
  BBmean <- 0
  BBvar <- 0
  if(n2==1) BBmean <- nncounts[,gr2]
  else if(n2>1)
  {
    BBmean <- apply(nncounts[,gr2],1,mean)
    BBvar <- apply(nncounts[,gr2],1,var)
    if(all(BBvar==0)) message("Counts in samples identical to each other in group:", ngr2)
  }
 
  if(n1<=1 &&n2>1) {
    tmean <- BBmean[ind]
    tvar <- BBvar[ind]
  }else if(n1>1 &&n2<=1)
  {
    tmean <- AAmean[ind]
    tvar <- AAvar[ind]
  }else if(n1>1 && n2>1)
  {
    tmean <- AAmean[ind]
    tmean <- c(tmean,BBmean[ind])
    tvar <- AAvar[ind]
    tvar <- c(tvar,BBvar[ind])
  }
  if(n1==1&&n2==1)
  {
    tmean <- c(AAmean[ind],BBmean[ind])
    tvar <- tmean
  }

  ##observed uncertainty
  varunca <- 0
  varuncb <- 0
  varuncc <- 0
  AAvar<-pmax(AAvar,AAmean)
  asiz <- sqrt(AAvar)/AAmean
  asiz[AAmean<=0]<-0
  varunca <- sqrt(AAvar)*(1+asiz)
  BBvar<-pmax(BBvar,BBmean)
  bsiz <- sqrt(BBvar)/BBmean
  bsiz[BBmean<=0]<-0
  varuncb <- sqrt(BBvar)*(1+bsiz)
  varund <- 0
  if(n1<=1&&n2>1) 
  {
    varund <- varuncb[ind]
  }
  else if(n2<=1&&n1>1){
    varund <- varunca[ind]
  }else
  {
    varund <- varunca+varuncb
    varund <- c(varund[ind],varund[ind])
  }
  ##baseline
  hideunc <- preAFold(tvar,tmean,AAmean,BBmean,varund,max(n1,n2),preval,qforkappa)
  precut <- hideunc[[2]]
  varunc <- varunca+varuncb+hideunc[[1]]
  totunc <- pmax(varunc,1.0e-9)
  scounts <- nncounts+totunc
  scounts <- log2(scounts)

  genSDA <- 0
  genSDB <- 0
  logAmean <- 0
  logAsd <- NULL
  if(n1==1) logAmean <- scounts[,gr1]
  else if(n1>1)
  {
    logAsd <- apply(scounts[,gr1],1,sd)
    logAmean <- apply(scounts[,gr1],1,mean)
    
  }
  
  logBmean <- 0
  logBsd <- NULL

  if(n2==1) logBmean <- scounts[,gr2]
  else if(n2>1)
  {
    logBsd <- apply(scounts[,gr2],1,sd)
    logBmean <- apply(scounts[,gr2],1,mean)

  }
  
  #final fold-change
  fold <- 0
  sdf <- 0
  if(!(n1==0||n2==0))
  {
    fold <- (logBmean-logAmean)
    sdf <- sd(fold)
    ##large DE will balance the CVs and thus less baseline
    precut <- precut/2^(max(sdf-preval*sqrt(2),0))
    if(n1>1)
    {
      ind <- logAsd>=precut
      qnum <- sum(ind)
      if(qnum>length(ind)*0.1||qnum>50) logAsd <- logAsd[ind]
      else logAsd <- rep(precut,length(logAmean))
      genSDA <- mean(logAsd)
    }
    if(n2>1)
    {
      ind <- logBsd>=precut
      qnum <- sum(ind)
      if(qnum>length(ind)*0.1||qnum>50) logBsd <- logBsd[ind]
      else logBsd <- rep(precut,length(logBmean))
      genSDB <- mean(logBsd)
    }
  }
  # estimate general SD

  genesd <- preval*sqrt(2)#as Poisson distribution
  if(!(genSDA==0 &&genSDB==0))
  {
    ns1 <- max(1,n1-1)+hideunc[[3]]
    ns2 <- max(1,n2-1)+hideunc[[3]]
    genSDA <- max(genSDA,genSDB)
    genesd <- genSDA*sqrt(1/ns1+1/ns2)
  }

  if(!missing(priorgenesd)) geneSD <- priorgenesd

  return(list(aFold=fold,geneSD=genesd,nncounts+totunc))
}
#' Calculate parameters for each gene (the moderating basemean, dispersions, moderated fold-change and general sd)
#' 
#' 
#' shifted and calculate a set of parameters from normalized counts table before \code{\link{callDEs}}
#'
#' @param object a \code{\link{ABSDataSet}} object.
#' @param replaceOutliers switch for outlier replacement, default is TRUE.
#' @param ... parameters past to \code{\link{ReplaceOutliersByMAD}}
#'
#' @return A ABSDataSet object with absolute differences, basemean, mean of each group, variance, 
#' log2 of foldchange, named as 'absD', 'baseMean', 'Amean', 'Bmean', 
#'  'Variance' and 'foldChange', respectively. Use the \code{\link{results}} to get access it and \code{\link{plotDifftoBase}} to plot it.
#'
#' @title Calculate parameters for differential expression test base on absolute counts differences
#' @note This function should run after \code{\link{normalFactors}} or providing size factors.
#' @examples
#'
#' data(simuN5)
#' obj <- ABSDataSet(counts=simuN5$counts, groups=factor(simuN5$groups))
#' obj <- normalFactors(obj)
#' obj <- callParameter(obj)
#' head(results(obj,c("foldChange","absD","baseMean")))
#' plotDifftoBase(obj)
#'
#' @export
callParameter <- function(object,replaceOutliers=TRUE,...) {
  if(!is(object,"ABSDataSet"))
  {
    stop("input is not an ABSDataSet object!")
  }
   if(is.null(sFactors(object)))
  {
    stop("Please run normalized function before 'callParameter'!")
  }
  if(length(sFactors(object)) !=length(groups(object)) || sum(sFactors(object)<=0)>0) 
  {
   stop("Size factors are wrong, not equal with sample size or with none positive number!")
  }
  
  igroups <- groups(object)
  ngr1 <- igroups[1]
  ngr2 <- igroups[igroups!=igroups[1]][1]
  gr1 <- igroups==ngr1
  gr2 <- igroups==ngr2
  n1 <- sum(igroups==ngr1)
  n2 <- sum(igroups==ngr2)
  
  object <- ReplaceOutliersByMAD(object,replaceOutlier=replaceOutliers,...)
  ex=excounts(object)
  mults <- object[["mults"]]
  m1 <- object[["m1"]]
  nncounts <- ex

  AAmean <- 0
  AAvar <- 0
  if(n1==1) AAmean <- nncounts[,gr1]
  else
  {
    AAvar <- apply(nncounts[,gr1],1,var)
    AAmean <- apply(nncounts[,gr1],1,mean)
  }
  BBmean <- 0
  BBvar <- 0
  if(n2==1) BBmean <- nncounts[,gr2]
  else
  {
    BBvar <- apply(nncounts[,gr2],1,var)
    BBmean <- apply(nncounts[,gr2],1,mean)
  }
  Mmax <- pmax(AAmean,BBmean)
  
  totalvar <- AAvar+BBvar
  ##paired
  if(paired(object))
  {
    totalvar <- apply(nncounts[,gr2]-nncounts[,gr1],1,var)
  }
  ##check variance
  if(all(totalvar==0)) 
  {
    stop("Variance across samples is 0! Check it! The program stops!")
  }
  
  sx <- Mmax
  sy <- totalvar#AAvar+BBvar
  sy <- (sy-Mmax)/Mmax^2
  LevelstoNormFC(object) <- min(sqrt(mean(totalvar)/2),LevelstoNormFC(object))#min(sqrt(mean(AAvar+BBvar)/2),LevelstoNormFC(object))
  ##in case 0 counts for all smaples
  sy[is.na(sy)] <- 0
  if(all(sy<=0))
  {
    stop("All dispersions is <= 0! The program stops!")
  }
  ## estimate penalty of dispersion
  if(length(object@minDispersion)==0)
  {
    mcut <- 10/max((max(n1,n2)-1),1)
    scut <- (sum(sy[sx>mcut]<0)/length(sy[sx>mcut]))
    
    mults1 <- quantile(sy[sy>0.005&sx>0],sqrt(scut))
    if(scut<0.05 || scut>0.9) mults1 <- scut*0.2
    minimalDispersion(object) <- mults1
  }
  mults1 <- minimalDispersion(object)
  ## in case abnormal of scut
  if(is.na(mults1)) mults1 <- 0
  allzero <- sx <=1
  if(any(allzero)) {
    sx <- sx[!allzero]
    sy <- sy[!allzero]
  }
  
  ab1 <- locfit(sqrt(sy[sy>0])~lp(log2(sx[sy>0]),deg=1,scale=F,nn=.5),family="gamma")
  
  a <- Mmax
  a[a==0] <- 1
  preabs <- predict(ab1,log2(a))
  
  baseadd <- sqrt((a+a^2*preabs^2)*LevelstoNormFC(object)/a/max(1,(max(n1,n2)-1)))
  
  
  ###shift counts and recalculate the mean and variances
  ncounts <- (nncounts+baseadd)
  Mmax <- Mmax+baseadd
  Mmax[Mmax==0] <- 1
  #excounts(object)=nncounts+ggad
  ###moderate fold-change
  
  afoldpara <- genAFold(nncounts,igroups,...)
  object[["genesd"]] <- afoldpara[[2]]
  object[["foldChange"]] <- afoldpara[[1]]
  ###paired
 
  abuf <- (totalvar-Mmax)/Mmax^2#(AAvar+BBvar-Mmax)/Mmax^2
  
  predis <- predict(ab1,log2(Mmax))^2/(max(n1,n2))

  mults1 <- max(0,mults1-min(predis))
  minimalDispersion(object) <- mults1
  object[["priors"]] <- 1/(pmax(predis+mults1+pmax(abuf,0),0))
  object[["absD"]] <- round(max(n1,n2)*abs(AAmean-BBmean)+0.1)
 
  
  design <- model.matrix(~0+igroups)
  colnames(design)  <- levels(igroups)
  
  #in case zero counts
  ncounts[ncounts==0] <- 1
  
  fit <- lmFit(log2(ncounts),design)
  #object[["foldChange"]] <- fit$coefficients[,ngr2]-fit$coefficients[,ngr1]
  object[["lowFC"]] <- fit$coefficients[,ngr2]-fit$coefficients[,ngr1]
  tvar <- mean((fit$sigma)*sqrt(fit$stdev.unscaled[,ngr1]^2+fit$stdev.unscaled[,ngr2]^2))
  Mmax <- pmax(Mmax,1)

  rats <- max(minRates(object),min(m1/mults/sqrt(8)/2+tvar/2,maxRates(object)))
  object[["mts"]] <- rats
  object[["baseMean"]] <- ((Mmax)*max(n1,n2)+sqrt(max(n1,n2)*pmin(Mmax,(totalvar))))*rats#((Mmax)*max(n1,n2)+sqrt(max(n1,n2)*pmin(Mmax,(AAvar+BBvar))))*rats
  object[["Variance"]] <- totalvar #AAvar+BBvar
  object[["preab"]] <- Mmax*max(n1,n2)
  return(object)
}
#' Using NB distribution to calculate p-value for each gene as well as adjust p-value
#'
#' This function firstly calls p-value used \code{\link{pnbinom}} to call pvalue based on sum of counts difference between two
#' groups or used \code{\link{pnorm}} to call pvalue via log2 fold-change, then adjusts the pvalues via \code{\link{p.adjust}} method. In addition, it also shrink the log2 fold-change towards a common dispersion
#' after pvalue calling.
#'
#' @title Testing the differential expression by counts difference
#' @param object an \code{\link{ABSDataSet}} object.
#' @param adjmethod the method for adjusting p-value, default is 'BH'. For details, see \code{\link{p.adjust.methods}}.
#' @param useaFold switch for DE detection through fold-change, which will use a normal distribution (N(0,sd)) to test the significance of log2 fold-change. The sd is estimated through a quantile function of gamma distribution at \code{\link{callParameter}}.
#' 
#' @return an \code{\link{ABSDataSet}} object with additional elements: shrinked log2 fold-change, pvalue and adjusted p-value,
#' denoted by foldChange pvalue and adj-pvalue, respectively. Use the \code{\link{results}} method to get access it. 
#' @note this function should run after \code{\link{callParameter}}
#' @examples
#'
#' data(simuN5)
#' obj <- ABSDataSet(counts=simuN5$counts, groups=factor(simuN5$groups))
#' obj <- normalFactors(obj)
#' obj <- callParameter(obj)
#' obj <- callDEs(obj)
#' head(results(obj))
#' @export
callDEs <- function(object,adjmethod="BH",useaFold=FALSE) {
  if(!is(object,"ABSDataSet"))
  {
    stop("input is not an ABSDataSet object!")
  }
  if(is.null(object[["absD"]]) | is.null(object[["baseMean"]]) | is.null(object[["Variance"]]))
  {
    stop("Please run 'callParameter' function before 'callDEs'!")
  }
  if(sum(adjmethod %in% p.adjust.methods)==0)
  {
    stop("Please provide a correct p-value adjust method according p.adjust.methods.")
  }
  if(useaFold==TRUE && is.null(object[["genesd"]]))
  {
    stop("To detect DE on fold-change, please run 'callParameter' function before 'callDEs'!")
  }
  if(useaFold==TRUE && object[["genesd"]]==0)
  {
    stop("To detect DE on fold-change, general sd should not be 0! Please check it!")
  }
  absD <- object[["absD"]]
  absF <- abs(object[["foldChange"]])
  if(useaFold)
  {
    message("Generating pvalues via aFold....")
    object[["pvalue"]] <- pnorm(absF/object[["genesd"]],lower.tail=F)*2
  }
  else object[["pvalue"]] <- pnbinom(absD,mu=object[["baseMean"]],size=object[["priors"]],lower.tail=F)
  object[["adj.pvalue"]] <- p.adjust(object[["pvalue"]],method=adjmethod[1]);
  ##get the end time
  object@ends <- Sys.time()
  return(object);
}

