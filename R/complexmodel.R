#' This function performs a default analysis by calling, in order, the functions:
#' \code{\link{normalFactors}},
#' \code{\link{aFoldcomplexDesign}},
#'
#' This function uses a linear model (limma-lmFit) to infer DE under complex design.
#' 
#' @title Differential expression analysis for complex desgin.
#' @param object a \code{\link{ABSDataSet}} object (not need 'groups' information).
#' @param design a numeric matrix for expriment, with samples and factors in rows and colnums, respectively. Design respresents the satuarated model.
#' @param condA a vector of factors for DE analysis, which could be redundant, see \code{\link{aFoldcomplexDesign}}.
#' @param condB a vector of factors for DE analysis, which could be redundant, default is null, if not provide, the DE analysis will switch to 
#' assess difference across factors in condA (analysis of variance). If provide, DE analysis will focus on contrast between condB and condA (condB-condA).
#' See \code{\link{aFoldcomplexDesign}}. The unique factors in condA+condB represents the reduced model.
#' @param lmodel switch of fit linear model from limma-lmFit under design, default is TRUE. If TRUE, a gene-specific residual varaince will
#' be estimated from (satuarated model - reduced model). Satuarated model includes all factors in design matrix and reduced model includes factors in condA+condB.
#' if satuarated model == reduced model, the DE analysis performs pairwise comparison or one-way analysis of variance. See \code{\link{aFoldcomplexDesign}}.
#' @param adjmethod defualt is 'BH', method for p-value adjusted, see \code{\link{p.adjust.methods}} for details
#' @param preval parameter for \code{\link{aFoldcomplexDesign}}, prior value for controlling of variance scale in case over-scaled, default is 0.05,
#' @param qforkappa parameter for \code{\link{aFoldcomplexDesign}}, quantile for estimating kappa(>=qforkappa), default is 0 (no trimming of data).
#' @param scale switch for scaling fold change according to common SD under log2 transformation, default is FALSE.
#' @param quiet default is FALSE, whether to print messages at each step
#' @param ... parameters passed to lmFit in limma
#' @return a result table with additional elements, including:
#' basemean, log of basemean, 
#' foldChange, shrinked (expression level and gene-specific) log2 of fold-change, B - A, or (SDs under log2 for analysis of variance)
#' pvalue, pvalue from NB distribution model, 
#' p.adj, adjuested p-value used p.adjust method.
#' scaledlogFC, scaled logFC if scale=TRUE.
#' @author Wentao Yang
#' @references Wentao Yang, Philip Rosenstiel & Hinrich Schulenburg: ABSSeq: a new RNA-Seq analysis method based on modelling absolute expression differences
#' @examples
#' 
#' data(simuN5)
#' groups=factor(simuN5$groups)
#' obj <- ABSDataSet(counts=simuN5$counts)
#' design <- model.matrix(~0+groups)
#' res <- ABSSeqlm(obj,design,condA=c("groups0"),condB=c("groups1"))
#' head(res)
#' @export
ABSSeqlm <- function(object,design,condA,condB=NULL,lmodel=TRUE,preval=0.05,qforkappa=0,adjmethod="BH",scale=FALSE,quiet=FALSE,...) {
  if(!is.matrix(design)) stop("design is supposed to be a numeric matrix!")
  if(ncol(design)<2) stop("design matrix contains only one factor!")
  if(any(!condA%in%colnames(design))||((!is.null(condB))&&any(!condB%in%colnames(design)))) stop("condA or condB not match with design matrix!")
  if (!quiet) message("eistimating size factors....")
  object <- normalFactors(object)
  ncounts <- counts(object,TRUE)
  res <- aFoldcomplexDesign(ncounts,design,condA,condB,lmodel,preval,qforkappa,...)
  logFC <- res[[1]]
  genesd <- res[[2]]
  basemean <- res[[4]]
  scls <- genesd/res[[5]]
  scaledlogFC <- logFC/scls
  pvalue <- pnorm(abs(logFC)/genesd,lower.tail=F)*2
  p.adj <- p.adjust(pvalue,method=adjmethod[1])
  restab <- cbind(basemean,logFC,pvalue,p.adj)
  if(scale) restab <- cbind(restab,scaledlogFC)
  rownames(restab) <- rownames(ncounts)
  return(restab)
}
#' Calculate mean, variance and uncertainty for a group of samples
#' 
#' 
#' Called by \code{\link{aFoldcomplexDesign}}
#'
#' @param ncounts normalized counts
#' @param cond name of group
#' @return A list object with uncertainty, mean and variance for each gene
#' 
#' @title Calculate mean, variance and uncertainty for a group of samples
#' @note Not available for user.
#' @examples
#'
callPergroup <- function(ncounts,cond)
{
  n1 <- ncol(ncounts)
  if(is.null(n1)) n1 <-1

  AAmean <- 0
  AAvar <- 0
  
  if(n1==1) AAmean <- ncounts
  else if(n1>1)
  {
    AAmean <- apply(ncounts,1,mean)
    AAvar <- apply(ncounts,1,var)
    if(all(AAvar==0)) message("Counts in sample identical to each other in group:", cond)
  }
  varunca <- 0
  AAvar<-pmax(AAvar,AAmean)
  asiz <- sqrt(AAvar)/AAmean
  asiz[AAmean<=0]<-0
  varunca <- sqrt(AAvar)*(1+asiz)

  return(list(varunca,AAmean,AAvar))
}
#' Calculate baseline and penalty for uncertainty
#' 
#' 
#' Called by \code{\link{aFoldcomplexDesign}}
#'
#' @param svar pooled variance of read count.
#' @param smean pooled mean of read count.
#' @param basem base difference across all factors.
#' @param mmean base mean of all factors.
#' @param sqrtmean average value of sqrt mean across factors
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
preAFoldComplex <- function(svar,smean,basem,mmean,sqrtmean,sunc,size,preval=0.05,qforkappa=0) {
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
  cpra <- 0
  absd <- pmax(basem,amin)
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
  apra <- (cpra/preval)^2/max(qper,ka)*sqrtmean
  thres <- min(sqrt(max(size-1,0))*preval,max(tmax/2,tmin))
  
  return(list(apra,thres,min(1,max(ka-size+1,0))))
}
#' Calculate aFold for each gene and general sd
#' 
#' 
#' shifted and calculate a set of parameters from normalized counts table
#'
#' @param nncounts matrix for read count.
#' @param design a numeric matrix for expriment, with samples and factors in rows and colnums, respectively.
#' @param condA a vector of factors for DE analysis, which could be redundant.
#' @param condB a vector of factors for DE analysis, which could be redundant, default is null. If not provide, the DE analysis will switch to 
#' assess difference across factors in condA (analysis of variance). If provide, DE analysis will focus on contrast between condB and condA (condB-condA). 
#' @param lmodel switch of fit linear model from limma-lmFit under design, default is TRUE. If TRUE, a gene-specific residual varaince will
#' be estimated from (satuarated model - reduced model). Satuarated model includes all factors in design matrix and reduced model includes factors in condA+condB.
#' @param preval pre-defined scale control for variance normalization, default is 0.05, a large value generally increases the fold-changes (decreases penalty of variances) under low expression.
#' @param qforkappa quantile for estimating kappa(>=qforkappa), default is 0 (without trimming of data). Please set up a value in [0,1) if you want to trim the low expressed data.
#' @param priorgenesd prior value  for general SD of fold change, if provided, the estimation of general SD will be replaced by this value.
#' @param ... parameters passed to lmFit in limma
#' @return A list with log2 foldchange, general SD (gene-specific SD if lmodel is TRUE) for calculating pvalue, variance stablized counts and basemean
#' 
#' @title Calculate parameters for differential expression test base on absolute counts differences
#' @note This function should run after \code{\link{normalFactors}}.
#' @examples
#'
#' data(simuN5)
#' groups=factor(simuN5$groups)
#' obj <- ABSDataSet(counts=simuN5$counts, groups=factor(simuN5$groups))
#' mtx <- counts(obj,TRUE)
#' design <- model.matrix(~0+groups)
#' aFold <- aFoldcomplexDesign(mtx,design,condA=c("groups0"),condB=c("groups1"))
#' hist(aFold[[1]])
#'
#' @export
aFoldcomplexDesign <- function(nncounts,design,condA,condB=NULL,lmodel=TRUE,preval=0.05,qforkappa=0, priorgenesd,...) {
  if(is.null(ncol(nncounts)))
  {
    stop("Please input a matrix as nncounts!")
  }
  if(ncol(nncounts)!=nrow(design))
  {
    stop("The columns of nncounts and length of design not match!")
  }
  if(!is.numeric(qforkappa)||qforkappa<0||qforkappa>1)
    stop("qforkappa should be in [0,1)!")
  if(!is.numeric(preval)||preval<0)
    stop("preval should be positive!")
  if(length(condA)<2&&is.null(condB)) stop("If not provide 'condB', require at least two factors for condA!")
  cond <- unique(condA)
  if(!is.null(condB)) cond <- unique(c(cond,unique(condB)))
  if(any(!cond%in%colnames(design))) stop("Factor is not present in desgin matrix")
  colnum <- length(cond)
  if(all(colSums(design[,cond])<2)) stop("No replicates for factors in condA and condB!")
  if(is.null(condB)) message("Run analysis of variance as \n\tcondA: ",paste(condA,collapse=", "))
  if(!is.null(condB)) message("Run analysis of contrast between \n\t-condA: ",paste(condA,collapse=", "),"\n\t+condB: ",paste(condB,collapse=", "))
  if(colnum!=length(condA)+length(condB)) message("Noticely: redundant factors in condA/condB!")
  gsize <- nrow(nncounts)
  tvar <- c()
  tmean <- c()
  mmean <- matrix(c(0),nrow=gsize,ncol=colnum)
  smean <- 0
  asize <- 1
  aveunc <- 0
  tnum <- 0
  rcond <- c()
  for(i in 1:colnum)
  {
    con <- design[,cond[i]]>0
    csize <- sum(con)
    paras <- callPergroup(nncounts[,con],cond[i])
    asize <- max(asize,csize)
    if(csize>1)
    {
      tmean <- c(tmean,paras[[2]])
      tvar <- c(tvar,paras[[3]])
      tnum <- tnum+1
      rcond <- c(rcond,cond[i])
    }
    aveunc <- aveunc+paras[[1]]
    mmean[,i] <- paras[[2]]
    smean <- smean+sqrt(paras[[2]])
  }
  aveunc <-aveunc/max(colnum-1,1)
  smean <- smean/max(colnum-1,1)
  tunc <- c()
  for(i in 1:colnum) tunc <-c(tunc,aveunc)
  basem <- 0
  if(colnum>2) basem <- apply(mmean,1,sd)
  else if(colnum==2) basem <- abs(mmean[,2]-mmean[,1])
  else basem <- mmean
  mmean <- apply(mmean,1,max)
 
  ##baseline
  hideunc <- preAFoldComplex(tvar,tmean,basem,mmean,smean,tunc,asize,preval,qforkappa)
  precut <- hideunc[[2]]
  varunc <- aveunc+hideunc[[1]]
  totunc <- pmax(varunc,1.0e-9)
  scounts <- nncounts+totunc
  scounts <- log2(scounts)
  
  gfold <- matrix(0,nrow=gsize,ncol=colnum)
  colnames(gfold) <- cond
  gsd <- matrix(0,nrow=gsize,ncol=tnum)
  colnames(gsd) <- rcond
  tnum <-0
  genesd <- 0
  for(i in 1:colnum)
  {
    con <- design[,cond[i]]>0
    csize <- sum(con)
    bufm <- scounts[,con]
    if(csize>1)
    {
      bufm <- apply(scounts[,con],1,mean)
      bufv <- apply(scounts[,con],1,sd)
      gsd[,cond[i]] <- bufv
    }
    gfold[,cond[i]] <- bufm
    ns1 <- max(1,csize-1)+hideunc[[3]]
    tnum <- tnum+1/ns1
  }
  
  #final fold-change
  fold <- 0
  scsd <- 0
  if(is.null(condB)) fold <- apply(gfold[,condA],1,sd)
  if(!is.null(condB))
  {
    mA <- gfold[,condA]
    mB <- gfold[,condB]
    vA <- 0
    vB <- 0
    if(length(condA)>1)
    {
      vA <- apply(mA,1,var)
      mA <- apply(mA,1,mean)
    }
    if(length(condB)>1)
    {
      vB <- apply(mB,1,var)
      mB <- apply(mB,1,mean)
    }
    fold <- mB-mA
    scsd <- vA/max(1,length(condA)-1)+vB/max(1,length(condB)-1)
  }
  sdf <- sd(fold)
  ##large DE will balance the CVs and thus less baseline
  precut <- precut/2^(max(sdf-preval*sqrt(2),0))
  for(i in 1:length(rcond))
  {
      bufv <- gsd[,i]
      logAsd <- bufv
      ind <- logAsd>=precut
      qnum <- sum(ind)
      if(qnum>length(ind)*0.1||qnum>50) logAsd <- logAsd[ind]
      else logAsd <- precut
      genesd <- max(genesd,mean(logAsd))
  }

  if(genesd==0) genesd <- preval*sqrt(2)
  
  ##linear model to estimate residual of variance per gene (satuarated model - reduced model)
  prior <- 0
  if(colnum<ncol(design)&&missing(priorgenesd)&&lmodel)
  {
    designr <- design[,cond]
    inds <- rowSums(designr)>0
    designr <- designr[inds,]
    fulfit <- lmFit(scounts,design,...)
    redfit <- lmFit(scounts[,inds],designr)
    prior <- (fulfit$sigma*fulfit$stdev.unscaled[,cond])^2-(redfit$sigma*redfit$stdev.unscaled[,cond])^2
    prior[prior<0] <-0
    prior <- rowSums(prior)
  }
  genesd <- genesd*sqrt(tnum/max(colnum-1,1))
  scls <- genesd
  genesd <- sqrt(genesd^2+prior/max(colnum-1,1)+scsd)
  if(!missing(priorgenesd)){
    genesd <- priorgenesd
    scls <- genesd
  }
  return(list(aFold=fold,geneSD=genesd,nncounts+totunc,log2(mmean+1),scls))
}
