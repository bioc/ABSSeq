
// includes from the plugin
#include <Rcpp.h>


#ifndef BEGIN_RCPP
#define BEGIN_RCPP
#endif

#ifndef END_RCPP
#define END_RCPP
#endif

using namespace Rcpp;


// user includes
#include <R.h>
#include <Rmath.h>
#include <stdio.h>
#include <time.h>
#include <vector>
#include <R_ext/Utils.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>
#include <R_ext/Visibility.h>

inline 
	double logGenPossion(int y,double mu,double alpha)
{
	double result,logyfac=lgamma(y+1);
	result = log(mu) - log(1+alpha*mu);
	result = y*result + (y-1)*log(1+alpha*y)-logyfac;
	result = result -mu*(1+alpha*y)/(1+alpha*mu);
	return result;
}
double logLhood(double a,double b,Rcpp::IntegerVector sp, Rcpp::IntegerVector tranData)
{
	double loglike=0;
	int lens=tranData.size();
	for(int i=0;i<lens;i++)
	{
		loglike+=logGenPossion(tranData[i],a,b)*sp[tranData[i]];
	}
	return loglike;
}
double unran(int *na, int *nb, int *nc)
{
	double random;
	*na=(171*(*na))%30269;
	*nb=(172*(*nb))%30307;
	*nc=(170*(*nc))%30323;
	random=(double) *na/30269.0+(double) *nb/30307.0+(double) *nc/30323.0;
	random=random-floor(random);
	return random;
}
void parameterEstimate(double &apra,double &bpra,Rcpp::IntegerVector sp, Rcpp::IntegerVector tranData)
	{// use Metropolis to get MLE of GP distribution parameters from BayesPeak. 		
		srand(time(NULL));
		int j;
		int nn;
		double aa[1000],bb[1000],sout[1000];
		double stepsize;
		int na,nb,nc;
		double rando,test,ratio;
		double aatry,bbtry;
		double aamax, bbmax;
		double oldlike, newlike;

		na=rand() +1;
		nb=rand() -1;
		nc=rand() ;
		nn = 1000;
		stepsize = 0.1;
		aa[0] = apra;
		bb[0] = bpra; 
		sout[0] = logLhood(aa[0],bb[0],sp,tranData);
		oldlike = sout[0];
		aamax = apra;
		bbmax = bpra;

		for(j=1;j<nn;j++)
		{
			rando = unran(&na, &nb, &nc);	
			aatry = aa[j-1] + 5*stepsize *(2.0*rando-1.0);
			bbtry = bb[j-1] + stepsize * (2.0*rando-1.0);
			newlike = logLhood(aatry,bbtry,sp,tranData);
			ratio = exp(newlike-oldlike);
			if(ratio>1)
			{
				aa[j] = aatry;
				bb[j] = bbtry;
				sout[j] = newlike;;
				aamax = aatry;
				bbmax = bbtry;
				oldlike = newlike;
			}
			else
			{
				test = unran(&na, &nb, &nc);
				if(test<ratio)
				{
		     		aa[j] = aatry;
     				bb[j] = bbtry;
					sout[j] = newlike;
					oldlike = newlike;
				}
				else
				{
					aa[j] = aa[j-1];
     				bb[j] = bb[j-1];
					sout[j] = sout[j-1];
				}
			}//end of else
		}//end of j
		//exit(0);
		apra=aamax;
		bpra=bbmax;
	}//end of parameterEstimate



double calPvalue(int ccs, double gam,double alpha,Rcpp::NumericVector nbm)
{
	
	if(ccs<1) return 1.0;
	double p0=logGenPossion(ccs,gam,alpha),p1,sums=0;
	//cout<<alps<<"\t"<<bes<<"\t"<<p0<<"\t"<<gammaln(0.1)<<endl;

	//cout<<p0<<endl;
	for(int j=0;j<=ccs;j++)
	{
		p1=logGenPossion(j,gam,alpha)+nbm[ccs-j]-p0;
		if(p1>500)
		{
			return 0.0;
		}
		sums+=exp(p1);
	}

	return 1.0/(1.0+sums);
}


// declarations
extern "C" {
SEXP fitGP( SEXP dMean);
SEXP getPvalue(SEXP dmean, SEXP logm, SEXP gams, SEXP alphas);
static const R_CallMethodDef callMethods[]  = {
	{"getPvalue", (DL_FUNC) &getPvalue, 4},
	{"fitGP", (DL_FUNC) &fitGP, 1},
   {NULL, NULL, 0}
};
void attribute_visible R_init_ABSSeq(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, callMethods, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
    R_forceSymbols(dll, TRUE);
}
}
SEXP fitGP(SEXP dMean)
{
	BEGIN_RCPP
	//Rcpp::IntegerVector y(dMean);
	Rcpp::NumericVector dy(dMean);
	Rcpp::NumericVector::iterator it=std::max_element(dy.begin(),dy.end());
	int maxc=ceil(*it);
	Rcpp::IntegerVector sp=rep(0,maxc+1),trans;
	for(it = dy.begin();it!=dy.end();++it)
	{
		//if(int(*it)>0) sp[int(*it)]++;
		sp[round(*it)]++;
	}
	for(int i=0;i<=maxc;i++)
	{
		if(sp[i]>0) trans.push_back(i);
	}
	double dmean=Rcpp::mean(dy), dvar=Rcpp::var(dy);
	dvar=(sqrt(dvar/dmean)-1)/dmean;
	if( dmean==0 || dvar<0) dvar=0;
	parameterEstimate(dmean,dvar,sp,trans);
	Rcpp::NumericVector logps(maxc+1);
	for(int i=0;i<=maxc;i++)
	{
		logps[i]=logGenPossion(i,dmean,dvar);
	}
	return Rcpp::List::create(Rcpp::Named("alp",dmean),
			  Rcpp::Named("betas",dvar),
			  Rcpp::Named("logp",logps));
END_RCPP
}

SEXP getPvalue(SEXP dmean, SEXP logm, SEXP gams, SEXP alphas)
{
BEGIN_RCPP
	Rcpp::NumericVector y(dmean);
	Rcpp::NumericVector pvec(logm),gvec(gams),avec(alphas),pvalues(y.size());

	for(int i=0;i<y.size();i++)
	{
		pvalues[i]=calPvalue(ceil(y[i]),gvec[i],avec[i],pvec);
	}

	return pvalues;
END_RCPP
}
