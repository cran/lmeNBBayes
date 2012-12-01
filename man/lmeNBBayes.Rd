\name{lmeNBBayes}
\alias{lmeNBBayes}
\alias{pointsgamma}
\alias{plotgamma}
\alias{plotbeta}
\alias{plotnbinom}
\alias{slim}
\alias{UsedPara}
\alias{plotGs}
\alias{paralnorm}
\alias{newCat}
\alias{getSample}
\alias{useSamp}
\alias{inpheatmap}
\alias{colmeansd}
\alias{getM}
\alias{Nuniq}
\alias{hm}
\alias{gsLabelnp}
\alias{repeatAsID}

\title{
Fitting a flexible mixed effect negative binomial regression.
}


\description{
  This function generates posterior samples from the negative binomial regression mixed effect model, where the patient-specific random effect is assumed to be from Dirichlet process mixture of beta distributions.
  The hierarchy of the model is as follows:
  
  Yij | gi, beta ~ NB(Yij; size=exp(Xij*beta),prob=gi) 
  beta ~ N(mu_beta,Sigma_beta)
  gi| aG, rG ~ Beta(gi; shape1=aG,shape2=rG)
  aG, rG | G ~ G
  G ~ DP(v,G0)
  where G0(aG,rG) ~ Unif(aG; min=0.5, max=max_aG)Unif(rG; min=0.5, max=max_aG) and v ~ Unif(v;a_D,ib_D)
  
  The default values of the hyperparameters are:
  
  \code{mu_beta} = \code{rep(0,p)}
  \code{Sigma_beta} = \code{diag(5,p)}
  \code{max_aG} = \code{30}
  \code{a_D} = \code{0.01}
  \code{ib_D} = \code{3}
  
  These selections of hyperparameters could be used as uninformative ones.
  
}


\usage{
lmeNBBayes(formula,data,ID, B = 105000, burnin = 5000,  
           printFreq = B, M = NULL,probIndex = FALSE,
           thin =1,labelnp=NULL, epsilonM = 0.001,
           para = list(mu_beta = NULL,Sigma_beta = NULL,a_D = 0.01, ib_D = 3,max_aG=30),
           DP=TRUE,seed=1,Reduce=1)
}

\arguments{
  \item{formula}{
    An object of class "formula" (or one that can be coerced to that class): a symbolic description of the model to be fitted.
    %
    The formula must contain an intercept term.
  }
 \item{data}{
   A data frame, list or environment (or object coercible by \code{as.data.frame} to a data frame) containing the variables in the model.
   Each row must contain the data corresponding to the repeated measure j of a subject  and the rows \emph{(i,j)}s must be ordered as 
   \emph{(1,1),...,(1,n_1),(2,1),...,(2,n_2),...,(N,n_N)}. Missing values are not accepted.
 }
 \item{ID}{
   A vector of length \emph{Ntot}, containing the patient index for each
   repeated measure.
 }
 \item{B}{
   A scalar, the number of McMC iterations.
 }
 \item{burnin}{
   A scalar for a burn-in period.
   The proposal variance of the Metoropolice-Hasting rate is adjusted during the burn-in preiod.
 }
 \item{printFreq}{
   An integer value, indicating the frequency of iterations to print
   during the McMC run.
 }
 \item{M}{
   A truncation value for the infinite mixture of beta distribution. \code{M} determines the number of mass points used in the approximation of the DP mixture. If it is NULL, M is selected so that the amount of probability assigned to the final mass point is expected to be \code{epsilonM}. The specification of \code{M} is unnecessary if \code{DP}=0.
 }
 \item{probIndex}{
   Logical, if it is \code{TRUE} then the conditional probability index
  is computed for each patient at every thinning after discarding
  burn-in.
}
\item{thin}{
  Thinning frequency. Necessary if \code{probIndex} is \code{TRUE}.
}
\item{labelnp}{
  A vector of length \emph{Ntot}, containing 0 or 1. Zero indicates that
  the corresponding repeated measure should be treated as pre-scan and 1
  indicates that it is a new scan. \code{labelnp} is necessary only if
  \code{probIndex} is \code{TRUE}.
}
\item{epsilonM}{
  A scalar. See the description of \code{M}.
}
 \item{para}{
A list containing hyperparameter values.  
If \code{DP}=0 then the followings must be specified: \code{mu_beta} (a vector of length \emph{p}), \code{Sigma_beta} (a \emph{p} by \emph{p} covariance matrix) and \code{max_aG} (a positive scaler).
If \code{DP}=1 then in addition to the above parameters,\code{a_D} (positive scaler) \code{ib_D} (positive scaler) must be specified. 
If some of these are not specified then the default values discussed in description are used. 
}
\item{DP}{
  If \code{DP}=1 then the flexible mixed effect negative binomial regression is fit to the dataset. 
If \code{DP}=0 then the random effect distribution is assumed to be a single beta distribution.
}
\item{seed}{
  Seed to set.
}
 \item{Reduce}{
If \code{Reduce}=0 then the target distribution of McMC is full posterior distribution.
If \code{Reduce}=1 then the target distribution of McMC is a partially marginalized posterior distribution, which integrate out the random effects.
We recommend the use of the partially marginalized posterior distribution because it is computationally faster and the McMC samples are less correlated.
}
}


\references{
Kondo, Y. and Zhao, Y., A Flexible Mixed Effect Negative Binomial Regression for Detecting
Abrupt Increase in MRI Lesion Counts in Multiple Sclerosis Patients.
}
\author{
Kondo, Y.
}

\seealso{
%\code{\link{lmeNBBayes}}
\code{\link{getDIC}}
\code{\link{dqmix}}
\code{\link{index.batch.Bayes}}
}

\examples{

\dontrun{
## Generate samples as in manuscript at review 1
d <- getSample(rev=1,iseed=1,mod=0,dist="YZ")
burnin <- 5000	
thin <- 10
B <- 105000
## Create a vector indicating which sample to keep and discard
useSample <- useSamp(thin,burnin,B)
## nonparametric Bayesian procedure
re <- lmeNBBayes(formula=Y~1+X2,             
             data=d,         
             ID=d$ID,        
             B = B, 
	     printFreq=20000,
             burnin = burnin,
	     probIndex=0,
             labelnp=d$labelnp)
## compute DIC to assess the fit of nonparametric Bayesian procedure
getDIC(re,data=d,ID=d$ID,useSample=useSample)

## Parametric Bayesian procedure                
re.p <- lmeNBBayes(formula=Y~1+X2,             
             data=d,         
               ID=d$ID,        
               B = B, 
	       printFreq=20000,
               burnin = burnin,
	       probIndex=0,
               labelnp=d$labelnp,
               DP=0)
## compute DIC to assess the fit of parametric Bayesian procedure
getDIC(re.p,data=d,ID=d$ID,useSample=useSample)

alphasNP.inc <- 0.001
alphasNP <- seq(0.001,0.99,alphasNP.inc)
cex.axis <- 1.76
lwd <- 2.7
ymax <- 6
## the density of the sampled random effect distribution
qM <- dqmix(
            weightH1 = re$weightH1[useSample,],
            aGs = re$aGs[useSample,],
            rGs = re$rGs[useSample,],
            alphas=alphasNP,dens=TRUE
            )[[1]]

## plot of sampled random effect distribution
colM <- colMeans(qM)
q0025 <- apply(qM,2,quantile,prob=0.025)
q0975 <- apply(qM,2,quantile,prob=0.975)
hist(slim(d$gsBASE,d$ID),probability=TRUE,
     breaks=seq(0,1,0.05),col="gray",main="",xlab="beta",
     axes=FALSE,cex.lab=1.3,ylim=c(0,ymax))
axis(1,seq(0,1,0.2),seq(0,1,0.2),lwd=3,lwd.tick=3,cex.axis=cex.axis)
axis(2,seq(0,ymax,2),seq(0,ymax,2),lwd=3,lwd.tick=3,las=2,cex.axis=cex.axis)
## estimate from Bayesian parametric
densP <- matrix(NA,nrow=length(useSample),ncol=length(alphasNP))
for (ialpha in 1 : length(alphasNP))
   densP[,ialpha] <- dbeta(alphasNP[ialpha],shape1=re.p$aG[useSample],shape2=re.p$rG[useSample])
colM.p <- colMeans(densP)
points(alphasNP,colM.p,type="l",ylim=c(0,ymax),lwd=lwd,col="blue")
points(alphasNP,colM,type="l",xlab="random effect",ylab="density",ylim=c(0,ymax),lwd=lwd)
points(alphasNP,q0025,type="l",col="gray30",lty=2,lwd=lwd*0.5)
## True random effect distribution: for gamma and lognormal r.e., dist'n, the probability parameter of
## the negative binomial is parametrized as p = 1/(g+1).transformed
## the following plot plots the p.d.f. of transformed p (not g).
dmix <- getSample(mod=2,dist="YZ",ts=alphasNP) 
points(alphasNP,dmix,type="l",col="red",lwd=lwd) ## true

}

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
