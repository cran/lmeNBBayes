\name{lmeNBBayes}
\alias{lmeNBBayes}
\alias{pointsgamma}
\alias{plotgamma}
\alias{plotbeta}
\alias{plotnbinom}
\alias{slim}
\alias{plotGs}
\alias{newCat}
\alias{useSamp}
\alias{adjustPosDef}
\alias{colmeansd}
\alias{getM}
\alias{Nuniq}
\alias{getSNormUnblind}
\alias{repeatAsID}
\alias{piM}
\alias{lnpara}
\alias{int.M}
\alias{E.KN}



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

lmeNBBayes(formula,data, ID, B = 105000, burnin = 5000,  
           printFreq = B, M = NULL, probIndex = FALSE,
           thin =1,labelnp=NULL, epsilonM = 1e-4, 
	   para = list(mu_beta = NULL,Sigma_beta = NULL,
	   	  max_aG=30,mu_lnD=NULL,sd_lnD=NULL),
           DP=TRUE,thinned.sample=FALSE, proposalSD = NULL)
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
If \code{DP}=1 then in addition to the above parameters,\code{mu_lnD} (positive scaler) \code{sd_lnD} (positive scaler) must be specified. 
If some of these are not specified then the default values discussed in description are used. 
}
\item{DP}{
  If \code{DP}=1 then the flexible mixed effect negative binomial regression is fit to the dataset. 
If \code{DP}=0 then the random effect distribution is assumed to be a single beta distribution.
}

\item{thinned.sample}{
Logical. 
If true then return only the thinned samples, else returns the entire MCMC sample of size B.
}
\item{proposalSD}{
List object containing two list objects \code{min} and \code{max}, which contain minimum and maximum values of the proposal standard deviations.

If \code{DP}=0 then a list object \code{min} (\code{max}) must contains 3 elements corresponding to minimum (maximum) values of the proposal standard deviation of aG, rG and beta. See details for beta.

If \code{DP}=1 then a list object \code{min} (\code{max}) must contains 4 elements corresponding to minimum (maximum) values of the proposal standard deviation of aG, rG, beta and ln D. See details for beta.

        
}

}

\details{
For the parameters with non-conjugate priors beta, D, aG and bG, the
Metropolis Hasting (MH) algorithm is employed to sample from their full conditional distributions.
For D, aG, bG, MH algorithm is
performed separately with a normal proposal distribution, where its proposal variance
is tuned during the burn-in to have the acceptance rates range between 0.2 and 0.6.
One can adjust the minimum and maximum of the proposal sd via the proposalSD arguments.
For beta, we found that updating
each regression coefficient with separate MH algorithm resulted in
poor mixing in the Markov chain when high correlation is assumed in some
of beta in our prior.
Therefore, MH algorithm is performed simultaneously for beta and a MVN proposal distribution is employed with
\emph{a}\emph{Sigma} as its proposal covariance matrix, where \emph{Sigma} is the
covariance of a prior for \emph{beta} and \emph{a} is a tuning scaler
adjusted during the burn-in period to have the acceptance rates range between 0.2 and 0.6.
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
## coming soon
}

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
