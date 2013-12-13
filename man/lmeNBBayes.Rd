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
\alias{repeatAsID}
\alias{piM}
\alias{lnpara}
\alias{int.M}
\alias{E.KN}



\title{
Generate posterior samples from a flexible mixed effect negative binomial regression model.
}

\description{

Let \eqn{Y_{ij}} be the response count at \eqn{j}th repeated measure from the \eqn{i}th patient (\eqn{i=1,\cdots,N} and \eqn{j=1,\cdots,n_i}).
The negative binomial mixed-effect independent model assumes that given the random effect \eqn{ G_i=g_i}{G[i]=g[i]}, 
the count response from the same subjects i.e., \eqn{Y_{ij}} and \eqn{Y_{ij'}} are conditionally independent and follow the negative binomial distribution:

\eqn{
Y_{ij} | G_i=g_i, \boldsymbol{\beta} \; i.i.d.\sim 
NB(Y_{ij}; \textrm{size}=\exp(\boldsymbol{X}_{ij}^T\boldsymbol{\beta}),\textrm{prob}=g_i) 
}{
Y[ij] | G[i]=g[i], beta i.i.d.~ NB(Y[ij]; size=exp(X[ij]*beta),prob=g[i])    
}

where \eqn{\boldsymbol{X}_{ij}} is the covariates for mean counts. This formulation results in \eqn{\log E(Y_{i,j}) = \log(\mu_{\frac{1}{G}} - 1 ) + \boldsymbol{X}_{ij}^T\boldsymbol{\beta}}.
To allow flexible form of a random effect distribution, we assume that the patient-specific random effect is assumed to be from Dirichlet process mixture of beta distributions. This essentially means that random effect \eqn{G_i}{G[i]} is from an infinite mixture of Beta distributions:

\eqn{
G_i| \{a_{G_h} r_{G_h}, \pi_h \}_{h=1}^{\infty} 
\sim \sum_{h=1}^{\infty} \pi_h Beta(G_i; \textrm{shape1}=a_{G_h},\textrm{shape2}=r_{G_h})
}{
G[i]| { a[G[h]], r[G[h]], pi[h] }[h=1]^{infty} ~ sum pi[h] Beta(G[i]; shape1=a[G[h]],shape2=r[G[h]])	  
},

where \eqn{\pi_h}{pi[h]} is modelled with the stick-breaking prior. Introducing latent variable 
\eqn{V_h,h=1,2,\cdots}{V[h],h=1,2,...}, this prior is defined as \eqn{\pi_1=V_1}{pi[1]=V[1]} and
\eqn{
 pi_h=V_h \prod_{l < h} (1-V_h) 
}{
 pi[h]=V[h]prod[l < h] (1-V[h])
}
for \eqn{h>1} \eqn{V_h \; i.i.d.\sim Beta(1,D)}{V[h] i.i.d. ~ Beta(1,D)}.


The rest of priors are specified as: 
\eqn{
\boldsymbol{\beta} \sim N(\boldsymbol{\mu},\Sigma)
}{
beta ~ N(mu,Sigma)
},


\eqn{
(a_G,r_G) \sim Unif(a_G; \textrm{min}=0.5, \textrm{max}=max_{a_G})Unif(r_G; \textrm{min}=0.5, \textrm{max}=max_{r_G}) 
}{
(a[G],r[G]) ~ Unif(a[G];\textrm{min}=0.5,\textrm{max}=max[a[G]])Unif(r[G];\textrm{min}=0.5, \textrm{max}=max[r[G]])  
},

\eqn{
D \sim Unif(v;\textrm{min}=a_D,\textrm{max}=ib_D)
}{
D ~ Unif(v;min=a_D,max=ib_D)
}. 

The default values of the hyperparameters are \eqn{\mu_{\boldsymbol{\beta}}}{mu[beta]} = \code{rep(0,p)}, \eqn{\Sigma_{\boldsymbol{\beta} }}{Sigma[beta]} = \code{diag(5,p)}, \eqn{\textrm{max}_{a_G}}{max[a[G]]} = \code{30}, \eqn{a_D}{a[D]} = \code{0.01} and \eqn{ib_D}{ib[D]} = \code{3}. These selections of hyperparameters could be used as uninformative ones.
  
The function \code{lmeNBBayes} also allows generating posterior samples from the parametric version of the model above which simply assumes that the random effect is from the single beta distribution. (The rest of the prior specifications are the same).

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
The formula must contain an intercept term.
  }
 \item{data}{
A data frame, list or environment (or object coercible by as.data.frame to a data frame) containing the variables in the model.
The each row must contains the data corresponding to the repeated measure \eqn{j}{j} of subjects and the rows \eqn{(i,j)}{(i,j)}s must be ordered in a way that measurements from a subject is clustered together as \eqn{(1,1)}{(1,1)},...,\eqn{(1,n_1)}{(1,n[1])},\eqn{(2,1)}{(2,1)},...,\eqn{(2,n_2)}{(2,n[2])},...,\eqn{(N,n_N)}{(N,n[N])}.
 }
 \item{ID}{
A vector of length \eqn{ \sum_{i=1}^N n_i }{ \sum[i=1]^N n[i] }, containing the patient IDs that corresponds to \code{data}. i.e., \code{c(rep(ID_1,n_1),}\code{rep(ID_2,n_2)}\code{,...,rep(ID_N,n_N))}. The length must be the same as the number of rows of \code{data}. Missing ID values are NOT accepted.
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
Necessary only if \code{DP=1}.
Our Gibbs sampler approximates the infinite mixture of beta distributions by truncating it with \code{M} components by setting \eqn{V_M=1}{V[M]=1} so that \eqn{pi_M = 1 - \sum_{h=1}^{M-1} \pi_h}. If \code{M} is NULL, \code{M} is selected so that the amount of probability assigned to the final mass point is expected to be \code{epsilonM}. i.e., \eqn{E(\pi_M)=E\{ (\pi_M | D )\}  = E (\{ 1- 1/(D+1)\}^{M-1}) < \epsilon}.
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
  A vector of length \eqn{ \sum_{i=1}^N n_i }{ \sum[i=1]^N n[i] }, containing 0 or 1. Zero indicates that
  the corresponding repeated measure should be treated as pre-scan and 1
  indicates that it is a new scan. \code{labelnp} is necessary only if
  \code{probIndex} is \code{TRUE}.
}
\item{epsilonM}{
  A scalar. See the description of \code{M}.
}
 \item{para}{
A list containing hyperparameter values.  
If \code{DP}=0 then the followings must be specified: \code{mu_beta} (a vector of length \eqn{p}{p}), \code{Sigma_beta} 
(a \eqn{p}{p} by \eqn{p}{p}  covariance matrix) and \code{max_aG} (a positive scaler).
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
For the parameters with non-conjugate priors \eqn{\boldsymbol{\beta}, D, a_G, b_G}{beta,D,a[G],b[G]}, 
the Metropolis Hasting (MH) algorithm is employed to sample from their full conditional distributions.
For \eqn{D, a_G, b_G}{D,a[G],b[G]}, the MH algorithm is
performed separately with a normal proposal distribution, where its proposal variance
is tuned during the burn-in to have the acceptance rates range between 0.2 and 0.6.
One can adjust the minimum and maximum of the proposal sd via the \code{proposalSD} arguments.
For each element of \eqn{\boldsymbol{\beta}}{beta}, we found that updating
each regression coefficient with separate MH algorithm resulted in
poor mixing in the Markov chain when high correlation is assumed in some
of \eqn{\boldsymbol{\beta}}{beta} in the prior.
Therefore, the MH algorithm is performed simultaneously for all \eqn{\boldsymbol{\beta}}{beta} and 
a MVN proposal distribution is employed with
\eqn{a}{a}\eqn{\Sigma}{Sigma} as its proposal covariance matrix, where \eqn{\Sigma}{Sigma} is the
covariance of a prior for \eqn{\boldsymbol{\beta}}{beta} and \eqn{a}{a} is a tuning scaler
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

## generate samples from DSMSB review 2
d <- getS.StatInMed(rev=2,iseed=1,dist="YZ",Scenario="full")
formula.fit <- Y ~ timeInt1:trtAss + timeInt2:trtAss

B <- 10000
burnin <- 1000
thin <- 2
fit <- lmeNBBayes(formula=formula.fit,data=d, ID=d$ID, 
                  B = B, burnin = burnin,  thin=thin)
## The output can be printed out:
fit 


## Now, compute the conditional probability index using the mean function of placebo patients.
## We need to modify two things in output of lmeNBBayes.
## 1st, change the formula so that it does not distinguish between treatment and placebo
fit$para$formula <- Y ~ timeInt1 + timeInt2
## 2nd, disregard the coefficient that corresponds to the treated patients
fit$beta <- fit$beta[,-c(3,5)]
cpi <- index.batch.Bayes(data=d,labelnp=d$labelnp,ID=d$ID,
                         olmeNBB=fit,printFreq=10^7)
cpi 

## finally access the accuracy of the CPI estimates in terms of RMSE
Npat <- length(unique(d$ID))
est <- cpi$condProbSummary[,1]
true <- d$probIndex[1:Npat]
sqrt( mean( ( est - true )^2 ,na.rm=TRUE) )
              

}



}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line

