\name{metaMultGibb}
\alias{metaMultGibb}
\alias{metaUnivGibb}
\alias{listMean}

\title{
Fit the random effect meta-analysis model.
}


\description{
The random effect meta-analysis model has the following hierarchy:
hat.beta_s| beta_s ~ N_p(beta_s,Phi_s), s = 1,...,\emph{S}
beta_s|mu_beta,Sigma_beta ~ N_p(mu_beta,Sigma_beta)
mu_beta ~ Unif(-Inf,Inf)  
Sigma_beta ~ inverse-Wishart(v0,S0)
}


\usage{
metaMultGibb(hat.betas, sigmas2, v0 = "uninfo", S0 = NULL, B = 10000) 
}

\arguments{
  \item{hat.betas}{A vector of length \emph{S}, containing the samples}
  
  \item{sigmas2}{A list of length \emph{S}, each containing the \emph{p}
    by \emph{p} sample variance covariance matrix.}

  \item{v0}{A scalar parameter for the inverse-Wishart distribution. 
    If \code{v0}="info" then \code{v0}=10 and \code{S0}=\code{diag(p) * (v0 - p - 1)} is specified.
    If \code{v0}="uninfo" then \code{v0} =\code{p + 2} and \code{S0} = \code{diag(c(0.5, 0.1))} is specified.
    These selections of \code{v0} is suggested by Hoff[2009]. 
  }

  \item{S0}{A positive definite matrix.}

  \item{B}{The McMC sample size.}
}


\references{
P.D. Hoff. A First Course in Bayesian Inference in Statist. Springer
Texts in Statistics, 2009.

Kondo, Y. and Zhao, Y., A Flexible Mixed Effect Negative Binomial Regression for Detecting
Abrupt Increase in MRI Lesion Counts in Multiple Sclerosis Patients.
}
\author{
Kondo, Y.
}

\seealso{
\code{\link{lmeNBBayes}}
}

\examples{
## See the examples of function nbinDP

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
