\name{getDIC}
\alias{getDIC}
\alias{llk.FG_i}
\title{
  Compute the DIC given the output from \code{\link{lmeNBBayes}}.
}


\description{
  If partially marginalized posterior distribution (i.e. \code{Reduce}=1 in the computation of \code{\link{lmeNBBayes}}) is a target distribution, the DIC is computed using the focused likelihood 

\eqn{
Pr(\boldsymbol{Y}_{i}=\boldsymbol{y}_{i}| \{ a_{G_h}, r_{G_h} \}_{h=1}^{M}, \pi_h \}_{h=1}^{\infty}, \boldsymbol{\beta})
}{
Pr(Y[i]=y[i]| { a[G[h]] }[h=1]^{M}, { r[G_h] }[h=1]^{M}, { H[h] }[h=1]^{M}, beta)
}. 

If not, then the DIC is computed using the focused likelihood
\eqn{  
Pr(\boldsymbol{Y}_{i}| a_{G_{h_i}},r_{G_{h_i}}, \boldsymbol{\beta})
}{
Pr(Y[i] |a[G[h[i]]],r[G[h[i]]], ,beta)
}.
}


\usage{
getDIC(olmeNBB, data,
       ID, useSample=NULL,focus = c("FG","G","aGh.rGh","para"), 
       lower.alpha=0.0001,upper.alpha=0.99999,inc.alpha=0.0005)
llk.FG_i(ys,rs,aGs,bGs,ps)
}

\arguments{
  \item{olmeNBB}{
    The output of \code{lmeNBBayes}.
  }
  \item{data}{
    See \code{\link{lmeNBBayes}}.
  }
  \item{ID}{
    See the description in \code{\link{lmeNBBayes}}.
  }
  
\item{useSample}{
  A vector of length the total number of repeated measures (i.e. the
  same as the length of \code{Y}), containing TRUE or FALSE, indicating
  which samples to be kept or discarded.
}
\item{focus}{
Focused likelihood
}
\item{lower.alpha}{
Used only when \code{focus} is \code{FG}. See details.
}
\item{upper.alpha}{
Used only when \code{focus} is \code{FG}. See details.
}
\item{inc.alpha}{
Used only when \code{focus} is \code{FG}. See details.
}

\item{ys}{
A vector containing the response values
}
\item{rs}{
A vector containing the size parameters of negative binomial.
The length of rs must be the same as the length of ys.
}
\item{aGs}{
A vector containing the shape1 parameters of the approximated infinite mixture of betas.
}
\item{bGs}{
A vector containing the shape2 parameters of the approximated infinite mixture of betas.
The length must be the same as the length of \code{bGs}.
}
\item{ps}{
A vector containing the probability parameters of the approximated infinite mixture of betas.
The length must be the same as the length of \code{aGs}.
}

}

\details{

Denote \eqn{P}{P} be a vector of "focused" parameters.

Using Spiegelhalter et. al.(2002)'s notation, the effective number of parmeters can be computed as:

\eqn{p_D=\bar{D} - D(\bar{ P  })}{p[D]=D.bar - D( P.bar)}

where \eqn{D}{D} is the deviance and the \eqn{ \bar{ P } }{P.bar} is the expectation of \eqn{P}{P}.


When \code{focus} = \code{FG} then the focused parameters, denoted as \eqn{P}, are the random effect distribution
(i.e.,infinite mixture of beta distribution) and the regression coefficients. In the computation, the expected regression coefficients are obtained by simply computing the mean of the posterior samples of coefficients. The expected infinite mixture of beta distribution is obtained in the following steps:

STEP 1: Provide a fine grids of points between [0,1]. We chose the grid of points to be 

\code{alphas <- seq(lower.alpha,upper.alpha,inc.alpha)}. 

STEP 2: For each sampled infinite mixture of betas, 
Evaluate its value at every grid provided from STEP 1 for each sample.
Obtain \code{B} by \code{length(alphas)} matrix.

STEP 3: Given the matrix from STEP2, at each grid of points, we compute the average value of density.
Obtain a vector of length \code{length(alpha)}, that contains the estimated expected random effect density at fine grid of points.

STEP 4: Given the estimated expected coefficients and the estimated expected random effect density, 
evaluate \eqn{D(\bar{P})}{D(P.bar)} by integrating the conditional likelihood given random effects with respect to the estimated expected random effect density from STEP 3.
}
\references{
Kondo, Y. and Zhao, Y., A Flexible Mixed Effect Negative Binomial Regression for Detecting
Abrupt Increase in MRI Lesion Counts in Multiple Sclerosis Patients.

Spiegelhalter, David J.; Best, Nicola G.; Carlin, Bradley P.; van der Linde, Angelika (2002). "Bayesian measures of model complexity and fit (with discussion)". Journal of the Royal Statistical Society, Series B 64 (4): 583-639.
}
\author{
Kondo, Y.
}

\seealso{
\code{\link{lmeNBBayes}}
%\code{\link{getDIC}}
\code{\link{dqmix}}
\code{\link{index.batch.Bayes}}
}

\examples{
## See the examples of function lmeNBBayes

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
