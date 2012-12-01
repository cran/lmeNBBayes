\name{dqmix}
\alias{dqmix}


\title{
Estimate the random effect distribution of the flexible mixed effect negative binomial regression.
}


\description{
Given the output of the \code{\link{lmeNBBayes}}, this function estimate the posterior density of the random effect at grids.
}


\usage{
dqmix(weightH1, aGs, rGs, alphas = seq(0, 1, 0.01), dens = TRUE) 
}

\arguments{

\item{weightH1}{
A \emph{B} by \emph{M} matrix, containing the probability components of the infinite mixture of beta distribution.
The output of \code{\link{lmeNBBayes}}.
}

\item{aGs}{
A \emph{B} by \emph{M} matrix, containing the shape1 parameters of the infinite mixture of beta distribution.
The output of \code{\link{lmeNBBayes}}.
}

\item{rGs}{
A \emph{B} by \emph{M} matrix, containing the shape2 parameters of the infinite mixture of beta distribution.
The output of \code{\link{lmeNBBayes}}.
}

\item{alphas}{
The grids of points at which the density is evaluated. Must be the points within 0 and 1.
}

\item{dens}{
If TRUE, the density estimate of the random effect is return, if FALSE then the quantiles of the random effects are estimated.
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
\code{\link{lmeNBBayes}}
\code{\link{getDIC}}
%\code{\link{dqmix}}
}

\examples{
## See the examples of function nbinDP

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
