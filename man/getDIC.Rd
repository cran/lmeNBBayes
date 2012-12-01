\name{getDIC}
\alias{getDIC}

\title{
  Compute the DIC given the output from \code{\link{lmeNBBayes}}.
}


\description{
  If partially marginalized posterior distribution (i.e. \code{Reduce}=1 in the computation of \code{\link{lmeNBBayes}}) is a target distribution, the DIC is computed using the focused likelihood P(Y| aGs, rGs, h1s, beta). 
  If not, then the DIC is computed using the focused likelihood
  P(Y|g,beta).
}


\usage{
getDIC(olmeNBBayes, data, ID, useSample)
}

\arguments{
  \item{olmeNBBayes}{
    The output of \code{nbinDP}.
  }
  \item{data}{
    A data frame, list or environment (or object coercible by \code{as.data.frame} to a data frame) containing the variables in the model.
    Each row must contain the data corresponding to the repeated measure j of a subject  and the rows \emph{(i,j)}s must be ordered as 
    \emph{(1,1),...,(1,n_1),(2,1),...,(2,n_2),...,(N,n_N)}. Missing values are not accepted.
  }
  \item{ID}{
    See the description in \code{nbinDP}.
  }
  
\item{useSample}{
  A vector of length the total number of repeated measures (i.e. the
  same as the length of \code{Y}), containing TRUE or FALSE, indicating
  which samples to be kept or discarded.
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
%\code{\link{getDIC}}
\code{\link{dqmix}}
\code{\link{index.batch.Bayes}}
}

\examples{
## See the examples of function nbinDP

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
