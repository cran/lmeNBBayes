\name{index.batch.Bayes}
\alias{index.batch.Bayes}
\alias{index.YZ}
\alias{index.b.each}
\alias{condProbCI}
\alias{prIndex}

\title{
Compute the conditional probability index.
}


\description{
Compute the posterior interval estimate of the conditional probability index for each patient.
}


\usage{

index.batch.Bayes(data,labelnp,ID,olmeNBBayes,thin=NULL,printFreq=10^5,unExpIncrease=TRUE)

}

\arguments{
 \item{data}{
   A data frame, list or environment (or object coercible by \code{as.data.frame} to a data frame) containing the variables in the model.
   Each row must contain the data corresponding to the repeated measure j of a subject  and the rows \emph{(i,j)}s must be ordered as 
   \emph{(1,1),...,(1,n_1),(2,1),...,(2,n_2),...,(N,n_N)}. Missing values are not accepted.
 }

\item{labelnp}{See the description in \code{\link{lmeNBBayes}}.
The \code{labelnp} used in \code{index.batch.Bayes} does not have to be the same as the \code{labelnp} used in \code{nbinDPmix}.
}
\item{ID}{See the description in \code{\link{lmeNBBayes}}.
The \code{ID} used in \code{index.batch.Bayes} does not have to be the same as the \code{ID} used in \code{nbinDPmix}.
}
\item{olmeNBBayes}{The output of the function \code{\link{lmeNBBayes}}.}

\item{thin}{The frequency of thinning}

\item{printFreq}{See the description in \code{\link{lmeNBBayes}}.}

\item{unExpIncrease}{Internal use only. Should be always TRUE.}
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
\code{\link{dqmix}}

}

\examples{
## See the examples of function nbinDP

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
