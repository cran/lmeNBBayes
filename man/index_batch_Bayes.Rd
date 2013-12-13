\name{index.batch.Bayes}
\alias{index.batch.Bayes}
\alias{index.YZ}
\alias{index.b.each}
\alias{prIndex}
\alias{condProbCI}



\title{

The main function to compute the point estimates and 95\% credible intervals of the conditional probabilities 
\eqn{ 
Pr(Y_{i,new+} \ge y_{i,new+}| \boldsymbol{Y}_{i,pre}=\boldsymbol{y}_{i,pre})
}{
Pr(Y[i,new+]>= y[i,new+]| Y[i,pre]=y[i,pre])
} for multiple subjects.

}



\description{
Let \eqn{m_i}{m[i]} be the number of pre-measurements and \eqn{n_i}{n[i]} be the total number of repeated measures. 
Then the repeated measure of a subject can be divided into a pre-measurement set and a new measurement set as 
\eqn{
\boldsymbol{Y}_{i}=(\boldsymbol{Y}_{i,pre},\boldsymbol{Y}_{i,new})
}{
Y[i]=(Y[i,pre],Y[i,new])
}, where
\eqn{
\boldsymbol{Y}_{i,pre}=(Y_{i,1},\cdots,Y_{i,m_i})
}{
Y[i,pre]=(y[i,1],\cdots,Y[i,m[i]])
} 
and 
\eqn{
\boldsymbol{Y}_{i,new}=(Y_{i,m_i+1},\cdots,Y_{i,n_i})
}{
Y[i,new]=(Y[i,m[i]+1],...,Y[i,n[i]])
}
.
Given an output of \code{\link{lmeNBBayes}},
this function computes the probability of observing the response counts as large as those new observations of subject \eqn{i}{i},
\eqn{
\boldsymbol{y}_{i,new}
}{
y[i,new]
} 
conditional on the subject's previous observations 
\eqn{
\boldsymbol{y}_{i,pre}
}{
y[i,pre]
} for subject \eqn{i}{i}.
That is, this function returns a point estimate and its asymptotic 95\% confidence interval (for a parametric model) of the conditional probability for each subject:
%

\eqn{ 
Pr(Y_{i,new+} \ge y_{i,new+}  | \boldsymbol{Y}_{i,pre}=\boldsymbol{y}_{i,pre})
}{ 
Pr( Y[i,new+] \ge y[i,new] | Y[i,pre]=y[i,pre]) 
},
where 
\eqn{
Y_{i,new+}=\sum_{j=m_i+1}^{n_i} Y_{ij} 
}{ 
Y[i,new+]=\sum[j=m[i]+1]^{n[i]} Y[ij]  
}.

}


%\description{
%Compute the posterior interval estimate of the conditional probability index for each patient.
%}


\usage{

index.batch.Bayes(data,labelnp,ID,olmeNBB,thin=NULL,printFreq=10^5,unExpIncrease=TRUE)

}

\arguments{
 \item{data}{
  See \code{\link{lmeNBBayes}}.
  This \code{data} does not have to be the same as the one used in the computations of negative binomial mixed effect regression (\code{\link{lmeNBBayes}}).
 }

\item{labelnp}{
See \code{\link{lmeNBBayes}}.
\code{nrow(data)} == \code{length(labelnp)} must be satisfied.
}

\item{ID}{
See the description in \code{\link{lmeNBBayes}}.
\code{nrow(data)} == \code{length(labelnp)} must be satisfied.
}
\item{olmeNBB}{The output of the function \code{\link{lmeNBBayes}}.}

\item{thin}{The frequency of thinning}

\item{printFreq}{See the description in \code{\link{lmeNBBayes}}.}

\item{unExpIncrease}{Internal use only. Should be always TRUE.}
}



\value{
\item{condProb}{
\code{(olmeNBB$para$B-olmeNBB$para$burnin)/thin} by the number of patients, \eqn{N} (=\code{length(unique(ID))}), matrix, containing the MCMC samples of the conditional probability index for each patient at every selected iteration (after discarding burn-in and thinning).
If some patients have 0 pre-scans or 0 new-scans, then \code{NA} is returned.
}
\item{condProbSummary}{
\code{4} by \eqn{N} matrix.
The first column contains the posterior estimates of the conditional probability index.
The second column contains the posterior SE.
The third column contains the lower bound of the 95\% credible interval.
The fourth column contains the upper bound of the 95\% credible interval. 
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
\code{\link{dqmix}}

}

\examples{
## See the examples of function lmeNBBayes

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
