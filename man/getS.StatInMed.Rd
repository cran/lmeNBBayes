\name{getS.StatInMed}
\alias{getS.StatInMed}


\title{
Generate samples from the flexible mixed-effect negative binomial distribution
}

\description{
This function yields samples from the simulation models specified in the paper by Kondo Y et al.
}







\usage{

getS.StatInMed(iseed = "random", rev = 4, dist = "b", 
               mod = 0, probs = seq(0, 0.99, 0.01), 
               ts = seq(0.001, 0.99, 0.001), trueCPI = FALSE, 
               full = FALSE, Scenario = "SPMS")
}


\arguments{
\item{iseed}{
Necessary only when \code{mod} = 0.
Integers to specify a seed. If iseed="random", seed is not specified.
}
\item{rev}{
 Necessary only when \code{mod} = 0.
At which DSMB reviews, data is generated.
}
\item{dist}{
 Necessary only when \code{mod} = 0.
 \code{dist} must be either \code{"b"} \code{"b2"} or \code{"YZ"}.
If \code{dist="b"} then random effect \eqn{G_{i}}{G[i]} is from a single beta, 

if \code{dist="b2"} then it is from a mixture of two betas and 

if \code{dist="YZ"} then it is transformed to range \eqn{[0,Inf)} and from a mixture of normal and gamma.

See \code{details} for more details.
}

 \item{mod}{
 If \code{mod} = 0 then \code{getS.StatInMed} generates a simulation sample.

 If \code{mod} = 1 then \code{getS.StatInMed} returns true quantiles of \eqn{G_{i}}{G[i]} at given \code{probs}.

 If \code{mod} = 2 then \code{getS.StatInMed} returns true densities of \eqn{G_{i}}{G[i]} at grids of points specified at \code{ts}.

 If \code{mod} = 3 then \code{getS.StatInMed} returns parameters of the simulation model.
 }
 \item{probs}{
 Necessary only when \code{mod} = 1. \code{probs} can be a vector of probabilities.
 }
 \item{ts}{
 Necessary only when \code{mod} = 3. \code{ts} can be a grid of points in (0,1).
 }



\item{trueCPI}{
 Necessary only when \code{mod} = 0.
 If \code{trueCPI=TRUE}, \code{getS.StatInMed} returns the true conditional probability indices of \eqn{N} patient, 
 computed using separate mean functions for the control and treated groups.
 }
 \item{full}{
 Necessary only when \code{mod} = 0.
 If \code{full=TRUE}, \code{rev} is ignored and \code{getS.StatInMed} returns complete dataset that contains 180 patients and all of them have 10 repeated measures.
}

\item{Scenario}{
Necessary only when \code{mod} = 0.
If \code{Scenario} specifies the prior for \eqn{\boldsymbol{\beta}}, and it must be either \code{"full"} or \code{"SPMS"}. 
See \code{details}.
} 


}











\value{
When \code{mod=0}, the it returns a dataframe that contains:

\code{Y}: the generated response counts

\code{Intercept}: all 1

\code{timeInt1}: 1 if the row corresponds to the count taken at the 1st 4-month followup interval, else 0.

\code{timeInt2}: 1 if the row corresponds to the count taken at the 2nd 4-month followup interval, else 0.

\code{ID}: the patient ID.

\code{gs}: generated random effect \eqn{G_i}.

\code{scan}: \eqn{-1,0,1,2,\cdots}. 
-1 indicates the screening scan, 0 indicates the baseline and 1,2,.. indicates 1,2-th followup scans.

\code{days}: \code{scan} + 2

\code{hs}: indicates which component the random effect is generated. 
If \code{dist="b"}, \code{hs} is all 1.
If \code{dist="b2"} or \code{dist="YZ"} then \code{hs} is either 1 or 2.

\code{trtAss}: indicates treatment assignments.

\code{labelnp}: 1 if scans correspond to new scans else 0.

\code{betPlcb}: The only first three elements are relevant. 
They are the generated intercept and time effects of placebo patients for mean counts.
i.e., \eqn{\alpha_0,\beta_{0,1},\beta_{0,2}}.

\code{betFull}: The only first five elements are relevant.  They are the generated intercept and time effects for mean counts.
i.e., \eqn{\alpha_0,\beta_{0,1},\beta_{1,1},\beta_{0,2},\beta_{1,2}}.

\code{probIndex}: The only first \eqn{N} elements are relevant. They are the true conditional probability indices of \eqn{N} patient, computed using true mean functions of the control patients (for both control and treated patients).

\code{probIndexTRUE}: This appear only if \code{trueCPI=TRUE}.
The only first \eqn{N} elements are relevant. They are the true conditional probability indices of \eqn{N} patient, computed using separate mean functions for the control and treated groups.
}






\details{

Simulation settings are as follows. 


Given the covariate vectors \eqn{\boldsymbol{X}_{ij}} for mean counts, response counts \eqn{Y_{ij}} of the \eqn{j}th repeated measure of \eqn{i}th patient are assumed to be from the mixed-effect negative binomial model: 

\eqn{
Y_{ij} | G_i=g_i, \boldsymbol{\beta} \; i.i.d.\sim 
NB(Y_{ij}; \textrm{size}=\exp(\boldsymbol{X}_{ij}^T\boldsymbol{\beta}),\textrm{prob}=g_i) 
}{
Y[ij] | G[i]=g[i], beta i.i.d.~ NB(Y[ij]; size=exp(X[ij]*beta),prob=g[i])    
}.

This formulation results in \eqn{\log E(Y_{i,j}) = \log(\mu_{\frac{1}{G}} - 1 ) + \boldsymbol{X}_{ij}^T\boldsymbol{\beta}}.
The mean count is modeled on the log scale as a constant over every four-month follow-up period, where the constants are allowed to dependent on the treatment assignment \eqn{A_i} (\eqn{A_i} = 1 for treatment, else 0). 

\eqn{
\log E(Y_{i,j}|A_i = a_i) = \log(\mu_{\frac{1}{G}} - 1 ) + \alpha_0 + 
\sum_{a=0}^{1} \sum_{t=1}^2 \beta_{a,t} I(j \in \boldsymbol{T}_t, a_i = a)
},

where \eqn{\mu_{\frac{1}{G}} = E(\frac{1}{G_i})} and \eqn{\boldsymbol{T}_1} and \eqn{\boldsymbol{T}_2} respectively contains indices corresponding to scans taken within the 1st and 2nd four-month interval during the follow-up.
 
The regression coefficients \eqn{\boldsymbol{\beta}=(\alpha_0,\beta_{0,1},\beta_{1,1},\beta_{0,2},\beta_{1,2})} 
are assumed to differ among studies and are generated from a multivariate normal distribution with 
\eqn{\boldsymbol{\mu}}{mu} and \eqn{\Sigma}{Sigma} replaced by the estimates from the full informative prior (Scenario A) or the SPMC informative prior (Scenario B) developed in Section 5.3 of the referred paper. To see their values set \code{mod=3}. (See example below). The proportion of treated patients is assumed to be 0.67.

\code{getS.StatInMed} allows three random effect model:

Setting 1: \eqn{G_i\sim Beta(3,0.8)}{G[i] ~ Beta(3,0.8)} which returns \eqn{(E(Y_{ij},SD(Y_{ij})))=(1.48,3.45)} and \eqn{(1.40,3.29)} at baseline under full and SPMC Scenarios respectively.

Setting 2: \eqn{G_i\sim 0.3 Beta(10,10) + 0.7 Beta(20,1)} which returns \eqn{(E(Y_{ij},SD(Y_{ij})))=(4.12,3.73)} and \eqn{(0.20,0.51)} at baseline for the patients whose REs are generated from the first and the second component of the mixture under full Scenario, and \eqn{(2.90,3.59)} and \eqn{(0.18,0.49)} under SPMS Scenario.

Setting 3: \eqn{G_i=\frac{1}{G_i^*+1}} where \eqn{G_i^*\sim 0.85 Gamma (0.176,2.226) + 0.15N(1.820,0.303)} where \eqn{Gamma(a,b)} represents the gamma CDF with variance \eqn{ab^2}. This mixture returns  \eqn{(E(Y_{ij},SD(Y_{ij})))=(1.46,4.16)} and \eqn{(6.75,4.56)} at baseline for the patients whose REs are drawn from the first and second component of the mixture under full Scenario, and \eqn{(1.38,3.97)} and \eqn{(6.39,4.41)} under SPMS Scenario. 

The choices of the random effect settings can be controlled by the input \code{dist}.

\code{getS.StatInMed} assumes MRI scans are taken monthly with a total of 10 MRI scans for each patients, a screening a baseline and eight follow-up scans. We assume that 15 patients are recruited every month so that 180 patients are recruited in 12 months, leading to a study duration of 21 month. DSMB reviews are assumed to occur every 4 months.
This means that by the DSMB review 1, 150 scans are available from 60 patients 
(15 patients have 4 scans each, 
next 15 patients have 3 scans each, 
next 15 patients have 2 scans each and the last 15 patient have 1 scan each.)
By the DSMB review 2, 540 scans are available from 120 patients
(15 patients have 8 scans each, the next 15 patients have 7 scans each,....the last 15 patients have only 1 scan each.)
The DSMB reviews can be specified by the input \code{rev}.
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
\code{\link{getS.StatInMed}}
\code{\link{getDIC}}
\code{\link{dqmix}}
\code{\link{index.batch.Bayes}}
}

\examples{

\dontrun{

## See the full informative prior for beta
temp <- getS.StatInMed(mod=3,Scenario="full")
temp$mu_beta
temp$Sigma_beta

## See the SPMS informative prior for beta
temp <- getS.StatInMed(mod=3,Scenario="SPMS")
temp$mu_beta
temp$Sigma_beta

## See also the examples in lmeNBBayes
}

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line

