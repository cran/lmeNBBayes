## This script contains the successful example of the code: 
## Data structure:
## no new scan, ni=10;sp1=0.5;sp2=3;Ntot=200 set.seed(1)-set.seed(10)
logL <- function(Y,ID,
                 X,vec.beta,vec.gPre
                 )
  {
    ## gPre, gNre is a vector of length N
    ## if (is.null(vec.gNew)) labelnp <- rep(0,length(Y))

    if (is.vector(X))
      X <- matrix(X,ncol=1)
    logL <- 0
    uniID <- unique(ID)
    for (ipat in 1 : length(uniID) )
      {
        patID <- uniID[ipat]
        iY <- Y[ID==patID]
        iX <- X[ID==patID,,drop=FALSE]
        i.gPre <- vec.gPre[ipat]
        logL <- logL + sum(dnbinom(iY,
                                   size=exp(iX%*%vec.beta),
                                   prob=i.gPre,
                                   log=TRUE)
                           )
      }
    return (logL)
  }

logL2 <- function(Y,ID,
                  X,vec.beta,
                  vec.aGs,vec.rGs ## LENGTH OF N
                  )
  {
    ## gPre, gNre is a vector of length N
    ## if (is.null(vec.gNew)) labelnp <- rep(0,length(Y))
    
    if (is.vector(X))
      X <- matrix(X,ncol=1)
    logL <- 0
    uniID <- unique(ID)
    for (ipat in 1 : length(uniID) )
      {
        patID <- uniID[ipat]
        iY <- Y[ID==patID]
        iX <- X[ID==patID,,drop=FALSE]
        aGh <- vec.aGs[ipat]
        rGh <- vec.rGs[ipat]
        logL <- logL + sum(log(choose(iY+exp(iX%*%vec.beta)-1,iY)))+
          log(beta(sum(exp(iX%*%vec.beta))+aGh,sum(iY)+rGh))-
            log(beta(aGh,rGh))
      }
    return (logL)
  }


index.batch.Bayes <- function(data,labelnp,ID,olmeNBBayes,thin=NULL,printFreq=10^5,unExpIncrease=TRUE)
  {
    burnin <- olmeNBBayes$para$burnin
    if (is.null(thin))
      {
        if (is.null(olmeNBBayes$para$thin))  thin <- 1
        else  thin <- olmeNBBayes$thin
      }
    formula <- olmeNBBayes$para$formula
    X <- model.matrix(object=formula,data=data)
    Y <- model.response(model.frame(formula=formula,data=data))

    if (is.vector(X)) X <- matrix(X,ncol=1)
    NtotAll <- length(Y)
    if (nrow(X)!= NtotAll) stop ("nrow(X) ! = length(Y)")
    if (length(ID)!= NtotAll)  stop ("length(ID)! = length(Y)")
    if (length(labelnp)!= NtotAll)  stop ("labelnp! = length(Y)")
  
    if (olmeNBBayes$para$DP==0)
      {
        olmeNBBayes$aGs <- matrix(olmeNBBayes$aG,ncol=1)
        olmeNBBayes$rGs <- matrix(olmeNBBayes$rG,ncol=1)
        olmeNBBayes$weightH1 <- matrix(1,ncol=1,nrow=length(olmeNBBayes$aGs))
      }
    
    useSample <- useSamp(thin=thin, burnin=burnin, B=nrow(olmeNBBayes$beta))
    maxni <- max(tapply(rep(1,length(ID)),ID,sum))

    nn <- length(unique(ID))
        ## change the index of ID to numeric from 1 to # patients
    temID <- ID  
    N <- length(unique(temID))
    uniID <- unique(temID)
    ID <- rep(NA,length(temID))
    for (i in 1 : length(uniID))
      {
        ID[temID == uniID[i]] <- i
      }
    ## ID starts from zero
    mID <- ID-1

    
    mat_betas <- olmeNBBayes$beta[useSample,,drop=FALSE]
    mat_aGs <- olmeNBBayes$aGs[useSample,,drop=FALSE]
    mat_rGs <- olmeNBBayes$rGs[useSample,,drop=FALSE]
    mat_weightH1 <- olmeNBBayes$weightH1[useSample,,drop=FALSE]
    B <- nrow(mat_betas) 
    M <- ncol(olmeNBBayes$aGs)
    
    if (unExpIncrease){
      
      re <- .Call("index_b_Bayes_UnexpIncrease",
                  as.numeric(Y), as.numeric(c(X)),
                  as.numeric(labelnp), as.integer(mID),
                  as.integer(B),as.integer(maxni),as.integer(M),
                  as.integer(nn),
                  as.numeric(c(t(mat_betas))), 
                  as.numeric(c(t(mat_aGs))),
                  as.numeric(c(t(mat_rGs))),
                  as.numeric(c(t(mat_weightH1))),
                  as.integer(printFreq),
                  package ="lmeNBBayes")
    }else{
      re <- .Call("index_b_Bayes_UnexpDecrease",
                  as.numeric(Y), as.numeric(c(X)),
                  as.numeric(labelnp), as.integer(mID),
                  as.integer(B),as.integer(maxni),as.integer(M),
                  as.integer(nn),
                  as.numeric(c(t(mat_betas))), 
                  as.numeric(c(t(mat_aGs))),
                  as.numeric(c(t(mat_rGs))),
                  as.numeric(c(t(mat_weightH1))),
                  as.integer(printFreq),
                  package ="lmeNBBayes")
    }
    
    re <- matrix(re[[1]],nrow=B,ncol=nn,byrow=TRUE)
    ## the labelnp of patients only with 1 (new scans) labels are replaced by all 0 (old scans)
    patwonew <- which(as.numeric(tapply((labelnp==0),ID,sum)==0)==1)
    for (i in 1 : length(patwonew)) labelnp[ID == patwonew[i]] <- 0
        
    patwoNorO <-  which(as.numeric(tapply((labelnp==1),ID,sum)==0)==1)
    if (length(patwoNorO)==0) patwoNorO <- NULL;

    re[,patwoNorO] <- NA

    res <- list()
    res$condProb <- re
    res$condProbSummary <- condProbCI(ID=ID,condProb=re)
    
    res$para$labelnp <- labelnp
    res$para$ID <-ID
    res$para$CEL <- Y
    res$para$thin <- thin
    res$para$B <- B
    res$para$burnin <- burnin
    class(res) <- "IndexBatch"
    return(res)
  }


condProbCI <- function(ID,condProb)
  {
    uniID <- unique(ID)
    
    condProb <- cbind(
                      colMeans(condProb),
                      apply(condProb,2,sd),
                      t(apply(condProb,2,
                              quantile,na.rm=TRUE,
                              prob=c(0.025,0.975)))
                      )
    colnames(condProb) <- c("CondProb","SE","2.5%","97.5%")
    rownames(condProb) <- uniID
    return(condProb)
  }


getDIC <- function(olmeNBBayes,data,ID,useSample)
{
  formula <- olmeNBBayes$para$formula
  X <- model.matrix(object=formula,data=data)
  Y <- model.response(model.frame(formula=formula,data=data))
  ## The effective number of parameters of the model is computed
  ## Dbar: posterior mean of the deviance bar(-2*logLlik)
  D.bar <- -2*mean(olmeNBBayes$logL[useSample]) ## + constant
  ## Dhat: -2*logLik(theta.bar)
  if (is.vector(olmeNBBayes$beta)) olmeNBBayes$beta <- matrix(olmeNBBayes$beta,ncol=1)
  ##if (is.vector(olmeNBBayes$gNew)) olmeNBBayes$gNew <- matrix(olmeNBBayes$gNew,ncol=1)
  if (olmeNBBayes$para$Reduce){
    if (olmeNBBayes$para$DP==0)
      {
        ## parametric procedure:
        olmeNBBayes$aGs_pat <- matrix(olmeNBBayes$aG,nrow=length(olmeNBBayes$aG),ncol=length(unique(ID)))
        olmeNBBayes$rGs_pat <- matrix(olmeNBBayes$rG,nrow=length(olmeNBBayes$rG),ncol=length(unique(ID)))
      }
    
    D.hat <- -2*logL2(Y=Y,ID=ID,X=X,
                      vec.beta = colMeans(olmeNBBayes$beta[useSample,,drop=FALSE]), 
                      vec.aGs = colMeans(olmeNBBayes$aGs_pat[useSample,,drop=FALSE]),
                      vec.rGs = colMeans(olmeNBBayes$rGs_pat[useSample,,drop=FALSE])
                     )
  }else{
    D.hat <- -2*logL(Y=Y,ID=ID,X=X,
                     vec.beta = colMeans(olmeNBBayes$beta[useSample,,drop=FALSE]),
                     vec.gPre = colMeans(olmeNBBayes$g1s[useSample,,drop=FALSE])
                     )
  }
  effect.para <- D.bar - D.hat
  DIC <- effect.para + D.bar
  return (list(DIC=DIC,effect.para=effect.para))
}



newCat <- function(label1,label2)
  {
    ## The length of label1 and label2 are the same
    newCat <- rep(0,length(label1))
    ## label1 and label2 must be at the same length
    for (ilabel1 in unique(label1))
      {
        for (ilabel2 in unique(label2))
          {
            newCat[label1==ilabel1 & label2 == ilabel2] <- paste(ilabel1,ilabel2,collapse=":",sep=":")
          }
      }
    return (as.factor(newCat))
  }

gsLabelnp <- function(d,olmeNBBayes,useSample){
  ## d: simulated dataset from getSample
  ## olmeNBBayes: must be the output of DPfit with model="nonpara"
  
  ## This function returns the estimated random effects of patients at each time point.
  ## the output is a vector of length sum ni
  ## for example if a patient ipat has 5 scans and the first two are treated as pre-scans and the rest as new-scans,
  ## then the d$labelnp[d$ID== opat] = 0,0,1,1,1
  ## the output is output[d$ID==ipat] = gPre,gPre,gNew,gNew,gNew
  ##            where gPre and gNew are the estimated random effects of patients 
  ## nonparametric changing random effects
  Npat <- length(unique(d$ID))
  nolnp0TF <- tapply(d$labelnp==0,d$ID,sum) == 0 ## TRUE if patients do not have prescans d$labelnp==0
  nolnp1TF <- tapply(d$labelnp==1,d$ID,sum) == 0 ## TRUE if patients do not have newscans d$labelnp==1
  patwG2 <- which(! (nolnp0TF|nolnp1TF)) ## patients with both new scans and pre scans
  
  gPres <- colMeans(olmeNBBayes$gPre[useSample,]); gPres_Index <- (1:Npat) - 0.5  ## everyone has g1
  gNews <- colMeans(olmeNBBayes$gNew[useSample, olmeNBBayes$gNew[1,]!= -1000,drop=FALSE]); gNews_Index <- patwG2
  ## combine gPres_gNews in the right order; in the order how they appear in the dataset d
  gPres_gNews <- c(gPres,gNews)[order(c(gPres_Index,gNews_Index))]
  g_long <- repeatAsID(values=gPres_gNews,ID=newCat(d$ID,d$labelnp))
  return (g_long)
}


repeatAsID <-function(values,ID)
  {
    ## ID does not have to be numerics
    ## values is a vector of length length(unique(ID))
    ivalue <- 1
    output <- rep(NA,length(ID))
    for (iID in unique(ID))
      {
        output[ID==iID] <- values[ivalue]
        ivalue <- ivalue + 1
      }
    return (output)
  }


slim <- function(vec,ID)
  {

    uniID <- unique(ID)
    re <- rep(NA,length(uniID))
    for (i in uniID) re[i] <- vec[ID==i][1]
    
    return(re)
  }


test <- function(x,mu,Sigma)
  {
    InvSigma <- solve(Sigma)
    evalue <- eigen(Sigma, symmetric = TRUE, only.values = TRUE)$values
    .Call("tempFun",
          as.numeric(x),
          as.numeric(mu),
          as.integer(length(mu)),
          as.numeric(evalue),
          as.numeric(c(InvSigma)),
          package ="lmeNBBayes")
  }



index.b.each <- function(Y,ID,labelnp,X,betas,aGs,rGs,pis)
  {
    M <- length(aGs)
    ## betas can be vector and X can be a matrix
    if (is.vector(X)) X <- matrix(X,ncol=1);
    ## compute the conditional probability of
    ## observing Y_{i,new+} >= y_{i,new+} given observing Y_{i,pre+}=y_{i,pew+}
    uniID <- unique(ID)
    Npat <- length(uniID)
    probs <- rep(NA,Npat);
    for (ipat in 1 : Npat)
      {
        Yi <- Y[ID==uniID[ipat]]
        Xi <- X[ID==uniID[ipat],,drop=FALSE]
        labelnpi <- labelnp[ID==uniID[ipat]]
        ## step 1: compute y_{i,pre+} and y_{i,new+} for all the patients:
        Ypresum <- sum(Yi[labelnpi==0])
        Ynewsum <- sum(Yi[labelnpi==1])
        if (sum(labelnpi==0)==0 || sum(labelnpi==1)==0 )
          {
            ## no new scans or no pre scans; no way to calculate the prob index!
            probs[ipat] <- NA
            next;
          }
        if (Ynewsum==0) ## If y_{i,new+} = 0 then the conditional prob of interest is one for that patient
          {
            probs[ipat] <- 1;
            next;
          }
        ## step 2: compute X_{ij}^T beta for all i and j
        rnewsum <- 0
        rpresum <- 0
        for (ivec in 1 : sum(ID==uniID[ipat]))
          {
            rij <- exp(sum(Xi[ivec,]*betas)) ## exp(Xij%*%beta)
            if (labelnpi[ivec]==0) rpresum <- rpresum + rij     ## sum_j exp(Xij%*%beta)
            else if (labelnpi[ivec]==1) rnewsum <- rnewsum + rij
          }
        num  <- 0
        ## step 3: compute the numerator sum_{k=0}^{y_{i,new+}-1} k
        for (k in 0:(Ynewsum-1))
          {
            BoverB <- 0
            for (ih in 1 : M )
              {
                if (pis[ih] == 0) next
                BoverB <- BoverB + pis[ih]*beta(rpresum+rnewsum+aGs[ih],k+Ypresum+rGs[ih])/beta(aGs[ih],rGs[ih])
              }
            num <- num + gamma(rnewsum+k)/(gamma(rnewsum)*gamma(k+1))*BoverB
          }
        ## step 4: compute denominator 
        den <- 0
        for (ih in 1 : M )
          {
            if (pis[ih] == 0) next
            den <- den + pis[ih]*beta(rpresum+aGs[ih],Ypresum+rGs[ih])/beta(aGs[ih],rGs[ih])
          }
        probs[ipat] <- 1- num/den
      }
    return(probs)
  }


numfun.norm <- function(g,t,Rnp,Rpp,Ypp,ph,ah,rh)
  {
    dnbinom(t,size=Rnp,prob=1/(1+g))*dnbinom(Ypp,size=Rpp,prob=1/(1+g))*ph*dgamma(g,shape=ah,scale=rh)
  }

denfun.norm <- function(g,Ypp,Rpp,ph,ah,rh)
  {
    dnbinom(Ypp,size=Rpp,prob=1/(1+g))*ph*dgamma(g,shape=ah,scale=rh)
  }

index.gammas <- function(Y,ID,labelnp,X,betas,
                         shapes, ## shape
                         scales, ## scale
                         pis)
  {
    M <- length(shapes)
    ## betas can be vector and X can be a matrix
    if (is.vector(X)) X <- matrix(X,ncol=1);
    ## compute the conditional probability of
    ## observing Y_{i,new+} >= y_{i,new+} given observing Y_{i,pre+}=y_{i,pew+}
    uniID <- unique(ID)
    Npat <- length(uniID)
    probs <- rep(NA,Npat);
    for (ipat in 1 : Npat)
      {
        Yi <- Y[ID==uniID[ipat]]
        Xi <- X[ID==uniID[ipat],,drop=FALSE]
        labelnpi <- labelnp[ID==uniID[ipat]]
        ## step 1: compute y_{i,pre+} and y_{i,new+} for all the patients:
        Ypresum <- sum(Yi[labelnpi==0])
        Ynewsum <- sum(Yi[labelnpi==1])
        if (sum(labelnpi==0)==0 || sum(labelnpi==1)==0 )
          {
            ## no new scans or no pre scans; no way to calculate the prob index!
            probs[ipat] <- NA
            next;
          }
        if (Ynewsum==0) ## If y_{i,new+} = 0 then the conditional prob of interest is one for that patient
          {
            probs[ipat] <- 1;
            next;
          }
        ## step 2: compute X_{ij}^T beta for all i and j
        rnewsum <- 0
        rpresum <- 0
        for (ivec in 1 : sum(ID==ipat))
          {
            rij <- exp(sum(Xi[ivec,]*betas))
            if (labelnpi[ivec]==0) rpresum <- rpresum + rij
            else if (labelnpi[ivec]==1) rnewsum <- rnewsum + rij
          }
        num  <- 0
        ## step 3: compute the numerator sum_{k=0}^{y_{i,new+}-1} k
        for (k in 0:(Ynewsum-1))
          {
            for (ih in 1 : M )
              {
                if (pis[ih] == 0) next
                num <- num + integrate(f=numfun.norm, lower=0, upper=Inf,
                                       t=k,Rnp=rnewsum,Rpp=rpresum,Ypp=Ypresum,ph=pis[ih],
                                       ah=shapes[ih],rh=scales[ih])$value


              }
          }
        ## step 4: compute denominator 
        den <- 0
        for (ih in 1 : M )
          {
            if (pis[ih] == 0) next
            den <- den +integrate(f=denfun.norm, lower=0, upper=Inf,
                                  Rpp=rpresum,Ypp=Ypresum,ph=pis[ih],
                                  ah=shapes[ih],rh=scales[ih])$value


          }
        probs[ipat] <- 1 - num/den
      }
    return(probs)
  }




numfun.YZ <- function(g,t,Rnp,Rpp,Ypp,ph,mu,sd)
  {
    dnbinom(t,size=Rnp,prob=1/(1+g))*
      dnbinom(Ypp,size=Rpp,prob=1/(1+g))*
        ph*dnorm(g,mean=mu,sd=sd)
  }

denfun.YZ <- function(g,Ypp,Rpp,ph,mu,sd)
  {
    dnbinom(Ypp,size=Rpp,prob=1/(1+g))*ph*dnorm(g,mean=mu,sd=sd)
  }

index.YZ <- function(Y,ID,labelnp,X,betas,
                     shape, ## shape
                     scale, ## scale
                     mu,sd,pi)
  {
    ## betas can be vector and X can be a matrix
    if (is.vector(X)) X <- matrix(X,ncol=1);
    ## compute the conditional probability of
    ## observing Y_{i,new+} >= y_{i,new+} given observing Y_{i,pre+}=y_{i,pew+}
    uniID <- unique(ID)
    Npat <- length(uniID)
    probs <- rep(NA,Npat);
    for (ipat in 1 : Npat)
      {
        Yi <- Y[ID==uniID[ipat]]
        Xi <- X[ID==uniID[ipat],,drop=FALSE]
        labelnpi <- labelnp[ID==uniID[ipat]]
        ## step 1: compute y_{i,pre+} and y_{i,new+} for all the patients:
        Ypresum <- sum(Yi[labelnpi==0])
        Ynewsum <- sum(Yi[labelnpi==1])
        if (sum(labelnpi==0)==0 || sum(labelnpi==1)==0 )
          {
            ## no new scans or no pre scans; no way to calculate the prob index!
            probs[ipat] <- NA
            next;
          }
        if (Ynewsum==0) ## If y_{i,new+} = 0 then the conditional prob of interest is one for that patient
          {
            probs[ipat] <- 1;
            next;
          }
        ## step 2: compute X_{ij}^T beta for all i and j
        rnewsum <- 0
        rpresum <- 0
        for (ivec in 1 : sum(ID==ipat))
          {
            rij <- exp(sum(Xi[ivec,]*betas))
            if (labelnpi[ivec]==0) rpresum <- rpresum + rij
            else if (labelnpi[ivec]==1) rnewsum <- rnewsum + rij
          }
        num  <- 0
        ## step 3: compute the numerator sum_{k=0}^{y_{i,new+}-1} k
        for (k in 0:(Ynewsum-1))
          {
            num <- num + integrate(f=numfun.norm, lower=0, upper=Inf,
                                   t=k,Rnp=rnewsum,Rpp=rpresum,Ypp=Ypresum,ph=pi,
                                   ah=shape,rh=scale)$value +
                                     integrate(f=numfun.YZ,lower=0,upper=Inf,
                                               t=k,Rnp=rnewsum,Rpp=rpresum,Ypp=Ypresum,ph=1-pi,
                                               mu=mu,sd=sd)$value
          }
        
        ## step 4: compute denominator 
        den <- 0
        
        den <- integrate(f=denfun.norm, lower=0, upper=Inf,
                         Rpp=rpresum,Ypp=Ypresum,ph=pi,
                         ah=shape,rh=scale)$value + integrate(f=denfun.YZ, lower=0, upper=Inf,
                                    Rpp=rpresum,Ypp=Ypresum,ph=1-pi,
                                    mu=mu,sd=sd)$value
        probs[ipat] <- 1 - num/den
      }
    return(probs)
  }


lmeNBBayes <- function(formula,          ##   A vector of length sum ni, containing responses
                       data,
                       ##   A sum ni by p matrix, containing covariate values. The frist column must be 1 (Intercept)
                       ID,         ##   A Vector of length sum ni, indicating patients
                       B = 105000, ##     A scalar, the number of Gibbs iteration 
                       burnin = 5000,  
                       printFreq = B,
                       M = NULL,
                       probIndex = FALSE,
                       thin =1, ## optional
                       labelnp=NULL, ## necessary if probIndex ==1
                       epsilonM = 0.001,## nonpara
                       para = list(mu_beta = NULL,Sigma_beta = NULL,a_D = 0.01, ib_D = 3,max_aG=30),
                       DP=TRUE,
                       seed=1,Reduce=1
                       )
  {
    set.seed(seed)
    
 
    X <- model.matrix(object=formula,data=data)
    covariatesNames<- colnames(X)
    Y <- model.response(model.frame(formula=formula,data=data))
    ## If ID is a character vector of length sum ni,
    ## it is modified to an integer vector, indicating the first appearing patient
    ## as 1, the second one as 2, and so on..

    ## This code generate samples from NBRE model with constant random effect ~ DP mixture of Beta
    if (is.vector(X)) X <- matrix(X,ncol=1)
    NtotAll <- length(Y)
    if (nrow(X)!= NtotAll) stop ("nrow(X) ! = length(Y)")
    if (length(ID)!= NtotAll)  stop ("length(ID)! = length(Y)")
    if (!is.null(labelnp) & length(labelnp)!= NtotAll)  stop ("labelnp! = length(Y)")

    
    dims <- dim(X)
    Ntot <- dims[1]
    pCov <- dims[2]
    ## cat("pCov",pCov)
    if (is.null(para$mu_beta))
      {
        para$mu_beta <- rep(0,pCov)
      }
    mu_beta <- para$mu_beta

    if (length(mu_beta)!=ncol(X))
      stop("The dimension of the fixed effect hyperparameter is wrong!")
    
    if (is.null(M)) M  <- round(1 + log(epsilonM)/log(para$ib_D/(1+para$ib_D)))
    if (is.null(para$Sigma_beta)) {
      para$Sigma_beta <-  diag(5,pCov)
    }
    
    evalue_sigma_beta <- eigen(para$Sigma_beta, symmetric = TRUE, only.values = TRUE)$values
    if (min(evalue_sigma_beta) <= 0) stop("Sigma_beta must be positive definite!")
    Inv_sigma_beta <- c( solve(para$Sigma_beta) )

    X <- c(X) ## {xij} = { x_{1,1},x_{2,1},..,x_{Ntot,1},x_{1,2},....,x_{Ntot,p} }

    ## change the index of ID to numeric from 1 to # patients
    temID <- ID  
    N <- length(unique(temID))
    uniID <- unique(temID)
    ID <- rep(NA,length(temID))
    for (i in 1 : length(uniID))
      {
        ID[temID == uniID[i]] <- i
      }
    
    mID <- ID-1
    ## the labelnp of patients only with 1 (new scans) labels are replaced by all 0 (old scans)
    
    maxni <- max(tapply(rep(1,length(ID)),ID,sum))
    Npat <- length(unique(ID))
    
    if (probIndex)
      {
        ## the labelnp of patients only with 1 (new scans) labels are replaced by all 0 (old scans)
        patwonew <- which(as.numeric(tapply((labelnp==0),ID,sum)==0)==1)
        for (i in 1 : length(patwonew)) labelnp[ID == patwonew[i]] <- 0
        
        patwoNorO <-  which(as.numeric(tapply((labelnp==1),ID,sum)==0)==1)
        if (length(patwoNorO)==0) patwoNorO <- -1000;
      }else{
        
        labelnp <- rep(0,length(Y))
      }
      

    if (B %% thin != 0 )
      stop("B %% thin !=0")
    if (burnin %% thin !=0)
      stop("burnin %% thin !=0")

    if (DP)
      {
        if (Reduce){
        re <- .Call("ReduceGibbs",
                    as.numeric(Y),           ## REAL
                    as.numeric(X),           ## REAL
                    as.integer(mID),         ## INTEGER
                    as.integer(B),           ## INTEGER
                    as.integer(maxni),       ## INTEGER
                    as.integer(Npat),        ## INTEGER
                    as.numeric(labelnp),     ## REAL
                    as.numeric(para$max_aG),
                    as.numeric(para$mu_beta),     ## REAL
                    as.numeric(evalue_sigma_beta),  ## REAL
                    as.numeric(Inv_sigma_beta),  ## REAL
                    as.numeric(para$a_D),
                    as.numeric(para$ib_D),
                    as.integer(M),
                    as.integer(burnin),      ## INTEGER
                    as.integer(printFreq),
                    as.integer(probIndex),
                    as.integer(thin),
                    package = "lmeNBBayes"
                    )

        for ( i in 13 : 14 ) re[[i]] <- matrix(re[[i]],B,Npat,byrow=TRUE)
        names(re) <- c("aGs","rGs","vs","weightH1",
                       "condProb","h1s","g1s",
                       "beta",
                       "D","logL",
                       "AR","prp","aGs_pat","rGs_pat")
      }else{
        re <- .Call("gibbs",
                    as.numeric(Y),           ## REAL
                    as.numeric(X),           ## REAL
                    as.integer(mID),         ## INTEGER
                    as.integer(B),           ## INTEGER
                    as.integer(maxni),       ## INTEGER
                    as.integer(Npat),        ## INTEGER
                    as.numeric(labelnp),     ## REAL
                    as.numeric(para$max_aG),
                    as.numeric(para$mu_beta),     ## REAL
                    as.numeric(evalue_sigma_beta),  ## REAL
                    as.numeric(Inv_sigma_beta),  ## REAL
                    as.numeric(para$a_D),
                    as.numeric(para$ib_D),
                    as.integer(M),
                    as.integer(burnin),      ## INTEGER
                    as.integer(printFreq),
                    as.integer(probIndex),
                    as.integer(thin),
                    package = "lmeNBBayes"
                    )
        names(re) <- c("aGs","rGs","vs","weightH1",
                       "condProb","h1s","g1s",
                       "beta",
                       "D","logL",
                       "AR","prp")
      }
        ## http://stackoverflow.com/questions/8720550/how-to-return-array-of-structs-from-call-to-c-shared-library-in-r
        ## for ( i in 1:4 ) re[[i]] <- matrix(re[[i]],B,M,byrow=TRUE)
        for ( i in 1 : 4 ) re[[i]] <- matrix(re[[i]],B,M,byrow=TRUE)
        re[[5]]  <- matrix(re[[5]],(B-burnin)/thin,Npat,byrow=TRUE)
        for ( i in 6 : 7 ) re[[i]] <- matrix(re[[i]],B,Npat,byrow=TRUE)
        re[[8]] <- matrix(re[[8]],B,pCov,byrow=TRUE)
 
        if (probIndex)
          {
            ## patients with no new scans
            if (sum(patwoNorO < 0) > 1) patwoNorO <- NULL
            re$condProb[,patwoNorO] <- NA
          }
        re$para$max_aG <- para$max_aG
        re$para$a_D <- para$a_D
        re$para$ib_D <- para$ib_D
        re$para$M <- M
        names(re$AR) <-c("aG", "rG","beta","D")
      }else{
        if (Reduce){
          re <- .Call("Beta1reduce",
                      as.numeric(Y),           ## REAL
                      as.numeric(X),           ## REAL
                      as.integer(mID),         ## INTEGER
                      as.integer(B),           ## INTEGER
                      as.integer(maxni),       ## INTEGER
                      as.integer(Npat),        ## INTEGER
                      as.numeric(labelnp),     ## REAL
                      as.numeric(para$max_aG),
                      as.numeric(para$mu_beta),     ## REAL
                      as.numeric(evalue_sigma_beta),  ## REAL
                      as.numeric(Inv_sigma_beta),  ## REAL
                      as.integer(burnin),      ## INTEGER
                      as.integer(printFreq),
                      as.integer(probIndex),
                      as.integer(thin),
                      package = "lmeNBBayes"
                    )

        }else{
        re <- .Call("Beta1",
                    as.numeric(Y),           ## REAL
                    as.numeric(X),           ## REAL
                    as.integer(mID),         ## INTEGER
                    as.integer(B),           ## INTEGER
                    as.integer(maxni),       ## INTEGER
                    as.integer(Npat),        ## INTEGER
                    as.numeric(labelnp),     ## REAL
                    as.numeric(para$max_aG),
                    as.numeric(para$mu_beta),     ## REAL
                    as.numeric(evalue_sigma_beta),  ## REAL
                    as.numeric(Inv_sigma_beta),  ## REAL
                    as.integer(burnin),      ## INTEGER
                    as.integer(printFreq),
                    as.integer(probIndex),
                    as.integer(thin),
                    package = "lmeNBBayes"
                    )
      }
        ## http://stackoverflow.com/questions/8720550/how-to-return-array-of-structs-from-call-to-c-shared-library-in-r
        ## for ( i in 1:4 ) re[[i]] <- matrix(re[[i]],B,M,byrow=TRUE)
        re[[3]] <- matrix(re[[3]],B,Npat,byrow=TRUE)
        re[[4]] <- matrix(re[[4]],B,pCov,byrow=TRUE)
        re[[7]]  <- matrix(re[[7]],(B-burnin)/thin,Npat,byrow=TRUE)
        names(re) <- c("aG","rG",
                       "g1s",
                       "beta",
                       "AR","prp",
                       "condProb",
                       "logL"
                       )
        if (probIndex)
          {
            ## patients with no new scans
            if (sum(patwoNorO < 0) > 1) patwoNorO <- NULL
            re$condProb[,patwoNorO] <- NA
          }
        re$para$max_aG <- para$max_aG
        names(re$AR) <-c("aG", "rG","beta")
        
      }
 

    if (probIndex)
      {
        re$para$labelnp <- labelnp
        re$para$CEL <- Y
        re$para$ID <- ID

        re$condProbSummary <- condProbCI(ID,re$condProb)
      }
    re$para$Sigma_beta <- para$Sigma_beta
    re$para$mu_beta <- para$mu_beta
    names(re$para$mu_beta) <- rownames(re$para$Sigma_beta) <-
      colnames(re$para$Sigma_beta) <- colnames(re$beta) <- covariatesNames
    re$para$B <- B
    re$para$burnin <- burnin
    re$para$thin <- thin
    re$para$probIndex <- probIndex
    re$para$Reduce <- Reduce
    re$para$burnin <- burnin
    re$para$DP <- DP
    re$para$formula <- formula
    class(re) <- "LinearMixedEffectNBBayes"
    return (re)
  }



print.LinearMixedEffectNBBayes <- function(x,...){

    para <- x$para
    if (is.vector(x$beta)) pCov <- 1 else pCov <- ncol(x$beta)
    cat("\n -----Negative binomial mixed effect regression----- ")
    if (para$DP){
      cat("\n ---random effect is from DP mixture of beta distributions---")
    }else{
      cat("\n ---random effect is from a single beta distribution---")
    }

    cat("\n====INPUT==== ")

    cat("\n formula: ");print(para$formula)
    if (is.null(para$thin))
      para$thin <- 1
    cat(paste(" MCMC parameters: B=",para$B,", burnin=",para$burnin,", thin=",para$thin,sep=""))
    if (para$Reduce)
      {
        cat("\n The target distribution in McMC is a partially marginalized posterior distribution")
        cat("\n (random effects are integrated out).")
      }
    cat("\n ------------------")
    cat("\n [[Hyperparameters]] ")
    cat("\n [Fixed Effect]       ")
    cat("\n beta ~ multivariate normal with mean and covariance matrix:\n  ")
    muSigma <- data.frame(para$mu_beta,"   |   ",para$Sigma_beta);
    colnames(muSigma) <- c("mean","   |   ","covariance",rep(" ",pCov-1))
    rownames(muSigma) <- names(para$mu_beta)
    print(muSigma)
    
    cat(paste(" [Random effect dist'n] aG,rG ~ Unif(min=",0.5,", max=",para$max_aG,")",sep=""))
    if (para$DP){
      cat(paste("\n [Precision parameter]      D ~ Unif(min=",para$a_D,", max=",para$ib_D,")",sep=""))
      cat(paste("\n [Truncation parameter]     M = ",para$M,sep=""))
    }

    
    if (!is.null(para$thin))
      {
        useSample <- useSamp(B=para$B,thin=para$thin,burnin=para$burnin)
        useSampleAll <- FALSE
      }
 
    cat("\n====OUTPUT====")

    cat("\n [Fixed effect]")
    if (para$DP){
       cat(" and [Precision parameter] \n")
       bb <- cbind(x$beta[useSample,],x$D[useSample])
       row.names =  c( names(para$mu_beta),"D")
    }else{
      cat("\n")
      bb <- cbind(x$beta[useSample,])
       row.names = names(para$mu_beta)
    }
    bsummary <- data.frame(colMeans(bb),apply(bb,2,sd),t(apply(bb,2,quantile,prob=c(0.025,0.975))),
                           row.names =  row.names)
    colnames(bsummary) <- c("Estimate","Std. Error","lower CI","upper CI")
    print(round(bsummary,4))
    cat("-----------------")
    cat("\n Reported estimates are posterior means computed after discarding burn-in")
    if (para$thin==1){
      cat(". No Thinning.")
    }else{
      cat(paste(" and thinning at every ",para$thin," samples.",sep=""))
    }
    cat("\n Posterior mean of the logLikelihood",mean(x$logL[useSample]))
    cat("\n Acceptance Rate in Metropolice Hasting rates")
    cat("\n ",paste(names(x$AR),round(x$AR,4),sep=":"))

    if (para$probIndex)
      {
        cat("\n-----------------\n")
        cat("\n ===== Conditional probability index ======\n")
        prIndex(x)
      }
  }




 print.IndexBatch <- prIndex <- function(x,...)
  {
    condProb <- x$condProbSummary
    ID <- x$para$ID
    uniID <- unique(ID)
    CEL <- x$para$CEL
    labelnp <- x$para$labelnp
    ## condProb:: A length(uniID) by 3 matrix, first column: a point estimate,
    ## the second column lower CI, the third column upper CI
    ## index of patient, reorder the patients in the increasing ways of probability index1 
    Suborderps <- order(condProb[,1])
    ## ps[order(ps)[1]] == min(ps)
    ## contains the index (1,2,...,Npat) of pat in decreasing way of index1
    orderps <- uniID[Suborderps]
    ## contains the ID (integer but the increment is not always 1) of pat in increasing way
    r <- NULL
    patwoNO <- patwNew0 <- rep(NA,length(uniID))
    ## m1G_1 must be > 0
    ## \hat{ E(1/G) } - 1
    ii <- rep(NA,length(uniID))
    for (i in 1 : length(uniID) )
      {
        i.orderp <- orderps[i]
        pickedIndex <- ID==i.orderp
        ## vector of length ID, containing TRUE if the position agree with the observations of ID orderps[i]
        pickedpt <- which(uniID==i.orderp)
        ## single element indicating the position of the orderps[i] in uniID
        CELpat <- CEL[pickedIndex]
        labelnppat <- labelnp[pickedIndex]
        CELpat0 <- CELpat[labelnppat==0]
        CELpat1 <- CELpat[labelnppat==1]
        if (length(CELpat0)==0) CELpat0 <- "--"
        if (length(CELpat1)==0)
          {
            CELpat1 <- "--"
            patwNew0[i] <- FALSE
          }else
        patwNew0[i] <- sum(CELpat1)>0 ## TRUE if patients have sum of CEL counts on new scan greater than zero
        ## Patients with no new/old scans are omitted at the end
        patwoNO[i] <- !is.na(condProb[pickedpt,1]) ## TRUE if patients have a new scan and old scan  
        r <- rbind(r,
                   c(
                     "pre-Scan"=paste(CELpat0,collapse="/"),
                     "new-Scan"=paste(CELpat1,collapse="/")
                     )
                   )
        ## ======== printable format =========;
        ii[i] <- paste(sprintf("%1.3f",condProb[pickedpt,1]),
                         " (",sprintf("%1.3f",condProb[pickedpt,3]),
                         ",",sprintf("%1.3f",condProb[pickedpt,4]),")",sep="")
      }
    
    outp <- cbind(ii,r)
    colnames(outp) <- c("conditional probability","pre-scan","new-scan")
    rownames(outp) <- orderps

    cat("\n Estimated conditional probability index and its 95 % CI.")
    if (!is.null(x$para$qsum))
      cat("\n The scalar function to summarize new scans:",x$para$qsum)
    cat("\n-----------------\n")
    print(data.frame(outp[patwoNO&patwNew0,]))
    cat("\n-----------------")
    cat("\n Patients are ordered by decreasing conditional probability.")
    cat("\n Patients with no pre-scans or no new scans are not reported.")
    cat("\n Patients whose new scans are all zero are not reported as conditional probability of such patients are always zero.")
  }
    





nbinDPREchange <- function(Y,          ##   A vector of length sum ni, containing responses 
                           X,          ##   A sum ni by p matrix, containing covariate values. The frist column must be 1 (Intercept)
                           ID,         ##   A Vector of length sum ni, indicating patients
                           B = 105000, ##     A scalar, the number of Gibbs iteration 
                           burnin = 5000,  
                           printFreq = B,
                           M = NULL,
                           labelnp, ## necessary if probIndex ==1
                           epsilonM = 0.01,## nonpara
                           para = list(mu_beta = NULL,Sigma_beta = NULL,a_D = 0.01, ib_D = 5,max_aG=30)
                           )
  {

    ## This code generate samples from NBRE model with constant random effect ~ DP mixture of Beta
    if (is.vector(X)) X <- matrix(X,ncol=1)
    NtotAll <- length(Y)
    if (nrow(X)!= NtotAll) stop("nrow(X) ! = length(Y)")
    if (length(ID)!= NtotAll)  stop("length(ID)! = length(Y)")
    if (length(labelnp)!= NtotAll)  stop ("labelnp! = length(Y)")

    
    dims <- dim(X)
    Ntot <- dims[1]
    pCov <- dims[2]
    ## cat("pCov",pCov)
    if (is.null(para$mu_beta))
      {
        para$mu_beta <- rep(0,pCov)
      }
    mu_beta <- para$mu_beta

    if (is.null(M)) M  <- round(1 + log(epsilonM)/log(para$ib_D/(1+para$ib_D)))
    if (is.null(para$Sigma_beta)) {
      para$Sigma_beta <-  diag(5,pCov)
    }
    Sigma_beta = para$Sigma_beta
    
    evalue_sigma_beta <- eigen(Sigma_beta, symmetric = TRUE, only.values = TRUE)$values
    if (min(evalue_sigma_beta) <= 0) stop("Sigma_beta must be positive definite!")
    Inv_sigma_beta <- c( solve(Sigma_beta) )

    X <- c(X) ## {xij} = { x_{1,1},x_{2,1},..,x_{Ntot,1},x_{1,2},....,x_{Ntot,p} }

    ## change the index of ID to numeric from 1 to # patients
    temID <- ID  
    N <- length(unique(temID))
    uniID <- unique(temID)
    ID <- rep(NA,length(temID))
    for (i in 1 : length(uniID))
      {
        ID[temID == uniID[i]] <- i
      }
    
    mID <- ID-1
    ## the labelnp of patients only with 1 (new scans) labels are replaced by all 0 (old scans)
    
    maxni <- max(tapply(rep(1,length(ID)),ID,sum))
    Npat <- length(unique(ID))


    ## the labelnp of patients only with 1 (new scans) labels are replaced by all 0 (old scans)
    patwonew <- which(as.numeric(tapply((labelnp==0),ID,sum)==0)==1)
    for (i in 1 : length(patwonew)) labelnp[ID == patwonew[i]] <- 0

    re <- .Call("gibbsREchange",
                as.numeric(Y),           ## REAL
                as.numeric(X),           ## REAL
                as.integer(mID),         ## INTEGER
                as.integer(B),           ## INTEGER
                as.integer(maxni),       ## INTEGER
                as.integer(Npat),        ## INTEGER
                as.numeric(labelnp),     ## REAL
                as.numeric(para$max_aG),
                as.numeric(para$mu_beta),     ## REAL
                as.numeric(evalue_sigma_beta),  ## REAL
                as.numeric(Inv_sigma_beta),  ## REAL
                as.numeric(para$a_D),
                as.numeric(para$ib_D),
                as.integer(M),
                as.integer(burnin),      ## INTEGER
                as.integer(printFreq),
                package = "lmeNBBayes"
                )


    ## http://stackoverflow.com/questions/8720550/how-to-return-array-of-structs-from-call-to-c-shared-library-in-r
    ## for ( i in 1:4 ) re[[i]] <- matrix(re[[i]],B,M,byrow=TRUE)
    for ( i in 1 : 4 ) re[[i]] <- matrix(re[[i]],B,M,byrow=TRUE)
    for ( i in 5 : 8 ) re[[i]] <- matrix(re[[i]],B,Npat,byrow=TRUE)
    re[[9]] <- matrix(re[[9]],B,pCov,byrow=TRUE)
    names(re) <- c("aGs","rGs","vs","weightH1",
                   "h1s","h2s","g1s","g2s",
                   "beta",
                   "logL","D",
                   "AR","prp")

    re$para <- para
    re$para$M <- M
    names(re$AR) <-names(re$prp) <- c("aG", "rG","beta","D")
    
    
    return (re)
  }


DPfit <- function(Y,          ##   A vector of length sum ni, containing responses 
                  X,          ##   A sum ni by p matrix, containing covariate values. The frist column must be 1 (Intercept)
                  ID,         ##   A Vector of length sum ni, indicating patients
                  B = 105000, ##     A scalar, the number of Gibbs iteration 
                  burnin = 5000,  
                  printFreq = B,
                  model = "nonpara-unif",
                  M = NULL,
                  epsilonM = 0.01,## nonpara
                  prob=0.9,        ## nonpara
                  labelnp=NULL,
                  para = list(mu_beta = NULL,Sigma_beta = NULL,a_D = 1, r_D = 1, a_qG = 1, r_qG = 1,
                    mu_aG = 0,sd_aG = 2, mu_rG = 0,sd_rG = 2, Nclust=NULL),
                  initBeta=NULL
                  )
  {

    Ntot <- length(Y)


    ## === prior input check =====## 
    if (is.vector(X)) X <- matrix(X,ncol=1)
    


    if (nrow(X)!= Ntot) stop ("nrow(X) ! = length(Y)")
    if (length(ID)!= Ntot)  stop ("length(ID)! = length(Y)")
    
    if (is.null(para$a_qG)){
      para$a_qG <- 1
      cat("\n a_qG needs to be specified!! set to 1")
    }
    a_qG = para$a_qG
    
    if (is.null(para$r_qG) ) {
      para$r_qG <- 1
      cat("\n r_qG needs to be specified!! set to 1")
    }
    r_qG = para$r_qG

    if (is.vector(X)) X <- matrix(X,ncol=1)
    pCov <- ncol(X)

    if (is.null(para$mu_beta))
      {
        para$mu_beta <- rep(0,pCov)
      }
    
    mu_beta = para$mu_beta
    
    if (is.null(para$Sigma_beta)) {
      para$Sigma_beta <-  diag(5,pCov)
    }
    Sigma_beta = para$Sigma_beta

    if (is.null(labelnp)) labelnp <- rep(0,length(Y))
    
    
    

    KK <- para$Nclust
    if (!is.null(KK)){
      if (!is.null(para$a_D)) stop("Both KK and a_D are selected! must be one of them...")
      ## the total number of patients observed by irev^th review
      NtotID <- length(unique(ID))
      ## the number of patients whose new scans are observed at the irev^th review
      NnewID <- length(unique(ID[labelnp==1]))
      ## parameters
      para$a_D <- uniroot(f=K_NisK,interval=c(0.003,50),
                          rD=para$r_D,N=(NtotID+NnewID)*2,K=KK)$root
      
    }
    a_D <- para$a_D
    
    if (is.null(M)) M  <- round(1 + log(0.01)/log(5/(1+5)))

    ## change the index of ID to numeric from 1 to # patients
    temID <- ID  
    N <- length(unique(temID))
    uniID <- unique(temID)
    ID <- rep(NA,length(temID))
    for (i in 1 : length(uniID))
      {
        ID[temID == uniID[i]] <- i
      }
    
    p <- pCov
    if( is.vector(X)) p <- 1 
    X <- c(X) ## {xij} = { x_{1,1},x_{2,1},..,x_{Ntot,1},x_{1,2},....,x_{Ntot,p} }
    mID <- ID-1

    ## The patients with labelnp = 1 (no old scans) for all repeated measures
    ## are treated as lack of new scans and used to estimate beta only.
    ## skip the computation of H2
    ## All patients have old scans
    
    ## the labelnp of patients only with 1 (new scans) labels are replaced by all 0 (old scans)
    patwonew <- which(as.numeric(tapply((labelnp==0),ID,sum)==0)==1)
    for (i in 1 : length(patwonew)) labelnp[ID == patwonew[i]] <- 0
    
    patwoNS <-  which(as.numeric(tapply((labelnp==1),ID,sum)==0)==1)
    if (length(patwoNS)==0) patwoNS <- -999;
    patwoNS <- patwoNS - 1
    
    maxni <- max(tapply(rep(1,length(ID)),ID,sum))
    Npat <- length(unique(ID))

    
    evalue_sigma_beta <- eigen(Sigma_beta, symmetric = TRUE, only.values = TRUE)$values
    if (min(evalue_sigma_beta) <= 0) stop("Sigma_beta must be positive definite!")
    Inv_sigma_beta <- c( solve(Sigma_beta) )
    
    if (model=="para-constantRE"| model==1)
      {
        
        
        if (is.null(para$mu_aG )){
          para$mu_aG <- 0.5
          cat("\n mu_aG needs to be specified!! set to 0.5")
        }
        mu_aG = para$mu_aG
        
        if (is.null(para$mu_rG)){
          para$mu_rG <- 0.5
          cat("\n mu_rG needs to be specified!! set to 0.5")
        }
        mu_rG = para$mu_rG
        
        if (is.null(para$sd_aG) ){
          para$sd_aG <- 2
          cat("\n mu_rG needs to be specified!! set to 2")
        }
        sd_aG = para$sd_aG
        
        if (is.null(para$sd_rG)) {
          para$sd_rG <- 2
          cat("\n mu_rG needs to be specified!! set to 2")
        }
        sd_rG = para$sd_rG
        
        ##parametric model: gi is constant over time 
        labelnp <- rep(0,length(Y))    
        re <- .Call("Beta1",
                    as.numeric(Y),           ## REAL
                    as.numeric(X),           ## REAL
                    as.integer(mID),         ## INTEGER
                    as.integer(B),           ## INTEGER
                    as.integer(maxni),       ## INTEGER
                    as.integer(Npat),        ## INTEGER
                    as.numeric(labelnp),     ## REAL
                    as.numeric(mu_aG),
                    as.numeric(sd_aG),
                    as.numeric(mu_rG),
                    as.numeric(sd_rG),
                    as.numeric(mu_beta),     ## REAL
                    as.numeric(evalue_sigma_beta),  ## REAL
                    as.numeric(Inv_sigma_beta),  ## REAL
                    as.integer(burnin),      ## INTEGER
                    as.integer(printFreq),
                    package = "lmeNBBayes"
                    )
        
        re[[3]] <- matrix(re[[3]],B,N,byrow=TRUE)
        re[[4]] <- matrix(re[[4]],B,p,byrow=TRUE)
        names(re) <- c("aG","rG","gPre","beta","AR","prp")
        para <- list(mu_beta = mu_beta,
                     Sigma_beta = Sigma_beta,
                     B=B,
                     model="para-constantRE",
                     burnin = burnin,
                     mu_aG = mu_aG,
                     sd_aG = sd_aG,
                     mu_rG = mu_rG,
                     sd_rG = sd_rG
                     )
        re$para <- para
        names(re$AR) <-c("aG", "rG", "beta")
        return(re)
        
      }else if (model=="para"){
        ## // Y_ij | Gij = gij ~ NB(size=exp(X_{ij}^T beta),prob=gij)
        ## // gij = g1 if j is in pre-scan 
        ## //     = g_new if j is in old-scan 
        ## // g_new = Ji * g1 + (1-Ji) * g2 
        ## // Ji ~ ber(qG)
        ## // qG ~ beta(a_qG,r_qG)
        ## // g1, g2 ~ beta(aG,rG)
        ## // beta ~ rnorm(mu_beta,sigma_beta)
        ## // aG, rG ~ lognorm(mu_aG,sd_aG),lognorm(mu_rG,sd_rG)
        
        
        if (is.null(para$mu_aG )){
          para$mu_aG <- 0.5
          cat("\n mu_aG needs to be specified!! set to 0.5")
        }
        mu_aG = para$mu_aG
        
        if (is.null(para$mu_rG)){
          para$mu_rG <- 0.5
          cat("\n mu_rG needs to be specified!! set to 0.5")
        }
        mu_rG = para$mu_rG
        
        if (is.null(para$sd_aG) ){
          para$sd_aG <- 2
          cat("\n mu_rG needs to be specified!! set to 2")
        }
        sd_aG = para$sd_aG
        
        if (is.null(para$sd_rG) ) {
          para$sd_rG <- 2
          cat("\n mu_rG needs to be specified!! set to 2")
        }
        sd_rG = para$sd_rG
        
        re <- .Call("Beta14",
                    as.numeric(Y),           ## REAL
                    as.numeric(X),           ## REAL
                    as.integer(mID),         ## INTEGER
                    as.integer(B),           ## INTEGER
                    as.integer(maxni),       ## INTEGER
                    as.integer(Npat),        ## INTEGER
                    as.numeric(labelnp),     ## REAL
                    as.numeric(a_qG),
                    as.numeric(r_qG),
                    as.numeric(mu_aG),
                    as.numeric(sd_aG),
                    as.numeric(mu_rG),
                    as.numeric(sd_rG),
                    as.numeric(mu_beta),     ## REAL
                    as.numeric(evalue_sigma_beta),  ## REAL
                    as.numeric(Inv_sigma_beta),  ## REAL
                    as.integer(burnin),      ## INTEGER
                    as.integer(printFreq),
                    as.integer(patwoNS),
                    package = "lmeNBBayes"
                    )
        ## http://stackoverflow.com/questions/8720550/how-to-return-array-of-structs-from-call-to-c-shared-library-in-r
        ## for ( i in 1:4 ) re[[i]] <- matrix(re[[i]],B,M,byrow=TRUE)
        for ( i in 4 : 7 ) re[[i]] <- matrix(re[[i]],B,N,byrow=TRUE)
        re[[8]] <- matrix(re[[8]],B,p,byrow=TRUE)
        names(re) <- c("aG","rG","qG",
                       "gPre","g2s","gNew","js",
                       "beta","AR","prp")
        para <- list(mu_beta = mu_beta,
                     Sigma_beta = Sigma_beta,
                     B=B,
                     model="para",
                     burnin = burnin,
                     mu_aG = mu_aG,
                     sd_aG = sd_aG,
                     mu_rG = mu_rG,
                     sd_rG = sd_rG,
                     a_qG=a_qG,
                     r_qG=r_qG
                     )
        re$para <- para
        names(re$AR) <-c("aG", "rG", "beta")
        return(re)
        
        
      }else if (model=="para-unif"){
        if (is.null(para$max_aG) )  para$max_aG <- 30
        max_aG <- para$max_aG
        ## // Y_ij | Gij = gij ~ NB(size=exp(X_{ij}^T beta),prob=gij)
        ## // gij = g1 if j is in pre-scan 
        ## //     = g_new if j is in old-scan 
        ## // g_new = Ji * g1 + (1-Ji) * g2 
        ## // Ji ~ ber(qG)
        ## // qG ~ beta(a_qG,r_qG)
        ## // g1, g2 ~ beta(aG,rG)
        ## // beta ~ rnorm(mu_beta,sigma_beta)
        ## // aG, rG ~ unif(0,max_aG)

        re <- .Call("Beta14Unif",
                    as.numeric(Y),           ## REAL
                    as.numeric(X),           ## REAL
                    as.integer(mID),         ## INTEGER
                    as.integer(B),           ## INTEGER
                    as.integer(maxni),       ## INTEGER
                    as.integer(Npat),        ## INTEGER
                    as.numeric(labelnp),     ## REAL
                    as.numeric(a_qG),
                    as.numeric(r_qG),
                    as.numeric(max_aG),
                    ## as.numeric(mu_aG),
                    ## as.numeric(sd_aG),
                    ## as.numeric(mu_rG),
                    ## as.numeric(sd_rG),
                    as.numeric(mu_beta),     ## REAL
                    as.numeric(evalue_sigma_beta),  ## REAL
                    as.numeric(Inv_sigma_beta),  ## REAL
                    as.integer(burnin),      ## INTEGER
                    as.integer(printFreq),
                    as.integer(patwoNS),
                    package = "lmeNBBayes"
                    )
        ## http://stackoverflow.com/questions/8720550/how-to-return-array-of-structs-from-call-to-c-shared-library-in-r
        ## for ( i in 1:4 ) re[[i]] <- matrix(re[[i]],B,M,byrow=TRUE)
        for ( i in 5 : 8 ) re[[i]] <- matrix(re[[i]],B,N,byrow=TRUE)
        re[[9]] <- matrix(re[[9]],B,p,byrow=TRUE)
        names(re) <- c("aG","rG","qG","logL",
                       "gPre","g2s","gNew","js",
                       "beta","AR","prp")
        para <- list(mu_beta = mu_beta,
                     Sigma_beta = Sigma_beta,
                     B=B,
                     model="para",
                     burnin = burnin,
                     max_aG = max_aG,
                     ## mu_aG = mu_aG,
                     ## sd_aG = sd_aG,
                     ## mu_rG = mu_rG,
                     ## sd_rG = sd_rG,
                     a_qG=a_qG,
                     r_qG=r_qG
                     )
        re$para <- para
        names(re$AR) <-c("aG", "rG", "beta")
        return(re)
      }else if (model== "nonpara" | model==8){

        ## nonparametric model
        ## // Y_ij | Gij = gij ~ NB(size=exp(X_{ij}^T beta),prob=gij)
        ## // gij = g1 if j is in pre-scan 
        ## //     = g_new if j is in old-scan 
        ## // g_new = Ji * g1 + (1-Ji) * g2 
        ## // Ji ~ ber(qG)
        ## // qG ~ beta(a_qG,r_qG)
        ## // g1, g2 ~ sum_{h=1}^M pi_h beta(aG_h,rG_h)
        ## // beta ~ mvrnorm(mu_beta,sigma_beta)
        ## // aG ~lognorm(mu_aG,sd_aG)
        ## // rG ~ lognorm(mu_rG,sd_rG)
        ##  nonparametric model,
        
        
        if (is.null(para$mu_aG )){
          para$mu_aG <- 0.5
          cat("\n mu_aG needs to be specified!! set to 0.5")
        }
        mu_aG = para$mu_aG
        
        if (is.null(para$mu_rG)){
          para$mu_rG <- 0.5
          cat("\n mu_rG needs to be specified!! set to 0.5")
        }
        mu_rG = para$mu_rG
        
        if (is.null(para$sd_aG) ){
          para$sd_aG <- 2
          cat("\n mu_rG needs to be specified!! set to 2")
        }
        sd_aG = para$sd_aG
        
        if (is.null(para$sd_rG) ) {
          para$sd_rG <- 2
          cat("\n mu_rG needs to be specified!! set to 2")
        }
        sd_rG = para$sd_rG
        if (is.null(para$r_D) ){
          cat("\n r_D needs to be specified!! set to 1")
          para$r_D <- 1
        }
        r_D = para$r_D
        re <- .Call("Beta24",
                    as.numeric(Y),           ## REAL
                    as.numeric(X),           ## REAL
                    as.integer(mID),         ## INTEGER
                    as.integer(B),           ## INTEGER
                    as.integer(maxni),       ## INTEGER
                    as.integer(Npat),        ## INTEGER
                    as.integer(M),           ## INTEGER
                    as.numeric(labelnp),     ## REAL
                    as.numeric(a_qG),        ## REAL
                    as.numeric(r_qG),        ## REAL
                    as.numeric(mu_aG),       ## REAL
                    as.numeric(sd_aG),       ## REAL
                    as.numeric(mu_rG),       ## REAL
                    as.numeric(sd_rG),       ## REAL
                    as.numeric(mu_beta),     ## REAL
                    as.numeric(evalue_sigma_beta),  ## REAL
                    as.numeric(Inv_sigma_beta),  ## REAL
                    as.numeric(a_D),         ## REAL
                    as.numeric(r_D),         ## REAL
                    as.integer(burnin),      ## INTEGER
                    as.integer(printFreq),
                    as.integer(patwoNS),
                    package = "lmeNBBayes"
                    )
        for ( i in 3 : 8 ) re[[i]] <- matrix(re[[i]],B,N,byrow=TRUE)
        for ( i in 9 : 12 ) re[[i]] <- matrix(re[[i]],B,M,byrow=TRUE)
        re[[13]] <- matrix(re[[13]],B,p,byrow=TRUE)
        names(re) <- c("qG","D",
                       "gPre","g2s","gNew","js","h1s","h2s",
                       "weightH1","vs","aGs","rGs",
                       "beta","AR","prp")
        re$para <- list(
                        burnin = burnin,
                        M=M,
                        B=B,
                        model=model,
                        Npat=N,
                        Ntot=length(Y),
                        para
                        )
        names(re$AR) <- names(re$prp) <- c(paste("aG",1:M,sep=""),
                                           paste("rG",1:M,sep=""),
                                           paste("beta",sep="")
                                           )
        re$h1s <- re$h1s + 1
        re$h2s <- re$h2s + 1
        re$MeanWH <- colMeans(re$weightH1);
        names(re$MeanWH) <- paste("cluster",1:M)
        return(re)
        
      }else if (model=="nonpara-unif"){
        
        ## nonparametric model D: is fixed
        ## // Y_ij | Gij = gij ~ NB(size=exp(X_{ij}^T beta),prob=gij)
        ## // gij = g1 if j is in pre-scan 
        ## //     = g_new if j is in old-scan 
        ## // g_new = Ji * g1 + (1-Ji) * g2 
        ## // Ji ~ ber(qG)
        ## // qG ~ beta(a_qG,r_qG)
        ## // g1, g2 ~ sum_{h=1}^M pi_h beta(aG_h,rG_h)
        ## // beta ~ mvrnorm(mu_beta,sigma_beta)
        ## // aG,rG ~ unif(0.0001,max_aG)
        ##  nonparametric model,
        if (is.null(para$max_aG) )  para$max_aG <- 30
        max_aG <- para$max_aG

        if (is.null(para$max_D))  para$max_D <- 5
        maxD <- para$max_D 

        if (is.null(initBeta)) initBeta <- rep(0, pCov)
        if (is.null(M)) M  <- 10 ## Need to fix this 
        
        re <- .Call("Beta24Unif",
                    as.numeric(Y),           ## REAL
                    as.numeric(X),           ## REAL
                    as.integer(mID),         ## INTEGER
                    as.integer(B),           ## INTEGER
                    as.integer(maxni),       ## INTEGER
                    as.integer(Npat),        ## INTEGER
                    as.integer(M),           ## INTEGER
                    as.numeric(labelnp),     ## REAL
                    as.numeric(a_qG),        ## REAL
                    as.numeric(r_qG),        ## REAL
                    as.numeric(max_aG),
                    ## as.numeric(mu_aG),    ## REAL
                    ## as.numeric(sd_aG),    ## REAL
                    ## as.numeric(mu_rG),    ## REAL
                    ## as.numeric(sd_rG),    ## REAL
                    as.numeric(mu_beta),     ## REAL
                    as.numeric(evalue_sigma_beta),  ## REAL
                    as.numeric(Inv_sigma_beta),  ## REAL
                    as.numeric(maxD),         ## REAL
                    ## as.numeric(r_D),         ## REAL
                    as.integer(burnin),      ## INTEGER
                    as.integer(printFreq),
                    as.integer(patwoNS),
                    as.numeric(initBeta),
                    package = "lmeNBBayes"
                    )
        for ( i in 4 : 9 ) re[[i]] <- matrix(re[[i]],B,N,byrow=TRUE)
        for ( i in 10 : 13 ) re[[i]] <- matrix(re[[i]],B,M,byrow=TRUE)
        re[[14]] <- matrix(re[[14]],B,p,byrow=TRUE)
        names(re) <- c("qG","D","logL",
                       "gPre","g2s","gNew","js","h1s","h2s",
                       "weightH1","vs","aGs","rGs",
                       "beta","AR","prp")
        
        para$mu_aG <- para$mu_rG <- para$sd_aG <- para$sd_rG <- para$a_D <- para$r_D <- NULL
        
        re$para <- list(
                        burnin = burnin,
                        M=M,
                        B=B,
                        model=model,
                        Npat=N,
                        Ntot=length(Y),
                        a_qG=a_qG,
                        r_qG=r_qG,
                        max_aG=max_aG,
                        mu_beta=mu_beta,
                        Sigma_beta=Sigma_beta,
                        maxD=maxD
                        )
        names(re$AR) <- names(re$prp) <- c(paste("aG",1:M,sep=""),
                                           paste("rG",1:M,sep=""),
                                           paste("beta",sep=""),
                                           "D"
                                           )
        re$h1s <- re$h1s + 1
        re$h2s <- re$h2s + 1
        re$MeanWH <- colMeans(re$weightH1);
        names(re$MeanWH) <- paste("cluster",1:M)
        return(re)

      }
  }



UsedPara <- function(olmeNBBayes,getSample=FALSE)
  {
    ## If getSample=TRUE, then outputDPFit must be the output of getSample with mod=4/3
    if (!getSample){
      o <- olmeNBBayes$para
    }else{
      o <- olmeNBBayes
      o$DP <- TRUE
    }

    UsedPara <- c(
                  "$B$"=o$B,burnin=o$burnin,
                  maxD=o$max_aG,
                  "$\\mu_{\\beta}$"=paste("(",paste(sprintf("%1.2f",c(o$mu_beta[1],o$mu_beta[2])),collapse=","),"...)",sep="",collapse=""),
                  "$\\Sigma_{\\beta}$"=paste("(",paste(sprintf("%1.2f",c(o$Sigma_beta[1,1],o$Sigma_beta[2,1])),collapse=","),
                    "...)")
                  )
    
      
    if (o$DP){
      UsedPara<-c(
                  UsedPara,"M"=o$M,
                  "$ib_D$"=sprintf("%1.2f",o$ib_D),"$a_D$"=sprintf("%1.2f",o$a_D),M=o$M
                  )
    }
    
    return (UsedPara)
    
  }




plotbeta <- function(shape1,shape2,poi=FALSE,col="red",lwd=1,lty=2)
  {
    xs <- seq(0,1,0.001)
    if (!poi)plot(xs,dbeta(xs,shape1,shape2),type="l",col=col,lwd=lwd,lty=lty)
    else points(xs,dbeta(xs,shape1,shape2),type="l",col=col,lwd=lwd,lty=lty)
  }
plotlnorm <- function(mulog,sdlog,poi=FALSE,col="red")
  {
    xs <- seq(0,10,0.001)
    if (!poi) plot(xs,dlnorm(xs,mulog,sdlog),type="l",col=col)
    else points(xs,dlnorm(xs,mulog,sdlog),type="l",col=col)
  }


useSamp <- function(thin,burnin,B)
  {
    temp <- (1:B)*thin
    unadj <- temp[burnin < temp & temp <= B]
    start <- unadj[1] - burnin -1
    return (unadj-start)
  }

dqmix <- function(weightH1,aGs,rGs,alphas=seq(0,1,0.01),dens=TRUE)
  {
    ## return a B by length(alphas) matrix, containing the quantile estimates at each corresponding alphas  
    ## weightH1s: B by M matrix
    ## aGs: B by M matrix
    ## rGs: B by M matrix
    B <- nrow(aGs)
    M <- ncol(aGs)
    pis <- c(t(weightH1))
    aGs <- c(t(aGs))
    rGs <- c(t(rGs))
    if (dens)
      {
        re <- .Call("mixdDist",
                    as.numeric(pis),
                    as.numeric(aGs),
                    as.numeric(rGs),as.numeric(alphas),
                    as.integer(B),as.integer(M),
                    package ="lmeNBBayes")
      }else{
        re <- .Call("mixQDist",
                    as.numeric(pis),
                    as.numeric(aGs),
                    as.numeric(rGs),as.numeric(alphas),
                    as.integer(B),as.integer(M),
                    package ="lmeNBBayes")
      }
    for (i in 1 : length(re)) re[[i]] <- matrix(re[[i]],nrow=B,ncol=length(alphas),byrow=TRUE)
    return (re)
  }


momentsBeta <- function(aG1,rG1,beta)
  {
    ## Y ~ NB(size=rij,prob=G)
    ## G ~ Beta(aG1,rG1)
    ## E(1/G) = ...some calculations... = (aG+rG-1)/(aG-1)
    ## E(Y)=exp(beta0)*(mu_{1/G}+1)=exp(1)*((aG1+rG1-1)/(aG1-1)+1)
    ## Var(1/G) = E(1/G^2) - E(1/G)^2 = (aG1+rG1-1)/(aG1-1)*( (aG1+rG1-2)/(aG1-2) - (aG1+rG1-1)/(aG1-1) ) = 0.06559506
    ## Var(Y) = rij*Var(1/G)*(rij+1) + rij*E(1/G)*(E(1/G)-1)
    ## calculate E(Y), Var(Y), E(1/G), and E(1/G^2) at initial time
    rij <- exp(beta) ## value of size parmeter at time zero
    E1G <- (aG1+rG1-1)/(aG1-1)
    V1G <- (aG1+rG1-1)/(aG1-1)*( (aG1+rG1-2)/(aG1-2) - (aG1+rG1-1)/(aG1-1) )
    EY <- rij*(E1G-1)
    VY <- rij*V1G*(rij+1)+rij*E1G*(E1G-1)
    return (c(EY=EY,SDY=sqrt(VY),
              E1G=E1G,V1G=V1G))
  }


momentsGL <- function(EG,VG,beta)
  {
    ## Y ~ NB(size=rij,prob=1/(G+1))
    ## G ~ dist(mean=EG,var=VG)
    rij <- exp(beta)
    ## E(Y|G) = rij*G
    ## Var(Y|G) = rij*G*(G+1)
    ## E(Y)=E(E(Y|G))=rij*E(G)
    EY <- rij*EG
    ## Var(Y) = Var(E(Y|G)) + E(Var(Y|G))
    ##        = Var(rij*G) + E(G*(G+1)*rij)
    ##        = rij^2*Var(G) + rij*(E(G^2)+E(G))
    ##        = rij^2*Var(G) + rij*(Var(G)+E(G)^2+E(G))
    VY <- rij^2*VG + rij*(VG+EG^2+EG)
    return (c(EY=EY,SDY=sqrt(VY),
              EG=EG,VG=VG))
  }


DevRatio <- function(gPre,gNew) ## gPre, gNew must be a matrix of length(useSample) by N
  {
    ## The deviation of CEL counts of single patient i from overall cohort at week j could be measured by:
    ## R_ij=E(Y_ij|G_ij)/E(Y_ij)= r_ij (1-g_ij)/gij * 1/(rij*(mu_{1/G}-1)) = (1-g_ij)/(gij*(mu_{1/G}-1))
    ## The change in the deviation R_ij at two time point could be measured by
    ## R_ij/R_ij' = (1-g_ij)/(gij*(mu_{1/G}-1))*(gij'*(mu_{1/G}-1))/(1-g_ij) = g_ij'*(1-g_ij)/(g_ij*(1-g_ij'))
    ## Hence the deviation ratio between the pre scan period and new scan period is:
    ## R_inew/R_ipre
    ## How much patient i increases the deviation from the overall trend between the pre-scan period and new-scan period

    ## the elemnt-wise operations
    m1G_1 <- mean(c(1/gPre,1/gNew)) - 1
    
    EYGpre <- 1/gPre - 1
    EYGnew <- 1/gNew - 1


    CI95_EYGpre <- apply(EYGpre,2,quantile,prob=c(0.5,0.025,0.975)) / m1G_1
    CI95_EYGnew <- apply(EYGnew,2,quantile,prob=c(0.5,0.025,0.975)) / m1G_1

    CI95_EYGpre[1,] <- colMeans(EYGpre)
    CI95_EYGnew[1,] <- colMeans(EYGnew)
    
    ##DevMat <- EYGnew/EYGpre ## B-Bburn by N
    ##CI95 <- apply(DevMat,2,quantile,prob=c(0.5,0.025,0.975))

    NewPre <-  CI95_EYGnew[1,]/CI95_EYGpre[1,]
    
    ##CI95[1,] <- colMeans(DevMat)
    
    ##rownames(CI95)[1] <- "mean"
    rownames(CI95_EYGpre)[1] <- "mean"
    rownames(CI95_EYGnew)[1] <- "mean"
    if (!is.null(rownames(gPre)))
      {
        names(CI95) <- colnames(gPre)
        colnames(CI95_EYGpre) <- colnames(gPre)
        colnames(CI95_EYGnew) <- colnames(gPre)
      }
    re <- list(Ratio=NewPre,
               Rpre=CI95_EYGpre,
               Rnew=CI95_EYGnew)
    
    return (re)
  }

paralnorm <- function(E,V=15^2){
  ## E: expectation of log-normally distributed r.v.
  ## V: variance of log-normally distributed r.v

  sd <- sqrt(log(V/E+1))
  mu <- log(E)-sd/2
  return(list(mu=mu,sd=sd))
}
F_inv_gam <- function(p,sp1,sc1,sp2,sc2,pi)
  {
    G = function (t)
      pi*pgamma(t,shape=sp1,scale=sc1)+(1-pi)*pgamma(t,shape=sp2,scale=sc2) - p
                return(uniroot(G,c(0,100))$root)
  }

getSampleS <- function(iseed=1,propTreat=2/3,RedFactor=0.5,detail=FALSE)
  {
    ## The subset of patients in treatment arm experience the treatment effect
    ## Those with treatment effect decrease the expectation of CEL counts by RedFactor at every time point.
    ## RedFactor is constant over time.

    ## propTreat is the proportion of the patients in treatment arm with treatment effect
    
    ## if full = TRUE then full dataset is returned and rev is ignored
    ni <- 8
    Npat <- 120;
    i.trt <- c(
               rep(1,round(Npat*propTreat*0.5)), ## 0.5 because half of the patients receive treatments
               rep(0,Npat-round(Npat*propTreat*0.5))
               )
    pi <- 0.3
    
    beta0 <- 1.5
    beta1 <- -0.05
    
    aG1 <- 10
    rG1 <- 10
    aG2 <- 20
    rG2 <- 1

    ## generate the initial random effect values of everyone
    hs <- sample(1:2,Npat,c(pi,1-pi),replace=TRUE)
    Npat_dist1 <- sum(hs==1)
    Npat_dist2 <- sum(hs==2)
    
    gsBASE <- rep(NA,Npat)
    gsBASE[hs==1] <- rbeta(Npat_dist1,aG1,rG1);
    gsBASE[hs==2] <- rbeta(Npat_dist2,aG2,rG2);

    ## generate samples for dist = "b","b3", "g", and "l" (except "b2")
    Y <- NULL
    timecov <- c(0,0,(1:(ni-2)))
    for ( ipat in 1 : Npat)
      {
        ## the number of repeated measures are the same
        ## we assume that the treatment effects occurs once after the screening and baseline
        got <- rnbinom(ni,size = exp(beta0+beta1*timecov)*RedFactor^(i.trt[ipat]*(timecov>0)), prob = gsBASE[ipat])         
        Y <- c(Y,got)
      }
    ## cat("beta0",beta0,"\n")
    d <- data.frame(Y=Y,
                    X1=rep(1,Npat*ni),
                    X2=rep(timecov,Npat),
                    labelnp=rep(c(rep(0,2),rep(1,ni-2)),Npat), ## the first two scans are pre-scans 
                    ID=rep(1:Npat,each=ni),
                    gsBASE = rep(gsBASE,each=ni),
                    Treat = rep(1:0,each=(Npat/2)*ni), ## the first half of the patients receive treatment
                    effectTreat = rep(i.trt,each=ni),
                    hs = rep(hs,each=ni)
                    )

    if (detail)
      {
        E1G.1 <- momentsBeta(aG1,rG1,beta0)[3] ## E1G
        E1G.2 <- momentsBeta(aG2,rG2,beta0)[3]
        
        E1G <- E1G.1*pi+E1G.2*(1-pi)

        E1s<-c(E1G,E1G.1,E1G.2)
        names(E1s) <- c("Overall","H=1","H=2")
        for (i in 1 : length(E1s))
          {
            noTreatEffect <- (E1s[i]-1)*exp(beta0+beta1*timecov)*RedFactor^(0*(timecov>0))
            TreatEffect <- (E1s[i]-1)*exp(beta0+beta1*timecov)*RedFactor^(1*(timecov>0))
            
            rr <- rbind(noTreatEffect,TreatEffect)
            colnames(rr) <- timecov
            cat("\n The expected CEL counts E(Yt)",names(E1s)[i],"\n")
            print(rr)
          }
      }
    
    return(d)
  }

getSample <- function(
                      iseed = 911,
                      rev = 4,
                      dist = "b",
                      mod = 0,
                      probs = seq(0,0.99,0.01),
                      ts = seq(0.001,0.99,0.001),
                      full = FALSE
                      )
  {
    ## mod = 0: generate sample
    ## mod = 1: quantiles of the true populations
    ## mod = 2: densities of the true populations
    ## mod = 3: parameters of the simulation model
    ## mod = 4: parameters for uninformative prior

    ## dist = "b","b2","g" 

    ## if full = TRUE then full dataset is returned and rev is ignored
    ni <- 11
    Npat <- 180; ## upto review 4, Npat=160 in total this number must be divisible by NpatEnterPerMonth
    NpatEnterPerMonth <- 15
    DSMBVisitInterval <- 4 ## months
    
    varInfoP <- c(0.1312,0.0148)
    ib_D <- 3
    ## === necessary for mod= 0, 3 and mod = 4
    ## d contains a full dataset
    days <- NULL
    for (ipatGroup in 1 : (Npat/NpatEnterPerMonth))
      {
        ScandaysForSingleGroup <- ipatGroup:(ipatGroup+ni-1)
        days <- c(days,rep(ScandaysForSingleGroup,NpatEnterPerMonth))
      }

    set.seed(iseed)
    if (dist=="g") ## need to went through by all mod=0,...,4
      {
        ## r.e. mixture of two gamma distributions
        pi <- 0.7
        sp1 <- 0.3; sc1 <- 0.5
        sp2 <- 4; sc2 <- 0.8
        beta0 <- 0.5
        beta1 <- 0
        hs <- sample(1:2,Npat,c(pi,1-pi),replace=TRUE)
        Npat_dist1 <- sum(hs==1)
        Npat_dist2 <- sum(hs==2)

        gsBASE <- rep(NA,Npat)
        gsBASE[hs==1] <- rgamma(Npat_dist1,shape=sp1,scale=sc1)
        gsBASE[hs==2] <- rgamma(Npat_dist2,shape=sp2,scale=sc2)
        gsBASE <- 1/(1+gsBASE)
        
        if (mod==1)
          {
            lq <- length(probs);
            quan <- rep(NA,lq);
            for (i in 1 : lq) quan[i] <- F_inv_gam(p=probs[i],sp1=sp1,sc1=sc1,sp2=sp2,sc2=sc2,pi=pi)
            quan = sort(1/(1+quan))
            
            return (rbind(probs=probs,quantile=quan))

          }else if (mod == 2){
            
            ts.trans <-1/ts-1
            return ((pi*dgamma(ts.trans,shape=sp1,scale=sc1)+(1-pi)*dgamma(ts.trans,shape=sp2,scale=sc2) )*(1/ts^2))
            
          }else if (mod == 3){
            
            c1 <- momentsGL(EG=sp1*sc1,VG=sp1*sc1^2,beta=beta0) ## expectation and variance of Y
            c2 <- momentsGL(EG=sp2*sc2,VG=sp2*sc2^2,beta=beta0)
            
            outputMod3 <- list(infoPara = list(mu_beta=c(beta0,beta1),
                                 Sigma_beta=diag(varInfoP),a_qG=100,r_qG=1),
                               beta0=beta0,beta1=beta1,K=2,
                               scales=c(sc1=sc1,sc2=sc2),
                               shapes=c(sp1=sp1,sp2=sp2),
                               c1=c1,c2=c2
                               )
          }
      }else if (dist == "b"){
       
        beta0 <- 1
        beta1 <- -0.05
        aG1 <- 3
        rG1 <- 0.8
        
        gsBASE <- rbeta(Npat,aG1,rG1)
        hs <- rep(1,Npat)
        if (mod==1){
          return( 
                 rbind(probs=probs,
                       quantile=qbeta(probs,shape1=aG1,shape2=rG1)
                       )
                 )
        }else if (mod==2){
          
          return (dbeta(ts,shape1=aG1,shape2=rG1))
          
        }else if (mod == 3){
          
          c1 <- momentsBeta(aG1,rG1,beta0)
          outputMod3 <- list(infoPara = list(mu_beta=c(beta0,beta1),
                               Sigma_beta=diag(varInfoP),
                               a_D=0.01,ib_D=ib_D,max_aG=30
                               ),##a_qG=10,r_qG=1),
                             beta0=beta0,beta1=beta1,K=1,c1=c1, aGs=c(aG1=aG1),
                             rGs=c(rG1=rG1))
        }

      }else if (dist == "b2")
        {
        ## mixture of two beta distributions
        ## cluster 1: E(Y)=exp(beta0)*mu_G=exp(2)*((aG1+rG1-1)/(aG1-1)+1)=3.098636 ## at initial time
        ## cluster 2: E(Y)=exp(beta0)*mu_G=exp(2)*((aG2+rG2-1)/(aG2-1)+1)=7.037196 ## at initial time
        pi <- 0.3
        
        beta0 <- 1.5
        beta1 <- -0.05
        
        aG1 <- 10
        rG1 <- 10
        aG2 <- 20
        rG2 <- 1

        ## generate the initial random effect values of everyone
        hs <- sample(1:2,Npat,c(pi,1-pi),replace=TRUE)
        Npat_dist1 <- sum(hs==1)
        Npat_dist2 <- sum(hs==2)
        
        gsBASE <- rep(NA,Npat)
        gsBASE[hs==1] <- rbeta(Npat_dist1,aG1,rG1);
        gsBASE[hs==2] <- rbeta(Npat_dist2,aG2,rG2);
        if (mod==1){
          F_inv_beta2 <- function(p,aG1,rG1,aG2,rG2,pi)
            {
              G = function (t) pi*pbeta(t,shape1=aG1,shape2=rG1)+(1-pi)*pbeta(t,shape1=aG2,shape2=rG2) - p
              return(uniroot(G,c(0,1000))$root)
            }
          lq <- length(probs);
          quant <- rep(NA,lq);
          for (i in 1 : lq) quant[i] <- F_inv_beta2(p=probs[i],aG1,rG1,aG2,rG2,pi=pi)
          return (rbind(probs=probs,quantile=quant))
        }
        else if (mod==2){

          return ((pi*dbeta(ts,shape1=aG1,shape2=rG1)+(1-pi)*dbeta(ts,shape1=aG2,shape2=rG2)) )
        }else if (mod == 3){

          c1 <- momentsBeta(aG1,rG1,beta0)
          c2 <- momentsBeta(aG2,rG2,beta0)
          outputMod3 <-  list(infoPara = list(mu_beta=c(beta0,beta1),
                                Sigma_beta=diag(varInfoP),a_D=0.01,ib_D=ib_D,max_aG=30),
                              beta0=beta0,beta1=beta1,K=2,c1=c1,c2=c2,
                              aGs=c(aG1=aG1,aG2=aG2),
                              rGs=c(rG1=rG1,rG2=rG2),
                              pi=pi
                              )
        }
      }
    if ( dist == "YZ")
      {
        pi <- 0.85

        alpha <- exp(-0.5)
        ## logalpha <- -0.5
        beta0 <- 0.905 ##0.405+0.5 beta - logalpha
        beta1 <- 0
        ## a bimodal distribution with 85 % of Gi from a gamma distribution with mean 0.647 and variance 2.374
        scale <- 2.374/0.647*alpha 
        shape <- 0.647^2/2.374
        mu <- 3*alpha
        sd <- sqrt(0.25)*alpha
       
        ## generate the initial random effect values of everyone
        hs <- sample(1:2,Npat,c(pi,1-pi),replace=TRUE)
        Npat_dist1 <- sum(hs==1)
        Npat_dist2 <- sum(hs==2)
        
        gsBASE <- rep(NA,Npat)
        gsBASE[hs==1] <- rgamma(Npat_dist1,shape=shape,scale=scale);
        gsBASE[hs==2] <- rnorm(Npat_dist2,mean=mu,sd=sd);
        gsBASE[gsBASE < 0 ] <- 0
        gsBASE <- 1/(1+gsBASE)
        
        if (mod==1)
          {
            return (NULL)

          }else if (mod == 2){
            
            ts.trans <-1/ts-1
            return ((pi*dgamma(ts.trans,shape=shape,scale=scale)+(1-pi)*dnorm(ts.trans,mean=mu,sd=sd) )*(1/ts^2))
            
          }else if (mod == 3){
           
            outputMod3 <- list(infoPara = list(mu_beta=c(beta0,beta1),
                                 Sigma_beta=diag(varInfoP),a_D=0.01,ib_D=ib_D,max_aG=30),
                               beta0=beta0,beta1=beta1,scale=scale,shape=shape,mu=mu,sd=sd,pi=pi,K=2
                               )
          }
      }
      
    ## generate samples for dist = "b","b3", "g", and "l" (except "b2")
    Y <- NULL
    timecov <- c(0,0,(1:(ni-2)))
    for ( ipat in 1 : Npat)
      {
        ## the number of repeated measures are the same
        ## we assume that the time effects occurs once after the treatments are in effect
        got <- rnbinom(ni,size = exp(beta0+beta1*timecov), prob = gsBASE[ipat])         
        Y <- c(Y,got)
      }
    ##cat("beta0",beta0,"\n")
    d <- data.frame(Y=Y,
                    X1=rep(1,Npat*ni),
                    X2=rep(timecov,Npat),
                    ID=rep(1:Npat,each=ni),
                    gsBASE = rep(gsBASE,each=ni),
                    scan = rep(1:ni,Npat),
                    ## day contains the day when the scan was taken
                    ## 10 patients enter a trial every month
                    days = days,
                    hs = rep(hs,each=ni)
                    )

    ## dSMB visit is assumed to be every 4 months
    if (full) return (d) 
    d <- subset(d,subset= days <= DSMBVisitInterval*rev)
    d$labelnp <- rep(0,nrow(d))
    d$labelnp[ DSMBVisitInterval*(rev-1) < d$days ] <- 1
    ## The first two scans (screening and base-line scans are treated as pre-scans)
    d$labelnp[ d$X2==0 ] <- 0

    if (dist=="b2")
      {
        temp <- index.b.each(Y=d$Y,ID=d$ID,labelnp=d$labelnp,X=cbind(d$X1,d$X2),
                             betas=c(beta0,beta1),aGs=c(aG1,aG2),rGs=c(rG1,rG2),pis=c(pi,1-pi))
        d$probIndex <- c(temp,rep(NA,nrow(d)-length(temp) ))
      }else if (dist=="b")
        {
          temp <- index.b.each(Y=d$Y,ID=d$ID,labelnp=d$labelnp,X=cbind(d$X1,d$X2),
                               betas=c(beta0,beta1),aGs=c(aG1),rGs=c(rG1),pis=c(1))
          d$probIndex <- c(temp,rep(NA,nrow(d)-length(temp)) )
        }else if (dist=="g"){

          temp <- index.gammas(Y=d$Y,ID=d$ID,labelnp=d$labelnp,X=cbind(d$X1,d$X2),
                              betas=c(beta0,beta1),
                              shapes=c(sc1,sc2), ## shape
                              scales=c(sp1,sp2), ## scale
                              pis=c(pi,1-pi))
          d$probIndex <- c(temp,rep(NA,nrow(d)-length(temp)) )
        }else if (dist == "YZ"){
          
          temp <- index.YZ(Y=d$Y,ID=d$ID,labelnp=d$labelnp,X=cbind(d$X1,d$X2),
                           betas=c(beta0,beta1),
                           shape=shape, ## shape
                           scale=scale, ## scale
                           mu=mu,
                           sd=sd,
                           pi=pi)
          d$probIndex <- c(temp,rep(NA,nrow(d)-length(temp)) )
            
          }

    if (mod==0) return(d)

    if (mod==3){
      
      Npat <- length(unique(d$ID))
      outputMod3$infoPara$maxD <- 5
      return(outputMod3)

    }else if (mod == 4){
      
      Npat <- length(unique(d$ID))
      mu_beta <- c(0,0)
      Sigma_beta <- diag(c(5,5))
      a_qG <- 1
      r_qG <- 1
      
      outputMod4 <- list(
                         mu_beta=mu_beta,Sigma_beta=Sigma_beta,
                         ##a_qG=a_qG,r_qG=r_qG,
                         a_D = 0.01, ib_D=ib_D,
                         ##maxD=5,
                         max_aG = 30
                         )
      return(list(UinfoPara=outputMod4))
    }
  }



## getSampleOld <- function(
##                       iseed=911,
##                       rev=4,
##                       dist="b",
##                       mod=0,
##                       probs=seq(0,0.99,0.01),
##                       ts = seq(0,0.99,0.01),
##                       full= FALSE
##                       )
##   {
##     ## mod = 0: generate sample
##     ## mod = 1: quantiles of the true populations
##     ## mod = 2: densities of the true populations
##     ## mod = 3: parameters of the simulation model
##     ## mod = 4: parameters for uninformative prior

##     ## dist = "b","b2","b3","l","g" 

##     ## if full = TRUE then full dataset is returned and rev is ignored
##     ni <- 11
##     Npat <- 180; ## upto review 4, Npat=160 in total this number must be divisible by NpatEnterPerMonth
##     NpatEnterPerMonth <- 15
##     DSMBVisitInterval <- 4 ## months
    

##     ## === necessary for mod= 0, 3 and mod = 4
##     ## d contains a full dataset
##     days <- NULL
##     for (ipatGroup in 1 : (Npat/NpatEnterPerMonth))
##       {
##         ScandaysForSingleGroup <- ipatGroup:(ipatGroup+ni-1)
##         days <- c(days,rep(ScandaysForSingleGroup,NpatEnterPerMonth))
##       }

##     set.seed(iseed)
##     if (dist=="g") ## need to went through by all mod=0,...,4
##       {
##         ## r.e. mixture of two gamma distributions
##         pi <- 0.7
##         sp1 <- 0.3; sc1 <- 0.5
##         sp2 <- 4; sc2 <- 0.8
##         beta0 <- 0.5
##         beta1 <- 0
##         hs <- sample(1:2,Npat,c(pi,1-pi),replace=TRUE)
##         Npat_dist1 <- sum(hs==1)
##         Npat_dist2 <- sum(hs==2)

##         gsBASE <- rep(NA,Npat)
##         gsBASE[hs==1] <- rgamma(Npat_dist1,shape=sp1,scale=sc1)
##         gsBASE[hs==2] <- rgamma(Npat_dist2,shape=sp2,scale=sc2)
##         gsBASE <- 1/(1+gsBASE)
        
##         if (mod==1)
##           {
##             F_inv_gam <- function(p,sp1,sc1,sp2,sc2,pi)
##               {
##                 G = function (t)
##                   pi*pgamma(t,shape=sp1,scale=sc1)+(1-pi)*pgamma(t,shape=sp2,scale=sc2) - p
##                 return(uniroot(G,c(0,100))$root)
##               }
##             lq <- length(probs);
##             quan <- rep(NA,lq);
##             for (i in 1 : lq) quan[i] <- F_inv_gam(p=probs[i],sp1=sp1,sc1=sc1,sp2=sp2,sc2=sc2,pi=pi)
##             quan = sort(1/(1+quan))
            
##             return (rbind(probs=probs,quantile=quan))

##           }else if (mod == 2){
            
##             ts.trans <-1/ts-1
##             return ((pi*dgamma(ts.trans,shape=sp1,scale=sc1)+(1-pi)*dgamma(ts.trans,shape=sp2,scale=sc2) )*(1/ts^2))
            
##           }else if (mod == 3){
            
##             c1 <- momentsGL(EG=sp1*sc1,VG=sp1*sc1^2,beta=beta0) ## expectation and variance of Y
##             c2 <- momentsGL(EG=sp2*sc2,VG=sp2*sc2^2,beta=beta0)
            
##             outputMod3 <- list(infoPara = list(mu_beta=c(beta0,beta1),
##                                  Sigma_beta=diag(c(0.05,0.05)),a_qG=100,r_qG=1),
##                                beta0=beta0,beta1=beta1,K=2,
##                                scales=c(sc1=sc1,sc2=sc2),
##                                shapes=c(sp1=sp1,sp2=sp2),
##                                c1=c1,c2=c2
##                                )
##           }
##       }else if (dist== "l"){
        
##         ## single log-normal distribution
##         ## E(Y)=exp(beta)*E(G)=exp(1)*exp(meanlog+sdlog^2/2)= 5.078419
##         beta0 <- 1.2
##         beta1 <- -0.1
##         meanlog <- -2
##         sdlog <- 1.5
##         gsBASE <- rlnorm(Npat,meanlog=meanlog,sdlog=sdlog)
##         gsBASE <- 1/(1+gsBASE)
##         quan <- qlnorm(probs,meanlog=meanlog,sdlog=sdlog)
##         hs <- rep(1,Npat)
##         if (mod==1){
##           return( 
##                  rbind(probs=probs,
##                        quantile=sort(1/(1+quan))
##                        )
##                  )
##         }else if (mod==2){
          
##           ts.trans <-1/ts-1
##           return(dlnorm(ts.trans,meanlog=meanlog,sdlog=sdlog)*(1/ts^2))
          
##         }else if (mod == 3){
          
##           c1 <- momentsGL(EG=exp(meanlog+sdlog^2/2),
##                           VG=(exp(sdlog^2)-1)*exp(2*meanlog+sdlog^2),
##                           beta=beta0)
          
##           outputMod3 <- list(infoPara = list(mu_beta=c(beta0,beta1),
##                                Sigma_beta=diag(c(0.05,0.05)),a_qG=100,r_qG=1),
##                              beta0=beta0,beta1=beta1,K=1,c1=c1,meanlog=meanlog,sdlog=sdlog
##                              )
##         }
        
##       }else if (dist == "b"){
        

##         beta0 <- 1
##         beta1 <- 0.1
##         aG1 <- 3
##         rG1 <- 0.8
        
##         gsBASE <- rbeta(Npat,aG1,rG1)
##         hs <- rep(1,Npat)
##         if (mod==1){
##           return( 
##                  rbind(probs=probs,
##                        quantile=qbeta(probs,shape1=aG1,shape2=rG1)
##                        )
##                  )
##         }else if (mod==2){
          
##           return (dbeta(ts,shape1=aG1,shape2=rG1))
          
##         }else if (mod == 3){
          
##           c1 <- momentsBeta(aG1,rG1,beta0)
##           outputMod3 <- list(infoPara = list(mu_beta=c(beta0,beta1),
##                                Sigma_beta=diag(c(0.05,0.05)),a_qG=10,r_qG=1),
##                              beta0=beta0,beta1=beta1,K=1,c1=c1, aGs=c(aG1=aG1),
##                              rGs=c(rG1=rG1))
##         }

##       }else if (dist == "b2"){
        
##         ## mixture of two beta distributions
##         ## cluster 1: E(Y)=exp(beta0)*mu_G=exp(2)*((aG1+rG1-1)/(aG1-1)+1)=3.098636 ## at initial time
##         ## cluster 2: E(Y)=exp(beta0)*mu_G=exp(2)*((aG2+rG2-1)/(aG2-1)+1)=7.037196 ## at initial time
##         pi <- 0.3
##         temp <- rmultinom(1,Npat,c(pi,1-pi))
##         Npat_dist1 <- temp[1,1]
##         Npat_dist2 <- temp[2,1]
        
##         pi <- 0.3
##         beta0 <- 1.5
##         beta1 <- 0.1
        
##         aG1 <- 10
##         rG1 <- 10
##         aG2 <- 20
##         rG2 <- 1

##         ## generate the initial random effect values of everyone
##         hs <- sample(1:2,Npat,c(pi,1-pi),replace=TRUE)
##         Npat_dist1 <- sum(hs==1)
##         Npat_dist2 <- sum(hs==2)
        
##         gsBASE <- rep(NA,Npat)
##         gsBASE[hs==1] <- rbeta(Npat_dist1,aG1,rG1);
##         gsBASE[hs==2] <- rbeta(Npat_dist2,aG2,rG2);
##         if (mod==1){
##           F_inv_beta2 <- function(p,aG1,rG1,aG2,rG2,pi)
##             {
##               G = function (t) pi*pbeta(t,shape1=aG1,shape2=rG1)+(1-pi)*pbeta(t,shape1=aG2,shape2=rG2) - p
##               return(uniroot(G,c(0,1000))$root)
##             }
##           lq <- length(probs);
##           quant <- rep(NA,lq);
##           for (i in 1 : lq) quant[i] <- F_inv_beta2(p=probs[i],aG1,rG1,aG2,rG2,pi=pi)
##           return (rbind(probs=probs,quantile=quant))
##         }
##         else if (mod==2){

##           return ((pi*dbeta(ts,shape1=aG1,shape2=rG1)+(1-pi)*dbeta(ts,shape1=aG2,shape2=rG2)) )
##         }else if (mod == 3){

##           c1 <- momentsBeta(aG1,rG1,beta0)
##           c2 <- momentsBeta(aG2,rG2,beta0)
##           outputMod3 <-  list(infoPara = list(mu_beta=c(beta0,beta1),
##                                 Sigma_beta=diag(c(0.05,0.05)),a_qG=100,r_qG=1),
##                               beta0=beta0,beta1=beta1,K=2,c1=c1,c2=c2,
##                               aGs=c(aG1=aG1,aG2=aG2),
##                               rGs=c(rG1=rG1,rG2=rG2)
##                               )
##         }
##       }
##       else if (dist == "b3")
##         {
##           ## mixture of three beta distributions
          
##           beta0 <- 1
##           beta1 <- 0.1
##           p1 <- 0.1
##           p2 <- 0.2
##           aG1 <- 5
##           rG1 <- 17
##           aG2 <- 25
##           rG2 <- 19
##           aG3 <- 20
##           rG3 <- 1
          
##           ## need modif
##           hs <- sample(1:3,Npat,c(p1,p2,1-p1-p2),replace=TRUE)
##           Npat_dist1 <- sum(hs==1)
##           Npat_dist2 <- sum(hs==2)
##           Npat_dist3 <- sum(hs==3)
          
##           gsBASE <- rep(NA,Npat)
##           gsBASE[hs==1] <- rbeta(Npat_dist1,aG1,rG1);
##           gsBASE[hs==2] <- rbeta(Npat_dist2,aG2,rG2);
##           gsBASE[hs==3] <- rbeta(Npat_dist3,aG3,rG3);
##           ## 

          
##           if (mod==1)
##             {
##               F_inv_beta3 <- function(p,aG1,rG1,aG2,rG2,aG3,rG3,p1,p2)
##                 {
##                   G = function (t)
##                     p1*pbeta(t,shape1=aG1,shape2=rG1)+p2*pbeta(t,shape1=aG2,shape2=rG2)+(1-p1-p2)*pbeta(t,shape1=aG3,shape2=rG3)-p
##                   return(uniroot(G,c(0,1000))$root)
##                 }
##               lq <- length(probs);
##               quant <- rep(NA,lq);
##               for (i in 1 : lq) quant[i] <- F_inv_beta3(p=probs[i],aG1,rG1,aG2,rG2,aG3,rG3,
##                                                         p1=p1,
##                                                         p2=p2)
##               return (rbind(probs=probs,quantile=quant))
##             }else if (mod==2)
##               {
##                 return (p1*dbeta(ts,shape1=aG1,shape2=rG1)+p2*dbeta(ts,shape1=aG2,shape2=rG2)+(1-p1-p2)*dbeta(ts,shape1=aG3,shape2=rG3 ))
##               }else if (mod == 3){
##                 c1 <- momentsBeta(aG1,rG1,beta0)
##                 c2 <- momentsBeta(aG2,rG2,beta0)
##                 c3 <- momentsBeta(aG3,rG3,beta0)
##                 outputMod3 <-  list(infoPara = list(mu_beta=c(beta0,beta1),
##                                       Sigma_beta=diag(c(0.05,0.05)),a_qG=100,r_qG=1),
##                                     beta0=beta0,beta1=beta1,K=3,c1=c1,c2=c2,c3=c3,
##                                     aGs=c(aG1=aG1,aG2=aG2,aG3=aG3),
##                                     rGs=c(rG1=rG1,rG2=rG2,rG3=rG3)
##                                     )
##               }
##         }

##     ## generate samples for dist = "b","b3", "g", and "l" (except "b2")
##     Y <- NULL
##     timecov <- c(0,0,(1:(ni-2)))
##     for ( ipat in 1 : Npat)
##       {
##         ## the number of repeated measures are the same
##         ## we assume that the time effects occurs once after the treatments are in effect
##         got <- rnbinom(ni,size = exp(beta0+beta1*timecov), prob = gsBASE[ipat])         
##         Y <- c(Y,got)
##       }
##     cat("beta0",beta0,"\n")


##     d <- data.frame(Y=Y,
##                     X1=rep(1,Npat*ni),
##                     X2=rep(timecov,Npat),
##                     ID=rep(1:Npat,each=ni),
##                     gsBASE = rep(gsBASE,each=ni),
##                     scan = rep(1:ni,Npat),
##                     ## day contains the day when the scan was taken
##                     ## 10 patients enter a trial every month
##                     days = days,
##                     hs = rep(hs,each=ni)
##                     )

##     ## dSMB visit is assumed to be every 4 months
##     if (full) return (d) 
##     d <- subset(d,subset= days <= DSMBVisitInterval*rev)
##     d$labelnp <- rep(0,nrow(d))
##     d$labelnp[ DSMBVisitInterval*(rev-1) < d$days ] <- 1
##     ## The first two scans (screening and base-line scans are treated as pre-scans)
##     d$labelnp[ d$X2==0 ] <- 0

##     if (dist=="b2")
##       {
##         temp <- index.b.each(Y=d$Y,ID=d$ID,labelnp=d$labelnp,X=cbind(d$X1,d$X2),
##                              betas=c(beta0,beta1),aGs=c(aG1,aG2),rGs=c(rG1,rG2),pis=c(pi,1-pi))
##         d$probIndex <- c(temp,rep(NA,nrow(d)-length(temp) ))
##       }else if (dist=="b")
##         {
##           temp <- index.b.each(Y=d$Y,ID=d$ID,labelnp=d$labelnp,X=cbind(d$X1,d$X2),
##                                betas=c(beta0,beta1),aGs=c(aG1),rGs=c(rG1),pis=c(1))
##           d$probIndex <- c(temp,rep(NA,nrow(d)-length(temp)) )
##         }else if (dist=="g"){

##           temp <- index.gammas(Y=d$Y,ID=d$ID,labelnp=d$labelnp,X=cbind(d$X1,d$X2),
##                               betas=c(beta0,beta1),
##                               shapes=c(sc1,sc2), ## shape
##                               scales=c(sp1,sp2), ## scale
##                               pis=c(pi,1-pi))
##           d$probIndex <- c(temp,rep(NA,nrow(d)-length(temp)) )
##         }

##     if (mod==0) return(d)

##     if (mod==3){
      
##       Npat <- length(unique(d$ID))
##       r_D <- 1
##       a_D <- 2.3 ## uniroot(f=K_NisK,
##       ##        interval=c(0.003,50),
##       ##        rD=r_D,N=100*2,
##       ##        K=10)$root
##       outputMod3$infoPara$a_D <- a_D
##       outputMod3$infoPara$r_D <- 1
##       outputMod3$infoPara$maxD <- 5
##       return(outputMod3)

##     }else if (mod == 4){
      
##       Npat <- length(unique(d$ID))
##       r_D <- 1
##       a_D <- 2.3 ## uniroot(f=K_NisK,
##       ##        interval=c(0.003,50),
##       ##        rD=r_D,N=100*2,
##       ##        K=10)$root
##       mu_beta <- c(0,0)
##       Sigma_beta <- diag(c(5^2,5^2))
##       a_qG <- 1
##       r_qG <- 1
      
##       outputMod4 <- list(
##                          mu_beta=mu_beta,Sigma_beta=Sigma_beta,
##                          a_qG=a_qG,r_qG=r_qG,
##                          a_D=a_D,r_D=r_D,
##                          M=M,maxD=5,max_aG = 30
##                          )
##       return(list(UinfoPara=outputMod4))
##     }
##   }


getTrueProb <- function(aG1=0.5,rG1=0.5,aG2=3,rG2=0.8,dist="gamma")
  {
    
    if (dist=="beta")
      {
        ## compute P(G1 < G2) = int_{y1=0}^1 (int_{y2=y1}^1 beta(y2;a2,b2) dy2 )*beta(y1;a1,b1) dy1
        insideIntBeta <- function(y1,a2,b2,a1,b1)
          integrate(dbeta,shape1=a2,shape2=b2,lower=y1, upper=1)$value*dbeta(y1,shape1=a1,shape2=b1)
        
        return (integrate(insideIntBeta,a2=aG2,b2=rG2,a1=aG1,b1=rG1,lower=0,upper=1)$value)
        
      }else if (dist=="gamma")
        {
          
          ## compute P(G1 < G2) = int_{y2=0}^Inf (int_{y1=0}^y1 gamma(y1;a1,b1) dy1 gamma(y2;a2,b2) dy2
          insideIntGam <- function(y2,a2,b2,a1,b1)
            {
              integrate(dgamma,shape=a1,scale=b1,lower=0, upper=y2)$value*dgamma(y2,shape=a2,scale=b2)
            }
          return  (integrate(insideIntGam,
                             a2=aG2,b2=rG2,a1=aG1,b1=rG1,lower=0,upper=Inf)$value)
        }
  }


mixgamma <- function(sp1,sc1,sp2,sc2,xs = seq(0,15,0.001),ylim=c(0,2),points=FALSE)
  {
    dg1 <- dgamma(xs,shape=sp1,scale=sc1)
    dg2 <- dgamma(xs,shape=sp2,scale=sc2)
    ys <- (1/2)*dg1+(1/2)*dg2
    if (!points){
      plot(xs,ys,ylim=ylim,type="l",
           main=paste("The mixture of gamma distributions:\n 1/2*dgamma(shape=",
             sp1,",scale=",sc1,")+1/2*dgamma(shape=",sp2,",scale=",sc2,")",sep=""))
      points(xs,dg1,type="l",col="red")
      points(xs,dg2,type="l",col="blue")
    }else{
      points(xs,ys,type="l",col="blue")
    }
  }


getM <- function(epsilonM,a_D,r_D,prob)
  {
    ## returns the appropriate truncation number M
    Dbig <- qgamma(prob,shape=a_D,scale=1/r_D)
    return(round(1+ log(epsilonM)/(log(Dbig/(1+Dbig))),0))
  }

Nuniq <-function(test,sampleindex){
  tt <- apply(test$h1[sampleindex,],1,unique)
  return( unlist(lapply(tt,length)))
}

colmeansd <- function (mat, name = NULL, sigDig = 2, sigDigSD = NULL) 
{
    if (is.vector(mat)) 
        mat <- matrix(mat, ncol = 1)
    if (length(sigDig) == 1) 
        sigDig <- rep(sigDig, ncol(mat))
    else if (length(sigDig) != ncol(mat)) 
        stop("length(sigDig) != ncol(mat)")
    
    if (length(sigDigSD) == 0)
      sigDigSD <- sigDig
    else if (length(sigDigSD) == 1) 
      sigDigSD <- rep(sigDigSD, ncol(mat))
    else if (length(sigDigSD) != ncol(mat)) 
      stop("length(sigDigSD) != ncol(mat)")
    
    mea <- colMeans(mat)
    sdd <- apply(mat, 2, sd)
    re <- NULL
    MysigDig <- paste("%1.", sigDig, "f", sep = "")
    MysigDigSD <- paste("%1.", sigDigSD, "f", sep = "")
    for (i in 1:ncol(mat)) {
        if (is.na(mea[i])) {
            re <- c(re, "------")
        }
        else {
            re <- c(re, paste(sprintf(MysigDig[i], mea[i]), " (", 
                sprintf(MysigDigSD[i], sdd[i], 3), ")", sep = ""))
        }
    }
    if (!is.null(name)) 
        names(re) <- name
    if (is.null(name) & !is.null(colnames(mat))) 
        names(re) <- colnames(mat)
    return(re)
}





## E(K|N,D) digamma(D+N)-digamma(D) ~= D*log(1+N/D) D*log(1+N/D)
K_N.D <- function(D,N=200) D*log(1+N/D)

## E(K|N) = E(E(K|N,D)) where D ~ gamma(shape=aD,scale1/rD)
## E(L|N) depends on aD and rD and N
K_N <- function(aD,rD,N)
  {
    f <- function(D,N,rD,aD) log(1+N/D)*exp(-rD*D)*D^aD
    k <- integrate(f=f,N=N,rD=rD,aD=aD,lower=0,upper=Inf)
    return(k$value*rD^aD/gamma(aD))
  }


K_NisK <- function(aD,rD,N,K)
  {
    return(K_N(aD,rD,N)-K)
  }


inpheatmap <- function(mat)
  {
    B = nrow(mat)
    mat <- .Call("map_c",
                 as.integer(c(mat)),
                 as.integer(ncol(mat)),
                 as.integer(B), package = "lmeNBBayes")/B
    diag(mat) = 1;
    return(mat);
  }



hm <- function(matBN,
               aveCEL,
               main="",
               minbin=0.5)
  {
    ## This function is designed to summarize a dataset such that
    ## the random effects of the first 100 pat are from dist1 and
    ## the random effects of second 100 pat are from dist2.

    ## Input is a B by N matrix, each row contains the cluster labels
    ## Based on this input, first, compute the similarity matrix:
    inphm <- inpheatmap(matBN)
    m.hmps <- 1 - inphm
    ord <- 1:ncol(matBN)

    ## reorder the patients within dist1/dist2 based on the dissimilarity measure 
    ## ord <- c(order.dendrogram(as.dendrogram(hclust(dist(m.hmps)))))

    spacedID <- rep(NA,ncol(matBN))
    for (i in 1 : ncol(matBN))
      {
        if (aveCEL[i]  < 0.5 )
          spacedID[i] <- ""
        if (aveCEL[i] >= 0.5 & aveCEL[i]  < 1 )
          spacedID[i] <- "-"
        else if (aveCEL[i] >= 1 & aveCEL[i]  < 2 )
          spacedID[i] <- "--"
        else if (aveCEL[i] >= 2 & aveCEL[i]  < 3 )
          spacedID[i] <- "---"
        else if (aveCEL[i] >= 3 & aveCEL[i]  < 4 )
          spacedID[i] <- "----"
        else if (aveCEL[i] >= 4 & aveCEL[i]  < 5 )
          spacedID[i] <- "-----"
        else if (aveCEL[i] >= 5 & aveCEL[i]  < 6 )
          spacedID[i] <- "------"
        else if (aveCEL[i] >= 6 & aveCEL[i]  < 7 )
          spacedID[i] <- "-------"
        else if (aveCEL[i] >= 7 & aveCEL[i]  < 8 )
          spacedID[i] <- "--------"
        else if (aveCEL[i] >= 8 & aveCEL[i]  < 9 )
          spacedID[i] <- "---------"
        else if (aveCEL[i] >= 9 & aveCEL[i]  < 10)
          spacedID[i] <- "----------"
        else if (aveCEL[i] >= 10)
          spacedID[i] <- "-----------"
      }
    rownames(inphm) <- colnames(inphm) <- spacedID
    
    ## Finally, plot the levelplot
    library(lattice)
    xmins <- c(0.01,0.5,0.01,0.5)
    xmaxs <- c(0.49,1,0.49,1)
    ymins <- c(0.5,0.5,0.01,0.01)
    ymaxs <- c(1,1,0.51,0.51)
    levelplot(inphm[ord,ord],
              labels = list(cex=0.01),
              at = do.breaks(c(minbin, 1.01), 20),
              scales = list(x = list(rot = 90)),
              aspect = "iso",
              colorkey=list(text=3,labels=list(cex=2)),
              xlab= "",ylab="",
              region = TRUE, 
              col.regions=gray.colors(100,start = 1, end = 0),
              main=paste(main,sep="")
              ##,par.settings=list(fontsize=list(text=15))
              )

    ## The posterior sample of label vector H1s contains the cluster labels.
    ## The distinct advantage of our DP mixture model is its ability to classify the patients into 
    ## the un-prespecified number of clusters.
    ## 
    ## Given B posterior samples of the label vectors of old scans H_1, for all combinations of the pairs of patients,
    ## the average times that paired patients belong to the same cluster are recorded to create single Npat by Npat
    ## similarity matrix. 
    ## The level map is created with this similarity matrix where 
    ## both row and column are reordered based on the result of hierarchical clustering on the similarity matrix (1-similarity matrix) 
    ##
  }

proG1LeG2 <- function(gsBbyN,g2sBbyN,beta=TRUE)
  {
    if (is.vector(gsBbyN))
      {
        gsBbyN <- matrix(gsBbyN,ncol=1)
        g2sBbyN <- matrix(g2sBbyN,ncol=1)
      }
    if(beta)
      {
        ## compare gPre >= gNew which is equivalent to 
        ## (1-gPre)/gPre <= (1-gNew)/gNew
        temp <- gsBbyN 
        gsBbyN <-g2sBbyN
        g2sBbyN <- temp
      }
    ##temp <- gsBbyN-g2sBbyN <= 0
    ## pG1G2 <- colMeans(temp)
    
    ##gsBbyN <- apply(gsBbyN,2,sort,decreasing=FALSE)
    ##g2sBbyN <- apply(g2sBbyN,2,sort,decreasing=FALSE)
    N <- ncol(gsBbyN)
    pG1G2 <- colMeans(gsBbyN < g2sBbyN)
    ## for (ipat in 1 : N)
    ##   {
    ##     g1 <- gsBbyN[,ipat]
    ##     g2 <- g2sBbyN[,ipat]
    ##     pG1G2[ipat] <-  .Call("pG1LeG2_c",
    ##                           as.numeric(g1),
    ##                           as.numeric(g2),
    ##                           package = "lmeNBBayes"
    ##                           )
    ##   }
    return(pG1G2)
  }

pointsgamma<-function(shape,scale,xmin=0,xmax=5,main="")
  {
    ## plotting the density of gamma dist'n with given shape and scale
    xs <-seq(xmin,xmax,length.out=1000)
    ys <- dgamma(xs,shape=shape,scale=scale)
    points(xs,ys,main=main,type="l",col="blue")
  }

plotgamma<-function(shape,scale,xmin=0,xmax=5,main="")
  {
    ## plotting the density of gamma dist'n with given shape and scale
    xs <-seq(xmin,xmax,length.out=1000)
    ys <- dgamma(xs,shape=shape,scale=scale)
    plot(xs,ys,main=main,type="l")
  }


plotnbinom <- function(size,prob,xmin=0,xmax=30,main="",size2=NULL,prob2=NULL,xlab="")
  {
    ## plotting the density of gamma dist'n with given size and prob
    xs <-xmin:xmax
    ys <- dnbinom(xs,size=size,prob=prob)
    plot(xs,ys,main=main,type="b",xlab=xlab)
    if (!is.null(size2))
      {
        ys <- dnbinom(xs,size=size2,prob=prob2)
        points(xs,ys,col="red",type="b")
      }
  }

plotHistAcf<- function(postsample,up.burnin=3000,mainplot,mainhist,xlim=NULL,col="black",everythin=1)
  {
    plot(postsample,type="l",main=mainplot,col=col)
    samp <- postsample[-(1:up.burnin)]
    samp <- samp[(1:(length(samp)/everythin))*everythin]
    med <- median(samp)
    if (is.null(xlim))
      {
        hist(samp,main=paste(mainhist," median:",round(med,3),sep=""))
      }else{
        hist(samp,main=paste(mainhist," median:",round(med,3),sep=""),
             xlim=xlim)##,breaks=0:1000)
      }
    ci <- round(quantile(samp,prob=c(0.025,0.975)),2)
    points(x=ci,y=c(0,0),type="b",col="blue",lwd=3)
    
    abline(v=med,col="blue")
    acf(samp,main="")
  }




plotGs <- function(vec,main,xlim=c(0,8),upto=10,length.out=1000,breaks=seq(0,15,0.001),ylim)
  {
    med <- round(median(vec),3)
    mea <- round(mean(vec),3)
    hist(vec,probability=TRUE,
         main=paste(main," mean: ",mea," median: ",med,sep=""),
         xlim=xlim,breaks=c(breaks,100000),ylim=ylim)
    ## density estimates are obtained at seq(0,upto,length.out=length.out)
    estden <- density(vec,kernel = "gaussian",from=0,to=upto,n=length.out)
    points(estden$x,estden$y,type="l")
  }



##=== Meta analysis function === ##
rigamma <- function (n, shape, scale) 
  {
    if (shape > 0 & scale > 0) 
      1/rgamma(n = n, shape=shape, scale=1/scale)
    else stop("rigamma: invalid parameters\n")
  }

metaUnivGibb <- function(ys,
### a vector of length Nstudy, containing the observed estimates
                         sigmas2,
### a vector of length Nstudy, containing the number of patients (or should it be the number of MS scan?) involved to estimate ys
                         a=1,
                         b=1,
                         B=100000)
  {
    Nstudy <- length(ys)
    
    ## Model:
    ## ys[istudy] ~ N(mean=theta[istudy],sd=sqrt(sigmas2[istudy]))
    ## theta_i ~ N(mu,tau2)
    ## mu ~ unif(-1000,1000)
    ## tau2 ~ inv-gamma(shape=a,scale=b)
    ## where i = 1,...,Nstudy and Nstudy is the number of dataset

    ## In our MS clinical analysis context,
    ## y_i = hat.beta_i (i=1,..,p) or hat.logaG_i
    
    ## informative prior construction
    ## Y ~ inv-gamma(shape=a,scale=b) then E(Y)=b/(a-1) for a > 1, Var(Y)=b^2/((a-1)^2(a-2)) for a > 2
    ## a <= 2 => infinite variance 
    postThetas <- matrix(NA,B,Nstudy)
    postMu <- rep(NA,B)
    postTau2 <- rep(NA,B)
    ## Initialization of unknown parameters
    thetas <- rep(0,Nstudy)
    mu <- mean(ys)
    tau2 <- mean(sigmas2)
    
    for (iB in 1 : B)
      {
        if (iB %% 5000 == 0) cat("\n",iB, "iterations are done" )
        ## Step 1: update theta_i
        for (istudy in 1 : Nstudy)
          {
            vari <- 1/(1/sigmas2[istudy]+1/tau2)
            thetas[istudy]<- rnorm(1,
                                   mean=(ys[istudy]/sigmas2[istudy] + mu/tau2)*vari,
                                   sd=sqrt(vari)
                                   )
          }
        
        mu <- rnorm(1,
                    mean=sum(thetas)/Nstudy,
                    sd=sqrt(tau2/Nstudy))
        
        tau2 <- rigamma(
                        1,
                        shape=(Nstudy/2+a),
                        scale=sum((thetas-mu)^2)/2+b
                        )
        
        postThetas[iB,] <- thetas
        postMu[iB] <- mu
        postTau2[iB] <- tau2
      }
    return(list(theta=postThetas,mu=postMu,tau2=postTau2))
  }



listSum<- function(listobj,burnin=0)
  {
    output <- listobj[[burnin+1]]
    for (i in (burnin+2) : length(listobj))
      {
        output <- output +listobj[[i]]
      }
    return (output)
  }


listMean <- function(listobj,useSample)
  {
    if (is.null(useSample))
      {
        useSample <- 1: length(listobj)
      }
    output <- listobj[[useSample[1]]]
    for (i in useSample[-1])
      {
        output <- output +listobj[[i]]
      }
    return (output/length(useSample))
  }



metaMultGibb <- function(hat.betas,
                         ## p by N.study matrix, where p = # covariates 
                         sigmas2,
                         ## a list of covariance matrix 
                         v0="uninfo",
                         S0=NULL,
                         B=10000)
  {
    ## TheWishart distribution is a multivariate analogue of the gamma distribution
    Nstudy <- ncol(hat.betas)
    p <- nrow(hat.betas)
    sigmas2Inv <- lapply(sigmas2,solve)
    ## Model:
    ## hat.betas[[istudy]] ~ mvnorm(mean=betas[,istudy], var=sigmas2[[istudy]])
    ## betas ~ mvnorm(mean=mu, var=tau2)
    ## :::: prior :::: 
    ## mu[i] ~ unif(-1000,1000), i = 1,...,p
    ## tau2 ~ inv-wishart(v0,S0) <=> solve(tau2) ~ wishart(v0,solve(S0))
    ## where i = 1,...,Nstudy and Nstudy is the number of dataset

    ## In our MS clinical analysis context,
    ## y_i = hat.beta_i (i=1,..,p) or hat.logaG_i


    
    ## The default values of hyperparameters
    ## inverse-wisrt is centered at diag(p) with large variance
    ## p109 Hogg
    if (v0=="info"){
      v0 <- 10
      S0 <- diag(p)*(v0-p-1) ## Identity matrix
      
    }else if (v0=="uninfo"){
      v0 <- p + 2
      S0 <- diag(c(0.5,0.10)) ## Identity matrix
    }

    
    ## storages
    postBetas <- postTau2 <- list()
    postMu <- matrix(NA,p,B)

    ## Initialization of unknown parameters
    betas <- matrix(0,nrow=p,ncol=Nstudy)
    mu <- rowMeans(hat.betas)
    tau2 <- listSum(sigmas2)/Nstudy
    tau2Inv <- solve( tau2 )
    
    for (iB in 1 : B)
      {
        if (iB %% 5000 == 0) cat("\n",iB, "iterations are done" )
        ## Step 1: update betas[,istudy]
        for (istudy in 1 : Nstudy)
          {
            vari <- solve(tau2Inv+sigmas2Inv[[istudy]] )
            betas[,istudy]<- mvrnorm(1,
                                     mu=vari%*%( tau2Inv%*%mu+sigmas2Inv[[istudy]]%*%hat.betas[,istudy]),
                                     Sigma=vari
                                     )
          }

        mu <- mvrnorm(1,
                      mu=rowSums(betas)/Nstudy, 
                      Sigma=tau2/Nstudy
                      )


        ## temp <- 0; for (istudy in 1 : Nstudy) temp <- temp + (betas[,istudy]-mu)%*%t(betas[,istudy]-mu)
        ## generate sample from Inverse-Wishart distribution
        Sn <- S0 + (betas-mu)%*%t(betas-mu) ## outer product sum_{istudy} (betas[,istudy]-mu)%*%t(betas[,istudy]-mu)
        tau2Inv <- matrix(rWishart(n=1,
                                   df=v0+Nstudy,
                                   Sigma=solve(Sn)
                                   ),p,p)

        postMu[,iB] <- mu
        postTau2[[iB]] <- tau2 <- solve(tau2Inv)
        postBetas[[iB]] <- betas
      }
    return(list(betas=postBetas,mu=postMu,tau2=postTau2))
  }











## nbinREDPmix <- function(Y,   ##     A vector of length sum ni, containing responses 
##                         X,   ##     A sum ni by p matrix, containing covariate values. The frist column must be 1 (Intercept)
##                         ID,  ##     A Vector of length sum ni, indicating patients
##                         B = 20000,   ##     A scalar, the number of Gibbs iteration 
##                         M = NULL,    ##     A scalar, the truncation value
##                         labelnp,
##                         dist = 1,      ##     dist == 1 means the base distribution of the random effects are assumed to be 
##                         max_aG = 3.0,    ##     from P0 = gamma(shape=aG,scale=1/rG)
##                         a_r = 3.0,       ##      where aG ~ Unif(0,max_aG ) and 
##                         r_r = 0.2,    ##     rG ~ gamma(shape=a_r,scale=1/r_r)
##                         a_D = 0.5,    ##     small shape a.D and small scale ib.D put weights on both small and large values of D
##                         r_D = 2,      ##      D ~ gamma(shape=a.D,scale=1/r.D)
##                         a_qG = 20.0,     ##     qG ~ beta(shape1=a_qG,scale=r_qG)
##                         r_qG = 1.0,
##                         mu_beta = 0,   ##     hyperpara of beta: beta ~ norm(mu_beta,sd=sigma_beta)
##                         sigma_beta = 5,##     hyperpara of beta: beta ~ norm(mu_beta,sd=sigma_beta) 
##                         burnin = round(B/3),
##                         epsilonM = 0.001,
##                         prob = 0.999,
##                         printfFreq = B
##                         )
##   {

##     ## ==== Description of the parameters in code ====
##     ## g1s:            A vector of length N, containing the random effect 1 values for each pat
##     ## g2s:            A vector of length N, containing the random effect 2 values for each pat
##     ## vs:             A vector of length M, containing the parameters to construct pis
##     ## weightH1(vs):   A vector of length M, containing the probabilities of the categorical distribution of Hi, function of vs
##     ## beta:           A vector of length p, containing the coefficients
##     ## h1s:            A vector of length N, containing the cluster labels of RE1 of each observation, cluster index is between 1 and M+1
##     ## h2s:            A vector of length N, containing the cluster labels of RE2 of each observation, cluster index is between 1 and M+1
##     ## aGs:            A vector of length M, containing the shape parameters of the gamma distribution
##     ## rGs:            A vector of length M, containing the scale parameters of the gamma distribution

##     ## ## === Model ===
##     ## // 1) y_ij | Gij = gij,beta, ~  NB(r_ij,p_i)
##     ## //    where r_ij = exp(t(x_ij)*beta), p_{ij} = 1/(g_ij+1)
##     ## //    gij  = g1i   if j is in pre scan
##     ## //         = g2i   if j is in new scan
##     ## //    and the distribution of g2i is specified as;
##     ## //    g2i  = g1i   with prob qG2
##     ## //         = g*    with prob 1-qG2   (The default value of p.g2 = 0)
##     ## //
##     ## // 2) Gi|ih ~ gamma(shape=aGs[ih],scale=1/rGs[ih])
##     ## //  aGs[ih],ih=1,...,M ~ Unif(0,max_aG )
##     ## //  rGs[ih],ih=1,...,M ~ gamma(shape=a_r,scale=1/r_r)
##     ## //
##     ## // 3) beta[i] ~ normal(0,sigma_beta) i=1,..,p
##     ## // 4) D ~ gamma(shape=a.D,scale=b.D)
##     if (is.null(M)) M = getM(epsilonM=epsilonM,a_D=a_D,r_D=r_D,prob=prob)

##     temID <- ID  
##     N <- length(unique(temID))
##     uniID <- unique(temID)
##     ID <- rep(NA,length(temID))
##     for (i in 1 : length(uniID))
##       {
##         ID[temID == uniID[i]] <- i
##       }
##     ## The patients with labelnp = 1 (no old scans) for all repeated measures
##     ## are treated as lack of new scans and used to estimate beta only.
##     ## skip the computation of H2
##     ## All patients have old scans

##     ## the labelnp of patients only with 1 (new scans) labels are replaced by all 0 (old scans)
##     patwonew <- which(as.numeric(tapply((labelnp==0),ID,sum)==0)==1)
##     for (i in 1 : length(patwonew)) labelnp[ID == patwonew[i]] <- 0

##     patwoONscan <-  which(as.numeric(tapply((labelnp==1),ID,sum)==0)==1) 
##     if (length(patwoONscan)==0) patwoONscan <- -999;
##     p <- ncol(X)
##     X <- c(X) ## {xij} = { x_{1,1},x_{2,1},..,x_{Ntot,1},x_{1,2},....,x_{Ntot,p} }
##     mID <- ID-1
##     mpatwoONscan <- patwoONscan - 1
##     maxni <- max(tapply(rep(1,length(ID)),ID,sum))
##     Npat <- length(unique(ID))

##     re <- .Call("res",
##                 as.numeric(Y),           ## REAL
##                 as.numeric(X),           ## REAL
##                 as.integer(mID),         ## INTEGER
##                 as.integer(B),           ## INTEGER
##                 as.integer(maxni),       ## INTEGER
##                 as.integer(M),           ## INTEGER
##                 as.integer(Npat),        ## INTEGER
##                 as.numeric(labelnp),     ## REAL
##                 as.numeric(max_aG),      ## REAL
##                 as.numeric(a_r),         ## REAL
##                 as.numeric(r_r),         ## REAL
##                 as.numeric(a_D),         ## REAL
##                 as.numeric(r_D),         ## REAL
##                 as.numeric(a_qG),        ## REAL 
##                 as.numeric(r_qG),        ## REAL
##                 as.numeric(mu_beta),     ## REAL
##                 as.numeric(sigma_beta),  ## REAL
##                 as.integer(mpatwoONscan),## INTEGER
##                 as.integer(dist),        ## INTEGER
##                 as.integer(burnin),      ## INTEGER
##                 as.integer(printfFreq),
##                 package = "lmeNBBayes"
##                 )

##     for ( i in 1 : 5) re[[i]] <- matrix(re[[i]],B,Npat,byrow=TRUE)

##     for ( i in 6 : 9) re[[i]] <- matrix(re[[i]],B,M,byrow=TRUE)
##     re[[10]] <- matrix(re[[10]],B,p,byrow=TRUE)

##     names(re) <- c("h1","h2","g1","g2","j",
##                    "aG","rG","v","weightH1",
##                    "beta","D","qG","AR","prp")

##     para <- list(mu_beta = mu_beta,
##                  sigma_beta = sigma_beta,
##                  M = M,
##                  a_D = a_D,
##                  r_D = r_D,
##                  dist = dist,
##                  burnin = burnin,
##                  a_r = a_r,
##                  r_r = r_r,
##                  a_qG = a_qG,
##                  r_qG = r_qG,
##                  max_aG = max_aG,
##                  patwoONscan=(patwoONscan+1))

##     re$para <- para
##     names(re$AR) <-c("gs", "aG", paste("beta",1:p))
##     re$hs <- re$h1 + 1  
##     re$h2s <- re$h2 + 1 
##     return(re)  
##   }



## nbinREDPmix2 <- function(Y,   ##     A vector of length sum ni, containing responses 
##                          X,   ##     A sum ni by p matrix, containing covariate values. The frist column must be 1 (Intercept)
##                          ID,  ##     A Vector of length sum ni, indicating patients
##                          B = 20000,   ##     A scalar, the number of Gibbs iteration 
##                          M = NULL,    ##     A scalar, the truncation value
##                          labelnp,
##                          dist = 1,    ##     dist == 1 means the base distribution of the random effects are assumed to be 
##                          max_aG = 3.0,    ##     from P0 = gamma(shape=aG,scale=1/rG)
##                          ##a_r = 3.0,       ##      where aG ~ Unif(0,max_aG ) and 
##                          ##r_r = 0.2,    ##     rG ~ gamma(shape=a_r,scale=1/r_r)
##                          a_D = 0.5,    ##     small shape a.D and small scale ib.D put weights on both small and large values of D
##                          r_D = 2,      ##      D ~ gamma(shape=a.D,scale=1/r.D)
##                          a_qG = 20.0,     ##     qG ~ beta(shape1=a_qG,scale=r_qG)
##                          r_qG = 1.0,
##                          mu_beta = 0,   ##     hyperpara of beta: beta ~ norm(mu_beta,sd=sigma_beta)
##                          sigma_beta = 5,##     hyperpara of beta: beta ~ norm(mu_beta,sd=sigma_beta) 
##                          burnin = round(B/3),
##                          epsilonM = 0.001,
##                          prob = 0.999,
##                          printfFreq = B
##                          )
##   {

##     ## ==== Description of the parameters in code ====
##     ## g1s:            A vector of length N, containing the random effect 1 values for each pat
##     ## g2s:            A vector of length N, containing the random effect 2 values for each pat
##     ## vs:             A vector of length M, containing the parameters to construct pis
##     ## weightH1(vs):   A vector of length M, containing the probabilities of the categorical distribution of Hi, function of vs
##     ## beta:           A vector of length p, containing the coefficients
##     ## h1s:            A vector of length N, containing the cluster labels of RE1 of each observation, cluster index is between 1 and M+1
##     ## h2s:            A vector of length N, containing the cluster labels of RE2 of each observation, cluster index is between 1 and M+1
##     ## aGs:            A vector of length M, containing the shape parameters of the gamma distribution
##     ## rGs:            A vector of length M, containing the scale parameters of the gamma distribution

##     ## ## === Model ===
##     ## // 1) y_ij | Gij = gij,beta, ~  NB(r_ij,p_i)
##     ## //    where r_ij = exp(t(x_ij)*beta), p_{ij} = 1/(g_ij+1)
##     ## //    gij  = g1i   if j is in pre scan
##     ## //         = g2i   if j is in new scan
##     ## //    and the distribution of g2i is specified as;
##     ## //    g2i  = g1i   with prob qG2
##     ## //         = g*    with prob 1-qG2   (The default value of p.g2 = 0)
##     ## //
##     ## // 2) Gi|ih ~ gamma(shape=aGs[ih],scale=1/rGs[ih])
##     ## //  aGs[ih],ih=1,...,M ~ Unif(0,max_aG )
##     ## //  rGs[ih],ih=1,...,M ~ gamma(shape=a_r,scale=1/r_r)
##     ## //
##     ## // 3) beta[i] ~ normal(0,sigma_beta) i=1,..,p
##     ## // 4) D ~ gamma(shape=a.D,scale=b.D)
##     if (is.null(M)) M = getM(epsilonM=epsilonM,a_D=a_D,r_D=r_D,prob=prob)

##     temID <- ID  
##     N <- length(unique(temID))
##     uniID <- unique(temID)
##     ID <- rep(NA,length(temID))
##     for (i in 1 : length(uniID))
##       {
##         ID[temID == uniID[i]] <- i
##       }
##     ## The patients with labelnp = 1 (no old scans) for all repeated measures
##     ## are treated as lack of new scans and used to estimate beta only.
##     ## skip the computation of H2
##     ## All patients have old scans

##     ## the labelnp of patients only with 1 (new scans) labels are replaced by all 0 (old scans)
##     patwonew <- which(as.numeric(tapply((labelnp==0),ID,sum)==0)==1)
##     for (i in 1 : length(patwonew)) labelnp[ID == patwonew[i]] <- 0

##     patwoONscan <-  which(as.numeric(tapply((labelnp==1),ID,sum)==0)==1) 
##     if (length(patwoONscan)==0) patwoONscan <- -999;
##     p <- ncol(X)
##     X <- c(X) ## {xij} = { x_{1,1},x_{2,1},..,x_{Ntot,1},x_{1,2},....,x_{Ntot,p} }
##     mID <- ID-1
##     mpatwoONscan <- patwoONscan - 1
##     maxni <- max(tapply(rep(1,length(ID)),ID,sum))
##     Npat <- length(unique(ID))

##     re <- .Call("resBeta",
##                 as.numeric(Y),           ## REAL
##                 as.numeric(X),           ## REAL
##                 as.integer(mID),         ## INTEGER
##                 as.integer(B),           ## INTEGER
##                 as.integer(maxni),       ## INTEGER
##                 as.integer(M),           ## INTEGER
##                 as.integer(Npat),        ## INTEGER
##                 as.numeric(labelnp),     ## REAL
##                 as.numeric(max_aG),      ## REAL
##                 ##as.numeric(a_r),         ## REAL
##                 ##as.numeric(r_r),         ## REAL
##                 as.numeric(a_D),         ## REAL
##                 as.numeric(r_D),         ## REAL
##                 as.numeric(a_qG),        ## REAL 
##                 as.numeric(r_qG),        ## REAL
##                 as.numeric(mu_beta),     ## REAL
##                 as.numeric(sigma_beta),  ## REAL
##                 as.integer(mpatwoONscan),## INTEGER
##                 as.integer(dist),        ## INTEGER
##                 as.integer(burnin),      ## INTEGER
##                 as.integer(printfFreq),
##                 package = "lmeNBBayes"
##                 )

##     for ( i in 1 : 5) re[[i]] <- matrix(re[[i]],B,Npat,byrow=TRUE)

##     for ( i in 6 : 9) re[[i]] <- matrix(re[[i]],B,M,byrow=TRUE)
##     re[[10]] <- matrix(re[[10]],B,p,byrow=TRUE)

##     names(re) <- c("h1","h2","g1","g2","j",
##                    "aG","rG","v","weightH1",
##                    "beta","D","qG","AR","prp")

##     para <- list(mu_beta = mu_beta,
##                  sigma_beta = sigma_beta,
##                  M = M,
##                  a_D = a_D,
##                  r_D = r_D,
##                  dist = dist,
##                  burnin = burnin,
##                  a_qG = a_qG,
##                  r_qG = r_qG,
##                  max_aG = max_aG,
##                  patwoONscan=(patwoONscan+1))

##     re$para <- para
##     names(re$AR) <-c("gs", "aG", paste("beta",1:p))
##     re$hs <- re$h1 + 1  
##     re$h2s <- re$h2 + 1 
##     return(re)  
##   }





## getSampleOld <- function(
##                       iseed=911,
##                       rev=4,
##                       dist="b",
##                       mod=0,
##                       beta1=0,
##                       probs=seq(0,0.99,0.01),
##                       ts = seq(0,0.99,0.01)
##                       )
## {
##   ## mod = 0: generate sample
##   ## mod = 1: quantiles of the true populations
##   ## mod = 2: densities of the true populations
##   ni <- 11
##   Npat <- 200;
##   ni_rev <- c(5,8,11)[rev-1]
##   set.seed(iseed)
##   ## generate unobservable latent variables gs:
##   if (dist=="g")
##     {
##       ## r.e. mixture of two gamma distributions
##       ## cluster1: E(Y)=exp(0.5)*E(G) = exp(0.5)*sp1*sc1 = 0.4121803
##       ## cluster2: E(Y)=exp(0.5)*E(G) = exp(0.5)*sp2*sc2 = 3.956931
##       pi <- 0.7
##       temp <- rmultinom(1,Npat,c(pi,1-pi))
##       Npat_dist1 <- temp[1,1]
##       Npat_dist2 <- temp[2,1]
##       sp1 <- 0.5
##       sc1 <- 0.5
##       sp2 <- 3
##       sc2 <- 0.8
##       beta0 <- 0.5
##       sampldist1 <- rgamma(Npat_dist1,shape=sp1,scale=sc1)
##       sampldist2 <- rgamma(Npat_dist2,shape=sp2,scale=sc2)
##       gsBASE <- c(sampldist1,sampldist2)
##       gsBASE = 1/(1+gsBASE)
##       if (mod==1)
##         {
##           F_inv_gam <- function(p,sp1,sc1,sp2,sc2,pi)
##             {
##               G = function (t)
##                 pi*pgamma(t,shape=sp1,scale=sc1)+(1-pi)*pgamma(t,shape=sp2,scale=sc2) - p
##               return(uniroot(G,c(0,100))$root)
##             }
##           lq <- length(probs);
##           quan <- rep(NA,lq);
##           for (i in 1 : lq) quan[i] <- F_inv_gam(p=probs[i],sp1=sp1,sc1=sc1,sp2=sp2,sc2=sc2,pi=pi)
##           quan = sort(1/(1+quan))
##           return (rbind(probs=probs,quantile=quan))
##         }else if (mod == 2){

##             ts.trans <-1/ts-1
##             return ((pi*dgamma(ts.trans,shape=sp1,scale=sc1)+(1-pi)*dgamma(ts.trans,shape=sp2,scale=sc2) )*(1/ts^2))
##           }else if (mod == 3){
##             return (list(beta0=beta0,K=2))
##           }
##     }else if (dist== "l")
##       {
##         ## single log-normal distribution
##         ## E(Y)=exp(beta)*E(G)=exp(1)*exp(meanlog+sdlog^2/2)= 5.078419
##         beta0 <- 1
##         meanlog <- -2.5
##         sdlog <- 2.5
##         gsBASE <- rlnorm(Npat,meanlog=meanlog,sdlog=sdlog)
##         gsBASE = 1/(1+gsBASE)
##         quan <- qlnorm(probs,meanlog=meanlog,sdlog=sdlog)
##         if (mod==1){
##           return( 
##                  rbind(probs=probs,
##                        quantile=sort(1/(1+quan))
##                        )
##                  )
##         }else if (mod==2){

##             ts.trans <-1/ts-1
##             return(dlnorm(ts.trans,meanlog=meanlog,sdlog=sdlog)*(1/ts^2))
##           }else if (mod == 3){
##             return (list(beta0=beta0,K=1))
##           }
##       }
##     else if (dist == "b")
##       {
##         ## single beta distribution
##         ## E(1/G) = ...some calculations... = (aG+rG-1)/(aG-1)
##         ## E(Y)=exp(beta0)*mu_{1/G}=exp(1)*((aG1+rG1-1)/(aG1-1)+1)= 4.238999
##         beta0 <- 1
##         aG1 <- 15.3
##         rG1 <- 8
##         gsBASE <- rbeta(Npat,aG1,rG1)
##         if (mod==1){
##           return( 
##                  rbind(probs=probs,
##                        quantile=qbeta(probs,shape1=aG1,shape2=rG1)
##                        )
##                  )
##         }else if (mod==2){

##             return (dbeta(ts,shape1=aG1,shape2=rG1))

##           }else if (mod == 3){

##             return (list(beta0=beta0,K=1))

##           }
##     }else if (dist == "b2")
##       {
##         ## mixture of two beta distributions
##         ## cluster 1: E(Y)=exp(beta0)*mu_G=exp(2)*((aG1+rG1-1)/(aG1-1)+1)=3.098636 ## at initial time
##         ## cluster 2: E(Y)=exp(beta0)*mu_G=exp(2)*((aG2+rG2-1)/(aG2-1)+1)=7.037196 ## at initial time
##         pi <- 0.3
##         temp <- rmultinom(1,Npat,c(pi,1-pi))
##         Npat_dist1 <- temp[1,1]
##         Npat_dist2 <- temp[2,1]
##         beta0 <- 2
##         aG1 <- 13
##         rG1 <- 18
##         aG2 <- 20
##         rG2 <- 1
##         sampldist1 <- rbeta(Npat_dist1,aG1,rG1);
##         sampldist2 <- rbeta(Npat_dist2,aG2,rG2);
##         gsBASE <- c(sampldist1,sampldist2)
##         if (mod==1){
##           F_inv_beta2 <- function(p,aG1,rG1,aG2,rG2,pi)
##             {
##               G = function (t) pi*pbeta(t,shape1=aG1,shape2=rG1)+(1-pi)*pbeta(t,shape1=aG2,shape2=rG2) - p
##               return(uniroot(G,c(0,1000))$root)
##             }
##           lq <- length(probs);
##           quant <- rep(NA,lq);
##           for (i in 1 : lq) quant[i] <- F_inv_beta2(p=probs[i],aG1,rG1,aG2,rG2,pi=pi)
##           return (rbind(probs=probs,quantile=quant))
##         }
##         else if (mod==2){

##           return ((pi*dbeta(ts,shape1=aG1,shape2=rG1)+(1-pi)*dbeta(ts,shape1=aG2,shape2=rG2)) )
##         }else if (mod == 3){
##           return(list(beta0=beta0,K=2))
##         }
##       }
##   else if (dist == "b3")
##     {
##       ## mixture of two beta distributions
##       ## cluster 1: E(Y)=exp(beta)*mu_G=exp(0.5)*((aG1+rG1-1)/(aG1-1)+1)=11.33496
##       ## cluster 2: E(Y)=exp(beta)*mu_G=exp(0.5)*((aG2+rG2-1)/(aG2-1)+1)=5.035284
##       ## cluster 3: E(Y)=exp(beta)*mu_G=exp(0.5)*((aG3+rG3-1)/(aG3-1)+1)=2.000417
##       beta0 <- 1
##       p1 <- 0.4
##       p2 <- 0.1
##       temp <- rmultinom(1,200,c(p1,p2,1-p1-p2))
##       Npat_dist1 <- temp[1,1]
##       Npat_dist2 <- temp[2,1]
##       Npat_dist3 <- temp[3,1]
##       aG1 <- 5
##       rG1 <- 19.5
##       aG2 <- 19.5
##       rG2 <- 19.5
##       aG3 <- 25
##       rG3 <- 0.01
##       sampldist1 <- rbeta(Npat_dist1,aG1,rG1);
##       sampldist2 <- rbeta(Npat_dist2,aG2,rG2);
##       sampldist3 <- rbeta(Npat_dist3,aG3,rG3);
##       gsBASE <- c(sampldist1,sampldist2,sampldist3)
##       if (mod==1)
##         {
##           F_inv_beta3 <- function(p,aG1,rG1,aG2,rG2,aG3,rG3,p1,p2)
##             {
##               G = function (t)
##                 p1*pbeta(t,shape1=aG1,shape2=rG1)+p2*pbeta(t,shape1=aG2,shape2=rG2)+(1-p1-p2)*pbeta(t,shape1=aG3,shape2=rG3)-p
##               return(uniroot(G,c(0,1000))$root)
##             }
##           lq <- length(probs);
##           quant <- rep(NA,lq);
##         for (i in 1 : lq) quant[i] <- F_inv_beta3(p=probs[i],aG1,rG1,aG2,rG2,aG3,rG3,
##                                                   p1=p1,
##                                                   p2=p2)
##           return (rbind(probs=probs,quantile=quant))
##         }else if (mod==2)
##           {
##             return (p1*dbeta(ts,shape1=aG1,shape2=rG1)+p2*dbeta(ts,shape1=aG2,shape2=rG2)+(1-p1-p2)*dbeta(ts,shape1=aG3,shape2=rG3 ))
##           }else if (mod == 3){
##             return( list(beta0=beta0,K=3))
##           }
##     }

##   Y <- NULL
##   for ( ipat in 1 : Npat)
##     {
##       ## the number of repeated measures are the same
##       got <- rnbinom(ni,size = exp(beta0+beta1*(1:ni)), prob = gsBASE[ipat])         
##       Y <- c(Y,got)
##     }
##   cat("beta0",beta0,"\n")

##   d <- data.frame(Y=Y,
##                   X=rep(1,Npat*ni),
##                   ID=rep(1:Npat,each=ni),
##                   gsBASE1 = rep(gsBASE,each=ni),
##                   scan = rep(1:ni,Npat)
##                   )

##   d <- subset(d,subset=scan <= ni_rev)
##   d$labelnp <- rep(c(rep(0,ni_rev-3),rep(1,3)),Npat)

##   return(d)
## }




## DPfitold <- function(Y,   ##   A vector of length sum ni, containing responses 
##                   X,   ##   A sum ni by p matrix, containing covariate values. The frist column must be 1 (Intercept)
##                   ID,  ##   A Vector of length sum ni, indicating patients
##                   B = 105000,     ##     A scalar, the number of Gibbs iteration 
##                   max_aG = 20.0,  ##     from P0 = gamma(shape=aG,scale=1/rG)
##                   mu_beta = NULL,    ##     hyperpara of beta: beta[ibeta] ~ norm(mu_beta[ibeta],sd=sigma_beta[ibeta])
##                   Sigma_beta = NULL, ## Beta24
##                   sigma_beta = NULL, ##    
##                   burnin = 5000,  
##                   printFreq = B,
##                   model=8,
##                   M=NULL,
##                   a_D=1,
##                   r_D=1,
##                   epsilonM=0.01,  ## nonpara
##                   prob=0.9,       ## nonpara
##                   labelnp=NULL, 
##                   sd_eta=0.3,      ## Beta13 and Beta23
##                   initial_beta=0.0,## Beta23
##                   a_qG=1,          ## Beta14 and Beta24
##                   r_qG=1,          ## Beta14 and Beta24
##                   mu_aG = 1.5,        ## Beta24
##                   sd_aG = 0.75,      ## Beta24
##                   mu_rG = 1.5,        ## Beta24
##                   sd_rG = 0.75      ## Beta24
##                   )
##   {
##     Ntot <- length(Y)

##     if (is.vector(X)) X <- matrix(X,ncol=1)
##     pCov <- ncol(X)
##     if (nrow(X)!= Ntot) stop ("nrow(X) ! = length(Y)")
##     if (length(ID)!= Ntot)  stop ("length(ID)! = length(Y)")

##     if (is.null(sigma_beta)) sigma_beta <- rep(5,pCov)
##     if (is.null(mu_beta)) mu_beta <- rep(0,pCov)
##     if (is.null(labelnp)) labelnp <- rep(0,length(Y))

##     ## if (is.null(M)) M = getM(epsilonM=epsilonM,a_D=a_D,r_D=r_D,prob=prob)
##     ## change the index of ID to numeric from 1 to # patients
##     temID <- ID  
##     N <- length(unique(temID))
##     uniID <- unique(temID)
##     ID <- rep(NA,length(temID))
##     for (i in 1 : length(uniID))
##       {
##         ID[temID == uniID[i]] <- i
##       }

##     p <- pCov
##     if( is.vector(X)) p <- 1 
##     X <- c(X) ## {xij} = { x_{1,1},x_{2,1},..,x_{Ntot,1},x_{1,2},....,x_{Ntot,p} }
##     mID <- ID-1

##     ## The patients with labelnp = 1 (no old scans) for all repeated measures
##     ## are treated as lack of new scans and used to estimate beta only.
##     ## skip the computation of H2
##     ## All patients have old scans

##     ## the labelnp of patients only with 1 (new scans) labels are replaced by all 0 (old scans)
##     patwonew <- which(as.numeric(tapply((labelnp==0),ID,sum)==0)==1)
##     for (i in 1 : length(patwonew)) labelnp[ID == patwonew[i]] <- 0

##     patwoNS <-  which(as.numeric(tapply((labelnp==1),ID,sum)==0)==1)
##     if (length(patwoNS)==0) patwoNS <- -999;
##     patwoNS <- patwoNS - 1

##     maxni <- max(tapply(rep(1,length(ID)),ID,sum))
##     Npat <- length(unique(ID))
##     if (model==1)
##       {
##         ##parametric model: gi is constant over time 
##         labelnp <- rep(0,length(Y))    
##         re <- .Call("Beta1",
##                     as.numeric(Y),           ## REAL
##                     as.numeric(X),           ## REAL
##                     as.integer(mID),         ## INTEGER
##                     as.integer(B),           ## INTEGER
##                     as.integer(maxni),       ## INTEGER
##                     as.integer(Npat),        ## INTEGER
##                     as.numeric(labelnp),     ## REAL
##                     as.numeric(mu_aG),
##                     as.numeric(sd_aG),
##                     as.numeric(mu_rG),
##                     as.numeric(sd_rG),
##                     as.numeric(mu_beta),     ## REAL
##                     as.numeric(sigma_beta),  ## REAL
##                     as.integer(burnin),      ## INTEGER
##                     as.integer(printFreq),
##                     package = "lmeNBBayes"
##                     )

##         re[[3]] <- matrix(re[[3]],B,N,byrow=TRUE)
##         re[[4]] <- matrix(re[[4]],B,p,byrow=TRUE)
##         names(re) <- c("aG","rG","gPre","beta","AR","prp")
##         para <- list(mu_beta = mu_beta,
##                      sigma_beta = sigma_beta,
##                      burnin = burnin,
##                      max_aG = max_aG
##                      )
##         re$para <- para
##         names(re$AR) <-c("aG", "rG", paste("beta",1:p))
##         return(re)
##       }
##     else if (model==2)
##       {
##         ## parametric model: gi is allowed to change between new scan g1 and old scan g2 period
##         ##                   g1 and g2 are assumed to be independent
##         re <- .Call("Beta2",
##                     as.numeric(Y),           ## REAL
##                     as.numeric(X),           ## REAL
##                     as.integer(mID),         ## INTEGER
##                     as.integer(B),           ## INTEGER
##                     as.integer(maxni),       ## INTEGER
##                     as.integer(Npat),        ## INTEGER
##                     as.numeric(labelnp),     ## REAL
##                     as.numeric(max_aG),      ## REAL
##                     as.numeric(mu_beta),     ## REAL
##                     as.numeric(sigma_beta),  ## REAL
##                     as.integer(burnin),      ## INTEGER
##                     as.integer(printFreq),
##                     package = "lmeNBBayes"
##                     )

##         for ( i in 3:4 ) re[[i]] <- matrix(re[[i]],B,N,byrow=TRUE)
##         re[[5]] <- matrix(re[[5]],B,p,byrow=TRUE)
##         names(re) <- c("aG","rG","g1s","g2s","beta","AR","prp")
##         para <- list(mu_beta = mu_beta,
##                      sigma_beta = sigma_beta,
##                      burnin = burnin,
##                      max_aG = max_aG
##                      )
##         re$para <- para
##         names(re$AR) <-c("aG", "rG", paste("beta",1:p))
##         return(re)
##       }
##     else if (model==3)
##       {
##         ## Nonparametric model: gi is constant over time 
##         labelnp <- rep(0,length(Y))
##         if (is.null(M)) M  <- getM(epsilonM=epsilonM,a_D=a_D,r_D=r_D,prob=prob)
##         re <- .Call("Beta12",
##                     as.numeric(Y),           ## REAL
##                     as.numeric(X),           ## REAL
##                     as.integer(mID),         ## INTEGER
##                     as.integer(B),           ## INTEGER
##                     as.integer(maxni),       ## INTEGER
##                     as.integer(Npat),        ## INTEGER
##                     as.numeric(labelnp),     ## REAL
##                     as.numeric(max_aG),      ## REAL
##                     as.numeric(mu_beta),     ## REAL
##                     as.numeric(sigma_beta),  ## REAL
##                     as.numeric(a_D),  ## REAL
##                     as.numeric(r_D),  ## REAL
##                     as.integer(M),
##                     as.integer(burnin),      ## INTEGER
##                     as.integer(printFreq),
##                     package = "lmeNBBayes"
##                     )
##         ## http://stackoverflow.com/questions/8720550/how-to-return-array-of-structs-from-call-to-c-shared-library-in-r
##         for ( i in 1:4 ) re[[i]] <- matrix(re[[i]],B,M,byrow=TRUE)
##         for ( i in 5:6 ) re[[i]] <- matrix(re[[i]],B,N,byrow=TRUE)
##         re[[7]] <- matrix(re[[7]],B,p,byrow=TRUE)
##         names(re) <- c("aG","rG","vs","weightH1",
##                        "h1s","g1s",
##                        "beta","D","AR","prp")
##         para <- list(mu_beta = mu_beta,
##                      sigma_beta = sigma_beta,
##                      burnin = burnin,
##                      max_aG = max_aG
##                      )
##         re$para <- para
##         names(re$AR) <-c("aG", "rG",paste("beta",1:p))
##         return(re)
##       }
##     else if (model==4)
##       {
##         ## nonparametric model: gi is allowed to change between new scan g1 and old scan g2 period
##         ##                      g1 and g2 are assumed to be independent 
##         if (is.null(M)) M  <- getM(epsilonM=epsilonM,a_D=a_D,r_D=r_D,prob=prob)
##         re <- .Call("Beta22",
##                     as.numeric(Y),           ## REAL
##                     as.numeric(X),           ## REAL
##                     as.integer(mID),         ## INTEGER
##                     as.integer(B),           ## INTEGER
##                     as.integer(maxni),       ## INTEGER
##                     as.integer(Npat),        ## INTEGER
##                     as.numeric(labelnp),     ## REAL
##                     as.numeric(max_aG),      ## REAL
##                     as.numeric(mu_beta),     ## REAL
##                     as.numeric(sigma_beta),  ## REAL
##                     as.numeric(a_D),         ## REAL
##                     as.numeric(r_D),         ## REAL
##                     as.integer(M),           ## INTEGER
##                     as.integer(burnin),      ## INTEGER
##                     as.integer(printFreq),
##                     package = "lmeNBBayes"
##                     )
##         ## http://stackoverflow.com/questions/8720550/how-to-return-array-of-structs-from-call-to-c-shared-library-in-r
##         for ( i in 1:4 ) re[[i]] <- matrix(re[[i]],B,M,byrow=TRUE)
##         for ( i in 5:8 ) re[[i]] <- matrix(re[[i]],B,N,byrow=TRUE)
##         re[[9]] <- matrix(re[[9]],B,p,byrow=TRUE)
##         names(re) <- c("aG","rG","vs","weightH1",
##                        "h1s","h2s","g1s","g2s",
##                        "beta","D","AR","prp")
##         para <- list(mu_beta = mu_beta,
##                      sigma_beta = sigma_beta,
##                      burnin = burnin,
##                      max_aG = max_aG,
##                      M=M,
##                      )
##         re$para <- para
##         names(re$AR) <-c("aG", "rG", paste("beta",1:p))
##         return(re)
##       }
##     else if (model == 5)
##       {
##         min_prec_eta = 500
##         max_prec_eta = 501
##         re <- .Call("Beta13",
##                     as.numeric(Y),           ## 1REAL
##                     as.numeric(X),           ## 2REAL
##                     as.integer(mID),         ## 3INTEGER
##                     as.integer(B),           ## 4INTEGER
##                     as.integer(maxni),       ## 5INTEGER
##                     as.integer(Npat),        ## 6INTEGER
##                     as.numeric(labelnp),     ## 7REAL
##                     as.numeric(min_prec_eta),   ## 8REAL
##                     as.numeric(max_prec_eta),   ## 9REAL
##                     as.numeric(max_aG),      ## 10REAL
##                     as.numeric(mu_beta),     ## 11REAL
##                     as.numeric(sigma_beta),  ## 12REAL
##                     as.integer(burnin),      ## 13INTEGER
##                     as.integer(printFreq),   ## 14INTEGER
##                     package = "lmeNBBayes"
##                     )
##         ## http://stackoverflow.com/questions/8720550/how-to-return-array-of-structs-from-call-to-c-shared-library-in-r
##         ## for ( i in 1:4 ) re[[i]] <- matrix(re[[i]],B,M,byrow=TRUE)
##         for ( i in 4 : 6 ) re[[i]] <- matrix(re[[i]],B,N,byrow=TRUE)
##         re[[7]] <- matrix(re[[7]],B,p,byrow=TRUE)
##         names(re) <- c("aG","rG","prec_eta",
##                        "g1s","etas","g2s",
##                        "beta","AR","prp")
##         para <- list(mu_beta = mu_beta,
##                      sigma_beta = sigma_beta,
##                      burnin = burnin,
##                      max_aG = max_aG
##                      )
##         re$para <- para
##         names(re$AR) <-c("g1s","etas","aG", "rG", "prec_eta",paste("beta",1:p))
##         return(re)

##       }
##     else if (model==6)
##       {
##         ##  nonparametric model,  
##         if (is.null(M)) M  <- getM(epsilonM=epsilonM,a_D=a_D,r_D=r_D,prob=prob)

##         re <- .Call("Beta23",
##                     as.numeric(Y),           ## REAL
##                     as.numeric(X),           ## REAL
##                     as.integer(mID),         ## INTEGER
##                     as.integer(B),           ## INTEGER
##                     as.integer(maxni),       ## INTEGER
##                     as.integer(Npat),        ## INTEGER
##                     as.integer(M),           ## INTEGER
##                     as.numeric(a_D),         ## REAL
##                     as.numeric(r_D),         ## REAL
##                     as.numeric(labelnp),     ## REAL
##                     as.numeric(sd_eta),      ## REAL
##                     as.numeric(max_aG),      ## REAL
##                     as.numeric(mu_beta),     ## REAL
##                     as.numeric(sigma_beta),  ## REAL
##                     as.integer(burnin),      ## INTEGER
##                     as.integer(printFreq),   ## INTEGER
##                     as.numeric(initial_beta),
##                     package = "lmeNBBayes"
##                     )
##         ## http://stackoverflow.com/questions/8720550/how-to-return-array-of-structs-from-call-to-c-shared-library-in-r
##         for ( i in 1 : 4 ) re[[i]] <- matrix(re[[i]],B,M,byrow=TRUE)
##         for ( i in 5 : 8 ) re[[i]] <- matrix(re[[i]],B,N,byrow=TRUE)
##         re[[5]] <- re[[5]]+1 ## "h1s"
##         re[[9]] <- matrix(re[[9]],B,p,byrow=TRUE)
##         names(re) <- c("aGs","rGs","vs","weightH1",
##                        "h1s","g1s","g2s","etas"
##                        ,"betas","D",##"prec_eta",
##                        "AR","prp")
##         para <- list(mu_beta = mu_beta,
##                      sigma_beta = sigma_beta,
##                      burnin = burnin,
##                      max_aG = max_aG,
##                      M=M,
##                      ID=ID,
##                      labelnp=labelnp,
##                      patWOnew=patwonew
##                      ) 
##         re$para <- para
##         names(re$AR) <- names(re$prp) <- c("g1s","etas",paste("aG",1:M,sep=""),
##                                            paste("rG",1:M,sep=""),
##                                            paste("beta",1:p,sep=""))##"prec_eta")
##         re$MeanWH <- colMeans(re$weightH1);
##         names(re$MeanWH) <- paste("cluster",1:M)
##         return(re)                
##       }
##     else if (model==7)
##       {

##         ## // Y_ij | Gij = gij ~ NB(size=exp(X_{ij}^T beta),prob=gij)
##         ## // gij = g1 if j is in pre-scan 
##         ## //     = g_new if j is in old-scan 
##         ## // g_new = Ji * g1 + (1-Ji) * g2 
##         ## // Ji ~ ber(qG)
##         ## // qG ~ beta(a_qG,r_qG)
##         ## // g1, g2 ~ beta(aG,rG)
##         ## // beta ~ rnorm(mu_beta,sigma_beta)
##         ## // aG, rG ~ unif(0,max_aG)

##         re <- .Call("Beta14",
##                     as.numeric(Y),           ## REAL
##                     as.numeric(X),           ## REAL
##                     as.integer(mID),         ## INTEGER
##                     as.integer(B),           ## INTEGER
##                     as.integer(maxni),       ## INTEGER
##                     as.integer(Npat),        ## INTEGER
##                     as.numeric(labelnp),     ## REAL
##                     as.numeric(a_qG),
##                     as.numeric(r_qG),
##                     as.numeric(max_aG),      ## REAL
##                     as.numeric(mu_beta),     ## REAL
##                     as.numeric(sigma_beta),  ## REAL
##                     ##as.numeric(a_D),         ## REAL
##                     ##as.numeric(r_D),         ## REAL
##                     ##as.integer(M),           ## INTEGER
##                     as.integer(burnin),      ## INTEGER
##                     as.integer(printFreq),
##                     as.integer(patwoNS),
##                     package = "lmeNBBayes"
##                     )
##         ## http://stackoverflow.com/questions/8720550/how-to-return-array-of-structs-from-call-to-c-shared-library-in-r
##         ## for ( i in 1:4 ) re[[i]] <- matrix(re[[i]],B,M,byrow=TRUE)
##         for ( i in 4 : 7 ) re[[i]] <- matrix(re[[i]],B,N,byrow=TRUE)
##         re[[8]] <- matrix(re[[8]],B,p,byrow=TRUE)
##         names(re) <- c("aG","rG","qG",
##                        "gPre","g2s","gNew","js",
##                        "beta","AR","prp")
##         para <- list(mu_beta = mu_beta,
##                      sigma_beta = sigma_beta,
##                      burnin = burnin,
##                      B=B,
##                      max_aG = max_aG,
##                      patwoNS= patwoNS+1
##                      )
##         re$para <- para
##         names(re$AR) <-c("aG", "rG", paste("beta",1:p))
##         return(re)
##       }
##     else if (model==8)
##       {
##         ## nonparametric model
##         ## // Y_ij | Gij = gij ~ NB(size=exp(X_{ij}^T beta),prob=gij)
##         ## // gij = g1 if j is in pre-scan 
##         ## //     = g_new if j is in old-scan 
##         ## // g_new = Ji * g1 + (1-Ji) * g2 
##         ## // Ji ~ ber(qG)
##         ## // qG ~ beta(a_qG,r_qG)
##         ## // g1, g2 ~ sum_{h=1}^M pi_h beta(aG_h,rG_h)
##         ## // beta ~ rnorm(mu_beta,sigma_beta)
##         ## // aG, rG ~ unif(0,max_aG)
##         ##  nonparametric model,

##         if (is.null(M)) M  <- getM(epsilonM=epsilonM,a_D=a_D,r_D=r_D,prob=prob)
##         if (is.null(Sigma_beta)) Sigma_beta <-  diag(5,pCov)
##         evalue_sigma_beta <- eigen(Sigma_beta, symmetric = TRUE, only.values = TRUE)$values
##         if (min(evalue_sigma_beta) <= 0) stop("Sigma_beta must be positive definite!")
##         Inv_sigma_beta <- c( solve(Sigma_beta) )

##         re <- .Call("Beta24",
##                     as.numeric(Y),           ## REAL
##                     as.numeric(X),           ## REAL
##                     as.integer(mID),         ## INTEGER
##                     as.integer(B),           ## INTEGER
##                     as.integer(maxni),       ## INTEGER
##                     as.integer(Npat),        ## INTEGER
##                     as.integer(M),           ## INTEGER
##                     as.numeric(labelnp),     ## REAL
##                     as.numeric(a_qG),        ## REAL
##                     as.numeric(r_qG),        ## REAL
##                     as.numeric(mu_aG),       ## REAL
##                     as.numeric(sd_aG),       ## REAL
##                     as.numeric(mu_rG),       ## REAL
##                     as.numeric(sd_rG),       ## REAL
##                     as.numeric(mu_beta),     ## REAL
##                     as.numeric(evalue_sigma_beta),  ## REAL
##                     as.numeric(Inv_sigma_beta),  ## REAL
##                     as.numeric(a_D),         ## REAL
##                     as.numeric(r_D),         ## REAL
##                     as.integer(burnin),      ## INTEGER
##                     as.integer(printFreq),
##                     as.integer(patwoNS),
##                     package = "lmeNBBayes"
##                     )
##         for ( i in 2 : 7 ) re[[i]] <- matrix(re[[i]],B,N,byrow=TRUE)
##         for ( i in 8 : 11 ) re[[i]] <- matrix(re[[i]],B,M,byrow=TRUE)
##         re[[12]] <- matrix(re[[12]],B,p,byrow=TRUE)
##         names(re) <- c("qG",
##                        "gPre","g2s","gNew","js","h1s","h2s",
##                        "weightH1","vs","aGs","rGs",
##                        "beta","AR","prp")
##         re$para <- list(mu_beta = mu_beta,
##                         sigma_beta = (Sigma_beta),
##                         burnin = burnin,
##                         M=M,
##                         B=B,
##                         model=model,
##                         mu_aG = mu_aG,
##                         sd_aG = sd_aG,
##                         mu_rG = mu_rG,
##                         sd_rG = sd_rG,
##                         a_D=a_D,
##                         r_D=r_D,
##                         a_qG=a_qG,
##                         r_qG=r_qG,
##                         Npat=N,
##                         Ntot=length(Y),
##                         patwoNS = patwoNS+1,
##                         CEL=Y,
##                         ID=ID
##                         )
##         names(re$AR) <- names(re$prp) <- c(paste("aG",1:M,sep=""),
##                                            paste("rG",1:M,sep=""),
##                                            paste("beta",sep="")
##                                            )
##         re$h1s <- re$h1s + 1
##         re$h2s <- re$h2s + 1
##         re$MeanWH <- colMeans(re$weightH1);
##         names(re$MeanWH) <- paste("cluster",1:M)
##         return(re)

##       }
##   }


## getSample2 <- function(iseed=1,
##                        mod=1, ## mod==2 returns the simulation parameter 
##                        aG1=10,rG1=10,aG2=20,rG2=1)
## {
##   ## The random effect changes its value over time
##   ## 20 % of patients switch their random effect values between pre and new scan period
##   ## By review 2, 80 patients are recruited to the study 
##   Npat <- 100; ## Npat must be divisible by NpatEnterPerMonth

##   beta0 <- 1.5
##   beta1 <- 0.1
##   set.seed(iseed)

##   cat("\n\n beta0",beta0,"beta1",beta1,"\n\n")
##   ## The number of patients entering the study per month
##   DSMBVisitInterval <- 4 ## months

##   ##
##   ## I_{ij} follows a Markov Chain with two states 0/1 with transition probability p_{i0} = 0.9 and p_{i1}=0.1 for i = 0 or 1
##   ## That is I_{ij} is more likely to stay at state 0.
##   ##
##   ## mixture of two beta distributions
##   ## cluster 1: E(Y)=exp(beta0)*mu_G=exp(1)*((aG1+rG1-1)/(aG1-1)+1)=3.098636 ## at initial time
##   ## cluster 2: E(Y)=exp(beta0)*mu_G=exp(1)*((aG2+rG2-1)/(aG2-1)+1)=7.037196 ## at initial time
##   Ys <- Gs <- Ispat <- ID <- Is <- NULL


##     cat("\n\n Random effect changes at time 2, at time 3 (DSMB review) and at time 4 \n\n")

##     Ys <- Gs <-Is <- NULL
##     ni <- 6
##     timecov <- c(0,0,1:(ni-2))
##     NpatEach <- 5
##     NpatEach1 <- 2
##     Ispat <- matrix(c(
##                       rep(c(0,0,0,1,1,1),NpatEach),
##                       rep(c(1,1,1,0,0,0),NpatEach),
##                       rep(c(0,0,1,1,1,1),NpatEach1), ## NpatEach 
##                       rep(c(0,0,0,0,1,1),NpatEach1),
##                       rep(c(1,1,0,0,0,0),NpatEach1),
##                       rep(c(1,1,1,1,0,0),NpatEach1),
##                       rep(rep(1,ni),40), ## 40 
##                       rep(rep(0,ni),130) ## 130 reasonabl
##                       ),ncol=ni,byrow=TRUE)
##     Npat <-nrow(Ispat)
##     for (ipat in 1 : Npat)
##       {
##         gs <- ys <- NULL
##         ## random effect corresponding to Ispat[ipat,ivec]==0
##         g0 <- rbeta(1,shape1=aG2,shape2=rG2);
##         g1 <- rbeta(1,shape1=aG1,shape2=rG1);
##         for (ivec in 1 : ni)
##           {
##             if (Ispat[ipat,ivec] == 0)
##               {
##                 ys[ivec] <- rnbinom(1,size = exp(beta0+beta1*timecov[ivec]),prob = g0)
##                 gs <- c(gs,g0)
##               }else if (Ispat[ipat,ivec] == 1)
##                 {
##                   ys[ivec] <- rnbinom(1,size = exp(beta0+beta1*timecov[ivec]),prob = g1)
##                   gs <- c(gs,g1)
##                 }
##           }
##         Ys <- c(Ys,ys)
##         Gs <- c(Gs,gs)
##         Is <- c(Is,Ispat[ipat,])
##       }
##     d <- data.frame(Y=Ys,gs=Gs,Is=Is,
##                     ID=rep(1:Npat,each=ni),
##                     days=rep(1:ni,Npat),
##                     X1=rep(1,Npat*ni),
##                     X2=rep(timecov,Npat))

##     d$labelnp <- rep(0,nrow(d))
##     d$labelnp[ 3 < d$days ] <- 1 

##     if (mod== 1)return(d)
##     else if (mod == 2){
##       ## return the parameter values:
##       c1 <- momentsBeta(aG1=aG1,rG1=rG1,beta=beta0)
##       c2 <- momentsBeta(aG1=aG2,rG1=rG2,beta=beta0)
##       REclass <- REclass(d)
##       return (list(c1=c1,c2=c2,beta0=beta0,beta1=beta1,aG1=aG1,rG1=rG1,aG2=aG2,rG2=rG2,REclass=REclass))
##     }

## }



## REclass <- function(d)
##   {
##     d <- data.frame(d)
##     Npat <- length(unique(d$ID))
##     REclass <- rep(NA,Npat)
##     for (ipat in 1 : Npat)
##       {
##         Is_pat <- d$Is[d$ID==ipat]
##         if (identical(Is_pat,c(1,1,1,1,1,1))){REclass[ipat] <- 1}
##         if (identical(Is_pat,c(0,0,0,0,0,0))){REclass[ipat] <- 2}
##         if (identical(Is_pat,c(1,1,0,0,0,0))){REclass[ipat] <- 3}
##         if (identical(Is_pat,c(0,0,0,0,1,1))){REclass[ipat] <- 4}
##         if (identical(Is_pat,c(1,1,1,1,0,0))){REclass[ipat] <- 5}
##         if (identical(Is_pat,c(0,0,1,1,1,1))){REclass[ipat] <- 6}
##         if (identical(Is_pat,c(0,0,0,1,1,1))){REclass[ipat] <- 7}
##         if (identical(Is_pat,c(1,1,1,0,0,0))){REclass[ipat] <- 8}
##       }
##     return(REclass)
##   }


## getSample2 <- function(iseed=1,
##                        mod=1, ## mod==2 returns the simulation parameter 
##                        aG1=10,rG1=10,aG2=20,rG2=1)
## {
##   ## The random effect changes its value over time
##   ## 20n % of patients switch their random effect values between pre and new scan period
##   ## By review 2, 80 patients are recruited to the study 
##   Npat <- 100; ## Npat must be divisible by NpatEnterPerMonth

##   beta0 <- 1.5
##   beta1 <- 0.1
##   set.seed(iseed)

##   cat("\n\n beta0",beta0,"beta1",beta1,"\n\n")
##   ## The number of patients entering the study per month
##   DSMBVisitInterval <- 4 ## months

##   ##
##   ## I_{ij} follows a Markov Chain with two states 0/1 with transition probability p_{i0} = 0.9 and p_{i1}=0.1 for i = 0 or 1
##   ## That is I_{ij} is more likely to stay at state 0.
##   ##
##   ## mixture of two beta distributions
##   ## cluster 1: E(Y)=exp(beta0)*mu_G=exp(1)*((aG1+rG1-1)/(aG1-1)+1)=3.098636 ## at initial time
##   ## cluster 2: E(Y)=exp(beta0)*mu_G=exp(1)*((aG2+rG2-1)/(aG2-1)+1)=7.037196 ## at initial time
##   Ys <- Gs <- Ispat <- ID <- Is <- NULL


##     cat("\n\n Random effect changes at time 2, at time 3 (DSMB review) and at time 4 \n\n")

##     Ys <- Gs <-Is <- NULL
##     ni <- 6
##     timecov <- c(0,0,1:(ni-2))
##     NpatEach <- 5
##     NpatEach1 <- 2
##     Ispat <- matrix(c(
##                       rep(c(0,0,0,1,1,1),NpatEach),
##                       rep(c(1,1,1,0,0,0),NpatEach),
##                       rep(c(0,0,1,1,1,1),NpatEach1), ## NpatEach 
##                       rep(c(0,0,0,0,1,1),NpatEach1),
##                       rep(c(1,1,0,0,0,0),NpatEach1),
##                       rep(c(1,1,1,1,0,0),NpatEach1),
##                       rep(rep(1,ni),40), ## 40 
##                       rep(rep(0,ni),130) ## 130 reasonabl
##                       ),ncol=ni,byrow=TRUE)
##     Npat <-nrow(Ispat)
##     for (ipat in 1 : Npat)
##       {
##         gs <- ys <- NULL
##         ## random effect corresponding to Ispat[ipat,ivec]==0
##         g0 <- rbeta(1,shape1=aG2,shape2=rG2);
##         g1 <- rbeta(1,shape1=aG1,shape2=rG1);
##         for (ivec in 1 : ni)
##           {
##             if (Ispat[ipat,ivec] == 0)
##               {
##                 ys[ivec] <- rnbinom(1,size = exp(beta0+beta1*timecov[ivec]),prob = g0)
##                 gs <- c(gs,g0)
##               }else if (Ispat[ipat,ivec] == 1)
##                 {
##                   ys[ivec] <- rnbinom(1,size = exp(beta0+beta1*timecov[ivec]),prob = g1)
##                   gs <- c(gs,g1)
##                 }
##           }
##         Ys <- c(Ys,ys)
##         Gs <- c(Gs,gs)
##         Is <- c(Is,Ispat[ipat,])
##       }
##     d <- data.frame(Y=Ys,gs=Gs,Is=Is,
##                     ID=rep(1:Npat,each=ni),
##                     days=rep(1:ni,Npat),
##                     X1=rep(1,Npat*ni),
##                     X2=rep(timecov,Npat))

##     d$labelnp <- rep(0,nrow(d))
##     d$labelnp[ 3 < d$days ] <- 1 

##     if (mod== 1)return(d)
##     else if (mod == 2){
##       ## return the parameter values:
##       c1 <- momentsBeta(aG1=aG1,rG1=rG1,beta=beta0)
##       c2 <- momentsBeta(aG1=aG2,rG1=rG2,beta=beta0)
##       REclass <- REclass(d)
##       return (list(c1=c1,c2=c2,beta0=beta0,beta1=beta1,aG1=aG1,rG1=rG1,aG2=aG2,rG2=rG2,REclass=REclass))
##     }

## }






 mvrnorm <- function (n = 1, mu, Sigma, tol = 1e-06, empirical = FALSE) 
{
   ## borrowed from library MASS
    p <- length(mu)
    if (!all(dim(Sigma) == c(p, p))) 
        stop("incompatible arguments")
    eS <- eigen(Sigma, symmetric = TRUE, EISPACK = TRUE)
    ev <- eS$values
    if (!all(ev >= -tol * abs(ev[1L]))) 
        stop("'Sigma' is not positive definite")
    X <- matrix(rnorm(p * n), n)
    if (empirical) {
        X <- scale(X, TRUE, FALSE)
        X <- X %*% svd(X, nu = 0)$v
        X <- scale(X, FALSE, TRUE)
    }
    X <- drop(mu) + eS$vectors %*% diag(sqrt(pmax(ev, 0)), p) %*% 
        t(X)
    nm <- names(mu)
    if (is.null(nm) && !is.null(dn <- dimnames(Sigma))) 
        nm <- dn[[1L]]
    dimnames(X) <- list(nm, NULL)
    if (n == 1) 
        drop(X)
    else t(X)
}
