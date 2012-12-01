# include <R.h>
# include <Rinternals.h>
# include <Rmath.h>
# include <R_ext/Linpack.h>
   



void matProd(const double A[],const int Arow,const int Acol,const double B[],const int Bcol,double temp[])
{
  // This function compute Amat %*% Bmat and return its output as a vector temp 
  // tempMat is a Arow by Bcol matrix so temp is a vector of length Arow*Bcol
  // Input 
  // A: a vector of length Arow*Acol
  // The original matrix Amat (Arow by Acol) must be transformed to A as:
  // A = c(Amat[,1],Amat[,2],...) i.e., A=c(Amat) in R
  // Similarly B is a vector of length Acol*Bcol

  for (int irow=0;irow < Arow;irow++)
    {
      for (int icolB=0;icolB < Bcol; icolB++)
	{
	  temp[irow+Arow*icolB]  = 0.0;
	  for (int icol=0;icol < Acol;icol++)
	    {
	      
	      temp[irow+Arow*icolB] += A[irow+Arow*icol]*B[icol+Acol*icolB];
	    }
	}
    }
}

double dmvnorm(const double x[], const double mu[], const int p, const double evalue[],const double InvSigma[],const int logScale)
{
  // density of mvnorm at point x
  // the code is based on dmvnorm of library mvtnorm 
  // Input 
  // x: vector of length p
  // mu: vector of length p
  // Let sigma be the variance covariance matrix from betas
  // evalue = eigen(sigma, symmetric = TRUE, only.values = TRUE)$values (evalue is all greater than zero as sigma is positive definite )
  // InvSigma = solve(sigma)
  double vec[p],vecInvSigma[p],dens=0.0;
  for (int ip=0;ip <p;ip++) vec[ip] =x[ip] - mu[ip];
  
  // compute the maharanis distance: (x-mu)^T InvSigma (x-mu) 
  // distval <- mahalanobis(x, center = mean, cov = sigma)
  // logdet <- sum(log(eigen(sigma, symmetric = TRUE, only.values = TRUE)$values))
  matProd(vec,1,p,InvSigma,p,vecInvSigma);
  for (int ip = 0; ip < p ;ip++)
    {
      dens += vecInvSigma[ip]*vec[ip] + log(evalue[ip]);
    }
  // logretval <- -(ncol(x) * log(2 * pi) + logdet + distval)/2
  dens += p*log(2*M_PI);
  dens /= -2;
  // exp(logretval)
  if (logScale) return (dens);
  else return (exp(dens));
}



void mvrnorm(const double mu[], const int p, const double evalue[], const double evec[],double Sample[])
{
  // generate a single (vector) sample from MVRnorm
  // the code is based on mvrnorm of library MASS 
  // evalue = sqrt(pmax(ev, 0)) where ev is eigenvalue of the covariance matrix i.e.,
  // eS <- eigen(Sigma, symmetric = TRUE, EISPACK = TRUE)
  // ev <- eS$values;evec <- eS$vector

  // evec is a vector of length p*p, the first p entries form an eigen vector corresponding to evalue[0].
  //                                 the second p entries form an eigen vector corresponding to evalue[1] and so on...
  double rv[p];
  
  for (int irow=0;irow < p ; irow++) rv[irow] = norm_rand();
 
  //  X <- matrix(rnorm(p * n), n)
  //  X <- mu + evec %*% diag(sqrt(pmax(ev, 0)), p) %*% t(X)
  //  Note that diag(evalue, p) %*% t(X) = evalue*X if X is a vector

  //  That is, using R calculation notation,
  //  1: compute temp = evalue*X  in R^p (element wise calculation)
  //  2: compute Sample = evec%*%temp i.e. Sample[irow] =sum( evec[irow][-]* temp[irow] )
  // for (int i=0;i<(p*p);i++) Rprintf("\v evec[%d] %f",i,evec[i]);

  for (int irow=0; irow < p; irow++)
    {
      Sample[irow] = mu[irow];
      // 
      Rprintf("\v\v");
      for (int icol=0; icol < p; icol++)
	{

	  Sample[irow] += evec[irow+p*icol]*evalue[icol]*rv[icol];

	  //Rprintf("\v rv[%d] %f",icol,rv[icol]);
	  //Rprintf("\v evalue[%d] %f",icol,evalue[icol]);
	  //Rprintf("\v evec[%d] %f",irow+p*icol,evec[irow+p*icol]);
	  //Rprintf("\v Sample[%d] %f",irow,Sample[irow]);
	  //Rprintf("\v\v ");
	}
    }
  //Rprintf("\n\n");
}

// extern "C" {
//   SEXP tempFun(SEXP x_,SEXP mu_,SEXP p_, SEXP evalue_,SEXP InvSigma_);
//     }
SEXP tempFun(SEXP x_,SEXP mu_,SEXP p_, SEXP evalue_,SEXP InvSigma_)
{
  double *x = REAL(x_);
  double *mu = REAL(mu_);
  int p = INTEGER(p_)[0];
  double *evalue = REAL(evalue_);
  double *InvSigma = REAL(InvSigma_);

  GetRNGstate();
  
  double re = dmvnorm(x, mu, p, evalue, InvSigma,1);
  
  Rprintf("\v %f", re);
  PutRNGstate();
  
 return(R_NilValue);
}

double mixtured(const double pis[],const double ti, const double aGs[], const double rGs[],const int iB, const int M)
{
  // cdf of mixture 
  double temp = 0;
  int idh=0;
  for (int ih =0; ih < M ; ih++) 
    //rN[0]/L*dWEI2(ts[it],Klambdas[M*iB+0],betas[M*iB+0],0))+...+rN[M-1]/L*dWEI2(ts[it],Klambdas[M*iB+M-1],betas[M*iB+M-1],0))
    {
      idh = ih+M*iB;
      if (pis[idh]>0)
	temp += pis[idh]*dbeta(ti,aGs[idh],rGs[idh],
			       0 // log
			       );
    }
  return temp;
}


// // declarations (necessary for C++ codes)
// extern "C" {
//   SEXP mixdDist(SEXP pis_, SEXP aGs_, SEXP rGs_,SEXP alphas_,SEXP B_,SEXP M_);
// }



SEXP mixdDist(SEXP pis_, SEXP aGs_, SEXP rGs_,SEXP alphas_,SEXP B_,SEXP M_)
{
  // alphas: a vector of length lalphas containing values between 0 and 1 (support of beta density). the density is evaluated at each alpha
  // aGs: a vector of length B*M containing aGs 
  // rGs: a vector of length B*M containing rGs
  // pis: a vector of length B*M containing weightH1s
  const int Lalphas = length(alphas_);
  const double *pis = REAL(pis_);
  const double *aGs = REAL(aGs_);
  const double *rGs = REAL(rGs_);
  const double *alphas = REAL(alphas_);
  const int B = INTEGER(B_)[0];
  const int M = INTEGER(M_)[0];

  // http://stackoverflow.com/questions/8720550/how-to-return-array-of-structs-from-call-to-c-shared-library-in-r
  SEXP res = PROTECT(allocVector(VECSXP, 1)); // result is stored in res
  SEXP quant = allocVector(REALSXP, B*Lalphas); 
  SET_VECTOR_ELT(res, 0, quant); 

  GetRNGstate();

  for (int iB=0; iB < B; iB++)
    {
      for (int ialpha = 0 ; ialpha < Lalphas; ialpha++)
	{
	  R_CheckUserInterrupt(); 
	  //R_ProcessEvents();
	  REAL(quant)[Lalphas*iB+ialpha]=mixtured(pis,alphas[ialpha],aGs,rGs,iB, M);
	}
    }
    
  PutRNGstate();
  UNPROTECT(1);
  return res;
}




double mixtureP(const double pis[],const double ti, const double aGs[], const double rGs[],const int iB, const int M,
		const double alpha)
{
  // cdf of mixture 
  double temp = 0;
  int idh=0;
  for (int ih =0; ih < M ; ih++) 
    //rN[0]/L*dWEI2(ts[it],Klambdas[M*iB+0],betas[M*iB+0],0))+...+rN[M-1]/L*dWEI2(ts[it],Klambdas[M*iB+M-1],betas[M*iB+M-1],0))
    {
      idh = ih+M*iB;
      if (pis[idh]>0)
	temp += pis[idh]*pbeta(ti,aGs[idh],rGs[idh],
			       1, //lower tail
			       0 // log
			       );
    }
  return temp-alpha;
}


// // declarations (necessary for C++ codes)
// extern "C" {
//   SEXP mixQDist(SEXP pis_, SEXP aGs_, SEXP rGs_,SEXP alphas_,SEXP B_,SEXP M_);
// }

SEXP mixQDist(SEXP pis_, SEXP aGs_, SEXP rGs_,SEXP alphas_,SEXP B_,SEXP M_)
{
  // alphas: a vector of length lalphas containing values between 0 and 1. the alpha^th quantile is evaluated
  // aGs: a vector of length B*M containing aGs 
  // rGs: a vector of length B*M containing rGs
  // pis: a vector of length B*M containing weightH1s
  const int Lalphas = length(alphas_);
  const double *pis = REAL(pis_);
  const double *aGs = REAL(aGs_);
  const double *rGs = REAL(rGs_);
  const double *alphas = REAL(alphas_);
  const int B = INTEGER(B_)[0];
  const int M = INTEGER(M_)[0];
  double X,Y,X1,X2,Y1,Y2;
  int OK,idh;

  // http://stackoverflow.com/questions/8720550/how-to-return-array-of-structs-from-call-to-c-shared-library-in-r
  SEXP res = PROTECT(allocVector(VECSXP, 2)); // result is stored in res
  SEXP quant = allocVector(REALSXP, B*Lalphas); 
  SET_VECTOR_ELT(res, 0, quant); 
  SEXP NconvItr = allocVector(INTSXP, B*Lalphas); 
  SET_VECTOR_ELT(res, 1, NconvItr); 


  GetRNGstate();

  for (int iB=0; iB < B; iB++)
    {
      for (int ialpha = 0 ; ialpha < Lalphas; ialpha++)
	{
	  // bisection algorithm
	  X1=0.0001; // Left hand side initial value
	  X2=1.0;   // Right hand side initial value 
	  Y1=mixtureP(pis,X1,aGs,rGs,iB, M,alphas[ialpha]);
	  Y2=mixtureP(pis,X2,aGs,rGs,iB, M,alphas[ialpha]);
	  OK = (Y1 > 0) != (Y2 > 0);
	  if (OK) 
	    {
	      idh = 0;
	      while ((fabs (X1 - X2) > 1.0e-10) && (idh < 1000))
		{
		  //Rprintf("         Y1 %f Y2 %f X1 %f X2 %f ",Y1,Y2,X1,X2);
		  X = (X1 + X2) / 2.0;
		  Y = mixtureP(pis,X,aGs,rGs,iB, M,alphas[ialpha]);
		  if ((Y1 < 0) && (Y < 0)) {
		    X1 = X;
		    Y1 = Y;
		  } else {
		    X2 = X;
		    Y2 = Y;
		  }
		  ++idh;
		}
	    }
	  REAL(quant)[Lalphas*iB+ialpha]=X1;
	  INTEGER(NconvItr)[Lalphas*iB+ialpha]=idh;
	}
    }
    
  PutRNGstate();
  UNPROTECT(1);
  return res;
}








double FisherBeta(const double shape1, const double shape2,const int N)
{
  // shape1=aG, shape2=rG
  // This function returns:
  // asymptotic SD(MLE of shape1)
  // = 1/sqrt( - second derivative of the log likelihood of N beta densities  with respect to shape1 ) 
  // dd loglik(shape1,shape2)/(d shape1)^2 = - N dd ln B(shape1,shape2)/(d shape1)^2
  //                                       = - N ( trigamma(shape1)-trigamma(shape1+shape2) )
  double secDervBeta = trigamma(shape1)-trigamma(shape2+shape1);
  return ( 1/(sqrt(secDervBeta*N)) );
  // ::NOTE::
  // SD(MLE of shape2) = 1/(sqrt(N*( trigamma(shape1)-trigamma(shape1+shape2) )))
  // One can obtain the MLE of shape2 by replacing the input of shape1 and shape2
}

double adjustcan(const double w, double can, double minV, double maxV)
{
  // w is acceptance rate
  if ( minV < can && can < maxV)
    {
      if (w > 0.6)  can *=1.1;
      else if (w <= 0.2) can *=0.9;
      // no change between 0.3-0.5
    }
  if (can < minV) can = minV;
  else if (can > maxV) can = maxV;
  return can;
}

void getwID(int *wIDsize, 
	    int wID[], // wID is treated as pointer even w/o 
	    SEXP ID, 
	    const int ipat)
{
  // This function modifies the values of wIDsize and wID:
  
  // ID:      A vector of length sum ni, patient labeling starts from ZERO 
  // ipat:    A scalar the patient index that you are interested in (0 to N-1)
  *wIDsize = 0; 
  for (int ivec=0; ivec < length(ID) ;ivec++)
    {
      if( INTEGER(ID)[ivec]==ipat) 
	{
	  wID[*wIDsize] = ivec; 
	  (*wIDsize)++; 
	}
    }
  //Rprintf("\v ipat %d wIDsize %d  :",ipat,*wIDsize);
  //for (int ii = 0; ii < *wIDsize; ii++) Rprintf(" %d ",wID[ii]);
}


double getsizeij(
		 const int ivec, 
		 SEXP X,  // reference even w/o &
		 const int wID[],  // reference even w/o &
		 const double beta[],
		 const int Ntot, // = sum_i ni
		 const int p
		 )// reference even w/o &
{
  // This function modifies the values of sizeij to be sizeij = exp(beta0+Xi1*beta1+Xi2*beta2+...):

  // wID:     A vector of length max(ni), containing the position index to indicate the locations of repeated measures 
  //          of the selected patient in the vector ID in the first ni entries. The rest of max(ni) - ni + 1 entries are junks
  // ivec:    A scalar, indicating the ivec^th repeated measure

  double sizeij = 0.0;
  
  for (int ip = 0 ; ip < p ; ip ++ ) 
    sizeij +=  REAL(X)[ wID[ivec]+Ntot*ip ]*beta[ip];
  
  sizeij = exp(sizeij);

  return(sizeij);

  // sizeij = exp(beta0+Xi1*beta1+Xi2*beta2+...):
}


double getpij(const double g1, 
	      const double g2, 
	      const double labelnp)
{
  // returns pij
  return( (1.0/(g1+1.0))*(1.0-labelnp)+(1.0/(g2+1.0))*labelnp) ;
}


double getpij2(const double g1, 
	       const double g2, 
	       const double labelnp)
{
  // returns pij
  return( g1*(1.0-labelnp)+g2*labelnp ) ;
}


double getpijEta(const double g1, 
		 const double eta, 
		 const double labelnp)
{
  // returns pij
  // logit(g2)=logit(g1)+eta
  if (labelnp==0.0) return(g1);
  else if (labelnp==1.0) return( g1/( (1.0-g1)*exp(-eta)+g1 ) ) ;
  Rprintf("should not show this! getpijEta");
}



int sample(double prob[], const int size)
{
  // randomly select the number between 0 to size-1 
  double p_tot=0.0;
  int ih,rN[size];
  for (ih = 0; ih < size; ih ++) p_tot += prob[ih];
  if (p_tot <0) error("At least one prob must be positive!");
  for (ih = 0; ih < size; ih ++) prob[ih] /= p_tot;

  rmultinom(1, prob, size, rN);

  for ( ih=0; ih < size; ih++)
    {
      if(rN[ih]==1) return(ih);
    }
}


SEXP index_b_Bayes_UnexpIncrease(SEXP Y_, SEXP X, SEXP labelnp_, SEXP ID, SEXP B_,SEXP maxni, SEXP M_, SEXP nn, SEXP mat_betas_, 
				 SEXP mat_aGs_, SEXP mat_rGs_, SEXP mat_weightH1_,SEXP FreqEvery)
{

  const double *Y = REAL(Y_); 
  const double *labelnp = REAL(labelnp_);
  const int Ntot = length(ID); // # patients ## ID must start from zero
  const int p = length(X)/Ntot; 
  const int N = INTEGER(nn)[0];
  const int M = INTEGER(M_)[0];
  const int B = INTEGER(B_)[0];
  const double *mat_aGs = REAL(mat_aGs_); // The first M entries of the vectors corresponds to the iB=1 mat_aGs = c(t(mat_aGs))
  const double *mat_rGs = REAL(mat_rGs_); // The first M entries of the vectors corresponds to the iB=1
  const double *mat_weightH1 = REAL(mat_weightH1_); // The first M entries of the vectors corresponds to the iB=1
  const double *mat_betas = REAL(mat_betas_); // The first M entries of the vectors corresponds to the iB=1
  // Note: the vector structure of X is different from others:
  // The first B entries of the vector X must corresponds to the first covariate. X =c(X) 
  int wID[INTEGER(maxni)[0]],wIDsize = 0, ih = 0, ip = 0, ipat = 0, ivec=0, idh=0, k=0;
  // wIDs[ipat][ivec] contains the position of the ivec^th repeated measure of the ipat^th patient at labelnp and Y
  //                  If ivec >= # repeated measures, wIDs[ipat][ivec] = -1000
  // wIDsizes[ipat] contains the number of repeated measure of ipat^th patient
  // sumYs[ipat][0] contains the sum of the CEL counts of ipat^th patients in the pre-scan period
  // sumYs[ipat][1] contains the sum of the CEL counts of ipat^th patients in the new-scan period
  int wIDs[N][INTEGER(maxni)[0]], wIDsizes[N];
  double sumYs[N][2];
  for (ipat = 0 ; ipat < N ; ipat++)
    {
      getwID( &wIDsize,wID,ID,ipat);
      wIDsizes[ipat] = wIDsize;
      sumYs[ipat][0] = 0;
      sumYs[ipat][1] = 0;
      for (ivec = 0 ; ivec < wIDsize; ivec++)
	{
	  wIDs[ipat][ivec]=wID[ivec];
	  if (labelnp[wIDs[ipat][ivec]]==0) sumYs[ipat][0] += Y[wIDs[ipat][ivec]];
	  if (labelnp[wIDs[ipat][ivec]]==1) sumYs[ipat][1] += Y[wIDs[ipat][ivec]];
	}
      //Rprintf("\n ipat %d sumYs[ipat][0] %f sumYs[ipat][1] %f ",ipat,sumYs[ipat][0],sumYs[ipat][1]);
      for (ivec = wIDsize ; ivec < INTEGER(maxni)[0]; ivec++)
	wIDs[ipat][ivec]=-1000;
    }

  double postprob[M],betas[p];
  double sizePre=0.0,sizeNew=0.0,w=0.0,
    candi=0.0,yij=0.0,num=0.0,den=0.0;

  SEXP res = PROTECT(allocVector(VECSXP, 1)); // result is stored in res
  SEXP probI = allocVector(REALSXP, N*B);
  SET_VECTOR_ELT(res, 0, probI);

  for (int iB = 0; iB < B ; iB ++)
    {
      R_CheckUserInterrupt(); 
      for (ip = 0 ; ip < p ; ip++ )
	{
	  betas[ip] = mat_betas[iB*p+ip]; 
	  //Rprintf("iB %d betas[ip] %f ",iB, betas[ip]);
	}
      for (ih = 0 ; ih < M ; ih++) postprob[ih] = beta(mat_aGs[iB*M+ih],mat_rGs[iB*M+ih]);
  
      for (ipat = 0 ; ipat < N; ipat++ )
	{
	  idh = ipat + iB*N;
	  if (sumYs[ipat][1]==0) REAL(probI)[idh]=1.0;
	  else {
	    sizePre = 0;
	    sizeNew = 0;
	    for (ivec = 0 ; ivec < wIDsizes[ipat] ;ivec++ )
	      {
		if(labelnp[wIDs[ipat][ivec]]==0) sizePre+= getsizeij(ivec, X, wIDs[ipat], betas, Ntot,p);//sum_{j in pre} exp(Xij%*%beta)
		if(labelnp[wIDs[ipat][ivec]]==1) sizeNew+=getsizeij(ivec, X, wIDs[ipat], betas, Ntot,p);//sum_{j in new} exp(Xij%*%beta)
	      }
	    den = 0;
	    for (ih = 0 ; ih < M ; ih++)
	      {
		if ( mat_weightH1[iB*M+ih]!= 0)
		  {
		    den += mat_weightH1[iB*M+ih]*beta(sizePre+mat_aGs[iB*M+ih],sumYs[ipat][0]+mat_rGs[iB*M+ih])/postprob[ih];
		  }
	      }
	    num = 0;
	    for (k = 0 ; k < sumYs[ipat][1]; k++)
	      {
		w = 0;
		for (ih = 0 ; ih < M ; ih++)
		  {
		    if ( mat_weightH1[iB*M+ih] != 0)
		      {
			w += mat_weightH1[iB*M+ih]*beta(sizeNew+sizePre+mat_aGs[iB*M+ih],
							k+sumYs[ipat][0]+mat_rGs[iB*M+ih])/postprob[ih];
		      }
		  }
		num += choose(sizeNew+k-1, k)*w;
	      }
	    REAL(probI)[idh]=1-num/den;
	    //if (ipat == 0) Rprintf("ipat %d idh %d prob %f sumYs[ipat][1] %f num %f den %f w %f ",ipat, idh, REAL(probI)[idh],sumYs[ipat][1],num,den,w);
	  }
	}
      if(iB > 0 && iB % INTEGER(FreqEvery)[0] == 0) 
	{
	  Rprintf("\v %d iterations are done...   ", iB);
	  R_CheckUserInterrupt(); 
	  //R_ProcessEvents();
	}
    }
  
  PutRNGstate();
  UNPROTECT(1);
  return res;
}




SEXP index_b_Bayes_UnexpDecrease(SEXP Y_, SEXP X, SEXP labelnp_, SEXP ID, SEXP B_,SEXP maxni, SEXP M_, SEXP nn, SEXP mat_betas_, 
				 SEXP mat_aGs_, SEXP mat_rGs_, SEXP mat_weightH1_,SEXP FreqEvery)
{

  const double *Y = REAL(Y_); 
  const double *labelnp = REAL(labelnp_);
  const int Ntot = length(ID); // # patients ## ID must start from zero
  const int p = length(X)/Ntot; 
  const int N = INTEGER(nn)[0];
  const int M = INTEGER(M_)[0];
  const int B = INTEGER(B_)[0];
  const double *mat_aGs = REAL(mat_aGs_); // The first M entries of the vectors corresponds to the iB=1 mat_aGs = c(t(mat_aGs))
  const double *mat_rGs = REAL(mat_rGs_); // The first M entries of the vectors corresponds to the iB=1
  const double *mat_weightH1 = REAL(mat_weightH1_); // The first M entries of the vectors corresponds to the iB=1
  const double *mat_betas = REAL(mat_betas_); // The first M entries of the vectors corresponds to the iB=1
  // Note: the vector structure of X is different from others:
  // The first B entries of the vector X must corresponds to the first covariate. X =c(X) 
  int wID[INTEGER(maxni)[0]],wIDsize = 0, ih = 0, ip = 0, ipat = 0, ivec=0, idh=0, k=0;
  // wIDs[ipat][ivec] contains the position of the ivec^th repeated measure of the ipat^th patient at labelnp and Y
  //                  If ivec >= # repeated measures, wIDs[ipat][ivec] = -1000
  // wIDsizes[ipat] contains the number of repeated measure of ipat^th patient
  // sumYs[ipat][0] contains the sum of the CEL counts of ipat^th patients in the pre-scan period
  // sumYs[ipat][1] contains the sum of the CEL counts of ipat^th patients in the new-scan period
  int wIDs[N][INTEGER(maxni)[0]], wIDsizes[N];
  double sumYs[N][2];
  for (ipat = 0 ; ipat < N ; ipat++)
    {
      getwID( &wIDsize,wID,ID,ipat);
      wIDsizes[ipat] = wIDsize;
      sumYs[ipat][0] = 0;
      sumYs[ipat][1] = 0;
      for (ivec = 0 ; ivec < wIDsize; ivec++)
	{
	  wIDs[ipat][ivec]=wID[ivec];
	  if (labelnp[wIDs[ipat][ivec]]==0) sumYs[ipat][0] += Y[wIDs[ipat][ivec]];
	  if (labelnp[wIDs[ipat][ivec]]==1) sumYs[ipat][1] += Y[wIDs[ipat][ivec]];
	}
      //Rprintf("\n ipat %d sumYs[ipat][0] %f sumYs[ipat][1] %f ",ipat,sumYs[ipat][0],sumYs[ipat][1]);
      for (ivec = wIDsize ; ivec < INTEGER(maxni)[0]; ivec++)
	wIDs[ipat][ivec]=-1000;
    }

  double postprob[M],betas[p];
  double sizePre=0.0,sizeNew=0.0,w=0.0,
    candi=0.0,yij=0.0,num=0.0,den=0.0;

  SEXP res = PROTECT(allocVector(VECSXP, 1)); // result is stored in res
  SEXP probI = allocVector(REALSXP, N*B);
  SET_VECTOR_ELT(res, 0, probI);

  for (int iB = 0; iB < B ; iB ++)
    {
      for (ip = 0 ; ip < p ; ip++ ) betas[ip] = mat_betas[iB*p+ip]; 
      for (ih = 0 ; ih < M ; ih++) postprob[ih] = beta(mat_aGs[iB*M+ih],mat_rGs[iB*M+ih]);
  
      for (ipat = 0 ; ipat < N; ipat++ )
	{
	  idh = ipat + iB*N;
	  sizePre = 0;
	  sizeNew = 0;
	  for (ivec = 0 ; ivec < wIDsizes[ipat] ;ivec++ )
	    {
	      if(labelnp[wIDs[ipat][ivec]]==0) sizePre += getsizeij(ivec, X, wIDs[ipat], betas, Ntot,p);//sum_{j in pre} exp(Xij%*%beta)
	      if(labelnp[wIDs[ipat][ivec]]==1) sizeNew+=getsizeij(ivec, X, wIDs[ipat], betas, Ntot,p);//sum_{j in new} exp(Xij%*%beta)
	    }
	  den = 0;
	  for (ih = 0 ; ih < M ; ih++)
	    {
	      if ( mat_weightH1[iB*M+ih]!= 0)
		{
		  den += mat_weightH1[iB*M+ih]*beta(sizePre+mat_aGs[iB*M+ih],sumYs[ipat][0]+mat_rGs[iB*M+ih])/postprob[ih];
		}
	    }
	  num = 0;
	  for (k = 0 ; k <= sumYs[ipat][1]; k++)
	    {
	      w = 0;
	      for (ih = 0 ; ih < M ; ih++)
		{
		  if ( mat_weightH1[iB*M+ih] != 0)
		    {
		      w += mat_weightH1[iB*M+ih]*beta(sizeNew+sizePre+mat_aGs[iB*M+ih],
						      k+sumYs[ipat][0]+mat_rGs[iB*M+ih])/postprob[ih];
		    }
		}
	      num += choose(sizeNew+k-1, k)*w;
	    }
	  REAL(probI)[idh]= num/den;
   //if (ipat == 0) Rprintf("ipat %d idh %d prob %f sumYs[ipat][1] %f num %f den %f w %f ",ipat, idh, REAL(probI)[idh],sumYs[ipat][1],num,den,w);

	}
      if(iB > 0 && iB % INTEGER(FreqEvery)[0] == 0) 
	{
	  Rprintf("\v %d iterations are done...   ", iB);
	  R_CheckUserInterrupt(); 
	  //R_ProcessEvents();
	}
    }
  
  PutRNGstate();
  UNPROTECT(1);
  return res;
}


SEXP Beta1(SEXP Y_,          // REAL
	   SEXP X,           // REAL 
	   SEXP ID,          // INTEGER
	   SEXP B_,          // INTEGER
	   SEXP maxni,       // INTEGER
	   SEXP nn,          // INTEGER
	   SEXP labelnp_,    // REAL
	   SEXP max_aG_,
	   SEXP mu_beta_,    // REAL
	   SEXP evalue_sigma_betas_, // REAL
	   SEXP Inv_sigma_beta_,
	   SEXP burnin_,     // INTEGER
	   SEXP FreqEvery,    // INTEGER
	   SEXP probInd_, // Int
	   SEXP thin_    // Int
	   )
{
  // === Model ===
  // Y_ij | Gi = gi ~ NB(size=exp(X_{ij}^T beta),prob=gi)
  // gi ~ beta(aG,rG)
  // beta ~ rnorm(mu_beta,sigma_beta)
  // aG, rG ~ norm
  const double *Y = REAL(Y_); 
  const double *labelnp = REAL(labelnp_);
  const int Ntot = length(ID); // # patients
  const int p = length(X)/Ntot; 
  const int N = INTEGER(nn)[0];
  const double max_aG = REAL(max_aG_)[0];
  const int B = INTEGER(B_)[0];
  const double *mu_betas = REAL(mu_beta_);
  const double *evalue_sigma_betas =REAL(evalue_sigma_betas_);
  const double *Inv_sigma_beta = REAL(Inv_sigma_beta_);
  const int burnin = INTEGER(burnin_)[0];
  const double min_aG = 0.5;
  const int thin = INTEGER(thin_)[0];
  const int probInd = INTEGER(probInd_)[0];
  int wID[INTEGER(maxni)[0]],wIDsize = 0, ih = 0, ip = 0, ipat = 0, ivec=0, idh=0, ibeta=0,iB=0;
  // wIDs[ipat][ivec] contains the position of the ivec^th repeated measure of the ipat^th patient at labelnp and Y
  //                  If ivec >= # repeated measures, wIDs[ipat][ivec] = -1000
  // wIDsizes[ipat] contains the number of repeated measure of ipat^th patient
  // sumYs[ipat][0] contains the sum of the CEL counts of ipat^th patients in the pre-scan period
  // sumYs[ipat][1] contains the sum of the CEL counts of ipat^th patients in the new-scan period
  int wIDs[N][INTEGER(maxni)[0]], wIDsizes[N];
  double sumYs[N][2];

  for (ipat = 0 ; ipat < N ; ipat++)
    {
      getwID( &wIDsize,wID,ID,ipat);
      wIDsizes[ipat] = wIDsize;
      sumYs[ipat][0] = 0;
      sumYs[ipat][1] = 0;
      for (ivec = 0 ; ivec < wIDsize; ivec++)
	{
	  wIDs[ipat][ivec]=wID[ivec];
	  if (labelnp[wIDs[ipat][ivec]]==0) sumYs[ipat][0] += Y[wIDs[ipat][ivec]];
	  if (labelnp[wIDs[ipat][ivec]]==1) sumYs[ipat][1] += Y[wIDs[ipat][ivec]];
	}
      for (ivec = wIDsize ; ivec < INTEGER(maxni)[0]; ivec++)
	wIDs[ipat][ivec]=-1000;
    }

  // int h1s[N],h2s[N],js[N];
  double betas[p],canbetas[p],aG,rG,g1s[N];
  // double vs[M],weightH1[M],postprob[M],g2s[N];
  double sizeij=0.0,cansizeij=0.0,w=0.0,pij=0.0,logLcan=0,logL=0,
    candi=0.0,MHrate=0.0,yij=0.0,num=0.0,den=0.0;
  // double  weipro=0.0,qG2=0.0,D=0.0;
  double att[2+1], acc[2+1],can[2+1];
  for (ih = 0 ; ih < 2+1; ih ++)
    {
      att[ih] = 0.0;
      acc[ih] = 0.0;
      can[ih] = 2.0;
    }

  // http://stackoverflow.com/questions/8720550/how-to-return-array-of-structs-from-call-to-c-shared-library-in-r
  SEXP res = PROTECT(allocVector(VECSXP, 8)); // result is stored in res
  SEXP postaGs = allocVector(REALSXP, B); 
  SET_VECTOR_ELT(res, 0, postaGs); 
  SEXP postrGs = allocVector(REALSXP, B); 
  SET_VECTOR_ELT(res, 1, postrGs); 
  SEXP postg1s = allocVector(REALSXP, N*B); 
  SET_VECTOR_ELT(res, 2, postg1s); 
  SEXP postbetas = allocVector(REALSXP, p*B); 
  SET_VECTOR_ELT(res, 3, postbetas); 
  SEXP AR = allocVector(REALSXP,2+1);
  SET_VECTOR_ELT(res, 4, AR);
  SEXP prp = allocVector(REALSXP,2+1);
  SET_VECTOR_ELT(res, 5, prp);
  SEXP probI = allocVector(REALSXP, N*(B-burnin)/thin );
  SET_VECTOR_ELT(res, 6, probI);
  SEXP logLlik = allocVector(REALSXP, B );
  SET_VECTOR_ELT(res, 7, logLlik);
  GetRNGstate();
  
  // ============ initialization ============ 
  for (ip=0;ip < p ; ip++) betas[ip]=0.0;
  aG = runif(min_aG,max_aG);
  rG = runif(min_aG,max_aG);
  for (ipat = 0 ; ipat < N ; ipat++)
    g1s[ipat]=0.5;

  // ================ MCMC =================
  for (iB = 0 ; iB < B; iB++)
    {
      R_CheckUserInterrupt(); 
      //R_ProcessEvents();
      // STEP 1 
      // Update g1s, a vector of length M+1. The first M terms are the random effects.
      // the gs[M]=-1000, place holder for the patients with no new scans
      for (ipat = 0 ; ipat < N; ipat++)
	{
	  w = 0.0;
	  for (ivec=0; ivec < wIDsizes[ipat]; ivec++ )
	    {
	      w += getsizeij(ivec, X, wIDs[ipat], betas, Ntot, p);
	      // sum_{j=-1}^{n_i} exp(Xij%*%beta)
	    }
	  //Rprintf("\v ipat  %d wIDsizes[ipat] %d sumYs[0][ipat] %f sum exp(xij beta) %f",ipat, wIDsizes[ipat],sumYs[0][ipat], w);
	  g1s[ipat] = rbeta(aG+w,rG+sumYs[ipat][0]+sumYs[ipat][1]);
	  if (g1s[ipat] > 1.0-1.0e-7) g1s[ipat]=1.0-1.0e-7;
	  else if (g1s[ipat] < 1.0e-7) g1s[ipat]=1.0e-7;
	}
      
      // ================= 
      // STEP 2: update beta, a vector of length p, containing the coefficients b0, b1,...
      // beta ~ MVR(mu_beta,Sigma_beta); 

      att[2]++;
      // proposal distribution
      for (ibeta = 0; ibeta <p ; ibeta ++)
	canbetas[ibeta] = rnorm(betas[ibeta], can[2]); // the mean and sd must be double (not NumericVector)
      
      // prior
      MHrate = dmvnorm(canbetas, mu_betas,p, evalue_sigma_betas,Inv_sigma_beta,1) 
	- dmvnorm(betas, mu_betas,p, evalue_sigma_betas,Inv_sigma_beta,1);
      logLcan = 0.0;
      logL = 0.0;
      for (ipat=0; ipat < N; ipat++)
	{
	  // if (printing)  Rprintf(" \v\v ipat: %d",ipat);
	  for ( ivec = 0; ivec < wIDsizes[ipat]; ivec++ )
	    {
	      yij = Y[ wIDs[ipat][ivec] ];
	      // yij is a vector b/c it is the first arg of dnbinom
	      pij = g1s[ipat]; 
	      //Rprintf("\v\v ipat %d js[ipat] %f ivec %d yij %f g1s[ipat] %f g2s[ipat] %f pij %f", 
	      //ipat,  js[ipat] , ivec,   yij,   g1s[ipat], g2s[ipat],  pij);
	      cansizeij = getsizeij(ivec, X, wIDs[ipat], canbetas, Ntot, p); // the value of cansizeij is changed
	      sizeij = getsizeij(ivec, X, wIDs[ipat], betas, Ntot, p); // the value of cansizeij is changed

	      logLcan += dnbinom(yij, cansizeij,  pij, 1);
	      logL += dnbinom(yij, sizeij,  pij, 1);
	    }
	  //Rprintf("\v");
	}    
      MHrate += logLcan-logL;

      if(runif(0.0,1.0) < exp(MHrate))
	{
	  for (ibeta = 0; ibeta < p;ibeta++) betas[ibeta] = canbetas[ibeta];
	  acc[2]++;
	  logL = logLcan;
	}
      
      // Update aG and rG, scalars for the base distribution of the random effect gi
      // gi ~ Beta(aG,rG)
      // ================
      candi=rnorm(aG,can[0]);
      if (candi > min_aG && candi < max_aG )
	{
	  att[0]++;
	  MHrate = 0;
	  for (ipat=0 ; ipat < N; ipat++ ) 
	    MHrate += dbeta(g1s[ipat],candi,rG,1)-dbeta(g1s[ipat],aG,rG,1); 
	  if (runif(0.0,1.0) < exp(MHrate) ) 
	    {
	      aG = candi;
	      acc[0]++;
	    }
	}
      
      candi=rnorm(rG,can[1]);
      if (candi > min_aG && candi < max_aG )
	{
	  att[1]++;
	  MHrate = 0;
	  for (ipat=0 ; ipat < N; ipat++ ) 
	    MHrate += dbeta(g1s[ipat],aG,candi,1)-dbeta(g1s[ipat],aG,rG,1);   
	  if (runif(0.0,1.0) < exp(MHrate) ) 
	    {
	      rG = candi;
	      acc[1]++;
	    }
	}
      // ======================
      REAL(logLlik)[iB] = logL;
      // SAVE the result from this iteration
      REAL(postaGs)[iB] = aG;
      REAL(postrGs)[iB] = rG;

      for (ibeta = 0; ibeta < p; ibeta++) 
	{
	  idh = ibeta + iB*p;
	  REAL(postbetas)[idh] = betas[ibeta];
	}
      
      for (ipat = 0; ipat < N; ipat++ ) 
	{
	  idh = ipat + iB*N;
	  REAL(postg1s)[idh] = g1s[ipat];
	}
      
      if (2*iB < burnin) //2*iB
	{
	  // use the first half of the burn-in period to adjust the proposal variance
	  for(ih=0; ih < (2+1) ; ih++)
	    {
	      // Keept the acceptance to range between 0.3 and 0.5
	      // The proposal variance is adjusted only during the burn-in period
	      if( ih==0 ) can[ih]=FisherBeta(aG,rG,N)*3;  // update can of aG
	      else if( ih==1 ) can[ih]=FisherBeta(rG,aG,N)*3; // update can of rG
	      else 
	      {
		if (att[ih]>50)  can[ih]=adjustcan(acc[ih]/att[ih],can[ih],0.001,0.08);
	      }
	    }
	  /* else if (iB*(2*5) == burnin ) */
	  /*   { */
	  /*     for(ih=0; ih < 2+1; ih++){ */
	  /* 	att[ih] = 0.0; */
	  /* 	acc[ih] = 0.0; */
	  /*     } */
	  /*   } */
	}
      else if (iB==burnin)
	{
	  for( ih=0; ih < (2+1); ih++){
	    att[ih] = 0.0;
	    acc[ih] = 0.0;
	  }
	}  
      
      // compute the probability index
      if (probInd==1 && iB >= burnin && iB % thin == 0) 
	{ 
	  // If B = 105000, thin =10, burnin = 5000 then the first iB accepting the condition is iB = 5009
	  // Rprintf("\v iB %d idh %d",iB, idh);
	  // beta(aGs[ih],rGs[ih]) does not change with patients
	  for (ipat = 0 ; ipat < N; ipat++ )
	    {
	      idh = ipat + (iB-burnin)*N/thin;
	      if (sumYs[ipat][1]==0) REAL(probI)[idh]=1.0;
	      else {
		sizeij = 0;
		cansizeij = 0;
		for (ivec = 0 ; ivec < wIDsizes[ipat] ;ivec++ )
		  {
		    if(labelnp[wIDs[ipat][ivec]]==0) sizeij += getsizeij(ivec, X, wIDs[ipat], betas, Ntot,p);//sum_{j in pre} exp(Xij%*%beta)
		    if(labelnp[wIDs[ipat][ivec]]==1) cansizeij+=getsizeij(ivec, X, wIDs[ipat], betas, Ntot,p);//sum_{j in new} exp(Xij%*%beta)
		  }
		den = beta(sizeij+aG,sumYs[ipat][0]+rG);
		num = 0;
		for (ibeta = 0 ; ibeta < sumYs[ipat][1]; ibeta++)
		  {
		    num += choose(cansizeij+ibeta-1, ibeta)*beta(cansizeij+sizeij+aG,ibeta+sumYs[ipat][0]+rG);
		  }
		REAL(probI)[idh]=1-num/den;
	     //if (ipat == 0) Rprintf("ipat %d idh %d prob %f sumYs[ipat][1] %f num %f den %f w %f ",ipat, idh, REAL(probI)[idh],sumYs[ipat][1],num,den,w);
	      }
	    }
	}

      if(iB > 0 && iB % INTEGER(FreqEvery)[0] == 0) Rprintf("\v %d iterations are done...   ", iB);
    } // iB
  for (ih=0;ih < 2+1; ih ++) 
    {
      REAL(AR)[ih] = acc[ih]/att[ih];
      REAL(prp)[ih] = can[ih];
    }
  
  PutRNGstate();
  UNPROTECT(1);
  return res;
}


SEXP Beta1_repara(SEXP Y_,          // REAL
		  SEXP X,           // REAL 
		  SEXP ID,          // INTEGER
		  SEXP B_,          // INTEGER
		  SEXP maxni,       // INTEGER
		  SEXP nn,          // INTEGER
		  SEXP labelnp_,    // REAL
		  SEXP max_aG_,
		  SEXP mu_beta_,    // REAL
		  SEXP evalue_sigma_betas_, // REAL
		  SEXP Inv_sigma_beta_,
		  SEXP burnin_,     // INTEGER
		  SEXP FreqEvery,    // INTEGER
		  SEXP probInd_, // Int
		  SEXP thin_    // Int
		  )
{
  // === Model ===
  // Y_ij | Gi = gi ~ NB(size=exp(X_{ij}^T beta),prob=gi)
  // gi ~ beta(aG,rG)
  // beta ~ rnorm(mu_beta,sigma_beta)
  // aG, rG ~ norm
  const double *Y = REAL(Y_); 
  const double *labelnp = REAL(labelnp_);
  const int Ntot = length(ID); // # patients
  const int p = length(X)/Ntot; 
  const int N = INTEGER(nn)[0];
  const double max_aG = REAL(max_aG_)[0];
  const int B = INTEGER(B_)[0];
  const double *mu_betas = REAL(mu_beta_);
  const double *evalue_sigma_betas =REAL(evalue_sigma_betas_);
  const double *Inv_sigma_beta = REAL(Inv_sigma_beta_);
  const int burnin = INTEGER(burnin_)[0];
  const double min_aG = 0.5;
  const int thin = INTEGER(thin_)[0];
  const int probInd = INTEGER(probInd_)[0];
  int wID[INTEGER(maxni)[0]],wIDsize = 0, ih = 0, ip = 0, ipat = 0, ivec=0, idh=0, ibeta=0,iB=0;
  // wIDs[ipat][ivec] contains the position of the ivec^th repeated measure of the ipat^th patient at labelnp and Y
  //                  If ivec >= # repeated measures, wIDs[ipat][ivec] = -1000
  // wIDsizes[ipat] contains the number of repeated measure of ipat^th patient
  // sumYs[ipat][0] contains the sum of the CEL counts of ipat^th patients in the pre-scan period
  // sumYs[ipat][1] contains the sum of the CEL counts of ipat^th patients in the new-scan period
  int wIDs[N][INTEGER(maxni)[0]], wIDsizes[N];
  double sumYs[N][2],meanBeta,precBeta;

  for (ipat = 0 ; ipat < N ; ipat++)
    {
      getwID( &wIDsize,wID,ID,ipat);
      wIDsizes[ipat] = wIDsize;
      sumYs[ipat][0] = 0;
      sumYs[ipat][1] = 0;
      for (ivec = 0 ; ivec < wIDsize; ivec++)
	{
	  wIDs[ipat][ivec]=wID[ivec];
	  if (labelnp[wIDs[ipat][ivec]]==0) sumYs[ipat][0] += Y[wIDs[ipat][ivec]];
	  if (labelnp[wIDs[ipat][ivec]]==1) sumYs[ipat][1] += Y[wIDs[ipat][ivec]];
	}
      for (ivec = wIDsize ; ivec < INTEGER(maxni)[0]; ivec++)
	wIDs[ipat][ivec]=-1000;
    }

  // int h1s[N],h2s[N],js[N];
  double betas[p],canbetas[p],aG,rG,g1s[N];
  // double vs[M],weightH1[M],postprob[M],g2s[N];
  double sizeij=0.0,cansizeij=0.0,w=0.0,pij=0.0,logLcan=0,logL=0,
    candi=0.0,MHrate=0.0,yij=0.0,num=0.0,den=0.0;
  // double  weipro=0.0,qG2=0.0,D=0.0;
  double att[2+1], acc[2+1],can[2+1];
  for (ih = 0 ; ih < 2+1; ih ++)
    {
      att[ih] = 0.0;
      acc[ih] = 0.0;
      can[ih] = 2.0;
    }

  // http://stackoverflow.com/questions/8720550/how-to-return-array-of-structs-from-call-to-c-shared-library-in-r
  SEXP res = PROTECT(allocVector(VECSXP, 8)); // result is stored in res
  SEXP postaGs = allocVector(REALSXP, B); 
  SET_VECTOR_ELT(res, 0, postaGs); 
  SEXP postrGs = allocVector(REALSXP, B); 
  SET_VECTOR_ELT(res, 1, postrGs); 
  SEXP postg1s = allocVector(REALSXP, N*B); 
  SET_VECTOR_ELT(res, 2, postg1s); 
  SEXP postbetas = allocVector(REALSXP, p*B); 
  SET_VECTOR_ELT(res, 3, postbetas); 
  SEXP AR = allocVector(REALSXP,2+1);
  SET_VECTOR_ELT(res, 4, AR);
  SEXP prp = allocVector(REALSXP,2+1);
  SET_VECTOR_ELT(res, 5, prp);
  SEXP probI = allocVector(REALSXP, N*(B-burnin)/thin );
  SET_VECTOR_ELT(res, 6, probI);
  SEXP logLlik = allocVector(REALSXP, B );
  SET_VECTOR_ELT(res, 7, logLlik);
  GetRNGstate();
  
  // ============ initialization ============ 
  for (ip=0;ip < p ; ip++) betas[ip]=0.0;
  meanBeta = runif(0.0,1.0);
  precBeta = runif(min_aG,max_aG);

  aG = meanBeta*precBeta;
  rG = precBeta*(1-meanBeta);
  for (ipat = 0 ; ipat < N ; ipat++)
    g1s[ipat]=0.5;

  // ================ MCMC =================
  for (iB = 0 ; iB < B; iB++)
    {
      R_CheckUserInterrupt(); 
      //R_ProcessEvents();
      // STEP 1 
      // Update g1s, a vector of length M+1. The first M terms are the random effects.
      // the gs[M]=-1000, place holder for the patients with no new scans
      for (ipat = 0 ; ipat < N; ipat++)
	{
	  w = 0.0;
	  for (ivec=0; ivec < wIDsizes[ipat]; ivec++ )
	    {
	      w += getsizeij(ivec, X, wIDs[ipat], betas, Ntot, p);
	      // sum_{j=-1}^{n_i} exp(Xij%*%beta)
	    }
	  //Rprintf("\v ipat  %d wIDsizes[ipat] %d sumYs[0][ipat] %f sum exp(xij beta) %f",ipat, wIDsizes[ipat],sumYs[0][ipat], w);
	  g1s[ipat] = rbeta(aG+w,rG+sumYs[ipat][0]+sumYs[ipat][1]);
	  if (g1s[ipat] > 1.0-1.0e-7) g1s[ipat]=1.0-1.0e-7;
	  else if (g1s[ipat] < 1.0e-7) g1s[ipat]=1.0e-7;
	}
      
      // ================= 
      // STEP 2: update beta, a vector of length p, containing the coefficients b0, b1,...
      // beta ~ MVR(mu_beta,Sigma_beta); 

      att[2]++;
      // proposal distribution
      for (ibeta = 0; ibeta <p ; ibeta ++)
	canbetas[ibeta] = rnorm(betas[ibeta], can[2]); // the mean and sd must be double (not NumericVector)
      
      // prior
      MHrate = dmvnorm(canbetas, mu_betas,p, evalue_sigma_betas,Inv_sigma_beta,1) 
	- dmvnorm(betas, mu_betas,p, evalue_sigma_betas,Inv_sigma_beta,1);
      logLcan = 0.0;
      logL = 0.0;
      for (ipat=0; ipat < N; ipat++)
	{
	  // if (printing)  Rprintf(" \v\v ipat: %d",ipat);
	  for ( ivec = 0; ivec < wIDsizes[ipat]; ivec++ )
	    {
	      yij = Y[ wIDs[ipat][ivec] ];
	      // yij is a vector b/c it is the first arg of dnbinom
	      pij = g1s[ipat]; 
	      //Rprintf("\v\v ipat %d js[ipat] %f ivec %d yij %f g1s[ipat] %f g2s[ipat] %f pij %f", 
	      //ipat,  js[ipat] , ivec,   yij,   g1s[ipat], g2s[ipat],  pij);
	      cansizeij = getsizeij(ivec, X, wIDs[ipat], canbetas, Ntot, p); // the value of cansizeij is changed
	      sizeij = getsizeij(ivec, X, wIDs[ipat], betas, Ntot, p); // the value of cansizeij is changed

	      logLcan += dnbinom(yij, cansizeij,  pij, 1);
	      logL += dnbinom(yij, sizeij,  pij, 1);
	    }
	  //Rprintf("\v");
	}    
      MHrate += logLcan-logL;

      if(runif(0.0,1.0) < exp(MHrate))
	{
	  for (ibeta = 0; ibeta < p;ibeta++) betas[ibeta] = canbetas[ibeta];
	  acc[2]++;
	  logL = logLcan;
	}
      
      // Update aG and rG, scalars for the base distribution of the random effect gi
      // gi ~ Beta(aG,rG)
      // ================
      candi=rnorm(meanBeta,can[0]);
      if (candi > 0.0 && candi < 1.0 )
	{
	  att[0]++;
	  MHrate = 0;
	  for (ipat=0 ; ipat < N; ipat++ ) 
	    MHrate += dbeta(g1s[ipat],candi*precBeta,precBeta*(1-candi),1)
	      -dbeta(g1s[ipat],meanBeta*precBeta,precBeta*(1-meanBeta),1); 
	  if (runif(0.0,1.0) < exp(MHrate) ) 
	    {
	      meanBeta = candi;
	      acc[0]++;
	    }
	}
      
      candi=rnorm(precBeta,can[1]);
      if (candi > min_aG && candi < max_aG )
	{
	  att[1]++;
	  MHrate = 0;
	  for (ipat=0 ; ipat < N; ipat++ ) 
	    MHrate += dbeta(g1s[ipat],meanBeta*candi,candi*(1-meanBeta),1)
	      -dbeta(g1s[ipat],meanBeta*precBeta,precBeta*(1-meanBeta),1); 
	  if (runif(0.0,1.0) < exp(MHrate) ) 
	    {
	      precBeta = candi;
	      acc[1]++;
	    }
	}
      aG=meanBeta*precBeta;
      rG=precBeta*(1-meanBeta);
      // ======================
      REAL(logLlik)[iB] = logL;
      // SAVE the result from this iteration
      REAL(postaGs)[iB] = aG;
      REAL(postrGs)[iB] = rG;

      for (ibeta = 0; ibeta < p; ibeta++) 
	{
	  idh = ibeta + iB*p;
	  REAL(postbetas)[idh] = betas[ibeta];
	}
      
      for (ipat = 0; ipat < N; ipat++ ) 
	{
	  idh = ipat + iB*N;
	  REAL(postg1s)[idh] = g1s[ipat];
	}
      
      if (2*iB < burnin) //2*iB
	{
	  // use the first half of the burn-in period to adjust the proposal variance
	  for(ih=0; ih < (2+1) ; ih++)
	    {
	      // Keept the acceptance to range between 0.3 and 0.5
	      // The proposal variance is adjusted only during the burn-in period
	      if( ih==0 ) can[ih]=FisherBeta(aG,rG,N)*3;  // update can of aG
	      else if( ih==1 ) can[ih]=FisherBeta(rG,aG,N)*3; // update can of rG
	      else 
	      {
		if (att[ih]>50)  can[ih]=adjustcan(acc[ih]/att[ih],can[ih],0.001,0.08);
	      }
	    }
	  /* else if (iB*(2*5) == burnin ) */
	  /*   { */
	  /*     for(ih=0; ih < 2+1; ih++){ */
	  /* 	att[ih] = 0.0; */
	  /* 	acc[ih] = 0.0; */
	  /*     } */
	  /*   } */
	}
      else if (iB==burnin)
	{
	  for( ih=0; ih < (2+1); ih++){
	    att[ih] = 0.0;
	    acc[ih] = 0.0;
	  }
	}  
      
      // compute the probability index
      if (probInd==1 && iB >= burnin && iB % thin == 0) 
	{ 
	  // If B = 105000, thin =10, burnin = 5000 then the first iB accepting the condition is iB = 5009
	  // Rprintf("\v iB %d idh %d",iB, idh);
	  // beta(aGs[ih],rGs[ih]) does not change with patients
	  for (ipat = 0 ; ipat < N; ipat++ )
	    {
	      idh = ipat + (iB-burnin)*N/thin;
	      if (sumYs[ipat][1]==0) REAL(probI)[idh]=1.0;
	      else {
		sizeij = 0;
		cansizeij = 0;
		for (ivec = 0 ; ivec < wIDsizes[ipat] ;ivec++ )
		  {
		    if(labelnp[wIDs[ipat][ivec]]==0) sizeij += getsizeij(ivec, X, wIDs[ipat], betas, Ntot,p);//sum_{j in pre} exp(Xij%*%beta)
		    if(labelnp[wIDs[ipat][ivec]]==1) cansizeij+=getsizeij(ivec, X, wIDs[ipat], betas, Ntot,p);//sum_{j in new} exp(Xij%*%beta)
		  }
		den = beta(sizeij+aG,sumYs[ipat][0]+rG);
		num = 0;
		for (ibeta = 0 ; ibeta < sumYs[ipat][1]; ibeta++)
		  {
		    num += choose(cansizeij+ibeta-1, ibeta)*beta(cansizeij+sizeij+aG,ibeta+sumYs[ipat][0]+rG);
		  }
		REAL(probI)[idh]=1-num/den;
	     //if (ipat == 0) Rprintf("ipat %d idh %d prob %f sumYs[ipat][1] %f num %f den %f w %f ",ipat, idh, REAL(probI)[idh],sumYs[ipat][1],num,den,w);
	      }
	    }
	}

      if(iB > 0 && iB % INTEGER(FreqEvery)[0] == 0) Rprintf("\v %d iterations are done...   ", iB);
    } // iB
  for (ih=0;ih < (2+1); ih ++) 
    {
      REAL(AR)[ih] = acc[ih]/att[ih];
      REAL(prp)[ih] = can[ih];
    }
  
  PutRNGstate();
  UNPROTECT(1);
  return res;
}




SEXP gibbs(SEXP Y_,          // REAL
	   SEXP X,           // REAL 
	   SEXP ID,          // INTEGER
	   SEXP B_,          // INTEGER
	   SEXP maxni,       // INTEGER
	   SEXP nn,          // INTEGER
	   SEXP labelnp_,    // REAL
	   SEXP max_aG_,     // REAL
	   SEXP mu_beta_,    // REAL
	   SEXP evalue_sigma_betas_,
	   SEXP Inv_sigma_betas_,
	   SEXP a_D_,         // REAL
	   SEXP ib_D_,        // REAL
	   SEXP M_,           // INTEGER
	   SEXP burnin_,     // INTEGER
	   SEXP FreqEvery,    // INTEGER
	   SEXP probInd_,
	   SEXP thin_
	   )
{
  // === Model ===
  // Nonparametric const random effect over time
  // Y_ij | Gi = gi ~ NB(size=exp(X_{ij}^T beta),prob=gi)
  // gi | Hi = h ~ beta(aG_h,rG_h)
  // beta ~ rnorm(mu_beta,sigma_beta)
  // aG_h, rG_h ~ unif(0,max_aG)
  // working well 
  const double *Y = REAL(Y_); 
  const double *labelnp = REAL(labelnp_); // must be all zero 
  const int Ntot = length(ID); // # patients
  const int p = length(X)/Ntot; 
  const int N = INTEGER(nn)[0];
  const double max_aG = REAL(max_aG_)[0];
  const int B = INTEGER(B_)[0];
  const double *mu_betas = REAL(mu_beta_);
  const double *evalue_sigma_betas=REAL(evalue_sigma_betas_);
  const double *Inv_sigma_beta=REAL(Inv_sigma_betas_);
  const int burnin = INTEGER(burnin_)[0];
  const double a_D = REAL(a_D_)[0];
  const double ib_D = REAL(ib_D_)[0];
  const int M = INTEGER(M_)[0];
  const int probInd = INTEGER(probInd_)[0];
  const int thin = INTEGER(thin_)[0];
  int h1s[N],wID[INTEGER(maxni)[0]],wIDsize = 0, ih = 0, ip = 0, ipat = 0, ivec=0, idh=0, ibeta=0,iB=0;
  const double min_aG = 0.5;

  // wIDs[ipat][ivec] contains the position of the ivec^th repeated measure of the ipat^th patient at labelnp and Y
  //                  If ivec >= # repeated measures, wIDs[ipat][ivec] = -1000
  // wIDsizes[ipat] contains the number of repeated measure of ipat^th patient
  // sumYs[ipat][0] contains the sum of the CEL counts of ipat^th patients in the pre-scan period
  // sumYs[ipat][1] contains the sum of the CEL counts of ipat^th patients in the new-scan period
  int wIDs[N][INTEGER(maxni)[0]], wIDsizes[N];
  double sumYs[N][2];
  for (ipat = 0 ; ipat < N ; ipat++)
    {
      getwID( &wIDsize,wID,ID,ipat);
      wIDsizes[ipat] = wIDsize;
      sumYs[ipat][0] = 0;
      sumYs[ipat][1] = 0;
      for (ivec = 0 ; ivec < wIDsize; ivec++)
	{
	  wIDs[ipat][ivec]=wID[ivec];
	  if (labelnp[wIDs[ipat][ivec]]==0) sumYs[ipat][0] += Y[wIDs[ipat][ivec]];
	  if (labelnp[wIDs[ipat][ivec]]==1) sumYs[ipat][1] += Y[wIDs[ipat][ivec]];
	}
      for (ivec = wIDsize ; ivec < INTEGER(maxni)[0]; ivec++)
	wIDs[ipat][ivec]=-1000;
    }

  // int h1s[N],h2s[N],js[N];
  double betas[p],canbetas[p],aGs[M],rGs[M],g1s[N];
  double vs[M],weightH1[M],postprob[M],D;
  double sizeij=0.0,cansizeij=0.0,w=0.0,pij=0.0,
    candi=0.0,MHrate=0.0,yij=0.0,logLcan=0.0,logL=0.0,num,den;
  int Nacc=4;
  double att[Nacc], acc[Nacc],can[Nacc];
  for (ih = 0 ; ih < Nacc; ih ++)
    {
      att[ih] = 0.0;
      acc[ih] = 0.0;
      can[ih] = 0.5;
    }

  // http://stackoverflow.com/questions/8720550/how-to-return-array-of-structs-from-call-to-c-shared-library-in-r
  SEXP res = PROTECT(allocVector(VECSXP, 12)); // result is stored in res
  
  SEXP postaGs = allocVector(REALSXP, M*B); 
  SET_VECTOR_ELT(res, 0, postaGs); 
  SEXP postrGs = allocVector(REALSXP, M*B); 
  SET_VECTOR_ELT(res, 1, postrGs); 
  SEXP postvs = allocVector(REALSXP, M*B); 
  SET_VECTOR_ELT(res, 2, postvs); 
  SEXP postweightH1 = allocVector(REALSXP, M*B);
  SET_VECTOR_ELT(res, 3, postweightH1); 

  SEXP probI = allocVector(REALSXP, N*(B-burnin)/thin );
  SET_VECTOR_ELT(res, 4, probI);
  SEXP posth1s = allocVector(INTSXP, N*B); 
  SET_VECTOR_ELT(res, 5, posth1s); 
  SEXP postg1s = allocVector(REALSXP, N*B); 
  SET_VECTOR_ELT(res, 6, postg1s); 

  SEXP postbetas = allocVector(REALSXP, p*B); 
  SET_VECTOR_ELT(res, 7, postbetas); 

  SEXP postD = allocVector(REALSXP, B); 
  SET_VECTOR_ELT(res, 8, postD); 
  SEXP logLlik = allocVector(REALSXP, B); 
  SET_VECTOR_ELT(res, 9, logLlik);
  SEXP AR = allocVector(REALSXP,Nacc);
  SET_VECTOR_ELT(res, 10, AR);
  SEXP prp = allocVector(REALSXP,Nacc);
  SET_VECTOR_ELT(res, 11, prp);


  
  GetRNGstate();
  
  // ============ initialization ============ 
  for (ip=0;ip < p ; ip++) betas[ip]=0.0;

  for (ih = 0 ; ih < M ; ih++)
    {
      aGs[ih] = runif(min_aG,max_aG);
      rGs[ih] = runif(min_aG,max_aG);
    }
  
  D = runif(a_D,ib_D);
  
  for (ih = 0 ; ih < (M-1) ; ih ++) vs[ih] = rbeta(1,D);
  vs[M-1]=1.0;
  
  weightH1[0] = vs[0];
  w = 1.0 ;
  for (ih=1; ih < M; ih++)
    {
      w *= (1.0-vs[ih-1]);
      weightH1[ih] = vs[ih]*w;
    }


  for ( ipat=0; ipat < N; ipat++) 
    {
      h1s[ipat] = sample(weightH1,M);
      g1s[ipat] = rbeta(aGs[h1s[ipat]],rGs[h1s[ipat]]);
      if (g1s[ipat] > 1.0-1.0e-7) g1s[ipat]=1.0-1.0e-7;
      else if (g1s[ipat] < 1.0e-7) g1s[ipat]=1.0e-7;
    }

  // ================ MCMC =================
  for (iB = 0 ; iB < B; iB++)
    {
      if (iB % 100 == 0)
	{
	  R_CheckUserInterrupt(); 
	  //R_ProcessEvents();
	}
      // STEP 1 
      // Update g1s, a vector of length M+1. The first M terms are the random effects.
      // the gs[M]=-1000, place holder for the patients with no new scans
      for (ipat = 0 ; ipat < N; ipat++)
	{
	  w = 0.0;
	  for (ivec=0; ivec < wIDsizes[ipat]; ivec++ )
	    {
	      w += getsizeij(ivec, X, wIDs[ipat], betas, Ntot, p);
	      // sum_{j=-1}^{n_i} exp(Xij%*%betas)
	    }
	  //Rprintf("\v ipat  %d wIDsizes[ipat] %d sumYs[0][ipat] %f sum exp(xij betas) %f",ipat, wIDsizes[ipat],sumYs[0][ipat], w);
	  g1s[ipat] = rbeta(aGs[h1s[ipat]]+w,rGs[h1s[ipat]]+sumYs[ipat][0]+sumYs[ipat][1]);
	  if (g1s[ipat] > 1.0-1.0e-7) g1s[ipat]=1.0-1.0e-7;
	  else if (g1s[ipat] < 1.0e-7) g1s[ipat]=1.0e-7;
	}

      // =======================
      // STEP 2: update beta, a vector of length p, containing the coefficients b0, b1,...
      // beta ~ MVR(mu_beta,Sigma_beta); 
 
      att[2]++;
      // proposal distribution
      for (ibeta = 0; ibeta <p ; ibeta ++)
	canbetas[ibeta] = rnorm(betas[ibeta], can[2]); // the mean and sd must be double (not NumericVector)
      
      // prior
      MHrate = dmvnorm(canbetas, mu_betas,p, evalue_sigma_betas,Inv_sigma_beta,1) 
	- dmvnorm(betas, mu_betas,p, evalue_sigma_betas,Inv_sigma_beta,1);
      logLcan = 0.0;
      logL = 0.0;
      for (ipat=0; ipat < N; ipat++)
	{
	  // if (printing)  Rprintf(" \v\v ipat: %d",ipat);
	  for ( ivec = 0; ivec < wIDsizes[ipat]; ivec++ )
	    {
	      yij = Y[ wIDs[ipat][ivec] ];
	      // yij is a vector b/c it is the first arg of dnbinom
	      pij = g1s[ipat]; 
	      //Rprintf("\v\v ipat %d js[ipat] %f ivec %d yij %f g1s[ipat] %f g2s[ipat] %f pij %f", 
	      //ipat,  js[ipat] , ivec,   yij,   g1s[ipat], g2s[ipat],  pij);
	      cansizeij = getsizeij(ivec, X, wIDs[ipat], canbetas, Ntot, p); // the value of cansizeij is changed
	      sizeij = getsizeij(ivec, X, wIDs[ipat], betas, Ntot, p); // the value of cansizeij is changed

	      logLcan += dnbinom(yij, cansizeij,  pij, 1);
	      logL += dnbinom(yij, sizeij,  pij, 1);
	    }
	  //Rprintf("\v");
	}    
      MHrate += logLcan-logL;

      if(runif(0.0,1.0) < exp(MHrate))
	{
	  for (ibeta = 0; ibeta < p;ibeta++) betas[ibeta] = canbetas[ibeta];
	  acc[2]++;
	  logL = logLcan;
	}


      // Step 3
      // Update aG and rG, scalars for the base distribution of the random effect gi
      // gi ~ Beta(aG,rG)
      // aG ~ unif(0,max_aG)
      // rG ~ unif(0,max_aG)
      for (ih =0; ih < M; ih ++)
	{
	  candi=rnorm(aGs[ih],can[0]);
	  if (candi > min_aG && candi < max_aG )
	    {
	      if (weightH1[ih] > 0.1) att[0]++;
	      MHrate = 0;
	      for (ipat=0 ; ipat < N; ipat++ ) 
		if (h1s[ipat]==ih) MHrate += dbeta(g1s[ipat],candi,rGs[ih],1)-dbeta(g1s[ipat],aGs[ih],rGs[ih],1); 
	      if (!R_finite(MHrate)) error("MHrate is in-finite");
	      if (runif(0.0,1.0) < exp(MHrate) ) 
		{
		  aGs[ih] = candi;
		  if (weightH1[ih] > 0.1) acc[0]++;
		}
	    }
	  
	  candi=rnorm(rGs[ih],can[1]);
	  if (candi > min_aG && candi < max_aG )
	    {
	      if (weightH1[ih] > 0.1) att[1]++;
	      MHrate = 0;
	      for (ipat=0 ; ipat < N; ipat++ ) 
		if (h1s[ipat]==ih) MHrate += dbeta(g1s[ipat],aGs[ih],candi,1)-dbeta(g1s[ipat],aGs[ih],rGs[ih],1);   
	      if (runif(0.0,1.0) < exp(MHrate) ) 
		{
		  rGs[ih] = candi;
		  if (weightH1[ih] > 0.1) acc[1]++;
		}
	    }
	}
      // Update vs, a vector of length M, containing the latent variable v, used to construct pis. 
      // Once vs[ih] is updated, update the probabilities of the categorical distribution of H1 based on the formula:
      // weightH1[0] = vs[0]
      // weightH1[ih]= vs[ih]*(1-vh[ih-1])**(1-vh[ih-2])*...**(1-vh[0]) for ih =1,2,...M-1
      w = 1.0;
      for (ih = 0 ; ih < (M-1) ; ih ++) // vs[M-1] = 1 always!
	{
	  ip = 0;    // # patients with h1s == ih
	  ibeta = 0; // # patients with h1s > ih
	  for (ipat = 0 ; ipat < N; ipat++ ) 
	    {
	      if ( h1s[ipat]==ih ) ip++ ;
	      if ( h1s[ipat] > ih ) ibeta++ ;
	    }
	  
	  vs[ih] = rbeta(
			 (double) 1.0 +ip, 
			 D + (double)ibeta 
			 );
	  //Rprintf("ih %d sp1 %f sp2 %f",ih, (double) 1+ip+idh ,  D + (double)ibeta + (double)ivec);
	  if (ih > 0) w *= (1.0-vs[ih-1]);
	  weightH1[ih] = vs[ih]*w;
	}
      weightH1[M-1] = w*(1-vs[M-2]);


      //  Update h1s, a vector of length N, containing the cluster index of random effect 1 of each patient.
      //  the cluster index is between 1 and M.
      //  conjugate: posterior distribution of h1s has a categorical distribution with M categories. 
      // if (printing && iB%printEvery==0)  Rprintf("\v\v Step 4");
      for (ipat = 0 ; ipat < N ; ipat++)
	{
	  
	  // Given patient ipat, compute P(H_1i=ih;-) for each ih 
	  // P(H_1i=ih;-) = prod_{i=1}^{# repeated measure} dnbinom(yij; size=exp(b0+b1*x1i+..), prob=1/(gij+1))
	  
	  for (ih=0; ih < M; ih++ ) postprob[ih] = dbeta(g1s[ipat],aGs[ih],rGs[ih],0)*weightH1[ih];
	  h1s[ipat] = sample(postprob,M); 
	  
	}

      // Update D
      // D ~ runif(0.01,ib_D) 
      candi = rnorm(D,can[3]);
      att[3]++;
      if (a_D < candi && candi < ib_D)
	{
	  MHrate = 0.0;
	  for (ih = 0 ; ih < (M-1);ih++)
	    {
	      MHrate += dbeta(vs[ih],1.0,candi,1)-dbeta(vs[ih],1.0,D,1);
	    }
	  if (runif(0.0,1.0) < exp(MHrate) ) 
	    {
	      D = candi;
	      acc[3]++;
	    }
	}
      
      // compute the probability index
      if (probInd==1 && iB >= burnin && iB % thin == 0) 
	{ 
	  // If B = 105000, thin =10, burnin = 5000 then the first iB accepting the condition is iB = 5009
	  // Rprintf("\v iB %d idh %d",iB, idh);
	  // beta(aGs[ih],rGs[ih]) does not change with patients
	  for (ih= 0 ; ih < M ; ih++) postprob[ih] = beta(aGs[ih],rGs[ih]);

	  for (ipat = 0 ; ipat < N; ipat++ )
	    {
	      idh = ipat + (iB-burnin)*N/thin;
	      if (sumYs[ipat][1]==0) REAL(probI)[idh]=1.0;
	      else {
		sizeij = 0;
		cansizeij = 0;
		for (ivec = 0 ; ivec < wIDsizes[ipat] ;ivec++ )
		  {
		    if(labelnp[wIDs[ipat][ivec]]==0) sizeij += getsizeij(ivec, X, wIDs[ipat], betas, Ntot,p);//sum_{j in pre} exp(Xij%*%beta)
		    if(labelnp[wIDs[ipat][ivec]]==1) cansizeij+=getsizeij(ivec, X, wIDs[ipat], betas, Ntot,p);//sum_{j in new} exp(Xij%*%beta)
		  }
		den = 0;
		for (ih = 0 ; ih < M ; ih++)
		  {
		    if ( weightH1[ih]!= 0)
		      {
			den += weightH1[ih]*beta(sizeij+aGs[ih],sumYs[ipat][0]+rGs[ih])/postprob[ih];
		      }
		  }
		num = 0;
		for (ibeta = 0 ; ibeta < sumYs[ipat][1]; ibeta++)
		  {
		    w = 0;
		    for (ih = 0 ; ih < M ; ih++)
		      {
			if ( weightH1[ih] != 0)
			  {
			    w += weightH1[ih]*beta(cansizeij+sizeij+aGs[ih],ibeta+sumYs[ipat][0]+rGs[ih])/postprob[ih];
			  }
		      }
		    num += choose(cansizeij+ibeta-1, ibeta)*w;
		  }
		REAL(probI)[idh]=1-num/den;
		//if (ipat == 0) Rprintf("ipat %d idh %d prob %f sumYs[ipat][1] %f num %f den %f w %f ",ipat, idh, REAL(probI)[idh],sumYs[ipat][1],num,den,w);
	      }
	     
	    }
	  
	}

      // SAVE the result from this iteration
      REAL(postD)[iB]=D;
      REAL(logLlik)[iB] = logL;
      for (ih = 0; ih < M ; ih++)
	{
	  idh = ih + iB*M;
	  REAL(postaGs)[idh] = aGs[ih];
	  REAL(postrGs)[idh] = rGs[ih];
	  REAL(postvs)[idh] = vs[ih];
	  REAL(postweightH1)[idh] = weightH1[ih];
	}
      for (ibeta = 0; ibeta < p; ibeta++) 
	{
	  idh = ibeta + iB*p;
	  REAL(postbetas)[idh] = betas[ibeta];
	}
      
      for (ipat = 0; ipat < N; ipat++ ) 
	{
	  idh = ipat + iB*N;
	  REAL(postg1s)[idh] = g1s[ipat];
	  INTEGER(posth1s)[idh] = h1s[ipat];
	}
      
      if (2*iB < burnin) //2*iB
	{
	  // use the first half of the burn-in period to adjust the proposal variance
	  for(ih=0; ih < Nacc ; ih++)
	    {
	      // Keept the acceptance to range between 0.3 and 0.5
	      // The proposal variance is adjusted only during the burn-in period
	      if( ih==0 ) 
		{
		  // aG
		  if (att[ih]>50) can[ih]=adjustcan(acc[ih]/att[ih],can[ih],0.05,5);
		}
	      else if( ih==1 )
		{
		  // rG
		  if (att[ih]>50) can[ih]=adjustcan(acc[ih]/att[ih],can[ih],0.05,5);
		}
	      else 
		{
		  // beta , D
		  if (att[ih]>50)  can[ih]=adjustcan(acc[ih]/att[ih],can[ih],0.02,2);
		}
	    }
	}
      else if (iB==burnin)
	{
	  for( ih=0; ih < Nacc; ih++){
	    att[ih] = 0.0;
	    acc[ih] = 0.0;
	  }
	}  
      
      if(iB > 0 && iB % INTEGER(FreqEvery)[0] == 0) 
	{
	  Rprintf("\v %d iterations are done...   ", iB);
	  R_FlushConsole(); 
	  //R_ProcessEvents();
	}
    }
  for (ih=0;ih < Nacc; ih ++) 
    {
      REAL(AR)[ih] = acc[ih]/att[ih];
      REAL(prp)[ih] = can[ih];
    }
  
  PutRNGstate();
  UNPROTECT(1);
  return res;
}

double betaYK (double a, double b)
{
  // No warning even for very large a and b 
  // such beta(a,b) is zero
  return (exp(lbeta(a,b)));
}
double IntegrateDensNB(const double sumYsipat[],
		       const int wIDsizesipat, SEXP X, const int wIDsipat[],const double betas[], 
		       const int Ntot, const int p,const double Y[], const double aGh, const double rGh,
		       const int LOG,double *sizeip,const double Beta_ar_ih)
{
  // Compute the density of Y_i (a vector of length ni+2) given aGs[pat], rGs[ipat] and beta
  // p(Y_{i,-1},...,Y_{i,ni}| beta, aGs[h1s[ipat]],rGs[h1s[ipat]] )
  //                = [ prod_{j=-1}^{ni} choose(sizeij+Yij-1,Yij)]
  //                   *beta(sizeip+aGs[h1s[ipat]],yip+rGs[h1s[ipat]])/beta(aGs[h1s[ipat]],rGs[h1s[ipat]])
  *sizeip = 0.0;
  double sizeij = 0.0; 
  double logdensity=0.0;
  for (int ivec=0; ivec < wIDsizesipat; ivec++ )
    {
      sizeij = getsizeij(ivec, X, wIDsipat, betas, Ntot, p);
      *sizeip += sizeij;
      logdensity += log(choose(sizeij+Y[wIDsipat[ivec]]-1,Y[wIDsipat[ivec]]));
    }
  //Rprintf("\n sizeip %f aGh %f sumYs %f rG %f beta %f",*sizeip,aGh,sumYsipat[0]+sumYsipat[1],rGh,
  //betaYK(*sizeip+aGh,sumYsipat[0]+sumYsipat[1]+rGh));
  logdensity += log(betaYK(*sizeip+aGh,sumYsipat[0]+sumYsipat[1]+rGh))-log(Beta_ar_ih);
  if (LOG) return(logdensity);
  else return(exp(logdensity));
}


SEXP ReduceGibbs(SEXP Y_,          // REAL
		 SEXP X,           // REAL 
		 SEXP ID,          // INTEGER
		 SEXP B_,          // INTEGER
		 SEXP maxni,       // INTEGER
		 SEXP nn,          // INTEGER
		 SEXP labelnp_,    // REAL
		 SEXP max_aG_,     // REAL
		 SEXP mu_beta_,    // REAL
		 SEXP evalue_sigma_betas_,
		 SEXP Inv_sigma_betas_,
		 SEXP a_D_,         // REAL
		 SEXP ib_D_,        // REAL
		 SEXP M_,           // INTEGER
		 SEXP burnin_,     // INTEGER
		 SEXP FreqEvery,    // INTEGER
		 SEXP probInd_,
		 SEXP thin_
		 )
{
  // === Model ===
  // Nonparametric const random effect over time
  // Y_ij | Gi = gi ~ NB(size=exp(X_{ij}^T beta),prob=gi)
  // gi | Hi = h ~ beta(aG_h,rG_h)
  // beta ~ rnorm(mu_beta,sigma_beta)
  // aG_h, rG_h ~ unif(0,max_aG)
  // working well 
  const double *Y = REAL(Y_); 
  const double *labelnp = REAL(labelnp_); // must be all zero 
  const int Ntot = length(ID); // # patients
  const int p = length(X)/Ntot; 
  const int N = INTEGER(nn)[0];
  const double max_aG = REAL(max_aG_)[0];
  const int B = INTEGER(B_)[0];
  const double *mu_betas = REAL(mu_beta_);
  const double *evalue_sigma_betas=REAL(evalue_sigma_betas_);
  const double *Inv_sigma_beta=REAL(Inv_sigma_betas_);
  const int burnin = INTEGER(burnin_)[0];
  const double a_D = REAL(a_D_)[0];
  const double ib_D = REAL(ib_D_)[0];
  const int M = INTEGER(M_)[0];
  const int probInd = INTEGER(probInd_)[0];
  const int thin = INTEGER(thin_)[0];
  int h1s[N],wID[INTEGER(maxni)[0]],wIDsize = 0, ih = 0, ip = 0, ipat = 0, ivec=0, idh=0, ibeta=0,iB=0;
  const double min_aG = 0.5;

  // wIDs[ipat][ivec] contains the position of the ivec^th repeated measure of the ipat^th patient at labelnp and Y
  //                  If ivec >= # repeated measures, wIDs[ipat][ivec] = -1000
  // wIDsizes[ipat] contains the number of repeated measure of ipat^th patient
  // sumYs[ipat][0] contains the sum of the CEL counts of ipat^th patients in the pre-scan period
  // sumYs[ipat][1] contains the sum of the CEL counts of ipat^th patients in the new-scan period
  int wIDs[N][INTEGER(maxni)[0]], wIDsizes[N];
  double sumYs[N][2];
  for (ipat = 0 ; ipat < N ; ipat++)
    {
      getwID( &wIDsize,wID,ID,ipat);
      wIDsizes[ipat] = wIDsize;
      sumYs[ipat][0] = 0;
      sumYs[ipat][1] = 0;
      for (ivec = 0 ; ivec < wIDsize; ivec++)
	{
	  wIDs[ipat][ivec]=wID[ivec];
	  if (labelnp[wIDs[ipat][ivec]]==0) sumYs[ipat][0] += Y[wIDs[ipat][ivec]];
	  if (labelnp[wIDs[ipat][ivec]]==1) sumYs[ipat][1] += Y[wIDs[ipat][ivec]];
	}
      for (ivec = wIDsize ; ivec < INTEGER(maxni)[0]; ivec++)
	wIDs[ipat][ivec]=-1000;
    }

  // int h1s[N],h2s[N],js[N];
  double betas[p],canbetas[p],aGs[M],rGs[M],Beta_ar[M]; 
  double vs[M],weightH1[M],postprob[M],D;
  double sizeij=0.0,cansizeij=0.0,w=0.0,pij=0.0,sizeip_s[N],cansizeip_s[N],
    candi=0.0,MHrate=0.0,yij=0.0,logLcan=0.0,logL=0.0,num,den,canlogDen=0,logDen=0;
  int Nacc=4;
  double att[Nacc], acc[Nacc],can[Nacc];
  for (ih = 0 ; ih < Nacc; ih ++)
    {
      att[ih] = 0.0;
      acc[ih] = 0.0;
      can[ih] = 0.5;
    }

  // http://stackoverflow.com/questions/8720550/how-to-return-array-of-structs-from-call-to-c-shared-library-in-r
  SEXP res = PROTECT(allocVector(VECSXP, 14)); // result is stored in res
  
  SEXP postaGs = allocVector(REALSXP, M*B); 
  SET_VECTOR_ELT(res, 0, postaGs); 
  SEXP postrGs = allocVector(REALSXP, M*B); 
  SET_VECTOR_ELT(res, 1, postrGs); 
  SEXP postvs = allocVector(REALSXP, M*B); 
  SET_VECTOR_ELT(res, 2, postvs); 
  SEXP postweightH1 = allocVector(REALSXP, M*B);
  SET_VECTOR_ELT(res, 3, postweightH1); 

  SEXP probI = allocVector(REALSXP, N*(B-burnin)/thin );
  SET_VECTOR_ELT(res, 4, probI);
  SEXP posth1s = allocVector(INTSXP, N*B); 
  SET_VECTOR_ELT(res, 5, posth1s); 
  SEXP postg1s = allocVector(REALSXP, N*B); 
  SET_VECTOR_ELT(res, 6, postg1s); 

  SEXP postbetas = allocVector(REALSXP, p*B); 
  SET_VECTOR_ELT(res, 7, postbetas); 

  SEXP postD = allocVector(REALSXP, B); 
  SET_VECTOR_ELT(res, 8, postD); 
  SEXP logLlik = allocVector(REALSXP, B); 
  SET_VECTOR_ELT(res, 9, logLlik);
  SEXP AR = allocVector(REALSXP,Nacc);
  SET_VECTOR_ELT(res, 10, AR);
  SEXP prp = allocVector(REALSXP,Nacc);
  SET_VECTOR_ELT(res, 11, prp);

  SEXP postaGs_pat = allocVector(REALSXP, N*B); 
  SET_VECTOR_ELT(res, 12, postaGs_pat); 
  SEXP postrGs_pat = allocVector(REALSXP, N*B); 
  SET_VECTOR_ELT(res, 13, postrGs_pat); 

  
  GetRNGstate();
  
  // ============ initialization ============ 
  for (ip=0;ip < p ; ip++) betas[ip]=0.0;

  for (ih = 0 ; ih < M ; ih++)
    {
      aGs[ih] = runif(min_aG,max_aG);
      rGs[ih] = runif(min_aG,max_aG);
      Beta_ar[ih] = beta(aGs[ih],rGs[ih]);
    }
  
  D = runif(a_D,ib_D);
  
  for (ih = 0 ; ih < (M-1) ; ih ++) vs[ih] = rbeta(1,D);
  vs[M-1]=1.0;
  
  weightH1[0] = vs[0];
  w = 1.0 ;
  for (ih=1; ih < M; ih++)
    {
      w *= (1.0-vs[ih-1]);
      weightH1[ih] = vs[ih]*w;
    }


  for ( ipat=0; ipat < N; ipat++) 
    {
      h1s[ipat] = sample(weightH1,M);
    }

  // ================ MCMC =================
  for (iB = 0 ; iB < B; iB++)
    {
      if (iB % 100 == 0)
	{
	  R_CheckUserInterrupt(); 
	}

      // =======================
      // STEP 1: update beta, a vector of length p, containing the coefficients b0, b1,...
      // beta ~ MVR(mu_beta,Sigma_beta); 
      att[2]++;
      // proposal distribution
      for (ibeta = 0; ibeta <p ; ibeta ++)
	canbetas[ibeta] = rnorm(betas[ibeta], can[2]); // the mean and sd must be double (not NumericVector)
      
      // prior
      MHrate = dmvnorm(canbetas, mu_betas,p, evalue_sigma_betas,Inv_sigma_beta,1) 
	- dmvnorm(betas, mu_betas,p, evalue_sigma_betas,Inv_sigma_beta,1);
      logLcan = 0.0;
      logL = 0.0;
      for (ipat=0; ipat < N; ipat++)
	{
	  logLcan += IntegrateDensNB(sumYs[ipat],wIDsizes[ipat],X,wIDs[ipat],canbetas, 
				     Ntot,p,Y,aGs[h1s[ipat]],rGs[h1s[ipat]],1,&cansizeip_s[ipat],Beta_ar[h1s[ipat]]);
	  logL += IntegrateDensNB(sumYs[ipat],wIDsizes[ipat],X,wIDs[ipat],betas, 
				  Ntot,p,Y,aGs[h1s[ipat]],rGs[h1s[ipat]],1,&sizeip_s[ipat],Beta_ar[h1s[ipat]]);
	}    
      MHrate += logLcan - logL;
      if(runif(0.0,1.0) < exp(MHrate))
	{
	  for (ibeta = 0; ibeta < p;ibeta++) betas[ibeta] = canbetas[ibeta];
	  acc[2]++;
	  logL = logLcan;
	  // sizeip_s contains the sum_j exp(Xij%*%beta) for each patient i with current beta
	  for (ipat = 0;ipat < N; ipat++) sizeip_s[ipat] = cansizeip_s[ipat];
	}

      // Step 2: Update aG and rG, scalars for the base distribution of the random effect gi
      // gi ~ Beta(aG,rG)
      // aG ~ unif(0,max_aG)
      // rG ~ unif(0,max_aG)
      for (ih =0; ih < M; ih ++)
	{
	  candi=rnorm(aGs[ih],can[0]);
	  logDen = log(Beta_ar[ih]);

	  if (weightH1[ih] > 0.1) att[0]++;

	  if (candi > min_aG && candi < max_aG )
	    {
	      canlogDen = log(beta(candi,rGs[ih]));

	      MHrate = 0;
	      for (ipat=0 ; ipat < N; ipat++ ) 
		if (h1s[ipat]==ih) 
		  MHrate += log(betaYK(sizeip_s[ipat]+candi,sumYs[ipat][0]+sumYs[ipat][1]+rGs[ih]))-canlogDen-
		    (log(betaYK(sizeip_s[ipat]+aGs[ih],sumYs[ipat][0]+sumYs[ipat][1]+rGs[ih]))-logDen);
	      if (!R_finite(MHrate)) error("MHrate is in-finite");
	      if (runif(0.0,1.0) < exp(MHrate) ) 
		{
		  aGs[ih] = candi;
		  logDen = canlogDen;
		  if (weightH1[ih] > 0.1) acc[0]++;
		}
	    }
	  // Step 3
	  candi=rnorm(rGs[ih],can[1]);
	  if (weightH1[ih] > 0.1) att[1]++;

	  if (candi > min_aG && candi < max_aG )
	    {
	      canlogDen = log(beta(aGs[ih],candi));
	      MHrate = 0;
	      for (ipat=0 ; ipat < N; ipat++ ) 
		if (h1s[ipat]==ih) 
		  MHrate += log(betaYK(sizeip_s[ipat]+aGs[ih],sumYs[ipat][0]+sumYs[ipat][1]+candi))-canlogDen-
		    (log(betaYK(sizeip_s[ipat]+aGs[ih],sumYs[ipat][0]+sumYs[ipat][1]+rGs[ih]))-logDen);
	  
	      if (runif(0.0,1.0) < exp(MHrate) ) 
		{
		  rGs[ih] = candi; 
		  logDen = canlogDen;
		  if (weightH1[ih] > 0.1) acc[1]++;
		}
	    }
	  Beta_ar[ih]=exp(logDen);
	}
      // Update vs, a vector of length M, containing the latent variable v, used to construct pis. 
      // Once vs[ih] is updated, update the probabilities of the categorical distribution of H1 based on the formula:
      // weightH1[0] = vs[0]
      // weightH1[ih]= vs[ih]*(1-vh[ih-1])**(1-vh[ih-2])*...**(1-vh[0]) for ih =1,2,...M-1
      w = 1.0;
      for (ih = 0 ; ih < (M-1) ; ih ++) // vs[M-1] = 1 always!
	{
	  ip = 0;    // # patients with h1s == ih
	  ibeta = 0; // # patients with h1s > ih
	  for (ipat = 0 ; ipat < N; ipat++ ) 
	    {
	      if ( h1s[ipat]==ih ) ip++ ;
	      if ( h1s[ipat] > ih ) ibeta++ ;
	    }
	  
	  vs[ih] = rbeta(
			 (double) 1.0 +ip, 
			 D + (double)ibeta 
			 );
	  //Rprintf("ih %d sp1 %f sp2 %f",ih, (double) 1+ip+idh ,  D + (double)ibeta + (double)ivec);
	  if (ih > 0) w *= (1.0-vs[ih-1]);
	  weightH1[ih] = vs[ih]*w;
	}
      weightH1[M-1] = w*(1-vs[M-2]);


      //  Update h1s, a vector of length N, containing the cluster index of random effect 1 of each patient.
      //  the cluster index is between 1 and M.
      //  conjugate: posterior distribution of h1s has a categorical distribution with M categories. 
      // if (printing && iB%printEvery==0)  Rprintf("\v\v Step 4");
      for (ipat = 0 ; ipat < N ; ipat++)
	{
	  // Given patient ipat, compute P(H_1i=ih;-) for each ih 
	  // P(H_1i=ih;-) = prod_{i=1}^{# repeated measure} dnbinom(yij; size=exp(b0+b1*x1i+..), prob=1/(gij+1))
	  for (ih=0; ih < M; ih++ ) 
	    postprob[ih]= (betaYK(sizeip_s[ipat]+aGs[ih],sumYs[ipat][0]+sumYs[ipat][1]+rGs[ih])/Beta_ar[ih])*weightH1[ih];
	  h1s[ipat] = sample(postprob,M); 
	  
	}

      // Update D
      // D ~ runif(0.01,ib_D) 
      candi = rnorm(D,can[3]);
      att[3]++;
      if (a_D < candi && candi < ib_D)
	{
	  MHrate = 0.0;
	  for (ih = 0 ; ih < (M-1);ih++)
	    {
	      MHrate += dbeta(vs[ih],1.0,candi,1)-dbeta(vs[ih],1.0,D,1);
	    }
	  if (runif(0.0,1.0) < exp(MHrate) ) 
	    {
	      D = candi;
	      acc[3]++;
	    }
	}
      
      // compute the probability index
      if (probInd==1 && iB >= burnin && iB % thin == 0) 
	{ 
	  // If B = 105000, thin =10, burnin = 5000 then the first iB accepting the condition is iB = 5009
	  // Rprintf("\v iB %d idh %d",iB, idh);
	  // beta(aGs[ih],rGs[ih]) does not change with patients

	  for (ipat = 0 ; ipat < N; ipat++ )
	    {
	      idh = ipat + (iB-burnin)*N/thin;
	      if (sumYs[ipat][1]==0) REAL(probI)[idh]=1.0;
	      else {
		sizeij = 0;
		cansizeij = 0;
		for (ivec = 0 ; ivec < wIDsizes[ipat] ;ivec++ )
		  {
		    if(labelnp[wIDs[ipat][ivec]]==0) sizeij += getsizeij(ivec, X, wIDs[ipat], betas, Ntot,p);//sum_{j in pre} exp(Xij%*%beta)
		    if(labelnp[wIDs[ipat][ivec]]==1) cansizeij+=getsizeij(ivec, X, wIDs[ipat], betas, Ntot,p);//sum_{j in new} exp(Xij%*%beta)
		  }
		den = 0;
		for (ih = 0 ; ih < M ; ih++)
		  {
		    if ( weightH1[ih]!= 0)
		      {
			den += weightH1[ih]*betaYK(sizeij+aGs[ih],sumYs[ipat][0]+rGs[ih])/Beta_ar[ih];
		      }
		  }
		num = 0;
		for (ibeta = 0 ; ibeta < sumYs[ipat][1]; ibeta++)
		  {
		    w = 0;
		    for (ih = 0 ; ih < M ; ih++)
		      {
			if ( weightH1[ih] != 0)
			  {
			    w += weightH1[ih]*betaYK(cansizeij+sizeij+aGs[ih],ibeta+sumYs[ipat][0]+rGs[ih])/Beta_ar[ih];
			  }
		      }
		    num += choose(cansizeij+ibeta-1, ibeta)*w;
		  }
		REAL(probI)[idh]=1-num/den;
		//if (ipat == 0) Rprintf("ipat %d idh %d prob %f sumYs[ipat][1] %f num %f den %f w %f ",ipat, idh, REAL(probI)[idh],sumYs[ipat][1],num,den,w);
	      }
	     
	    }
	  
	}

      // SAVE the result from this iteration
      REAL(postD)[iB]=D;
      REAL(logLlik)[iB] = logL;
      for (ih = 0; ih < M ; ih++)
	{
	  idh = ih + iB*M;
	  REAL(postaGs)[idh] = aGs[ih];
	  REAL(postrGs)[idh] = rGs[ih];
	  REAL(postvs)[idh] = vs[ih];
	  REAL(postweightH1)[idh] = weightH1[ih];
	}
      for (ibeta = 0; ibeta < p; ibeta++) 
	{
	  idh = ibeta + iB*p;
	  REAL(postbetas)[idh] = betas[ibeta];
	}
      
      for (ipat = 0; ipat < N; ipat++ ) 
	{
	  idh = ipat + iB*N;
	  //REAL(postg1s)[idh] = g1s[ipat];
	  INTEGER(posth1s)[idh] = h1s[ipat];
	  REAL(postaGs_pat)[idh] = aGs[h1s[ipat]];
	  REAL(postrGs_pat)[idh] = rGs[h1s[ipat]];
	}
      
      if (2*iB < burnin) //2*iB
	{
	  // use the first half of the burn-in period to adjust the proposal variance
	  for(ih=0; ih < Nacc ; ih++)
	    {
	      // Keept the acceptance to range between 0.3 and 0.5
	      // The proposal variance is adjusted only during the burn-in period
	      if( ih==0 ) 
		{
		  // aG
		  if (att[ih]>50) can[ih]=adjustcan(acc[ih]/att[ih],can[ih],0.05,5);
		}
	      else if( ih==1 )
		{
		  // rG
		  if (att[ih]>50) can[ih]=adjustcan(acc[ih]/att[ih],can[ih],0.05,5);
		}
	      else 
		{
		  // beta , D
		  if (att[ih]>50)  can[ih]=adjustcan(acc[ih]/att[ih],can[ih],0.02,2);
		}
	    }
	}
      else if (iB==burnin)
	{
	  for( ih=0; ih < Nacc; ih++){
	    att[ih] = 0.0;
	    acc[ih] = 0.0;
	  }
	}  
      
      if(iB > 0 && iB % INTEGER(FreqEvery)[0] == 0) 
	{
	  Rprintf("\v %d iterations are done...   ", iB);
	  R_FlushConsole(); 
	  //R_ProcessEvents();
	}
    }
  for (ih=0;ih < Nacc; ih ++) 
    {
      REAL(AR)[ih] = acc[ih]/att[ih];
      REAL(prp)[ih] = can[ih];
    }
  
  PutRNGstate();
  UNPROTECT(1);
  return res;
}





SEXP Beta1reduce(SEXP Y_,          // REAL
		 SEXP X,           // REAL 
		 SEXP ID,          // INTEGER
		 SEXP B_,          // INTEGER
		 SEXP maxni,       // INTEGER
		 SEXP nn,          // INTEGER
		 SEXP labelnp_,    // REAL
		 SEXP max_aG_,
		 SEXP mu_beta_,    // REAL
		 SEXP evalue_sigma_betas_, // REAL
		 SEXP Inv_sigma_beta_,
		 SEXP burnin_,     // INTEGER
		 SEXP FreqEvery,    // INTEGER
		 SEXP probInd_, // Int
		 SEXP thin_    // Int
	   )
{
  // === Model ===
  // Y_ij | Gi = gi ~ NB(size=exp(X_{ij}^T beta),prob=gi)
  // gi ~ beta(aG,rG)
  // beta ~ rnorm(mu_beta,sigma_beta)
  // aG, rG ~ norm
  const double *Y = REAL(Y_); 
  const double *labelnp = REAL(labelnp_);
  const int Ntot = length(ID); // # patients
  const int p = length(X)/Ntot; 
  const int N = INTEGER(nn)[0];
  const double max_aG = REAL(max_aG_)[0];
  const int B = INTEGER(B_)[0];
  const double *mu_betas = REAL(mu_beta_);
  const double *evalue_sigma_betas =REAL(evalue_sigma_betas_);
  const double *Inv_sigma_beta = REAL(Inv_sigma_beta_);
  const int burnin = INTEGER(burnin_)[0];
  const double min_aG = 0.5;
  const int thin = INTEGER(thin_)[0];
  const int probInd = INTEGER(probInd_)[0];
  int wID[INTEGER(maxni)[0]],wIDsize = 0, ih = 0, ip = 0, ipat = 0, ivec=0, idh=0, ibeta=0,iB=0;
  // wIDs[ipat][ivec] contains the position of the ivec^th repeated measure of the ipat^th patient at labelnp and Y
  //                  If ivec >= # repeated measures, wIDs[ipat][ivec] = -1000
  // wIDsizes[ipat] contains the number of repeated measure of ipat^th patient
  // sumYs[ipat][0] contains the sum of the CEL counts of ipat^th patients in the pre-scan period
  // sumYs[ipat][1] contains the sum of the CEL counts of ipat^th patients in the new-scan period
  int wIDs[N][INTEGER(maxni)[0]], wIDsizes[N];
  double sumYs[N][2];


  for (ipat = 0 ; ipat < N ; ipat++)
    {
      getwID( &wIDsize,wID,ID,ipat);
      wIDsizes[ipat] = wIDsize;
      sumYs[ipat][0] = 0;
      sumYs[ipat][1] = 0;
      for (ivec = 0 ; ivec < wIDsize; ivec++)
	{
	  wIDs[ipat][ivec]=wID[ivec];
	  if (labelnp[wIDs[ipat][ivec]]==0) sumYs[ipat][0] += Y[wIDs[ipat][ivec]];
	  if (labelnp[wIDs[ipat][ivec]]==1) sumYs[ipat][1] += Y[wIDs[ipat][ivec]];
	}
      for (ivec = wIDsize ; ivec < INTEGER(maxni)[0]; ivec++)
	wIDs[ipat][ivec]=-1000;
    }

  // int h1s[N],h2s[N],js[N];
  double betas[p],canbetas[p],aG,rG,g1s[N],sizeip_s[N],cansizeip_s[N],Beta_ar;
  // double vs[M],weightH1[M],postprob[M],g2s[N];
  double sizeij=0.0,cansizeij=0.0,w=0.0,pij=0.0,logLcan=0,logL=0,
    candi=0.0,MHrate=0.0,yij=0.0,num=0.0,den=0.0,logDen=0.0,canlogDen=0.0;
  // double  weipro=0.0,qG2=0.0,D=0.0;
  const int Nacc= 2+1;
  double att[Nacc], acc[Nacc],can[Nacc];
  for (ih = 0 ; ih < Nacc; ih ++)
    {
      att[ih] = 0.0;
      acc[ih] = 0.0;
      can[ih] = 2.0;
    }

  // http://stackoverflow.com/questions/8720550/how-to-return-array-of-structs-from-call-to-c-shared-library-in-r
  SEXP res = PROTECT(allocVector(VECSXP, 8)); // result is stored in res
  SEXP postaGs = allocVector(REALSXP, B); 
  SET_VECTOR_ELT(res, 0, postaGs); 
  SEXP postrGs = allocVector(REALSXP, B); 
  SET_VECTOR_ELT(res, 1, postrGs); 
  SEXP postg1s = allocVector(REALSXP, N*B); 
  SET_VECTOR_ELT(res, 2, postg1s); 
  SEXP postbetas = allocVector(REALSXP, p*B); 
  SET_VECTOR_ELT(res, 3, postbetas); 
  SEXP AR = allocVector(REALSXP,Nacc);
  SET_VECTOR_ELT(res, 4, AR);
  SEXP prp = allocVector(REALSXP,Nacc);
  SET_VECTOR_ELT(res, 5, prp);
  SEXP probI = allocVector(REALSXP, N*(B-burnin)/thin );
  SET_VECTOR_ELT(res, 6, probI);
  SEXP logLlik = allocVector(REALSXP, B );
  SET_VECTOR_ELT(res, 7, logLlik);
  GetRNGstate();
  
  // ============ initialization ============ 
  for (ip=0;ip < p ; ip++) betas[ip]=0.0;
  aG = runif(min_aG,max_aG);
  rG = runif(min_aG,max_aG);
  Beta_ar = beta(aG,rG);
  for (ipat = 0 ; ipat < N ; ipat++)
    g1s[ipat]=0.5;

  // ================ MCMC =================
  for (iB = 0 ; iB < B; iB++)
    {
      R_CheckUserInterrupt(); 
      // ================= 
      // STEP 2: update beta, a vector of length p, containing the coefficients b0, b1,...
      // beta ~ MVR(mu_beta,Sigma_beta); 

      att[2]++;
      // proposal distribution
      for (ibeta = 0; ibeta <p ; ibeta ++)
	canbetas[ibeta] = rnorm(betas[ibeta], can[2]); // the mean and sd must be double (not NumericVector)
      
      // prior
      MHrate = dmvnorm(canbetas, mu_betas,p, evalue_sigma_betas,Inv_sigma_beta,1) 
	- dmvnorm(betas, mu_betas,p, evalue_sigma_betas,Inv_sigma_beta,1);
      logLcan = 0.0;
      logL = 0.0;
      for (ipat=0; ipat < N; ipat++)
	{
	  logLcan += IntegrateDensNB(sumYs[ipat],wIDsizes[ipat],X,wIDs[ipat],canbetas, 
				     Ntot,p,Y,aG,rG,1,&cansizeip_s[ipat],Beta_ar);
	  logL += IntegrateDensNB(sumYs[ipat],wIDsizes[ipat],X,wIDs[ipat],betas, 
				  Ntot,p,Y,aG,rG,1,&sizeip_s[ipat],Beta_ar);
	  //Rprintf("\v");
	}    
      MHrate += logLcan-logL;

      if(runif(0.0,1.0) < exp(MHrate))
	{
	  for (ibeta = 0; ibeta < p;ibeta++) betas[ibeta] = canbetas[ibeta];
	  acc[2]++;
	  logL = logLcan;
	  for (ipat=0;ipat<N;ipat++)sizeip_s[ipat] = cansizeip_s[ipat];
	}
      
      // Update aG and rG, scalars for the base distribution of the random effect gi
      // gi ~ Beta(aG,rG)
      // ================

      logDen = log(Beta_ar);

      candi=rnorm(aG,can[0]);
      att[0]++;
      if (candi > min_aG && candi < max_aG )
	{
	  canlogDen = log(beta(candi,rG));
	  MHrate = 0;
	  
	  for (ipat=0 ; ipat < N; ipat++ ) 
	      MHrate += log(betaYK(sizeip_s[ipat]+candi,sumYs[ipat][0]+sumYs[ipat][1]+rG))-canlogDen-
		    (log(betaYK(sizeip_s[ipat]+aG,sumYs[ipat][0]+sumYs[ipat][1]+rG))-logDen);
	  if (!R_finite(MHrate)) error("MHrate is in-finite");
	  if (runif(0.0,1.0) < exp(MHrate) ) 
	    {
	      aG = candi;
	      logDen = canlogDen;
	      acc[0]++;
	    }
	}
	  
      candi=rnorm(rG,can[1]);
      att[1]++;
      if (candi > min_aG && candi < max_aG )
	{
	  canlogDen = log(beta(aG,candi));
	  MHrate = 0;
	  for (ipat=0 ; ipat < N; ipat++ ) 
	    MHrate += log(betaYK(sizeip_s[ipat]+aG,sumYs[ipat][0]+sumYs[ipat][1]+candi))-canlogDen-
	      (log(betaYK(sizeip_s[ipat]+aG,sumYs[ipat][0]+sumYs[ipat][1]+rG))-logDen);
	  
	  if (runif(0.0,1.0) < exp(MHrate) ) 
	    {
	      rG = candi;
	      logDen = canlogDen;
	      acc[1]++;
	    }
	}
      Beta_ar=exp(logDen);

      // ======================
      REAL(logLlik)[iB] = logL;
      // SAVE the result from this iteration
      REAL(postaGs)[iB] = aG;
      REAL(postrGs)[iB] = rG;

      for (ibeta = 0; ibeta < p; ibeta++) 
	{
	  idh = ibeta + iB*p;
	  REAL(postbetas)[idh] = betas[ibeta];
	}
      
      for (ipat = 0; ipat < N; ipat++ ) 
	{
	  idh = ipat + iB*N;
	  REAL(postg1s)[idh] = g1s[ipat];
	}
      
      if (2*iB < burnin) //2*iB
	{
	  // use the first half of the burn-in period to adjust the proposal variance
	  for(ih=0; ih < Nacc ; ih++)
	    {
	      // Keept the acceptance to range between 0.3 and 0.5
	      // The proposal variance is adjusted only during the burn-in period
	      if( ih==0 ) can[ih]=FisherBeta(aG,rG,N)*3;  // update can of aG
	      else if( ih==1 ) can[ih]=FisherBeta(rG,aG,N)*3; // update can of rG
	      else 
	      {
		if (att[ih]>50)  can[ih]=adjustcan(acc[ih]/att[ih],can[ih],0.001,0.08);
	      }
	    }
	  /* else if (iB*(2*5) == burnin ) */
	  /*   { */
	  /*     for(ih=0; ih < Nacc; ih++){ */
	  /* 	att[ih] = 0.0; */
	  /* 	acc[ih] = 0.0; */
	  /*     } */
	  /*   } */
	}
      else if (iB==burnin)
	{
	  for( ih=0; ih < Nacc; ih++){
	    att[ih] = 0.0;
	    acc[ih] = 0.0;
	  }
	}  
      
      // compute the probability index
      if (probInd==1 && iB >= burnin && iB % thin == 0) 
	{ 
	  // If B = 105000, thin =10, burnin = 5000 then the first iB accepting the condition is iB = 5009
	  // Rprintf("\v iB %d idh %d",iB, idh);
	  // beta(aGs[ih],rGs[ih]) does not change with patients
	  for (ipat = 0 ; ipat < N; ipat++ )
	    {
	      idh = ipat + (iB-burnin)*N/thin;
	      if (sumYs[ipat][1]==0) REAL(probI)[idh]=1.0;
	      else {
		sizeij = 0;
		cansizeij = 0;
		for (ivec = 0 ; ivec < wIDsizes[ipat] ;ivec++ )
		  {
		    if(labelnp[wIDs[ipat][ivec]]==0) sizeij += getsizeij(ivec, X, wIDs[ipat], betas, Ntot,p);//sum_{j in pre} exp(Xij%*%beta)
		    if(labelnp[wIDs[ipat][ivec]]==1) cansizeij+=getsizeij(ivec, X, wIDs[ipat], betas, Ntot,p);//sum_{j in new} exp(Xij%*%beta)
		  }
		den = beta(sizeij+aG,sumYs[ipat][0]+rG);
		num = 0;
		for (ibeta = 0 ; ibeta < sumYs[ipat][1]; ibeta++)
		  {
		    num += choose(cansizeij+ibeta-1, ibeta)*beta(cansizeij+sizeij+aG,ibeta+sumYs[ipat][0]+rG);
		  }
		REAL(probI)[idh]=1-num/den;
	     //if (ipat == 0) Rprintf("ipat %d idh %d prob %f sumYs[ipat][1] %f num %f den %f w %f ",ipat, idh, REAL(probI)[idh],sumYs[ipat][1],num,den,w);
	      }
	    }
	}

      if(iB > 0 && iB % INTEGER(FreqEvery)[0] == 0) Rprintf("\v %d iterations are done...   ", iB);
    } // iB
  for (ih=0;ih < Nacc; ih ++) 
    {
      REAL(AR)[ih] = acc[ih]/att[ih];
      REAL(prp)[ih] = can[ih];
    }
  
  PutRNGstate();
  UNPROTECT(1);
  return res;
}


// // declarations (necessary for C++ codes)
// extern "C" {
//   SEXP Beta12(SEXP Y_,          // REAL
// 	      SEXP X,           // REAL 
// 	      SEXP ID,          // INTEGER
// 	      SEXP B_,          // INTEGER
// 	      SEXP maxni,       // INTEGER
// 	      SEXP nn,          // INTEGER
// 	      SEXP labelnp_,    // REAL
// 	      SEXP max_aG_,     // REAL
// 	      SEXP mu_beta_,    // REAL
// 	      SEXP sigma_beta_, // REAL
// 	      SEXP a_D_,         // REAL
// 	      SEXP ib_D_,        // REAL
// 	      SEXP M_,           // INTEGER
// 	      SEXP burnin_,     // INTEGER
// 	      SEXP FreqEvery    // INTEGER
// 	      ) ;
// }


// SEXP Beta12(SEXP Y_,          // REAL
// 	    SEXP X,           // REAL 
// 	    SEXP ID,          // INTEGER
// 	    SEXP B_,          // INTEGER
// 	    SEXP maxni,       // INTEGER
// 	    SEXP nn,          // INTEGER
// 	    SEXP labelnp_,    // REAL
// 	    SEXP max_aG_,     // REAL
// 	    SEXP mu_beta_,    // REAL
// 	    SEXP sigma_beta_, // REAL
// 	    SEXP a_D_,         // REAL
// 	    SEXP ib_D_,        // REAL
// 	    SEXP M_,           // INTEGER
// 	    SEXP burnin_,     // INTEGER
// 	    SEXP FreqEvery    // INTEGER
// 	    )
// {
//   // === Model ===
//   // Nonparametric const random effect over time
//   // Y_ij | Gi = gi ~ NB(size=exp(X_{ij}^T beta),prob=gi)
//   // gi | Hi = h ~ beta(aG_h,rG_h)
//   // beta ~ rnorm(mu_beta,sigma_beta)
//   // aG_h, rG_h ~ unif(0,max_aG)
//   // working well 
//   const double *Y = REAL(Y_); 
//   const double *labelnp = REAL(labelnp_);
//   const int Ntot = length(ID); // # patients
//   const int p = length(X)/Ntot; 
//   const int N = INTEGER(nn)[0];
//   const double max_aG = REAL(max_aG_)[0];
//   const int B = INTEGER(B_)[0];
//   const double *mu_betas = REAL(mu_beta_);
//   const double *sigma_betas = REAL(sigma_beta_);
//   const int burnin = INTEGER(burnin_)[0];
//   const double a_D = REAL(a_D_)[0];
//   const double ib_D = REAL(ib_D_)[0];
//   const int M = INTEGER(M_)[0];
//   int h1s[N],wID[INTEGER(maxni)[0]],wIDsize = 0, ih = 0, ip = 0, ipat = 0, ivec=0, idh=0, ibeta=0,iB=0;
//   // wIDs[ipat][ivec] contains the position of the ivec^th repeated measure of the ipat^th patient at labelnp and Y
//   //                  If ivec >= # repeated measures, wIDs[ipat][ivec] = -1000
//   // wIDsizes[ipat] contains the number of repeated measure of ipat^th patient
//   // sumYs[ipat][0] contains the sum of the CEL counts of ipat^th patients in the pre-scan period
//   // sumYs[ipat][1] contains the sum of the CEL counts of ipat^th patients in the new-scan period
//   int wIDs[N][INTEGER(maxni)[0]], wIDsizes[N];
//   double sumYs[N][2];
//   for (ipat = 0 ; ipat < N ; ipat++)
//     {
//       getwID( &wIDsize,wID,ID,ipat);
//       wIDsizes[ipat] = wIDsize;
//       sumYs[ipat][0] = 0;
//       sumYs[ipat][1] = 0;
//       for (ivec = 0 ; ivec < wIDsize; ivec++)
// 	{
// 	  wIDs[ipat][ivec]=wID[ivec];
// 	  if (labelnp[wIDs[ipat][ivec]]==0) sumYs[ipat][0] += Y[wIDs[ipat][ivec]];
// 	  if (labelnp[wIDs[ipat][ivec]]==1) sumYs[ipat][1] += Y[wIDs[ipat][ivec]];
// 	}
//       for (ivec = wIDsize ; ivec < INTEGER(maxni)[0]; ivec++)
// 	wIDs[ipat][ivec]=-1000;
//     }

//   // int h1s[N],h2s[N],js[N];
//   double beta[p],canbeta[p],aGs[M],rGs[M],g1s[N];
//   double vs[M],weightH1[M],postprob[M],D;
//   double sizeij=0.0,cansizeij=0.0,w=0.0,pij=0.0,
//     candi=0.0,MHrate=0.0,yij=0.0;

//   double att[2+p], acc[2+p],can[2+p];
//   for (ih = 0 ; ih < 2+p; ih ++)
//     {
//       att[ih] = 0.0;
//       acc[ih] = 0.0;
//       can[ih] = 2.0;
//     }

//   // http://stackoverflow.com/questions/8720550/how-to-return-array-of-structs-from-call-to-c-shared-library-in-r
//   SEXP res = PROTECT(allocVector(VECSXP, 10)); // result is stored in res
//   SEXP postaGs = allocVector(REALSXP, M*B); 
//   SET_VECTOR_ELT(res, 0, postaGs); 
//   SEXP postrGs = allocVector(REALSXP, M*B); 
//   SET_VECTOR_ELT(res, 1, postrGs); 
//   SEXP postvs = allocVector(REALSXP, M*B); 
//   SET_VECTOR_ELT(res, 2, postvs); 
//   SEXP postweightH1 = allocVector(REALSXP, M*B); 
//   SET_VECTOR_ELT(res, 3, postweightH1); 
//   SEXP posth1s = allocVector(INTSXP, N*B); 
//   SET_VECTOR_ELT(res, 4, posth1s); 
//   SEXP postg1s = allocVector(REALSXP, N*B); 
//   SET_VECTOR_ELT(res, 5, postg1s); 
//   SEXP postbetas = allocVector(REALSXP, p*B); 
//   SET_VECTOR_ELT(res, 6, postbetas); 
//   SEXP postD = allocVector(REALSXP, B); 
//   SET_VECTOR_ELT(res, 7, postD); 
//   SEXP AR = allocVector(REALSXP,2+p);
//   SET_VECTOR_ELT(res, 8, AR);
//   SEXP prp = allocVector(REALSXP,2+p);
//   SET_VECTOR_ELT(res, 9, prp);


//   GetRNGstate();
  
//   // ============ initialization ============ 
//   for (ip=0;ip < p ; ip++) beta[ip]=0.0;

//   for (ih = 0 ; ih < M ; ih++)
//     {
//       aGs[ih] = runif(0.0,max_aG);
//       rGs[ih] = runif(0.0,max_aG);
//     }
  
//   D = rgamma(a_D,1/ib_D);
  
//   for (ih = 0 ; ih < (M-1) ; ih ++) vs[ih] = rbeta(1,D);
//   vs[M-1]=1.0;
  
//   weightH1[0] = vs[0];
//   w = 1.0 ;
//   for (ih=1; ih < M; ih++)
//     {
//       w *= (1.0-vs[ih-1]);
//       weightH1[ih] = vs[ih]*w;
//     }


//   for ( ipat=0; ipat < N; ipat++) 
//     {
//       h1s[ipat] = sample(weightH1,M);
//       g1s[ipat] = rbeta(aGs[h1s[ipat]],rGs[h1s[ipat]]);
//       if (g1s[ipat] > 1.0-1.0e-7) g1s[ipat]=1.0-1.0e-7;
//       else if (g1s[ipat] < 1.0e-7) g1s[ipat]=1.0e-7;
//     }

//   // ================ MCMC =================
//   for (iB = 0 ; iB < B; iB++)
//     {
//       // STEP 1 
//       // Update g1s, a vector of length M+1. The first M terms are the random effects.
//       // the gs[M]=-1000, place holder for the patients with no new scans
//       for (ipat = 0 ; ipat < N; ipat++)
// 	{
// 	  w = 0.0;
// 	  for (ivec=0; ivec < wIDsizes[ipat]; ivec++ )
// 	    {
// 	      sizeij = getsizeij(ivec, X, wIDs[ipat], beta, Ntot, p);
// 	      w += sizeij*(1.0-labelnp[wIDs[ipat][ivec]]); // sum_{j in pre} exp(Xij%*%beta)
// 	    }
// 	  //Rprintf("\v ipat  %d wIDsizes[ipat] %d sumYs[0][ipat] %f sum exp(xij beta) %f",ipat, wIDsizes[ipat],sumYs[0][ipat], w);
// 	  g1s[ipat] = rbeta(aGs[h1s[ipat]]+w,rGs[h1s[ipat]]+sumYs[ipat][0]);
// 	  if (g1s[ipat] > 1.0-1.0e-7) g1s[ipat]=1.0-1.0e-7;
// 	  else if (g1s[ipat] < 1.0e-7) g1s[ipat]=1.0e-7;
// 	}

//       // STEP 2: 
//       // Update beta, a vector of length p, containing the coefficients b0, b1,...
//       // M-H 
//       // beta will be updated till the posterior categorical distribution of H1 contains at least one non zero probability.

//       for ( ibeta = 0; ibeta < p; ibeta++)
// 	{
// 	  for (ip = 0 ; ip < p ;ip ++ ) canbeta[ip]=beta[ip]; // copy beta
// 	  att[ibeta+2]++;

// 	  // proposal distribution
// 	  canbeta[ibeta] = rnorm(beta[ibeta], can[ibeta+2]); // the mean and sd must be double (not NumericVector)
	  
// 	  // prior
// 	  MHrate = dnorm(canbeta[ibeta], mu_betas[ibeta], sigma_betas[ibeta], 1) - dnorm(beta[ibeta], mu_betas[ibeta], sigma_betas[ibeta], 1);
	  
// 	  for (ipat=0; ipat < N; ipat++)
// 	    {

// 	      // if (printing)  Rprintf(" \v\v ipat: %d",ipat);
// 	      for ( ivec = 0; ivec < wIDsizes[ipat]; ivec++ )
// 		{
// 		  yij = Y[ wIDs[ipat][ivec] ];
// 		  pij = getpij2(g1s[ipat],0.0,labelnp[wIDs[ipat][ivec]]); 
// 		  //Rprintf("\v\v ipat %d ivec %d yij %f g1s[ipat] %f pij %f", 
// 		  //ipat,  ivec,   yij,   g1s[ipat],   pij);

// 		  cansizeij = getsizeij(ivec, X, wIDs[ipat], canbeta, Ntot, p); // the value of cansizeij is changed

// 		  sizeij = getsizeij(ivec, X, wIDs[ipat], beta, Ntot, p); // the value of cansizeij is changed

// 		  //Rprintf("canbeta %f cansizeij %f beta %f sizeij %f", canbeta[0],   cansizeij,   beta[0],   sizeij);

// 		  MHrate += dnbinom(yij, cansizeij,  pij, 1)-dnbinom(yij, sizeij,  pij, 1); 
// 		  //Rprintf("\v dnbinom(yij, cansizeij,  pij, 1) %f",dnbinom(yij, cansizeij,  pij, 1));
// 		  //Rprintf("\v dnbinom(yij, sizeij,  pij, 1) %f",dnbinom(yij, sizeij,  pij, 1));
// 		  //Rprintf("\v summed MHrate %f", MHrate);
// 		}
// 	      //Rprintf("\v");
// 	    }    
// 	  //Rprintf("\v beta[ibeta] %f canbeta[ibeta] %f",beta[ibeta],canbeta[ibeta]);
// 	  if(runif(0.0,1.0) < exp(MHrate))
// 	    {
// 	      beta[ibeta] = canbeta[ibeta];
// 	      acc[ibeta+2]++;
// 	    }
// 	}
//       // Step 3
//       // Update aG and rG, scalars for the base distribution of the random effect gi
//       // gi ~ Beta(aG,rG)
//       // aG ~ unif(0,max_aG)
//       // rG ~ unif(0,max_aG)
//       for (ih =0; ih < M; ih ++)
// 	{
// 	  candi=rnorm(aGs[ih],can[0]);
// 	  if (candi > 1.0e-7 && candi < max_aG )
// 	    {
// 	      if (weightH1[ih] > 0.1) att[0]++;
// 	      MHrate = 0;
// 	      for (ipat=0 ; ipat < N; ipat++ ) 
// 		if (h1s[ipat]==ih) MHrate += dbeta(g1s[ipat],candi,rGs[ih],1)-dbeta(g1s[ipat],aGs[ih],rGs[ih],1); 
// 	      if (!R_finite(MHrate)) error("MHrate is in-finite");
// 	      if (runif(0.0,1.0) < exp(MHrate) ) 
// 		{
// 		  aGs[ih] = candi;
// 		  if (weightH1[ih] > 0.1) acc[0]++;
// 		}
// 	    }
	  
// 	  candi=rnorm(rGs[ih],can[1]);
// 	  if (candi > 1.0e-7 && candi < max_aG )
// 	    {
// 	      if (weightH1[ih] > 0.1) att[1]++;
// 	      MHrate = 0;
// 	      for (ipat=0 ; ipat < N; ipat++ ) 
// 		if (h1s[ipat]==ih) MHrate += dbeta(g1s[ipat],aGs[ih],candi,1)-dbeta(g1s[ipat],aGs[ih],rGs[ih],1);   
// 	      if (runif(0.0,1.0) < exp(MHrate) ) 
// 		{
// 		  rGs[ih] = candi;
// 		  if (weightH1[ih] > 0.1) acc[1]++;
// 		}
// 	    }
// 	}
//       // Update vs, a vector of length M, containing the latent variable v, used to construct pis. 
//       // Once vs[ih] is updated, update the probabilities of the categorical distribution of H1 based on the formula:
//       // weightH1[0] = vs[0]
//       // weightH1[ih]= vs[ih]*(1-vh[ih-1])**(1-vh[ih-2])*...**(1-vh[0]) for ih =1,2,...M-1
//       w = 1.0;
//       for (ih = 0 ; ih < (M-1) ; ih ++) // vs[M-1] = 1 always!
// 	{
// 	  ip = 0;    // # patients with h1s == ih
// 	  ibeta = 0; // # patients with h1s > ih
// 	  for (ipat = 0 ; ipat < N; ipat++ ) 
// 	    {
// 	      if ( h1s[ipat]==ih ) ip++ ;
// 	      if ( h1s[ipat] > ih ) ibeta++ ;
// 	    }
	  
// 	  vs[ih] = rbeta(
// 			 (double) 1.0 +ip, 
// 			 D + (double)ibeta 
// 			 );
// 	  //Rprintf("ih %d sp1 %f sp2 %f",ih, (double) 1+ip+idh ,  D + (double)ibeta + (double)ivec);
// 	  if (ih > 0) w *= (1.0-vs[ih-1]);
// 	  weightH1[ih] = vs[ih]*w;
// 	}
//       weightH1[M-1] = w*(1-vs[M-2]);


//       //  Update h1s, a vector of length N, containing the cluster index of random effect 1 of each patient.
//       //  the cluster index is between 1 and M.
//       //  conjugate: posterior distribution of h1s has a categorical distribution with M categories. 
//       // if (printing && iB%printEvery==0)  Rprintf("\v\v Step 4");
//       for (ipat = 0 ; ipat < N ; ipat++)
// 	{
	  
// 	  // Given patient ipat, compute P(H_1i=ih;-) for each ih 
// 	  // P(H_1i=ih;-) = prod_{i=1}^{# repeated measure} dnbinom(yij; size=exp(b0+b1*x1i+..), prob=1/(gij+1))
	  
// 	  for (ih=0; ih < M; ih++ ) postprob[ih] = dbeta(g1s[ipat],aGs[ih],rGs[ih],0)*weightH1[ih];
// 	  h1s[ipat] = sample(postprob,M); 
	  
// 	}

//       // Update D
//       // D ~ gamma(scale=a.D,shape=b.D) = gamma(scale=a.D,shape=1/ib.D)  
//       w = 0.0;
//       for (ih=0;ih<(M-1);ih++) w += log(1.0-vs[ih]);
//       D = rgamma( (a_D+ (double) M-1.0), (1.0/(-w+ib_D)) ); // N+w is the number of random var following DP given js
//       if (ISNAN(D)){ 
// 	D = 1.0; 
// 	Rprintf("\v ISNAN(D): a_D %f ib_D %f sum log(1-vs[ih]) %f ",a_D,ib_D,w);
//       }
//       if( D <0.01) D = 0.01; // if D is too small, rbeta(1,a,0) returns NaN for any a
      
//       // SAVE the result from this iteration
//       REAL(postD)[iB]=D;
//       for (ih = 0; ih < M ; ih++)
// 	{
// 	  idh = ih + iB*M;
// 	  REAL(postaGs)[idh] = aGs[ih];
// 	  REAL(postrGs)[idh] = rGs[ih];
// 	  REAL(postvs)[idh] = vs[ih];
// 	  REAL(postweightH1)[idh] = weightH1[ih];
// 	}
//       for (ibeta = 0; ibeta < p; ibeta++) 
// 	{
// 	  idh = ibeta + iB*p;
// 	  REAL(postbetas)[idh] = beta[ibeta];
// 	}
      
//       for (ipat = 0; ipat < N; ipat++ ) 
// 	{
// 	  idh = ipat + iB*N;
// 	  REAL(postg1s)[idh] = g1s[ipat];
// 	  INTEGER(posth1s)[idh] = h1s[ipat];
// 	}
      
//       if (2*iB < burnin) //2*iB
// 	{
// 	  // use the first half of the burn-in period to adjust the proposal variance
// 	  for(ih=0; ih < 2+p ; ih++)
// 	    {
// 	      // Keept the acceptance to range between 0.3 and 0.5
// 	      // The proposal variance is adjusted only during the burn-in period
// 	      if( ih==0 ) 
// 		{
// 		  if (att[ih]>50) can[ih]=adjustcan(acc[ih]/att[ih],can[ih],0.05);
// 		}
// 	      else if( ih==1 )
// 		{
// 		  if (att[ih]>50) can[ih]=adjustcan(acc[ih]/att[ih],can[ih],0.05);
// 		}
// 	      else 
// 		{
// 		  if (att[ih]>50)  can[ih]=adjustcan(acc[ih]/att[ih],can[ih],0.1);
// 		}
// 	    }
// 	}
//       else if (iB==burnin)
// 	{
// 	  for( ih=0; ih < 2+p; ih++){
// 	    att[ih] = 0.0;
// 	    acc[ih] = 0.0;
// 	  }
// 	}  
      
//       if(iB > 0 && iB % INTEGER(FreqEvery)[0] == 0) Rprintf("\v %d iterations are done...   ", iB);
//     }
//   for (ih=0;ih < 2+p; ih ++) 
//     {
//       REAL(AR)[ih] = acc[ih]/att[ih];
//       REAL(prp)[ih] = can[ih];
//     }
  
//   PutRNGstate();
//   UNPROTECT(1);
//   return res;
// }


// // declarations (necessary for C++ codes)
// extern "C" {
//   SEXP Beta2(SEXP Y_,          // REAL
// 	     SEXP X,           // REAL 
// 	     SEXP ID,          // INTEGER
// 	     SEXP B_,          // INTEGER
// 	     SEXP maxni,       // INTEGER
// 	     SEXP nn,          // INTEGER
// 	     SEXP labelnp_,    // REAL
// 	     SEXP max_aG_,     // REAL
// 	     SEXP mu_beta_,    // REAL
// 	     SEXP sigma_beta_, // REAL
// 	     SEXP burnin_,     // INTEGER
// 	     SEXP FreqEvery    // INTEGER
// 	     );
// }

// SEXP Beta2(SEXP Y_,          // REAL
// 	   SEXP X,           // REAL 
// 	   SEXP ID,          // INTEGER
// 	   SEXP B_,          // INTEGER
// 	   SEXP maxni,       // INTEGER
// 	   SEXP nn,          // INTEGER
// 	   SEXP labelnp_,    // REAL
// 	   SEXP max_aG_,     // REAL
// 	   SEXP mu_beta_,    // REAL
// 	   SEXP sigma_beta_, // REAL
// 	   SEXP burnin_,     // INTEGER
// 	   SEXP FreqEvery    // INTEGER
// 	   )
// {
//   // === Model ===
//   // Y_ij | Gij = gi ~ NB(size=exp(X_{ij}^T beta),prob=gij)
//   // gij = g1s if j is in pre-scan 
//   //     = g2s if j is in old-scan 
//   // qG ~ beta(a_qG,r_qG)
//   // g1, g2 ~ beta(aG,rG)
//   // beta ~ rnorm(mu_beta,sigma_beta)
//   // aG, rG ~ unif(0,max_aG)
//   const double *Y = REAL(Y_); 
//   const double *labelnp = REAL(labelnp_);
//   const int Ntot = length(ID); // # patients
//   const int p = length(X)/Ntot; 
//   const int N = INTEGER(nn)[0];
//   const double max_aG = REAL(max_aG_)[0];
//   const int B = INTEGER(B_)[0];
//   const double *mu_betas = REAL(mu_beta_);
//   const double *sigma_betas = REAL(sigma_beta_);
//   const int burnin = INTEGER(burnin_)[0];
//   int wID[INTEGER(maxni)[0]],wIDsize = 0, ih = 0, ip = 0, ipat = 0, ivec=0, idh=0, ibeta=0,iB=0;
//   // wIDs[ipat][ivec] contains the position of the ivec^th repeated measure of the ipat^th patient at labelnp and Y
//   //                  If ivec >= # repeated measures, wIDs[ipat][ivec] = -1000
//   // wIDsizes[ipat] contains the number of repeated measure of ipat^th patient
//   // sumYs[ipat][0] contains the sum of the CEL counts of ipat^th patients in the pre-scan period
//   // sumYs[ipat][1] contains the sum of the CEL counts of ipat^th patients in the new-scan period
//   int wIDs[N][INTEGER(maxni)[0]], wIDsizes[N];
//   double sumYs[N][2];
//   for (ipat = 0 ; ipat < N ; ipat++)
//     {
//       getwID( &wIDsize,wID,ID,ipat);
//       wIDsizes[ipat] = wIDsize;
//       sumYs[ipat][0] = 0;
//       sumYs[ipat][1] = 0;
//       for (ivec = 0 ; ivec < wIDsize; ivec++)
// 	{
// 	  wIDs[ipat][ivec]=wID[ivec];
// 	  if (labelnp[wIDs[ipat][ivec]]==0) sumYs[ipat][0] += Y[wIDs[ipat][ivec]];
// 	  if (labelnp[wIDs[ipat][ivec]]==1) sumYs[ipat][1] += Y[wIDs[ipat][ivec]];
// 	}
//       for (ivec = wIDsize ; ivec < INTEGER(maxni)[0]; ivec++)
// 	wIDs[ipat][ivec]=-1000;
//     }

//   //h1s[N],h2s[N]
//   double beta[p],canbeta[p],aG,rG,qG,g1s[N],g2s[N];
//   // double vs[M],weightH1[M],postprob[M];
//   double sizeij=0.0,cansizeij=0.0,w=0.0,pij=0.0,
//     candi=0.0,MHrate=0.0,yij=0.0,sumgs=0.0;
//   // double  weipro=0.0,qG2=0.0,D=0.0;
//   double att[2+p], acc[2+p],can[2+p];
//   for (ih = 0 ; ih < 2+p; ih ++)
//     {
//       att[ih] = 0.0;
//       acc[ih] = 0.0;
//       can[ih] = 5.0;
//     }

//   // http://stackoverflow.com/questions/8720550/how-to-return-array-of-structs-from-call-to-c-shared-library-in-r
//   SEXP res = PROTECT(allocVector(VECSXP, 7)); // result is stored in res
//   SEXP postaGs = allocVector(REALSXP, B); 
//   SET_VECTOR_ELT(res, 0, postaGs); 
//   SEXP postrGs = allocVector(REALSXP, B); 
//   SET_VECTOR_ELT(res, 1, postrGs); 
 
//   SEXP postg1s = allocVector(REALSXP, N*B); 
//   SET_VECTOR_ELT(res, 2, postg1s); 
//   SEXP postg2s = allocVector(REALSXP, N*B); 
//   SET_VECTOR_ELT(res, 3, postg2s); 

//   SEXP postbetas = allocVector(REALSXP, p*B); 
//   SET_VECTOR_ELT(res, 4, postbetas); 
//   SEXP AR = allocVector(REALSXP,2+p);
//   SET_VECTOR_ELT(res, 5, AR);
//   SEXP prp = allocVector(REALSXP,2+p);
//   SET_VECTOR_ELT(res, 6, prp);

//   GetRNGstate();
  
//   // ============ initialization ============ 
//   for (ip=0;ip < p ; ip++) beta[ip]=0.0;
//   aG = 0.05;
//   rG = 0.05;
//   qG = 0.5;
//   for (ipat = 0 ; ipat < N ; ipat++)
//     {
//       g1s[ipat]=0.5;
//       g2s[ipat]=0.25;
//     }
//   // ================ MCMC =================
//   for (iB = 0 ; iB < B; iB++)
//     {
//       // STEP 1 
//       // Update g1s and g2s, a vector of length N. 
//       // The g2s are updated ASSUMING THAT ji = 0 for all i to be used in updating ji
//       for (ipat = 0 ; ipat < N; ipat++)
// 	{
// 	  w = 0.0;
// 	  sumgs = 0.0;
// 	  for (ivec=0; ivec < wIDsizes[ipat]; ivec ++)
// 	    {
// 	      sizeij = getsizeij(ivec, X, wIDs[ipat], beta, Ntot, p);
// 	      w += sizeij*(1.0-labelnp[wIDs[ipat][ivec]]); // sum_{j in pre} exp(Xij%*%beta)
// 	      sumgs += sizeij*labelnp[wIDs[ipat][ivec]]; // sum_{j in new} exp(Xij%*%beta)
// 	    }
// 	  g1s[ipat] = rbeta(aG+w,rG+sumYs[ipat][0]);
// 	  if (g1s[ipat] > 1.0-1.0e-5) g1s[ipat]= 1.0-1.0e-5;
	   
// 	  g2s[ipat] = rbeta(aG+sumgs,rG+sumYs[ipat][1]);
// 	  if (g2s[ipat] > 1.0-1.0e-5) g2s[ipat]= 1.0-1.0e-5;
// 	}


//       // STEP 2: 
//       // Update beta, a vector of length p, containing the coefficients b0, b1,...
//       // M-H 
//       // beta will be updated till the posterior categorical distribution of H1 contains at least one non zero probability.
//       for ( ibeta = 0; ibeta < p; ibeta++)
// 	{
// 	  for (ip = 0 ; ip < p ;ip ++ ) canbeta[ip]=beta[ip]; // copy beta
// 	  att[ibeta+2]++;

// 	  // proposal distribution
// 	  canbeta[ibeta] = rnorm(beta[ibeta], can[ibeta+2]); // the mean and sd must be double (not NumericVector)
	  
// 	  // prior
// 	  MHrate = dnorm(canbeta[ibeta], mu_betas[ibeta], sigma_betas[ibeta], 1) - dnorm(beta[ibeta], mu_betas[ibeta], sigma_betas[ibeta], 1);
	  
// 	  for (ipat=0; ipat < N; ipat++)
// 	    {
// 	      // if (printing)  Rprintf(" \v\v ipat: %d",ipat);
// 	      for ( ivec = 0; ivec < wIDsizes[ipat]; ivec++ )
// 		{
// 		  yij = Y[ wIDs[ipat][ivec] ];
// 		  // yij is a vector b/c it is the first arg of dnbinom
// 		  pij = getpij2(g1s[ipat],
// 				g2s[ipat],
// 				labelnp[wIDs[ipat][ivec]]); 

// 		  cansizeij = getsizeij(ivec, X, wIDs[ipat], canbeta, Ntot, p); // the value of cansizeij is changed

// 		  sizeij = getsizeij(ivec, X, wIDs[ipat], beta, Ntot, p); // the value of cansizeij is changed

// 		  //Rprintf("canbeta %f cansizeij %f beta %f sizeij %f", canbeta[0],   cansizeij,   beta[0],   sizeij);

// 		  MHrate += dnbinom(yij, cansizeij,  pij, 1)-dnbinom(yij, sizeij,  pij, 1); 
// 		  //Rprintf("\v dnbinom(yij, cansizeij,  pij, 1) %f",dnbinom(yij, cansizeij,  pij, 1));
// 		  //Rprintf("\v dnbinom(yij, sizeij,  pij, 1) %f",dnbinom(yij, sizeij,  pij, 1));
// 		  //Rprintf("\v summed MHrate %f", MHrate);
// 		}
// 	      //Rprintf("\v");
// 	    }    

// 	  if(runif(0.0,1.0) < exp(MHrate))
// 	    {
// 	      beta[ibeta] = canbeta[ibeta];
// 	      acc[ibeta+2]++;
// 	    }
// 	}


      
//       // Update aG and rG, scalars for the base distribution of the random effect gi
//       // gi ~ Beta(aG,rG)
//       // aG ~ unif(0,max_aG)
//       // rG ~ unif(0,max_aG)
//       candi=rnorm(aG,can[0]);
//       if (candi > 1.0e-5 && candi < max_aG )
// 	{
// 	  att[0]++;
// 	  MHrate = 0;
// 	  for (ipat=0 ; ipat < N; ipat ++ ) 
// 	    {
// 	      MHrate += dbeta(g1s[ipat],candi,rG,1)-dbeta(g1s[ipat],aG,rG,1); 
// 	      MHrate += dbeta(g2s[ipat],candi,rG,1)-dbeta(g2s[ipat],aG,rG,1); 
// 	    }
// 	  if (runif(0.0,1.0) < exp(MHrate) ) 
// 	    {
// 	      aG = candi;
// 	      acc[0]++;
// 	    }
// 	}
	
      
//       candi=rnorm(rG,can[1]);
//       if (candi > 1.0e-5 && candi < max_aG )
// 	{
// 	  att[1]++;
// 	  MHrate = 0;
// 	  for (ipat=0 ; ipat < N; ipat ++ ) 
// 	    {
// 	      MHrate += dbeta(g1s[ipat],aG,candi,1)-dbeta(g1s[ipat],aG,rG,1);   
// 	      MHrate += dbeta(g2s[ipat],aG,candi,1)-dbeta(g2s[ipat],aG,rG,1); 
// 	    }
// 	  if (runif(0.0,1.0) < exp(MHrate) ) 
// 	    {
// 	      rG = candi;
// 	      acc[1]++;
// 	    }
// 	}

    
//       // http://stackoverflow.com/questions/8720550/how-to-return-array-of-structs-from-call-to-c-shared-library-in-r
//       // SAVE the result from this iteration
//       REAL(postaGs)[iB] = aG;
//       REAL(postrGs)[iB] = rG;

//       for (ibeta = 0; ibeta < p; ibeta++) 
// 	{
// 	  idh = ibeta + iB*p;
// 	  REAL(postbetas)[idh] = beta[ibeta];
// 	}
      
//       for (ipat = 0; ipat < N; ipat++ ) 
// 	{
// 	  idh = ipat + iB*N;
// 	  REAL(postg1s)[idh] = g1s[ipat];
// 	  REAL(postg2s)[idh] = g2s[ipat];
// 	}

//       if (iB < burnin ) //2*iB
// 	{
// 	  for(ih=0; ih < 2+p ; ih++)
// 	    {
// 	      // Keept the acceptance to range between 0.3 and 0.5
// 	      // The proposal variance is adjusted only during the burn-in period
// 	      if( ih==0 ) can[ih]=FisherBeta(aG,rG,N)*3; 
// 	      else if( ih==1 ) can[ih]=FisherBeta(rG,aG,N)*3; 
// 	      else 
// 		{
// 		  if (att[ih]>50)  can[ih]=adjustcan(acc[ih]/att[ih],can[ih],0.1);
// 		}
// 	    }
	 
// 	  /* if (iB % 100==0){ */
// 	  /*   for(ih=0; ih < 2+p; ih++){ */
// 	  /*     att[ih] = 0.0; */
// 	  /*     acc[ih] = 0.0; */
// 	  /*   } */
// 	}else if (iB==burnin){
// 	for( ih=0; ih < 2+p; ih++){
// 	  att[ih] = 0.0;
// 	  acc[ih] = 0.0;
// 	}
//       }  
      
//       if(iB > 0 && iB % INTEGER(FreqEvery)[0] == 0) Rprintf("\v %d iterations are done...   ", iB);
//     }
//   for (ih=0;ih < 2+p; ih ++) 
//     {
//       REAL(AR)[ih] = acc[ih]/att[ih];
//       REAL(prp)[ih] = can[ih];
//     }
  
//   PutRNGstate();
//   UNPROTECT(1);
//   return res;
// }



// extern "C" {
//   SEXP Beta22(SEXP Y_,          // REAL
// 	      SEXP X,           // REAL 
// 	      SEXP ID,          // INTEGER
// 	      SEXP B_,          // INTEGER
// 	      SEXP maxni,       // INTEGER
// 	      SEXP nn,          // INTEGER
// 	      SEXP labelnp_,    // REAL
// 	      SEXP max_aG_,     // REAL
// 	      SEXP mu_beta_,    // REAL
// 	      SEXP sigma_beta_, // REAL
// 	      SEXP a_D_,         // REAL
// 	      SEXP ib_D_,        // REAL
// 	      SEXP M_,           // INTEGER
// 	      SEXP burnin_,     // INTEGER
// 	      SEXP FreqEvery    // INTEGER
// 	      );
// }

SEXP gibbsREchange(SEXP Y_,          // REAL
		   SEXP X,           // REAL 
		   SEXP ID,          // INTEGER
		   SEXP B_,          // INTEGER
		   SEXP maxni,       // INTEGER
		   SEXP nn,          // INTEGER
		   SEXP labelnp_,    // REAL
		   SEXP max_aG_,     // REAL
		   SEXP mu_beta_,    // REAL
		   SEXP evalue_sigma_betas_,
		   SEXP Inv_sigma_betas_,
		   SEXP a_D_,         // REAL
		   SEXP ib_D_,        // REAL
		   SEXP M_,           // INTEGER
		   SEXP burnin_,     // INTEGER
		   SEXP FreqEvery    // INTEGER
		   )
{
  // === Model ===
  // Y_ij | Gi = gi ~ NB(size=exp(X_{ij}^T beta),prob=gi)
  // gi | Hi = h ~ beta(aG_h,rG_h)
  // beta ~ rnorm(mu_beta,sigma_beta)
  // aG_h, rG_h ~ unif(0,max_aG)
  const double *Y = REAL(Y_); 
  const double *labelnp = REAL(labelnp_);
  const int Ntot = length(ID); // # patients
  const int p = length(X)/Ntot; 
  const int N = INTEGER(nn)[0];
  const double max_aG = REAL(max_aG_)[0];
  const int B = INTEGER(B_)[0];
  const double *mu_betas = REAL(mu_beta_);
  const double *evalue_sigma_betas=REAL(evalue_sigma_betas_);
  const double *Inv_sigma_beta=REAL(Inv_sigma_betas_);
  const int burnin = INTEGER(burnin_)[0];
  const double a_D = REAL(a_D_)[0];
  const double ib_D = REAL(ib_D_)[0];
  const int M = INTEGER(M_)[0];
  int h1s[N],h2s[N],wID[INTEGER(maxni)[0]],wIDsize = 0, ih = 0, ip = 0, ipat = 0, ivec=0, idh=0, ibeta=0,iB=0;
  const double min_aG = 0.5;
  // wIDs[ipat][ivec] contains the position of the ivec^th repeated measure of the ipat^th patient at labelnp and Y
  //                  If ivec >= # repeated measures, wIDs[ipat][ivec] = -1000
  // wIDsizes[ipat] contains the number of repeated measure of ipat^th patient
  // sumYs[ipat][0] contains the sum of the CEL counts of ipat^th patients in the pre-scan period
  // sumYs[ipat][1] contains the sum of the CEL counts of ipat^th patients in the new-scan period
  int wIDs[N][INTEGER(maxni)[0]], wIDsizes[N];
  double sumYs[N][2];
  for (ipat = 0 ; ipat < N ; ipat++)
    {
      getwID( &wIDsize,wID,ID,ipat);
      wIDsizes[ipat] = wIDsize;
      sumYs[ipat][0] = 0;
      sumYs[ipat][1] = 0;
      for (ivec = 0 ; ivec < wIDsize; ivec++)
	{
	  wIDs[ipat][ivec]=wID[ivec];
	  if (labelnp[wIDs[ipat][ivec]]==0) sumYs[ipat][0] += Y[wIDs[ipat][ivec]];
	  if (labelnp[wIDs[ipat][ivec]]==1) sumYs[ipat][1] += Y[wIDs[ipat][ivec]];
	}
      for (ivec = wIDsize ; ivec < INTEGER(maxni)[0]; ivec++)
	wIDs[ipat][ivec]=-1000;
    }

  // int h1s[N],h2s[N],js[N];
  double betas[p],canbetas[p],aGs[M],rGs[M],g1s[N],g2s[N];
  double vs[M],weightH1[M],postprob[M],D;
  double sizeij=0.0,cansizeij=0.0,w=0.0,pij=0.0,
    candi=0.0,MHrate=0.0,yij=0.0, logLcan=0,logL=0;
  int Nacc = 4;
  double att[Nacc], acc[Nacc],can[Nacc];
  for (ih = 0 ; ih < Nacc; ih ++)
    {
      att[ih] = 0.0;
      acc[ih] = 0.0;
      can[ih] = 2.0;
    }

  // http://stackoverflow.com/questions/8720550/how-to-return-array-of-structs-from-call-to-c-shared-library-in-r
  SEXP res = PROTECT(allocVector(VECSXP, 13)); // result is stored in res
  SEXP postaGs = allocVector(REALSXP, M*B); 
  SET_VECTOR_ELT(res, 0, postaGs); 
  SEXP postrGs = allocVector(REALSXP, M*B); 
  SET_VECTOR_ELT(res, 1, postrGs); 

  SEXP postvs = allocVector(REALSXP, M*B); 
  SET_VECTOR_ELT(res, 2, postvs); 
  SEXP postweightH1 = allocVector(REALSXP, M*B); 
  SET_VECTOR_ELT(res, 3, postweightH1); 

  SEXP posth1s = allocVector(INTSXP, N*B); 
  SET_VECTOR_ELT(res, 4, posth1s); 
  SEXP posth2s = allocVector(INTSXP, N*B); 
  SET_VECTOR_ELT(res, 5, posth2s);
 
  SEXP postg1s = allocVector(REALSXP, N*B); 
  SET_VECTOR_ELT(res, 6, postg1s); 
  SEXP postg2s = allocVector(REALSXP, N*B); 
  SET_VECTOR_ELT(res, 7, postg2s); 

  SEXP postbetas = allocVector(REALSXP, p*B); 
  SET_VECTOR_ELT(res, 8, postbetas);
 
  SEXP logLlik = allocVector(REALSXP, B); 
  SET_VECTOR_ELT(res, 9, logLlik);
  SEXP postD = allocVector(REALSXP, B); 
  SET_VECTOR_ELT(res, 10, postD); 
  SEXP AR = allocVector(REALSXP,Nacc);
  SET_VECTOR_ELT(res, 11, AR);
  SEXP prp = allocVector(REALSXP,Nacc);
  SET_VECTOR_ELT(res, 12, prp);
  
  

  GetRNGstate();
  
  // ============ initialization ============ 
  for (ip=0;ip < p ; ip++) betas[ip]=0.0;

  for (ih = 0 ; ih < M ; ih++)
    {
      aGs[ih] = runif(min_aG,max_aG);
      rGs[ih] = runif(min_aG,max_aG);
    }
  
  D = rgamma(a_D,1/ib_D);
  
  for (ih = 0 ; ih < (M-1) ; ih ++) vs[ih] = rbeta(1,D);
  vs[M-1]=1.0;
  
  weightH1[0] = vs[0];
  w = 1.0 ;
  for (ih=1; ih < M; ih++)
    {
      w *= (1.0-vs[ih-1]);
      weightH1[ih] = vs[ih]*w;
    }


  for ( ipat=0; ipat < N; ipat++) 
    {
      h1s[ipat] = sample(weightH1,M);
      h2s[ipat] = h1s[ipat];
      g1s[ipat] = rbeta(aGs[h1s[ipat]],rGs[h1s[ipat]]);
      g2s[ipat] = g1s[ipat];
      if (g1s[ipat] > 1.0-1.0e-7) g1s[ipat]=1.0-1.0e-7;
      else if (g1s[ipat] < 1.0e-7) g1s[ipat]=1.0e-7;
    }

  // ================ MCMC =================
  for (iB = 0 ; iB < B; iB++)
    {
      R_CheckUserInterrupt(); 
      //R_ProcessEvents();
      // STEP 1 
      // Update g1s and g2s, two vectors of length M. The first M terms are the random effects.
      // the gs[M]=-1000, place holder for the patients with no new scans
      for (ipat = 0 ; ipat < N; ipat++)
	{
	  w = 0.0;
	  pij = 0.0;
	  for (ivec=0; ivec < wIDsizes[ipat]; ivec++ )
	    {
	      sizeij = getsizeij(ivec, X, wIDs[ipat], betas, Ntot, p);
	      w += sizeij*(1.0-labelnp[wIDs[ipat][ivec]]); // sum_{j in pre} exp(Xij%*%betas)
	      pij += sizeij*labelnp[wIDs[ipat][ivec]]; // sum_{j in new} exp(Xij%*%betas)
	    }
	  //Rprintf("\v ipat  %d wIDsizes[ipat] %d sumYs[0][ipat] %f sum exp(xij betas) %f",ipat, wIDsizes[ipat],sumYs[0][ipat], w);
	  g1s[ipat] = rbeta(aGs[h1s[ipat]]+w,rGs[h1s[ipat]]+sumYs[ipat][0]);
	  g2s[ipat] = rbeta(aGs[h2s[ipat]]+pij,rGs[h2s[ipat]]+sumYs[ipat][1]);
	  if (g1s[ipat] > 1.0-1.0e-7) g1s[ipat]=1.0-1.0e-7;
	  else if (g1s[ipat] < 1.0e-7) g1s[ipat]=1.0e-7;
	  if (g2s[ipat] > 1.0-1.0e-7) g2s[ipat]=1.0-1.0e-7;
	  else if (g2s[ipat] < 1.0e-7) g2s[ipat]=1.0e-7;
	}


      // ==================
      att[2]++;
      // proposal distribution
      for (ibeta = 0; ibeta <p ; ibeta ++)
	canbetas[ibeta] = rnorm(betas[ibeta], can[2]); // the mean and sd must be double (not NumericVector)
      
      // prior
      MHrate = dmvnorm(canbetas, mu_betas,p, evalue_sigma_betas,Inv_sigma_beta,1)
	- dmvnorm(betas, mu_betas,p, evalue_sigma_betas,Inv_sigma_beta,1);
      logLcan = 0.0;
      logL = 0.0;
      for (ipat=0; ipat < N; ipat++)
	{
	  // if (printing)  Rprintf(" \v\v ipat: %d",ipat);
	  for ( ivec = 0; ivec < wIDsizes[ipat]; ivec++ )
	    {
	      yij = Y[ wIDs[ipat][ivec] ];
	      // yij is a vector b/c it is the first arg of dnbinom
	      pij = getpij2(g1s[ipat],g2s[ipat],labelnp[wIDs[ipat][ivec]]);
	      //Rprintf("\v\v ipat %d js[ipat] %f ivec %d yij %f g1s[ipat] %f g2s[ipat] %f pij %f", 
	      //ipat,  js[ipat] , ivec,   yij,   g1s[ipat], g2s[ipat],  pij);
	      cansizeij = getsizeij(ivec, X, wIDs[ipat], canbetas, Ntot, p); // the value of cansizeij is changed
	      sizeij = getsizeij(ivec, X, wIDs[ipat], betas, Ntot, p); // the value of cansizeij is changed

	      logLcan += dnbinom(yij, cansizeij,  pij, 1);
	      logL += dnbinom(yij, sizeij,  pij, 1);
	    }
	  //Rprintf("\v");
	}    
      MHrate += logLcan-logL;
      if(runif(0.0,1.0) < exp(MHrate))
	{
	  for (ibeta = 0; ibeta < p;ibeta++) betas[ibeta] = canbetas[ibeta];
	  acc[2]++;
	  logL = logLcan;
	}

      // ==================

      // Step 3
      // Update aG and rG, scalars for the base distribution of the random effect gi
      // gi ~ Beta(aG,rG)
      // aG ~ unif(0,max_aG)
      // rG ~ unif(0,max_aG)
      for (ih =0; ih < M; ih ++)
	{
	  candi=rnorm(aGs[ih],can[0]);
	  if (candi > min_aG && candi < max_aG )
	    {
	      if (weightH1[ih] > 0.1) att[0]++;
	      MHrate = 0;
	      for (ipat=0 ; ipat < N; ipat++ ) 
		{
		  if (h1s[ipat]==ih) MHrate += dbeta(g1s[ipat],candi,rGs[ih],1)-dbeta(g1s[ipat],aGs[ih],rGs[ih],1);
		  if (h2s[ipat]==ih) MHrate += dbeta(g2s[ipat],candi,rGs[ih],1)-dbeta(g2s[ipat],aGs[ih],rGs[ih],1);
		}
	      if (!R_finite(MHrate)) error("MHrate is in-finite");
	      if (runif(0.0,1.0) < exp(MHrate) ) 
		{
		  aGs[ih] = candi;
		  if (weightH1[ih] > 0.1) acc[0]++;
		}
	    }
	  
	  candi=rnorm(rGs[ih],can[1]);
	  if (candi > min_aG && candi < max_aG )
	    {
	      if (weightH1[ih] > 0.1) att[1]++;
	      MHrate = 0;
	      for (ipat=0 ; ipat < N; ipat++ )  
		{
		  if (h1s[ipat]==ih) MHrate += dbeta(g1s[ipat],aGs[ih],candi,1)-dbeta(g1s[ipat],aGs[ih],rGs[ih],1);   
		  if (h2s[ipat]==ih) MHrate += dbeta(g2s[ipat],aGs[ih],candi,1)-dbeta(g2s[ipat],aGs[ih],rGs[ih],1);
		}
	      if (runif(0.0,1.0) < exp(MHrate) ) 
		{
		  rGs[ih] = candi;
		  if (weightH1[ih] > 0.1) acc[1]++;
		}
	    }
	}
      // Update vs, a vector of length M, containing the latent variable v, used to construct pis. 
      // Once vs[ih] is updated, update the probabilities of the categorical distribution of H1 based on the formula:
      // weightH1[0] = vs[0]
      // weightH1[ih]= vs[ih]*(1-vh[ih-1])**(1-vh[ih-2])*...**(1-vh[0]) for ih =1,2,...M-1
      w = 1.0;
      for (ih = 0 ; ih < (M-1) ; ih ++) // vs[M-1] = 1 always!
	{
	  ip = 0;    // # patients with hs == ih
	  idh = 0;   // # patients with h2s == ih and js==0
	  ibeta = 0; // # patients with hs > ih
	  ivec = 0;  // # patients with h2s > ih and js==0
	  for (ipat = 0 ; ipat < N; ipat++ ) 
	    {
	      if ( h1s[ipat]==ih ) ip++ ;
	      if ( h2s[ipat]==ih )  idh++; // h2s[ipat]=M+1 is not counted because they correspond to js[ipat]=-1000
	      if ( h1s[ipat] > ih ) ibeta++ ;
	      if ( h2s[ipat] > ih )  ivec++;
	    }
	  
	  vs[ih] = rbeta(
			 (double) 1+ip+idh, 
			 D + (double)ibeta + (double)ivec
			 );
	  //Rprintf("ih %d sp1 %f sp2 %f",ih, (double) 1+ip+idh ,  D + (double)ibeta + (double)ivec);
	  if (ih > 0) w *= (1.0-vs[ih-1]);
	  weightH1[ih] = vs[ih]*w;
	}
      weightH1[M-1] = w*(1-vs[M-2]);
      // weightH1[M-1]= vs[M-1]*(1-vs)[M-2]*(1-vs)[M-3]*...*(1-vs)[0]
      //              =       1*(1-vs)[M-2]*(1-vs)[M-3]*...*(1-vs)[0]

      //  Update h1s and h2s, two vectors of length N, containing the cluster index of random effect 1 of each patient.
      //  the cluster index is between 1 and M.
      //  conjugate: posterior distribution of h1s has a categorical distribution with M categories. 
      // if (printing && iB%printEvery==0)  Rprintf("\v\v Step 4");
      for (ipat = 0 ; ipat < N ; ipat++)
	{
	  
	  // Given patient ipat, compute P(H_1i=ih;-) for each ih 
	  // P(H_1i=ih;-) = prod_{i=1}^{# repeated measure} dnbinom(yij; size=exp(b0+b1*x1i+..), prob=1/(gij+1))
	  for (ih=0; ih < M; ih++ ) postprob[ih] = dbeta(g1s[ipat],aGs[ih],rGs[ih],0)*weightH1[ih];
	  h1s[ipat] = sample(postprob,M); 
	  
	  for (ih=0; ih < M; ih++ ) postprob[ih] = dbeta(g2s[ipat],aGs[ih],rGs[ih],0)*weightH1[ih];
	  h2s[ipat] = sample(postprob,M); 
	}

      // Update D
      // D ~ runif(0.01,maxD) //gamma(scale=a.D,shape=b.D) = gamma(scale=a.D,shape=1/ib.D)  
      candi = rnorm(D,can[3]);
      att[3]++;
      if (a_D < candi && candi < ib_D)
	{
	  MHrate = 0.0;
	  for (ih = 0 ; ih < (M-1);ih++)
	    {
	      MHrate += dbeta(vs[ih],1.0,candi,1)-dbeta(vs[ih],1.0,D,1);
	    }
	  if (runif(0.0,1.0) < exp(MHrate) ) 
	    {
	      D = candi;
	      acc[3]++;
	    }
	}
      
      // SAVE the result from this iteration
      REAL(postD)[iB]=D;
      REAL(logLlik)[iB] = logL;

      for (ih = 0; ih < M ; ih++)
	{
	  idh = ih + iB*M;
	  REAL(postaGs)[idh] = aGs[ih];
	  REAL(postrGs)[idh] = rGs[ih];
	  REAL(postvs)[idh] = vs[ih];
	  REAL(postweightH1)[idh] = weightH1[ih];
	}
      for (ibeta = 0; ibeta < p; ibeta++) 
	{
	  idh = ibeta + iB*p;
	  REAL(postbetas)[idh] = betas[ibeta];
	}
      
      for (ipat = 0; ipat < N; ipat++ ) 
	{
	  idh = ipat + iB*N;
	  REAL(postg1s)[idh] = g1s[ipat];
	  INTEGER(posth1s)[idh] = h1s[ipat];
	  REAL(postg2s)[idh] = g2s[ipat];
	  INTEGER(posth2s)[idh] = h2s[ipat];
	}
      
      if (2*iB < burnin) //2*iB
	{
	  // use the first half of the burn-in period to adjust the proposal variance
	  for(ih=0; ih < Nacc ; ih++)
	    {
	      // Keept the acceptance to range between 0.3 and 0.5
	      // The proposal variance is adjusted only during the burn-in period
	      if( ih==0 ) 
		{
		  if (att[ih]>50) can[ih]=adjustcan(acc[ih]/att[ih],can[ih],0.05,2);
		}
	      else if( ih==1 )
		{
		  if (att[ih]>50) can[ih]=adjustcan(acc[ih]/att[ih],can[ih],0.05,2);
		}
	      else 
		{
		  if (att[ih]>50)  can[ih]=adjustcan(acc[ih]/att[ih],can[ih],0.1,2);
		}
	    }
	}
      else if (iB==burnin)
	{
	  for( ih=0; ih < Nacc; ih++){
	    att[ih] = 0.0;
	    acc[ih] = 0.0;
	  }
	}  
      
      if(iB > 0 && iB % INTEGER(FreqEvery)[0] == 0) Rprintf("\v %d iterations are done...   ", iB);
    }
  for (ih=0;ih < Nacc; ih ++) 
    {
      REAL(AR)[ih] = acc[ih]/att[ih];
      REAL(prp)[ih] = can[ih];
    }
  
  PutRNGstate();
  UNPROTECT(1);
  return res;
}




// extern "C" {
//   SEXP Beta13(SEXP Y_,          // REAL
// 	      SEXP X,           // REAL 
// 	      SEXP ID,          // INTEGER
// 	      SEXP B_,          // INTEGER
// 	      SEXP maxni,       // INTEGER
// 	      SEXP nn,              // INTEGER
// 	      SEXP labelnp_,        // REAL
// 	      SEXP min_prec_eta_,  // REAL
// 	      SEXP max_prec_eta_,  // REAL
// 	      SEXP max_aG_,     // REAL
// 	      SEXP mu_beta_,    // REAL
// 	      SEXP sigma_beta_, // REAL
// 	      SEXP burnin_,     // INTEGER
// 	      SEXP FreqEvery    // INTEGER
// 	      );
// }


// SEXP Beta13(SEXP Y_,          // REAL
// 	    SEXP X,           // REAL 
// 	    SEXP ID,          // INTEGER
// 	    SEXP B_,          // INTEGER
// 	    SEXP maxni,       // INTEGER
// 	    SEXP nn,              // INTEGER
// 	    SEXP labelnp_,        // REAL
// 	    SEXP min_prec_eta_,  // REAL
// 	    SEXP max_prec_eta_,  // REAL
// 	    SEXP max_aG_,     // REAL
// 	    SEXP mu_beta_,    // REAL
// 	    SEXP sigma_beta_, // REAL
// 	    SEXP burnin_,     // INTEGER
// 	    SEXP FreqEvery    // INTEGER
// 	    )
// {
//   // === Model ===
//   // Y_ij | Gij = gi ~ NB(size=exp(X_{ij}^T beta),prob=gij)
//   // gij = g1                     if j is in pre-scan 
//   //     = g1/[(1-g1)exp(-eta)+g1] if j is in old-scan
//   // That is: logit(g2) = logit(g1)+eta where eta \in (-inf,inf)
//   // g1 < g2 <=> 0 < eta 
//   // g1 ~ beta(aG,rG)
//   // eta ~ norm(0.0,sigma_eta)
//   // beta ~ rnorm(mu_beta,sigma_beta)
//   // aG, rG ~ unif(0,max_aG)
//   // Currently, this model is very sensible to the choice of the min value that sigma_eta takes
//   const double *Y = REAL(Y_); 
//   const double *labelnp = REAL(labelnp_);
//   const double min_prec_eta=REAL(min_prec_eta_)[0];
//   const double max_prec_eta=REAL(max_prec_eta_)[0];
//   const int Ntot = length(ID); // # patients
//   const int p = length(X)/Ntot; 
//   const int N = INTEGER(nn)[0];
//   const double max_aG = REAL(max_aG_)[0];
//   const int B = INTEGER(B_)[0];
//   const double *mu_betas = REAL(mu_beta_);
//   const double *sigma_betas = REAL(sigma_beta_);
//   const int burnin = INTEGER(burnin_)[0];
//   int wID[INTEGER(maxni)[0]],wIDsize = 0, ih = 0, ip = 0, ipat = 0, ivec=0, idh=0, ibeta=0,iB=0;
//   // wIDs[ipat][ivec] contains the position of the ivec^th repeated measure of the ipat^th patient at labelnp and Y
//   //                  If ivec >= # repeated measures, wIDs[ipat][ivec] = -1000
//   // wIDsizes[ipat] contains the number of repeated measure of ipat^th patient
//   // sumYs[ipat][0] contains the sum of the CEL counts of ipat^th patients in the pre-scan period
//   // sumYs[ipat][1] contains the sum of the CEL counts of ipat^th patients in the new-scan period
//   int wIDs[N][INTEGER(maxni)[0]], wIDsizes[N];
//   double sumYs[N][2];
//   for (ipat = 0 ; ipat < N ; ipat++)
//     {
//       getwID( &wIDsize,wID,ID,ipat);
//       wIDsizes[ipat] = wIDsize;
//       sumYs[ipat][0] = 0;
//       sumYs[ipat][1] = 0;
//       for (ivec = 0 ; ivec < wIDsize; ivec++)
// 	{
// 	  wIDs[ipat][ivec]=wID[ivec];
// 	  if (labelnp[wIDs[ipat][ivec]]==0) sumYs[ipat][0] += Y[wIDs[ipat][ivec]];
// 	  if (labelnp[wIDs[ipat][ivec]]==1) sumYs[ipat][1] += Y[wIDs[ipat][ivec]];
// 	}
//       for (ivec = wIDsize ; ivec < INTEGER(maxni)[0]; ivec++)
// 	wIDs[ipat][ivec]=-1000;
//     }

//   //h1s[N],h2s[N]
//   double beta[p],canbeta[p],aG,rG,qG,g1s[N],etas[N];
//   double prec_eta = (min_prec_eta+max_prec_eta)/2;
//   // double vs[M],weightH1[M],postprob[M];
//   double sizeij=0.0,cansizeij=0.0,w=0.0,canpij=0.0,pij=0.0,
//     candi=0.0,MHrate=0.0,yij=0.0,sumgs=0.0;
//   // double  weipro=0.0,qG2=0.0,D=0.0;
//   double att[5+p], acc[5+p],can[5+p];
//   for (ih = 0 ; ih < 5+p; ih ++)
//     {
//       att[ih] = 0.0;
//       acc[ih] = 0.0;
//       can[ih] = 2.0;
//     }

  
//   // http://stackoverflow.com/questions/8720550/how-to-return-array-of-structs-from-call-to-c-shared-library-in-r
//   SEXP res = PROTECT(allocVector(VECSXP, 9)); // result is stored in res
//   SEXP postaGs = allocVector(REALSXP, B); 
//   SET_VECTOR_ELT(res, 0, postaGs); 
//   SEXP postrGs = allocVector(REALSXP, B); 
//   SET_VECTOR_ELT(res, 1, postrGs); 
//   SEXP post_prec_eta= allocVector(REALSXP, B); 
//   SET_VECTOR_ELT(res, 2, post_prec_eta);
//   SEXP postg1s = allocVector(REALSXP, N*B); 
//   SET_VECTOR_ELT(res, 3, postg1s); 
//   SEXP postetas = allocVector(REALSXP, N*B); 
//   SET_VECTOR_ELT(res, 4, postetas); 
//   SEXP postg2s = allocVector(REALSXP, N*B); 
//   SET_VECTOR_ELT(res, 5, postg2s); 

//   SEXP postbetas = allocVector(REALSXP, p*B); 
//   SET_VECTOR_ELT(res, 6, postbetas); 
//   SEXP AR = allocVector(REALSXP,5+p);
//   SET_VECTOR_ELT(res, 7, AR);
//   SEXP prp = allocVector(REALSXP,5+p);
//   SET_VECTOR_ELT(res, 8, prp);

//   GetRNGstate();
  
//   // ============ initialization ============ 
//   for (ip=0;ip < p ; ip++) beta[ip]=0.0;
//   aG = 0.05;
//   rG = 0.05;
//   for (ipat = 0 ; ipat < N ; ipat++)
//     {
//       g1s[ipat]=0.5;
//       etas[ipat]=0.0;
//     }
//   // ================ MCMC =================
//   for (iB = 0 ; iB < B; iB++)
//     {
//       // STEP 1 
//       // Update g1s, a vector of length N via M-H. 
//       for (ipat=0; ipat < N; ipat++)
// 	{
// 	  // proposal distribution
// 	  candi = rnorm(g1s[ipat],can[0]);
// 	  att[0]++;
// 	  if (candi > 0 && candi < 1)
// 	    {
// 	      MHrate = dbeta(candi, aG, rG, 1) - dbeta(g1s[ipat], aG, rG, 1);
// 	      for ( ivec = 0; ivec < wIDsizes[ipat]; ivec++ )
// 		{
// 		  yij = Y[ wIDs[ipat][ivec] ];
// 		  canpij = getpijEta(candi,etas[ipat], labelnp[wIDs[ipat][ivec]]);
// 		  pij = getpijEta(g1s[ipat],etas[ipat], labelnp[wIDs[ipat][ivec]]);
// 		  sizeij = getsizeij(ivec, X, wIDs[ipat], beta, Ntot, p); // the value of cansizeij is changed
// 		  MHrate += dnbinom(yij, sizeij,  canpij, 1)-dnbinom(yij, sizeij,  pij, 1); 
// 		}
// 	      if(runif(0.0,1.0) < exp(MHrate))
// 		{
// 		  g1s[ipat] = candi;
// 		  acc[0]++;
// 		}
// 	    }
// 	}    
//       // STEP 2
//       // Update etas, a vector of length N via M-H. 
//       for (ipat=0; ipat < N; ipat++)
// 	{
// 	  // proposal distribution
// 	  candi = rnorm(etas[ipat],can[1]);
// 	  att[1]++;
// 	  // prior of g is norm(mean=0,sd=1/sqrt(prec_eta)); prec_eta ~ gamma(shp_prec_eta,1/rate_prec_eta)
// 	  MHrate = dnorm(candi, 0.0, 1/sqrt(prec_eta), 1) - dnorm(g1s[ipat], 0.0, 1/sqrt(prec_eta), 1);
// 	  for ( ivec = 0; ivec < wIDsizes[ipat]; ivec++ )
// 	    {
// 	      if (labelnp[wIDs[ipat][ivec]] == 1) // only the new scans
// 		{
// 		  yij = Y[ wIDs[ipat][ivec] ];
// 		  canpij = getpijEta(g1s[ipat],candi, labelnp[wIDs[ipat][ivec]]);
// 		  pij = getpijEta(g1s[ipat],etas[ipat], labelnp[wIDs[ipat][ivec]]);
// 		  sizeij = getsizeij(ivec, X, wIDs[ipat], beta, Ntot, p); // the value of cansizeij is changed
// 		  MHrate += dnbinom(yij, sizeij,  canpij, 1)-dnbinom(yij, sizeij,  pij, 1); 
// 		}
// 	    }
// 	  if(runif(0.0,1.0) < exp(MHrate))
// 	    {
// 	      etas[ipat] = candi;
// 	      acc[1]++;
// 	    }
// 	}

//       // prec_eta
//       att[4]++;
//       candi = rnorm(prec_eta,can[4]/5);
//       //Rprintf("prec_eta %f candi %f can[4] %f",prec_eta,candi,can[4]);
//       if (min_prec_eta < candi && candi < max_prec_eta )
// 	{
// 	  MHrate = 0;
// 	  for (ipat = 0 ; ipat < N ; ipat++)
// 	    MHrate += dnorm(etas[ipat],0.0,1/sqrt(candi),1) - dnorm(etas[ipat],0.0,1/sqrt(prec_eta),1);
// 	  if (runif(0.0,1.0)< exp(MHrate) ) 
// 	    {
// 	      acc[4]++;
// 	      prec_eta=candi;
// 	    }
// 	}
      

//       // STEP 3: 
//       // Update beta, a vector of length p, containing the coefficients b0, b1,...
//       // M-H 
//       // beta will be updated till the posterior categorical distribution of H1 contains at least one non zero probability.
//       for ( ibeta = 0; ibeta < p; ibeta++)
// 	{
// 	  for (ip = 0 ; ip < p ;ip ++ ) canbeta[ip]=beta[ip]; // copy beta
// 	  att[ibeta+5]++;
// 	  // proposal distribution
// 	  canbeta[ibeta] = rnorm(beta[ibeta], can[ibeta+5]); // the mean and sd must be double (not NumericVector)
// 	  // prior
// 	  MHrate = dnorm(canbeta[ibeta], mu_betas[ibeta], sigma_betas[ibeta], 1) - dnorm(beta[ibeta], mu_betas[ibeta], sigma_betas[ibeta], 1);
	  
// 	  for (ipat=0; ipat < N; ipat++)
// 	    {
// 	      // if (printing)  Rprintf(" \v\v ipat: %d",ipat);
// 	      for ( ivec = 0; ivec < wIDsizes[ipat]; ivec++ )
// 		{
// 		  yij = Y[ wIDs[ipat][ivec] ];
// 		  pij = getpijEta(g1s[ipat],etas[ipat], labelnp[wIDs[ipat][ivec]]);
// 		  cansizeij = getsizeij(ivec, X, wIDs[ipat], canbeta, Ntot, p); // the value of cansizeij is changed
// 		  sizeij = getsizeij(ivec, X, wIDs[ipat], beta, Ntot, p); // the value of cansizeij is changed
// 		  MHrate += dnbinom(yij, cansizeij,  pij, 1)-dnbinom(yij, sizeij,  pij, 1); 
// 		}
// 	    }    

// 	  if(runif(0.0,1.0) < exp(MHrate))
// 	    {
// 	      beta[ibeta] = canbeta[ibeta];
// 	      acc[ibeta+5]++;
// 	    }
// 	}

//       // Update aG and rG, scalars for the base distribution of the random effect gi
//       // gi ~ Beta(aG,rG)
//       // aG ~ unif(0,max_aG)
//       // rG ~ unif(0,max_aG)
//       candi=rnorm(aG,can[2]);
//       if (candi > 1.0e-5 && candi < max_aG )
// 	{
// 	  att[2]++;
// 	  MHrate = 0;
// 	  for (ipat=0 ; ipat < N; ipat ++ )  MHrate += dbeta(g1s[ipat],candi,rG,1)-dbeta(g1s[ipat],aG,rG,1); 
// 	  if (runif(0.0,1.0) < exp(MHrate) ) 
// 	    {
// 	      aG = candi;
// 	      acc[2]++;
// 	    }
// 	}
	
      
//       candi=rnorm(rG,can[3]);
//       if (candi > 1.0e-5 && candi < max_aG )
// 	{
// 	  att[3]++;
// 	  MHrate = 0;
// 	  for (ipat=0 ; ipat < N; ipat ++ ) MHrate += dbeta(g1s[ipat],aG,candi,1)-dbeta(g1s[ipat],aG,rG,1);   
// 	  if (runif(0.0,1.0) < exp(MHrate) ) 
// 	    {
// 	      rG = candi;
// 	      acc[3]++;
// 	    }
// 	}
      
//       // http://stackoverflow.com/questions/8720550/how-to-return-array-of-structs-from-call-to-c-shared-library-in-r
//       // SAVE the result from this iteration
//       REAL(postaGs)[iB] = aG;
//       REAL(postrGs)[iB] = rG;
//       REAL(post_prec_eta)[iB] = prec_eta;
//       for (ibeta = 0; ibeta < p; ibeta++) 
// 	{
// 	  idh = ibeta + iB*p;
// 	  REAL(postbetas)[idh] = beta[ibeta];
// 	}
      
//       for (ipat = 0; ipat < N; ipat++ ) 
// 	{
// 	  idh = ipat + iB*N;
// 	  REAL(postg1s)[idh] = g1s[ipat];
// 	  REAL(postetas)[idh] = etas[ipat];
// 	  REAL(postg2s)[idh] = getpijEta(g1s[ipat],etas[ipat],1);
// 	}

//       if (iB < burnin ) //2*iB
// 	{
// 	  for(ih=0; ih < 5+p ; ih++)
// 	    {
// 	      // Keept the acceptance to range between 0.3 and 0.5
// 	      // The proposal variance is adjusted only during the burn-in period
// 	      if( ih==2 ) can[ih]=FisherBeta(aG,rG,N)*3; 
// 	      else if( ih==3 ) can[ih]=FisherBeta(rG,aG,N)*3; 
// 	      else 
// 		{
// 		  if (att[ih]>50)  can[ih]=adjustcan(acc[ih]/att[ih],can[ih],0.1);
// 		}
// 	    }
// 	}else if (iB==burnin){
// 	for( ih=0; ih < 5+p; ih++){
// 	  att[ih] = 0.0;
// 	  acc[ih] = 0.0;
// 	}
//       }  
      
//       if(iB > 0 && iB % INTEGER(FreqEvery)[0] == 0) Rprintf("\v %d iterations are done...   ", iB);
//     }
//   for (ih=0;ih < 5+p; ih ++) 
//     {
//       REAL(AR)[ih] = acc[ih]/att[ih];
//       REAL(prp)[ih] = can[ih];
//     }
  
//   PutRNGstate();
//   UNPROTECT(1);
//   return res;
// }





// extern "C" {
//   SEXP Beta23(SEXP Y_,          // REAL
// 	      SEXP X,           // REAL 
// 	      SEXP ID,          // INTEGER
// 	      SEXP B_,          // INTEGER
// 	      SEXP maxni,       // INTEGER
// 	      SEXP nn,          // INTEGER
// 	      SEXP M_,          // INTEGER
// 	      SEXP aD_,         // REAL
// 	      SEXP ibD_,        // REAL
// 	      SEXP labelnp_,    // REAL
// 	      SEXP sd_eta_,  // REAL
// 	      SEXP max_aG_,     // REAL
// 	      SEXP mu_beta_,    // REAL
// 	      SEXP sigma_beta_, // REAL
// 	      SEXP burnin_,     // INTEGER
// 	      SEXP FreqEvery,    // INTEGER
// 	      SEXP initial_beta_
// 	      );
// }




// SEXP Beta23(SEXP Y_,          // REAL
// 	    SEXP X,           // REAL 
// 	    SEXP ID,          // INTEGER
// 	    SEXP B_,          // INTEGER
// 	    SEXP maxni,       // INTEGER
// 	    SEXP nn,          // INTEGER
// 	    SEXP M_,          // INTEGER
// 	    SEXP aD_,         // REAL
// 	    SEXP ibD_,        // REAL
// 	    SEXP labelnp_,    // REAL
// 	    SEXP sd_eta_,  // REAL
// 	    SEXP max_aG_,     // REAL
// 	    SEXP mu_beta_,    // REAL
// 	    SEXP sigma_beta_, // REAL
// 	    SEXP burnin_,     // INTEGER
// 	    SEXP FreqEvery,    // INTEGER
// 	    SEXP initial_beta_
// 	    )
// {
//   // === Model ===
//   // Y_ij | Gij = gi ~ NB(size=exp(X_{ij}^T beta),prob=gij)
//   // gij = g1                     if j is in pre-scan 
//   //     = g1/[(1-g1)exp(-eta)+g1] if j is in old-scan
//   // That is: logit(g2) = logit(g1)+eta where eta \in (-inf,inf)
//   // g1 < g2 <=> 0 < eta 
//   // g1 ~ sum_{h=1}^{inf} beta(aG_h,rG_h) pi_h
//   // pi_h = v_h prod_{l < h} (1-v_l)
//   // v_h ~ beta(1,D)
//   // D ~ gamma(a_D,1/r_D)
//   // eta ~ norm(mean=0.0,sd=sd_eta) where sd_eta = 0.3
//   // beta ~ rnorm(mu_beta,sigma_beta)
//   // aG, rG ~ unif(0,max_aG)
//   // sd_eta is not random as if it is, the MSE of random effects are not as good.
//   const double *Y = REAL(Y_); 
//   const double *labelnp = REAL(labelnp_);
//   const int Ntot = length(ID); // # patients
//   const int p = length(X)/Ntot; 
//   const int N = INTEGER(nn)[0];
//   const int M = INTEGER(M_)[0];
//   const double a_D = REAL(aD_)[0];
//   const double ib_D = REAL(ibD_)[0];
//   const double max_aG = REAL(max_aG_)[0];
//   const int B = INTEGER(B_)[0];
//   const double *mu_betas = REAL(mu_beta_);
//   const double *sigma_betas = REAL(sigma_beta_);
//   const int burnin = INTEGER(burnin_)[0];
//   const int Nacc = 2+ M*2 + p; 
//   const double smallest = 1.0e-7;
//   const double sd_eta = REAL(sd_eta_)[0];
//   // initial_beta_ can be changed for Gelman-Rubin Diagnostics
//   const double initial_beta = REAL(initial_beta_)[0];
//   // acc[0] = g1s, acc[1] = eta, 
//   // acc[2]= aGs[0],...,acc[M+1]=aGs[M-1],
//   // acc[M+2]=rGs[0],...,acc[2M+1]=rGs[M-1],
//   // acc[2M+2]=beta[0],acc[2M+3]=beta[1],...,acc[2M+2+p-1]=beta[p-1],
//   int wID[INTEGER(maxni)[0]],wIDsize = 0, ih = 0, ip = 0, ipat = 0, ivec=0, idh=0, ibeta=0,iB=0;
//   // wIDs[ipat][ivec] contains the position of the ivec^th repeated measure of the ipat^th patient at labelnp and Y
//   //                  If ivec >= # repeated measures, wIDs[ipat][ivec] = -1000
//   // wIDsizes[ipat] contains the number of repeated measure of ipat^th patient
//   // sumYs[ipat][0] contains the sum of the CEL counts of ipat^th patients in the pre-scan period
//   // sumYs[ipat][1] contains the sum of the CEL counts of ipat^th patients in the new-scan period
//   int wIDs[N][INTEGER(maxni)[0]], wIDsizes[N];
//   double sumYs[N][2];
//   for (ipat = 0 ; ipat < N ; ipat++)
//     {
//       getwID( &wIDsize,wID,ID,ipat);
//       wIDsizes[ipat] = wIDsize;
//       sumYs[ipat][0] = 0;
//       sumYs[ipat][1] = 0;
//       for (ivec = 0 ; ivec < wIDsize; ivec++)
// 	{
// 	  wIDs[ipat][ivec]=wID[ivec];
// 	  if (labelnp[wIDs[ipat][ivec]]==0) sumYs[ipat][0] += Y[wIDs[ipat][ivec]];
// 	  if (labelnp[wIDs[ipat][ivec]]==1) sumYs[ipat][1] += Y[wIDs[ipat][ivec]];
// 	}
//       for (ivec = wIDsize ; ivec < INTEGER(maxni)[0]; ivec++)
// 	wIDs[ipat][ivec]=-1000;
//     }

//   int h1s[N];
//   double beta[p],canbeta[p],aGs[M],rGs[M],g1s[N],etas[N];
//   //double prec_eta = (min_prec_eta+max_prec_eta)/2;
//   double vs[M],weightH1[M],postprob[M];
//   double sizeij=0.0,cansizeij=0.0,w=0.0,canpij=0.0,pij=0.0,D=1,
//     candi=0.0,MHrate=0.0,yij=0.0,sumgs=0.0;
//   double att[Nacc], acc[Nacc],can[Nacc];
//   for (ih = 0 ; ih < Nacc; ih ++)
//     {
//       att[ih] = 0.0;
//       acc[ih] = 0.0;
//       can[ih] = 2.0;
//     }

//   // http://stackoverflow.com/questions/8720550/how-to-return-array-of-structs-from-call-to-c-shared-library-in-r
//   SEXP res = PROTECT(allocVector(VECSXP, 12)); // result is stored in res
//   SEXP postaGs = allocVector(REALSXP, M*B); 
//   SET_VECTOR_ELT(res, 0, postaGs); 
//   SEXP postrGs = allocVector(REALSXP, M*B); 
//   SET_VECTOR_ELT(res, 1, postrGs); 
//   SEXP postvs = allocVector(REALSXP, M*B); 
//   SET_VECTOR_ELT(res, 2, postvs); 
//   SEXP postweightH1 = allocVector(REALSXP, M*B); 
//   SET_VECTOR_ELT(res, 3, postweightH1); 
//   SEXP posth1s = allocVector(INTSXP, N*B); 
//   SET_VECTOR_ELT(res, 4, posth1s); 
//   SEXP postg1s = allocVector(REALSXP, N*B); 
//   SET_VECTOR_ELT(res, 5, postg1s); 
//   SEXP postg2s = allocVector(REALSXP, N*B); 
//   SET_VECTOR_ELT(res, 6, postg2s); 
//   SEXP postetas = allocVector(REALSXP, N*B); 
//   SET_VECTOR_ELT(res, 7, postetas); 
//   SEXP postbetas = allocVector(REALSXP, p*B); 
//   SET_VECTOR_ELT(res, 8, postbetas); 
//   SEXP postD = allocVector(REALSXP, B); 
//   SET_VECTOR_ELT(res, 9, postD); 
//   //SEXP post_prec_eta = allocVector(REALSXP, B); 
//   //SET_VECTOR_ELT(res, 10, post_prec_eta);
//   SEXP AR = allocVector(REALSXP,Nacc);
//   SET_VECTOR_ELT(res, 10, AR);
//   SEXP prp = allocVector(REALSXP,Nacc);
//   SET_VECTOR_ELT(res, 11, prp);


//   GetRNGstate();
  
//   // ============ initialization ============ 
//   for (ip=0;ip < p ; ip++) beta[ip]=initial_beta;
//   for ( ih = 0 ; ih < M ; ih++)
//     {
//       aGs[ih] = runif(smallest,max_aG);
//       rGs[ih] = runif(smallest,max_aG);
//     }

//   D = rgamma(a_D,1/ib_D);
  
//   for (ih = 0 ; ih < (M-1) ; ih ++) vs[ih] = rbeta(1,D);
//   vs[M-1]=1.0;
  
//   weightH1[0] = vs[0];
//   w = 1.0 ;
//   for (ih=1; ih < M; ih++)
//     {
//       w *= (1.0-vs[ih-1]);
//       weightH1[ih] = vs[ih]*w;
//     }

//   for ( ipat=0; ipat < N; ipat++) 
//     {
//       etas[ipat]=0.0;
//       h1s[ipat] = sample(weightH1,M);
//       g1s[ipat] = rbeta(aGs[h1s[ipat]],rGs[h1s[ipat]]);
//       if (g1s[ipat] > 1.0-smallest) g1s[ipat]=1.0-smallest;
//       else if (g1s[ipat] < smallest) g1s[ipat]=smallest;
//     }


//   // ================ MCMC =================
//   for (iB = 0 ; iB < B; iB++)
//     {
//       // STEP 1 
//       // Update g1s, a vector of length N via M-H. 
//       for (ipat=0; ipat < N; ipat++)
// 	{
// 	  // proposal distribution
// 	  candi = rnorm(g1s[ipat],can[0]);
// 	  att[0]++;
// 	  if (candi > smallest && candi < 1.0-smallest)
// 	    {
// 	      MHrate = dbeta(candi, aGs[h1s[ipat]], rGs[h1s[ipat]], 1) - dbeta(g1s[ipat], aGs[h1s[ipat]], rGs[h1s[ipat]], 1);
// 	      for ( ivec = 0; ivec < wIDsizes[ipat]; ivec++ )
// 		{
// 		  yij = Y[ wIDs[ipat][ivec] ];
// 		  canpij = getpijEta(candi,etas[ipat], labelnp[wIDs[ipat][ivec]]);
// 		  pij = getpijEta(g1s[ipat],etas[ipat], labelnp[wIDs[ipat][ivec]]);
// 		  sizeij = getsizeij(ivec, X, wIDs[ipat], beta, Ntot, p); // the value of cansizeij is changed

// 		  //Rprintf("\v ipat %d ivec %d labelnp %f",ipat,ivec,labelnp[wIDs[ipat][ivec]]);

// 		  MHrate += dnbinom(yij, sizeij,  canpij, 1)-dnbinom(yij, sizeij,  pij, 1); 
// 		}
// 	      if(runif(0.0,1.0) < exp(MHrate))
// 		{
// 		  g1s[ipat] = candi;
// 		  acc[0]++;
// 		}
// 	    }
// 	}    
//       // STEP 2
//       // Update etas, a vector of length N via M-H. 
//       for (ipat=0; ipat < N; ipat++)
// 	{
// 	  // proposal distribution
// 	  candi = rnorm(etas[ipat],can[1]);
// 	  att[1]++;
// 	  // prior of eta is norm(mean=0,sd=sd_eta); 
// 	  MHrate = dnorm(candi, 0.0, sd_eta , 1) - dnorm(g1s[ipat], 0.0, sd_eta, 1);
// 	  for ( ivec = 0; ivec < wIDsizes[ipat]; ivec++ )
// 	    {
// 	      if (labelnp[wIDs[ipat][ivec]] == 1.0) // only the new scans
// 		{
// 		  yij = Y[ wIDs[ipat][ivec] ];
// 		  canpij = getpijEta(g1s[ipat],candi, 1);
// 		  pij = getpijEta(g1s[ipat],etas[ipat], 1);
// 		  sizeij = getsizeij(ivec, X, wIDs[ipat], beta, Ntot, p); // the value of cansizeij is changed
// 		  MHrate += dnbinom(yij, sizeij,  canpij, 1)-dnbinom(yij, sizeij,  pij, 1); 
// 		}
// 	    }
// 	  if(runif(0.0,1.0) < exp(MHrate))
// 	    {
// 	      etas[ipat] = candi;
// 	      acc[1]++;
// 	    }
// 	}
      
//       // for (ipat = 0 ; ipat < N;ipat++) etas[ipat]=0.0;
//       /* // prec_eta */
//       /* att[Nacc-1]++; */
//       /* candi = rnorm(prec_eta,can[Nacc-1]/5); */
//       /* //Rprintf("prec_eta %f candi %f can[4] %f",prec_eta,candi,can[4]); */
//       /* if (min_prec_eta < candi && candi < max_prec_eta ) */
//       /* 	{ */
//       /* 	  MHrate = 0; */
//       /* 	  for (ipat = 0 ; ipat < N ; ipat++) */
//       /* 	    MHrate += dnorm(etas[ipat],0.0,1/sqrt(candi),1) - dnorm(etas[ipat],0.0,1/sqrt(prec_eta),1); */
//       /* 	  if (runif(0.0,1.0)< exp(MHrate) )  */
//       /* 	    { */
//       /* 	      acc[Nacc-1]++; */
//       /* 	      prec_eta=candi; */
//       /* 	    } */
//       /* 	} */

//       // STEP 3: 
//       // Update beta, a vector of length p, containing the coefficients b0, b1,...
//       // M-H 
//       // beta will be updated till the posterior categorical distribution of H1 contains at least one non zero probability.
//       for ( ibeta = 0; ibeta < p; ibeta++)
// 	{
// 	  for (ip = 0 ; ip < p ;ip ++ ) canbeta[ip]=beta[ip]; // copy beta
// 	  att[2*M+2+ibeta]++;
// 	  // proposal distribution
// 	  canbeta[ibeta] = rnorm(beta[ibeta], can[2*M+2+ibeta]); // the mean and sd must be double (not NumericVector)
// 	  // prior
// 	  MHrate = dnorm(canbeta[ibeta], mu_betas[ibeta], sigma_betas[ibeta], 1) - dnorm(beta[ibeta], mu_betas[ibeta], sigma_betas[ibeta], 1);
	  
// 	  for (ipat=0; ipat < N; ipat++)
// 	    {
// 	      // if (printing)  Rprintf(" \v\v ipat: %d",ipat);
// 	      for ( ivec = 0; ivec < wIDsizes[ipat]; ivec++ )
// 		{
// 		  yij = Y[ wIDs[ipat][ivec] ];
// 		  pij = getpijEta(g1s[ipat],etas[ipat], labelnp[wIDs[ipat][ivec]]);
// 		  cansizeij = getsizeij(ivec, X, wIDs[ipat], canbeta, Ntot, p); // the value of cansizeij is changed
// 		  sizeij = getsizeij(ivec, X, wIDs[ipat], beta, Ntot, p); // the value of cansizeij is changed
// 		  MHrate += dnbinom(yij, cansizeij,  pij, 1)-dnbinom(yij, sizeij,  pij, 1); 
// 		}
// 	    }    

// 	  if(runif(0.0,1.0) < exp(MHrate))
// 	    {
// 	      beta[ibeta] = canbeta[ibeta];
// 	      acc[2*M+2+ibeta]++;
// 	    }
// 	}
    

//       // Update aGs and rGs, scalars for the base distribution of the random effect gi
//       // gi ~ Beta(aG,rG)
//       // aG ~ unif(0,max_aG)
//       // rG ~ unif(0,max_aG)
//       for (ih = 0 ; ih < M ; ih ++)
// 	{
// 	  candi=rnorm(aGs[ih],can[2+ih]);
// 	  att[2+ih]++;
// 	  if (candi > smallest && candi < max_aG )
// 	    {
// 	      MHrate = 0;
// 	      for (ipat=0 ; ipat < N; ipat ++ )  
// 		{
// 		  if (h1s[ipat]==ih)
// 		    MHrate += dbeta(g1s[ipat],candi,rGs[ih],1)-dbeta(g1s[ipat],aGs[ih],rGs[ih],1); 
// 		}
// 	      if (runif(0.0,1.0) < exp(MHrate) ) 
// 		{
// 		  aGs[ih] = candi;
// 		  acc[2+ih]++;
// 		}
// 	    }
 	  
	  
// 	  candi=rnorm(rGs[ih],can[M+2+ih]);
// 	  att[M+2+ih]++;
// 	  if (candi > smallest && candi < max_aG )
// 	    {
// 	      MHrate = 0;
// 	      for (ipat=0 ; ipat < N; ipat ++ ) 
// 		{
// 		  if (h1s[ipat]==ih) 
// 		    MHrate += dbeta(g1s[ipat],aGs[ih],candi,1)-dbeta(g1s[ipat],aGs[ih],rGs[ih],1);  
// 		}
// 	      if (runif(0.0,1.0) < exp(MHrate) ) 
// 		{
// 		  rGs[ih] = candi;
// 		  acc[M+2+ih]++;
// 		}
// 	    }
// 	}
      
//       // Update vs, a vector of length M, containing the latent variable v, used to construct pis. 
//       // Once vs[ih] is updated, update the probabilities of the categorical distribution of H1 based on the formula:
//       // weightH1[0] = vs[0]
//       // weightH1[ih]= vs[ih]*(1-vh[ih-1])**(1-vh[ih-2])*...**(1-vh[0]) for ih =1,2,...M-1
//       w = 1.0;
//       for (ih = 0 ; ih < (M-1) ; ih ++) // vs[M-1] = 1 always!
// 	{
// 	  ip = 0;    // # patients with h1s == ih
// 	  ibeta = 0; // # patients with h1s > ih
// 	  for (ipat = 0 ; ipat < N; ipat++ ) 
// 	    {
// 	      if ( h1s[ipat]==ih ) ip++ ;
// 	      if ( h1s[ipat] > ih ) ibeta++ ;
// 	    }
	  
// 	  vs[ih] = rbeta(
// 			 (double) 1.0 +ip, 
// 			 D + (double)ibeta 
// 			 );
// 	  //Rprintf("ih %d sp1 %f sp2 %f",ih, (double) 1+ip+idh ,  D + (double)ibeta + (double)ivec);
// 	  if (ih > 0) w *= (1.0-vs[ih-1]);
// 	  weightH1[ih] = vs[ih]*w;
// 	}
//       weightH1[M-1] = w*(1-vs[M-2]);


//       //  Update h1s, a vector of length N, containing the cluster index of random effect 1 of each patient.
//       //  the cluster index is between 1 and M.
//       //  conjugate: posterior distribution of h1s has a categorical distribution with M categories. 
//       // if (printing && iB%printEvery==0)  Rprintf("\v\v Step 4");
//       for (ipat = 0 ; ipat < N ; ipat++)
// 	{
	  
// 	  // Given patient ipat, compute P(H_1i=ih;-) for each ih 
// 	  // P(H_1i=ih;-) = prod_{i=1}^{# repeated measure} dnbinom(yij; size=exp(b0+b1*x1i+..), prob=1/(gij+1))
	  
// 	  for (ih=0; ih < M; ih++ ) postprob[ih] = dbeta(g1s[ipat],aGs[ih],rGs[ih],0)*weightH1[ih];
// 	  h1s[ipat] = sample(postprob,M); 
	  
// 	}

//       // Update D
//       // D ~ gamma(scale=a.D,shape=b.D) = gamma(scale=a.D,shape=1/ib.D)  
//       w = 0.0;
//       for (ih=0;ih<(M-1);ih++) w += log(1.0-vs[ih]);
//       D = rgamma( (a_D+ (double) M-1.0), (1.0/(-w+ib_D)) ); // N+w is the number of random var following DP given js
//       if (ISNAN(D)){ 
// 	D = 1.0; 
// 	Rprintf("\v ISNAN(D): a_D %f ib_D %f sum log(1-vs[ih]) %f ",a_D,ib_D,w);
//       }
//       if( D <0.01) D = 0.01; // if D is too small, rbeta(1,a,0) returns NaN for any a
      



//       if (iB < burnin ) //2*iB
// 	{
// 	  for(ih=0; ih < Nacc ; ih++)
// 	    {
// 	      // Keept the acceptance to range between 0.3 and 0.5
// 	      // The proposal variance is adjusted only during the burn-in period
// 	      if( 1 < ih && ih < M+2 ) can[ih]=FisherBeta(aGs[ih-2],rGs[ih-2],N)*3;  // update can of aGs[ih]
// 	      else if( M+1 < ih && ih < 2*M+2 ) can[ih]=FisherBeta(rGs[ih-M-2],aGs[ih-M-2],N)*3;  // update can of rGs[ih]
// 	      else{ 
// 		if(att[ih]>50)  can[ih]=adjustcan(acc[ih]/att[ih],can[ih],0.1);
// 	      }
// 	      if (!R_FINITE(can[ih]))
// 		Rprintf("!");
// 	    }
// 	}else if (iB==burnin){
// 	for( ih=0; ih < Nacc; ih++){
// 	  att[ih] = 0.0;
// 	  acc[ih] = 0.0;
// 	}
//       }  


//       // http://stackoverflow.com/questions/8720550/how-to-return-array-of-structs-from-call-to-c-shared-library-in-r
//       // SAVE the result from this iteration
//       for (ih=0;ih < M ; ih++)
// 	{
// 	  idh = ih+iB*M;
// 	  REAL(postaGs)[idh] = aGs[ih];
// 	  REAL(postrGs)[idh] = rGs[ih];
// 	  REAL(postweightH1)[idh] = weightH1[ih];
// 	  REAL(postvs)[idh] = vs[ih];
// 	}

//       REAL(postD)[iB] = D;

//       for (ibeta = 0; ibeta < p; ibeta++) 
// 	{
// 	  idh = ibeta + iB*p;
// 	  REAL(postbetas)[idh] = beta[ibeta];
// 	}
      
//       for (ipat = 0; ipat < N; ipat++ ) 
// 	{
// 	  idh = ipat + iB*N;
// 	  REAL(postg1s)[idh] = g1s[ipat];
// 	  REAL(postg2s)[idh] = getpijEta(g1s[ipat],etas[ipat],1);
// 	  REAL(postetas)[idh] = etas[ipat];
// 	  INTEGER(posth1s)[idh] = h1s[ipat];
// 	}

      
//       if(iB > 0 && iB % INTEGER(FreqEvery)[0] == 0) Rprintf("\v %d iterations are done...   ", iB);
//     }
//   for (ih=0;ih < Nacc; ih ++) 
//     {
//       REAL(AR)[ih] = acc[ih]/att[ih];
//       REAL(prp)[ih] = can[ih];
//     }
  
//   PutRNGstate();
//   UNPROTECT(1);
//   return res;
// }



// extern "C" {
// SEXP Beta14Unif(SEXP Y_,          // REAL
// 		SEXP X,           // REAL 
// 		SEXP ID,          // INTEGER
// 		SEXP B_,          // INTEGER
// 		SEXP maxni,       // INTEGER
// 		SEXP nn,          // INTEGER
// 		SEXP labelnp_,    // REAL
// 		SEXP a_qG_,       // REAL
// 		SEXP r_qG_,       // REAL
// 		SEXP max_aG_,
// 		//SEXP mu_aG_,
// 		//SEXP sd_aG_,
// 		//SEXP mu_rG_,
// 		//SEXP sd_rG_,
// 		SEXP mu_beta_,    // REAL
// 		SEXP evalue_sigma_betas_,
// 		SEXP Inv_sigma_beta_,
// 		SEXP burnin_,     // INTEGER
// 		SEXP FreqEvery,   // INTEGER
// 		SEXP patwoNS_      // INTEGER
// 		);
// }

SEXP Beta14Unif(SEXP Y_,          // REAL
		SEXP X,           // REAL 
		SEXP ID,          // INTEGER
		SEXP B_,          // INTEGER
		SEXP maxni,       // INTEGER
		SEXP nn,          // INTEGER
		SEXP labelnp_,    // REAL
		SEXP a_qG_,       // REAL
		SEXP r_qG_,       // REAL
		SEXP max_aG_,
		//SEXP mu_aG_,
		//SEXP sd_aG_,
		//SEXP mu_rG_,
		//SEXP sd_rG_,
		SEXP mu_beta_,    // REAL
		SEXP evalue_sigma_betas_,
		SEXP Inv_sigma_beta_,
		SEXP burnin_,     // INTEGER
		SEXP FreqEvery,   // INTEGER
		SEXP patwoNS_      // INTEGER
		)
{
  // Last modified:: Oct 26, 2012
  // The patients without new scans can be inputed.
  // === Model ===
  // Y_ij | Gij = gij ~ NB(size=exp(X_{ij}^T beta),prob=gij)
  // gij = g1 if j is in pre-scan 
  //     = g_new if j is in old-scan 
  // g_new = Ji * g1 + (1-Ji) * g2 
  // Ji ~ ber(qG)
  // qG ~ beta(a_qG,r_qG)
  // g1, g2 ~ beta(aG,rG)
  // beta ~ rnorm(mu_beta,sigma_beta)
  // aG, rG ~ unif(0,max_aG)
  const double *Y = REAL(Y_); 
  const double *labelnp = REAL(labelnp_);
  const double a_qG = REAL(a_qG_)[0];
  const double r_qG = REAL(r_qG_)[0];
  const int Ntot = length(ID); // # patients
  const int p = length(X)/Ntot; 
  const int N = INTEGER(nn)[0];
  const double max_aG = REAL(max_aG_)[0];
  //const double mu_aG = REAL(mu_aG_)[0];
  //const double sd_aG = REAL(sd_aG_)[0];
  //const double mu_rG = REAL(mu_rG_)[0];
  //const double sd_rG = REAL(sd_rG_)[0];
  const int B = INTEGER(B_)[0];
  const double *mu_betas = REAL(mu_beta_);
  const double *evalue_sigma_betas = REAL(evalue_sigma_betas_);
  const double *Inv_sigma_beta = REAL(Inv_sigma_beta_);
  const int burnin = INTEGER(burnin_)[0];
  const int *patwoNS = INTEGER(patwoNS_);// If all patients have new scans, then patwoNS[0]=-1000.
  int NpatwoNS = length(patwoNS_);  // The number of patients without new scans
  const double smallest2 = 0.01; 
  if (patwoNS[0]==-1000) NpatwoNS = 0; 
 
  int wID[INTEGER(maxni)[0]],wIDsize = 0, ih = 0, ip = 0, ipat = 0, ivec=0, idh=0, ibeta=0,iB=0;
  // wIDs[ipat][ivec] contains the position of the ivec^th repeated measure of the ipat^th patient at labelnp and Y
  //                  If ivec >= # repeated measures, wIDs[ipat][ivec] = -1000
  // wIDsizes[ipat] contains the number of repeated measure of ipat^th patient
  // sumYs[ipat][0] contains the sum of the CEL counts of ipat^th patients in the pre-scan period
  // sumYs[ipat][1] contains the sum of the CEL counts of ipat^th patients in the new-scan period
  int wIDs[N][INTEGER(maxni)[0]], wIDsizes[N];
  double sumYs[N][2];
  for (ipat = 0 ; ipat < N ; ipat++)
    {
      getwID( &wIDsize,wID,ID,ipat);
      wIDsizes[ipat] = wIDsize;
      sumYs[ipat][0] = 0;
      sumYs[ipat][1] = 0;
      for (ivec = 0 ; ivec < wIDsize; ivec++)
	{
	  wIDs[ipat][ivec]=wID[ivec];
	  if (labelnp[wIDs[ipat][ivec]]==0) sumYs[ipat][0] += Y[wIDs[ipat][ivec]];
	  if (labelnp[wIDs[ipat][ivec]]==1) sumYs[ipat][1] += Y[wIDs[ipat][ivec]];
	}
      for (ivec = wIDsize ; ivec < INTEGER(maxni)[0]; ivec++)
	wIDs[ipat][ivec]=-1000;
    }

  //h1s[N],h2s[N]
  double beta[p],canbeta[p],aG,rG,qG,g1s[N],g2s[N],js[N];
  // double vs[M],weightH1[M],postprob[M];
  double sizeij=0.0,cansizeij=0.0,w=0.0,canpij=0.0,pij=0.0,
    candi=0.0,MHrate=0.0,yij=0.0,sumgs=0.0;
  // double  weipro=0.0,qG2=0.0,D=0.0;
  double att[2+1], acc[2+1],can[2+1];
  double logLcan, logL; // focused likelihood values
  for (ih = 0 ; ih < 3; ih ++)
    {
      att[ih] = 0.0;
      acc[ih] = 0.0;
      can[ih] = 5.0;
    }
  // http://stackoverflow.com/questions/8720550/how-to-return-array-of-structs-from-call-to-c-shared-library-in-r
  SEXP res = PROTECT(allocVector(VECSXP, 11)); // result is stored in res
  SEXP postaGs = allocVector(REALSXP, B); 
  SET_VECTOR_ELT(res, 0, postaGs); 
  SEXP postrGs = allocVector(REALSXP, B); 
  SET_VECTOR_ELT(res, 1, postrGs); 
  SEXP postqG = allocVector(REALSXP, B); 
  SET_VECTOR_ELT(res, 2, postqG);
  SEXP logLlik = allocVector(REALSXP, B); 
  SET_VECTOR_ELT(res, 3, logLlik);
  
  SEXP postgPre = allocVector(REALSXP, N*B); 
  SET_VECTOR_ELT(res, 4, postgPre); 
  SEXP postg2s = allocVector(REALSXP, N*B); 
  SET_VECTOR_ELT(res, 5, postg2s); 
  SEXP postgNew = allocVector(REALSXP, N*B); 
  SET_VECTOR_ELT(res, 6, postgNew); 
  SEXP postjs = allocVector(REALSXP, N*B); 
  SET_VECTOR_ELT(res, 7, postjs); 

  SEXP postbetas = allocVector(REALSXP, p*B); 
  SET_VECTOR_ELT(res, 8, postbetas); 
  SEXP AR = allocVector(REALSXP,3);
  SET_VECTOR_ELT(res, 9, AR);
  SEXP prp = allocVector(REALSXP,3);
  SET_VECTOR_ELT(res, 10, prp);

  GetRNGstate();
  
  // ============ initialization ============ 
  for (ip=0;ip < p ; ip++) beta[ip]=0.0;
  aG = 0.05;
  rG = 0.05;
  qG = 0.5;
  for (ipat = 0 ; ipat < N ; ipat++)
    {
      g1s[ipat]=0.5;
      g2s[ipat]=0.25;
      js[ipat]=1.0;
      for (ih = 0 ; ih < NpatwoNS; ih++) // If the patient does not have any new scans,g2s[ipat]= js[ipat]=-1000 during the whole algorithm.
	{
	  if (ipat==patwoNS[ih])
	    {
	      js[ipat] = -1000;
	      g2s[ipat] = -1000;
	    }
	}
    }
  // ================ MCMC =================
  for (iB = 0 ; iB < B; iB++)
    {
      // STEP 1 
      // Update g1s and g2s, a vector of length N. 
      // The g2s are updated ASSUMING THAT ji = 0 for all i to be used in updating ji
      for (ipat = 0 ; ipat < N; ipat++)
	{
	  w = 0.0;
	  sumgs = 0.0;
	  for (ivec=0; ivec < wIDsizes[ipat]; ivec ++)
	    {
	      sizeij = getsizeij(ivec, X, wIDs[ipat], beta, Ntot, p);
	      w += sizeij*(1.0-labelnp[wIDs[ipat][ivec]]); // sum_{j in pre} exp(Xij%*%beta)
	      sumgs += sizeij*labelnp[wIDs[ipat][ivec]]; // sum_{j in new} exp(Xij%*%beta)
	    }
	  if (js[ipat]==0) 
	    {
	      // update g1s
	      g1s[ipat] = rbeta(aG+w,rG+sumYs[ipat][0]);
	      // update g2s
	      g2s[ipat] = rbeta(aG+sumgs,rG+sumYs[ipat][1]);
	    }
	  if (js[ipat]==1) 
	    {
	      // update g1s
	      g1s[ipat] = rbeta(aG+w+sumgs,rG+sumYs[ipat][0]+sumYs[ipat][1]);
	      // update g2s
	      g2s[ipat] = rbeta(aG,rG);//prior

	    }
	  if (js[ipat]==-1000) // no new scans
	    {
	      g1s[ipat] = rbeta(aG+w,rG+sumYs[ipat][0]);
	      g2s[ipat] = -1000;
	    }
	  if (g1s[ipat] > 1.0-1.0e-2) g1s[ipat]= 1.0-1.0e-2;
	  if (g2s[ipat] > 1.0-1.0e-2) g2s[ipat]= 1.0-1.0e-2;
	}
      
      // Step 2: update js
      for (ipat = 0; ipat < N; ipat++)
        {
	  if (js[ipat]!=-1000){
	    w = qG;
	    sumgs = (1.0-qG) ;
	    for (ivec = 0 ; ivec < wIDsizes[ipat]; ivec++)
	      {
		if (labelnp[wIDs[ipat][ivec]]==1.0) // going through only the new scans
		  {
		    yij = Y[ wIDs[ipat][ivec] ];
		    sizeij = getsizeij(ivec, X, wIDs[ipat], beta, Ntot, p);
		    
		    // corresponds to Ji=1 gNew = gPre =g1s
		    canpij = getpij2(0.0,g1s[ipat],1.0);
		    w *= dnbinom(yij, sizeij, canpij, 0);
		    // corresponds to Ji=0 gNew = g2s != gPre =g1s
		    pij = getpij2(0.0,g2s[ipat],1.0);
		    sumgs *= dnbinom(yij, sizeij,  pij, 0);
		  }
	      }
	    js[ipat] = rbinom(1.0,w/(sumgs+w));
	    //if (iB == 5000) Rprintf(" \v ipat %d js[ipat] %f g1s[ipat] %f g2s[ipat] %f w %f sumgs %f",
	    //ipat, js[ipat], g1s[ipat], g2s[ipat], w, sumgs);
	  }
        }  

      //==

      // STEP 3: update beta, a vector of length p, containing the coefficients b0, b1,...
      // beta ~ MVR(mu_beta,Sigma_beta); 

      att[2]++;
      // proposal distribution
      for (ibeta = 0; ibeta <p ; ibeta ++)
	canbeta[ibeta] = rnorm(beta[ibeta], can[2]); // the mean and sd must be double (not NumericVector)
      
      // prior
      MHrate = dmvnorm(canbeta, mu_betas,p, evalue_sigma_betas,Inv_sigma_beta,1) 
	- dmvnorm(beta, mu_betas,p, evalue_sigma_betas,Inv_sigma_beta,1);
      logLcan = 0.0;
      logL = 0.0;
      for (ipat=0; ipat < N; ipat++)
	{
	  // if (printing)  Rprintf(" \v\v ipat: %d",ipat);
	  for ( ivec = 0; ivec < wIDsizes[ipat]; ivec++ )
	    {
	      yij = Y[ wIDs[ipat][ivec] ];
	      // yij is a vector b/c it is the first arg of dnbinom
	      pij = getpij2(g1s[ipat],
			    (g1s[ipat]*js[ipat]+g2s[ipat]*(1.0-js[ipat])),
			    labelnp[wIDs[ipat][ivec]]); 
	      //Rprintf("\v\v ipat %d js[ipat] %f ivec %d yij %f g1s[ipat] %f g2s[ipat] %f pij %f", 
	      //ipat,  js[ipat] , ivec,   yij,   g1s[ipat], g2s[ipat],  pij);
	      
	      cansizeij = getsizeij(ivec, X, wIDs[ipat], canbeta, Ntot, p); // the value of cansizeij is changed
	      sizeij = getsizeij(ivec, X, wIDs[ipat], beta, Ntot, p); // the value of cansizeij is changed

	      logLcan += dnbinom(yij, cansizeij,  pij, 1);
	      logL += dnbinom(yij, sizeij,  pij, 1);
	      
	      //MHrate += dnbinom(yij, cansizeij,  pij, 1)-dnbinom(yij, sizeij,  pij, 1); 
	    }
	  //Rprintf("\v");
	}    
      MHrate += logLcan-logL;

      if(runif(0.0,1.0) < exp(MHrate))
	{
	  for (ibeta = 0; ibeta < p;ibeta++) beta[ibeta] = canbeta[ibeta];
	  acc[2]++;
	  logL = logLcan;
	}
  
      // step 3: update qG
      w=0.0;
      for (ipat = 0 ; ipat < N ; ipat ++) 
	{
	  if (js[ipat]==1) w += js[ipat];
	}
      qG = rbeta(a_qG+w,r_qG+N-NpatwoNS-w);
      
      // step 4: Update aG and rG, scalars for the base distribution of the random effect gi
      // gi ~ Beta(aG,rG)
      // aG ~ unif(0,max_aG)
      // rG ~ unif(0,max_aG)

      // Update aG and rG, scalars for the base distribution of the random effect gi
      // gi ~ Beta(aG,rG)
      // aG,rG ~ runif(smallest2,max_aG)

       candi=rnorm(aG,can[0]);
       att[0]++;
       if (candi > smallest2 && candi < max_aG )
	 {
	   MHrate = 0;
	   for (ipat=0 ; ipat < N; ipat ++ )  
	     {
		 MHrate += dbeta(g1s[ipat],candi,rG,1)-dbeta(g1s[ipat],aG,rG,1); 
	     }
	   if (runif(0.0,1.0) < exp(MHrate) ) 
	     {
	       aG = candi;
	       acc[0]++;
	     }
	 }	  
       candi=rnorm(rG,can[1]);
       att[1]++;
       if (candi > smallest2 && candi < max_aG )
	 {
	   MHrate = 0;
	   for (ipat=0 ; ipat < N; ipat ++ ) 
	     {
		 MHrate += dbeta(g1s[ipat],aG,candi,1)-dbeta(g1s[ipat],aG,rG,1);  
	     }
	   if (runif(0.0,1.0) < exp(MHrate) ) 
	     {
	       rG = candi;
	       acc[1]++;
	     }
	 }

      // http://stackoverflow.com/questions/8720550/how-to-return-array-of-structs-from-call-to-c-shared-library-in-r
      // SAVE the result from this iteration
      REAL(postaGs)[iB] = aG;
      REAL(postrGs)[iB] = rG;
      REAL(postqG)[iB] = qG;
      REAL(logLlik)[iB] = logL;
      for (ibeta = 0; ibeta < p; ibeta++) 
	{
	  idh = ibeta + iB*p;
	  REAL(postbetas)[idh] = beta[ibeta];
	}
      
      for (ipat = 0; ipat < N; ipat++ ) 
	{
	  idh = ipat + iB*N;
	  REAL(postgPre)[idh] = g1s[ipat];
	  REAL(postg2s)[idh] = g2s[ipat];
	  if (js[ipat]!=-1000) REAL(postgNew)[idh] = g1s[ipat]*js[ipat]+g2s[ipat]*(1.0-js[ipat]);
	  else REAL(postgNew)[idh] = -1000;
	  REAL(postjs)[idh] = js[ipat];
	}

      if (iB < burnin ) //2*iB
	{
	  for(ih=0; ih < 3 ; ih++)
	    {
	      // Keept the acceptance to range between 0.3 and 0.5
	      // The proposal variance is adjusted only during the burn-in period
	      if( ih==0 ) can[ih]=FisherBeta(aG,rG,N)*3; 
	      else if( ih==1 ) can[ih]=FisherBeta(rG,aG,N)*3; 
	      else 
		{
		  if (att[ih]>50)  can[ih]=adjustcan(acc[ih]/att[ih],can[ih],0.01,10);
		}
	    }
	 
	  /* if (iB % 100==0){ */
	  /*   for(ih=0; ih < 3; ih++){ */
	  /*     att[ih] = 0.0; */
	  /*     acc[ih] = 0.0; */
	  /*   } */
	}else if (iB==burnin){
	for( ih=0; ih < 3; ih++){
	  att[ih] = 0.0;
	  acc[ih] = 0.0;
	}
      }  
      
      if(iB > 0 && iB % INTEGER(FreqEvery)[0] == 0) Rprintf("\v %d iterations are done...   ", iB);
    }
  for (ih=0;ih < 3; ih ++) 
    {
      REAL(AR)[ih] = acc[ih]/att[ih];
      REAL(prp)[ih] = can[ih];
    }
  
  PutRNGstate();
  UNPROTECT(1);
  return res;
}





// extern "C" {
// SEXP Beta14(SEXP Y_,          // REAL
// 	    SEXP X,           // REAL 
// 	    SEXP ID,          // INTEGER
// 	    SEXP B_,          // INTEGER
// 	    SEXP maxni,       // INTEGER
// 	    SEXP nn,          // INTEGER
// 	    SEXP labelnp_,    // REAL
// 	    SEXP a_qG_,       // REAL
// 	    SEXP r_qG_,       // REAL
// 	    SEXP mu_aG_,
// 	    SEXP sd_aG_,
// 	    SEXP mu_rG_,
// 	    SEXP sd_rG_,
// 	    SEXP mu_beta_,    // REAL
// 	    SEXP evalue_sigma_betas_,
// 	    SEXP Inv_sigma_beta_,
// 	    SEXP burnin_,     // INTEGER
// 	    SEXP FreqEvery,   // INTEGER
// 	    SEXP patwoNS_      // INTEGER
// 	    );
// }

SEXP Beta14(SEXP Y_,          // REAL
	    SEXP X,           // REAL 
	    SEXP ID,          // INTEGER
	    SEXP B_,          // INTEGER
	    SEXP maxni,       // INTEGER
	    SEXP nn,          // INTEGER
	    SEXP labelnp_,    // REAL
	    SEXP a_qG_,       // REAL
	    SEXP r_qG_,       // REAL
	    SEXP mu_aG_,
	    SEXP sd_aG_,
	    SEXP mu_rG_,
	    SEXP sd_rG_,
	    SEXP mu_beta_,    // REAL
	    SEXP evalue_sigma_betas_,
	    SEXP Inv_sigma_beta_,
	    SEXP burnin_,     // INTEGER
	    SEXP FreqEvery,   // INTEGER
	    SEXP patwoNS_      // INTEGER
	    )
{
  // Last modified:: Oct 26, 2012
  // The patients without new scans can be inputed.
  // === Model ===
  // Y_ij | Gij = gij ~ NB(size=exp(X_{ij}^T beta),prob=gij)
  // gij = g1 if j is in pre-scan 
  //     = g_new if j is in old-scan 
  // g_new = Ji * g1 + (1-Ji) * g2 
  // Ji ~ ber(qG)
  // qG ~ beta(a_qG,r_qG)
  // g1, g2 ~ beta(aG,rG)
  // beta ~ rnorm(mu_beta,sigma_beta)
  // aG, rG ~ unif(0,max_aG)
  const double *Y = REAL(Y_); 
  const double *labelnp = REAL(labelnp_);
  const double a_qG = REAL(a_qG_)[0];
  const double r_qG = REAL(r_qG_)[0];
  const int Ntot = length(ID); // # patients
  const int p = length(X)/Ntot; 
  const int N = INTEGER(nn)[0];
  const double mu_aG = REAL(mu_aG_)[0];
  const double sd_aG = REAL(sd_aG_)[0];
  const double mu_rG = REAL(mu_rG_)[0];
  const double sd_rG = REAL(sd_rG_)[0];
  const int B = INTEGER(B_)[0];
  const double *mu_betas = REAL(mu_beta_);
  const double *evalue_sigma_betas = REAL(evalue_sigma_betas_);
  const double *Inv_sigma_beta = REAL(Inv_sigma_beta_);
  const int burnin = INTEGER(burnin_)[0];
  const int *patwoNS = INTEGER(patwoNS_);// If all patients have new scans, then patwoNS[0]=-1000.
  int NpatwoNS = length(patwoNS_);  // The number of patients without new scans
  if (patwoNS[0]==-1000) NpatwoNS = 0; 
 
  int wID[INTEGER(maxni)[0]],wIDsize = 0, ih = 0, ip = 0, ipat = 0, ivec=0, idh=0, ibeta=0,iB=0;
  // wIDs[ipat][ivec] contains the position of the ivec^th repeated measure of the ipat^th patient at labelnp and Y
  //                  If ivec >= # repeated measures, wIDs[ipat][ivec] = -1000
  // wIDsizes[ipat] contains the number of repeated measure of ipat^th patient
  // sumYs[ipat][0] contains the sum of the CEL counts of ipat^th patients in the pre-scan period
  // sumYs[ipat][1] contains the sum of the CEL counts of ipat^th patients in the new-scan period
  int wIDs[N][INTEGER(maxni)[0]], wIDsizes[N];
  double sumYs[N][2];
  for (ipat = 0 ; ipat < N ; ipat++)
    {
      getwID( &wIDsize,wID,ID,ipat);
      wIDsizes[ipat] = wIDsize;
      sumYs[ipat][0] = 0;
      sumYs[ipat][1] = 0;
      for (ivec = 0 ; ivec < wIDsize; ivec++)
	{
	  wIDs[ipat][ivec]=wID[ivec];
	  if (labelnp[wIDs[ipat][ivec]]==0) sumYs[ipat][0] += Y[wIDs[ipat][ivec]];
	  if (labelnp[wIDs[ipat][ivec]]==1) sumYs[ipat][1] += Y[wIDs[ipat][ivec]];
	}
      for (ivec = wIDsize ; ivec < INTEGER(maxni)[0]; ivec++)
	wIDs[ipat][ivec]=-1000;
    }

  //h1s[N],h2s[N]
  double beta[p],canbeta[p],aG,rG,qG,g1s[N],g2s[N],js[N];
  // double vs[M],weightH1[M],postprob[M];
  double sizeij=0.0,cansizeij=0.0,w=0.0,canpij=0.0,pij=0.0,
    candi=0.0,MHrate=0.0,yij=0.0,sumgs=0.0;
  // double  weipro=0.0,qG2=0.0,D=0.0;
  double att[2+1], acc[2+1],can[2+1];
  for (ih = 0 ; ih < 3; ih ++)
    {
      att[ih] = 0.0;
      acc[ih] = 0.0;
      can[ih] = 5.0;
    }
  // http://stackoverflow.com/questions/8720550/how-to-return-array-of-structs-from-call-to-c-shared-library-in-r
  SEXP res = PROTECT(allocVector(VECSXP, 10)); // result is stored in res
  SEXP postaGs = allocVector(REALSXP, B); 
  SET_VECTOR_ELT(res, 0, postaGs); 
  SEXP postrGs = allocVector(REALSXP, B); 
  SET_VECTOR_ELT(res, 1, postrGs); 
  SEXP postqG = allocVector(REALSXP, B); 
  SET_VECTOR_ELT(res, 2, postqG);
 
  SEXP postgPre = allocVector(REALSXP, N*B); 
  SET_VECTOR_ELT(res, 3, postgPre); 
  SEXP postg2s = allocVector(REALSXP, N*B); 
  SET_VECTOR_ELT(res, 4, postg2s); 
  SEXP postgNew = allocVector(REALSXP, N*B); 
  SET_VECTOR_ELT(res, 5, postgNew); 
  SEXP postjs = allocVector(REALSXP, N*B); 
  SET_VECTOR_ELT(res, 6, postjs); 

  SEXP postbetas = allocVector(REALSXP, p*B); 
  SET_VECTOR_ELT(res, 7, postbetas); 
  SEXP AR = allocVector(REALSXP,3);
  SET_VECTOR_ELT(res, 8, AR);
  SEXP prp = allocVector(REALSXP,3);
  SET_VECTOR_ELT(res, 9, prp);

  GetRNGstate();
  
  // ============ initialization ============ 
  for (ip=0;ip < p ; ip++) beta[ip]=0.0;
  aG = 0.05;
  rG = 0.05;
  qG = 0.5;
  for (ipat = 0 ; ipat < N ; ipat++)
    {
      g1s[ipat]=0.5;
      g2s[ipat]=0.25;
      js[ipat]=1.0;
      for (ih = 0 ; ih < NpatwoNS; ih++) // If the patient does not have any new scans,g2s[ipat]= js[ipat]=-1000 during the whole algorithm.
	{
	  if (ipat==patwoNS[ih])
	    {
	      js[ipat] = -1000;
	      g2s[ipat] = -1000;
	    }
	}
    }
  // ================ MCMC =================
  for (iB = 0 ; iB < B; iB++)
    {
      // STEP 1 
      // Update g1s and g2s, a vector of length N. 
      // The g2s are updated ASSUMING THAT ji = 0 for all i to be used in updating ji
      for (ipat = 0 ; ipat < N; ipat++)
	{
	  w = 0.0;
	  sumgs = 0.0;
	  for (ivec=0; ivec < wIDsizes[ipat]; ivec ++)
	    {
	      sizeij = getsizeij(ivec, X, wIDs[ipat], beta, Ntot, p);
	      w += sizeij*(1.0-labelnp[wIDs[ipat][ivec]]); // sum_{j in pre} exp(Xij%*%beta)
	      sumgs += sizeij*labelnp[wIDs[ipat][ivec]]; // sum_{j in new} exp(Xij%*%beta)
	    }
	  if (js[ipat]==0) 
	    {
	      // update g1s
	      g1s[ipat] = rbeta(aG+w,rG+sumYs[ipat][0]);
	      // update g2s
	      g2s[ipat] = rbeta(aG+sumgs,rG+sumYs[ipat][1]);
	    }
	  if (js[ipat]==1) 
	    {
	      // update g1s
	      g1s[ipat] = rbeta(aG+w+sumgs,rG+sumYs[ipat][0]+sumYs[ipat][1]);
	      // update g2s
	      g2s[ipat] = rbeta(aG,rG);//prior

	    }
	  if (js[ipat]==-1000) // no new scans
	    {
	      g1s[ipat] = rbeta(aG+w,rG+sumYs[ipat][0]);
	      g2s[ipat] = -1000;
	    }
	  if (g1s[ipat] > 1.0-1.0e-2) g1s[ipat]= 1.0-1.0e-2;
	  if (g2s[ipat] > 1.0-1.0e-2) g2s[ipat]= 1.0-1.0e-2;
	}
      
      // Step 2: update js
      for (ipat = 0; ipat < N; ipat++)
        {
	  if (js[ipat]!=-1000){
	    w = qG;
	    sumgs = (1.0-qG) ;
	    for (ivec = 0 ; ivec < wIDsizes[ipat]; ivec++)
	      {
		if (labelnp[wIDs[ipat][ivec]]==1.0) // going through only the new scans
		  {
		    yij = Y[ wIDs[ipat][ivec] ];
		    sizeij = getsizeij(ivec, X, wIDs[ipat], beta, Ntot, p);
		    
		    // corresponds to Ji=1 gNew = gPre =g1s
		    canpij = getpij2(0.0,g1s[ipat],1.0);
		    w *= dnbinom(yij, sizeij, canpij, 0);
		    // corresponds to Ji=0 gNew = g2s != gPre =g1s
		    pij = getpij2(0.0,g2s[ipat],1.0);
		    sumgs *= dnbinom(yij, sizeij,  pij, 0);
		  }
	      }
	    js[ipat] = rbinom(1.0,w/(sumgs+w));
	    //if (iB == 5000) Rprintf(" \v ipat %d js[ipat] %f g1s[ipat] %f g2s[ipat] %f w %f sumgs %f",
	    //ipat, js[ipat], g1s[ipat], g2s[ipat], w, sumgs);
	  }
        }  

      //==


      // STEP 3: update beta, a vector of length p, containing the coefficients b0, b1,...
      // beta ~ MVR(mu_beta,Sigma_beta); 

      att[2]++;
      // proposal distribution
      for (ibeta = 0; ibeta <p ; ibeta ++)
	canbeta[ibeta] = rnorm(beta[ibeta], can[2]); // the mean and sd must be double (not NumericVector)
      
      // prior
      MHrate = dmvnorm(canbeta, mu_betas,p, evalue_sigma_betas,Inv_sigma_beta,1)
	- dmvnorm(beta, mu_betas,p, evalue_sigma_betas,Inv_sigma_beta,1);


      for (ipat=0; ipat < N; ipat++)
	{
	  // if (printing)  Rprintf(" \v\v ipat: %d",ipat);
	  for ( ivec = 0; ivec < wIDsizes[ipat]; ivec++ )
	    {
	      yij = Y[ wIDs[ipat][ivec] ];
	      // yij is a vector b/c it is the first arg of dnbinom
	      pij = getpij2(g1s[ipat],
			    (g1s[ipat]*js[ipat]+g2s[ipat]*(1.0-js[ipat])),
			    labelnp[wIDs[ipat][ivec]]); 
	      //Rprintf("\v\v ipat %d js[ipat] %f ivec %d yij %f g1s[ipat] %f g2s[ipat] %f pij %f", 
	      //ipat,  js[ipat] , ivec,   yij,   g1s[ipat], g2s[ipat],  pij);
	      
	      cansizeij = getsizeij(ivec, X, wIDs[ipat], canbeta, Ntot, p); // the value of cansizeij is changed
	      sizeij = getsizeij(ivec, X, wIDs[ipat], beta, Ntot, p); // the value of cansizeij is changed	      
	      MHrate += dnbinom(yij, cansizeij,  pij, 1)-dnbinom(yij, sizeij,  pij, 1); 
	    }
	  //Rprintf("\v");
	}    
      
      if(runif(0.0,1.0) < exp(MHrate))
	{
	  for (ibeta = 0; ibeta < p;ibeta++) beta[ibeta] = canbeta[ibeta];
	  acc[2]++;
	}
  
      // step 3: update qG
      w=0.0;
      for (ipat = 0 ; ipat < N ; ipat ++) 
	{
	  if (js[ipat]==1) w += js[ipat];
	}
      qG = rbeta(a_qG+w,r_qG+N-NpatwoNS-w);
      
      // step 4: Update aG and rG, scalars for the base distribution of the random effect gi
      // gi ~ Beta(aG,rG)
      // aG ~ unif(0,max_aG)
      // rG ~ unif(0,max_aG)

      // Update aG and rG, scalars for the base distribution of the random effect gi
      // gi ~ Beta(aG,rG)
      // log(aG) ~ norm(mu_aG,sigma_aG)
      // log(rG) ~ norm(mu_rG,sigma_rG)

      candi=rnorm(aG,can[0]);
      if (candi > 1.0e-7)
	{
	  att[0]++;
	  MHrate = dnorm(log(candi),mu_aG,sd_aG,1) - dnorm(log(aG),mu_aG,sd_aG,1); // prior
	  for (ipat=0 ; ipat < N; ipat++ ) 
	    MHrate += dbeta(g1s[ipat],candi,rG,1)-dbeta(g1s[ipat],aG,rG,1); 
	  if (!R_finite(MHrate)) error("MHrate is in-finite");
	  if (runif(0.0,1.0) < exp(MHrate) ) 
	    {
	      aG = candi;
	      acc[0]++;
	    }
	}
		      
      candi=rnorm(rG,can[1]);
      if (candi > 1.0e-7 )
	{
	  att[1]++;
	  MHrate = dnorm(log(candi),mu_rG,sd_rG,1) - dnorm(log(rG),mu_rG,sd_rG,1); // prior
	  for (ipat=0 ; ipat < N; ipat++ ) 
	    MHrate += dbeta(g1s[ipat],aG,candi,1)-dbeta(g1s[ipat],aG,rG,1);   
	  if (runif(0.0,1.0) < exp(MHrate) ) 
	    {
	      rG = candi;
	      acc[1]++;
	    }
	}
     
      // http://stackoverflow.com/questions/8720550/how-to-return-array-of-structs-from-call-to-c-shared-library-in-r
      // SAVE the result from this iteration
      REAL(postaGs)[iB] = aG;
      REAL(postrGs)[iB] = rG;
      REAL(postqG)[iB] = qG;

      for (ibeta = 0; ibeta < p; ibeta++) 
	{
	  idh = ibeta + iB*p;
	  REAL(postbetas)[idh] = beta[ibeta];
	}
      
      for (ipat = 0; ipat < N; ipat++ ) 
	{
	  idh = ipat + iB*N;
	  REAL(postgPre)[idh] = g1s[ipat];
	  REAL(postg2s)[idh] = g2s[ipat];
	  if (js[ipat]!=-1000) REAL(postgNew)[idh] = g1s[ipat]*js[ipat]+g2s[ipat]*(1.0-js[ipat]);
	  else REAL(postgNew)[idh] = -1000;
	  REAL(postjs)[idh] = js[ipat];
	}

      if (iB < burnin ) //2*iB
	{
	  for(ih=0; ih < 3 ; ih++)
	    {
	      // Keept the acceptance to range between 0.3 and 0.5
	      // The proposal variance is adjusted only during the burn-in period
	      if( ih==0 ) can[ih]=FisherBeta(aG,rG,N)*3; 
	      else if( ih==1 ) can[ih]=FisherBeta(rG,aG,N)*3; 
	      else 
		{
		  if (att[ih]>50)  can[ih]=adjustcan(acc[ih]/att[ih],can[ih],0.1,10);
		}
	    }
	 
	  /* if (iB % 100==0){ */
	  /*   for(ih=0; ih < 3; ih++){ */
	  /*     att[ih] = 0.0; */
	  /*     acc[ih] = 0.0; */
	  /*   } */
	}else if (iB==burnin){
	for( ih=0; ih < 3; ih++){
	  att[ih] = 0.0;
	  acc[ih] = 0.0;
	}
      }  
      
      if(iB > 0 && iB % INTEGER(FreqEvery)[0] == 0) Rprintf("\v %d iterations are done...   ", iB);
    }
  for (ih=0;ih < 3; ih ++) 
    {
      REAL(AR)[ih] = acc[ih]/att[ih];
      REAL(prp)[ih] = can[ih];
    }
  
  PutRNGstate();
  UNPROTECT(1);
  return res;
}







// extern "C" {
//   SEXP Beta24(SEXP Y_,          // REAL
// 	      SEXP X,           // REAL 
// 	      SEXP ID,          // INTEGER
// 	      SEXP B_,          // INTEGER
// 	      SEXP maxni,       // INTEGER
// 	      SEXP nn,          // INTEGER
// 	      SEXP M_,          // INTEGER
// 	      SEXP labelnp_,    // REAL
// 	      SEXP a_qG_,       // REAL
// 	      SEXP r_qG_,       // REAL
// 	      //SEXP max_aG_,   // REAL
// 	      SEXP mu_aG_,      // REAL
// 	      SEXP sd_aG_,      // REAL
// 	      SEXP mu_rG_,      // REAL
// 	      SEXP sd_rG_,      // REAL
// 	      SEXP mu_beta_,    // REAL
// 	      SEXP evalue_sigma_beta_, // REAL 
// 	      SEXP Inv_sigma_beta_,
// 	      SEXP a_D_,        // REAL
// 	      SEXP ib_D_,       // REAL
// 	      SEXP burnin_,     // INTEGER
// 	      SEXP FreqEvery,   // INTEGER
// 	      SEXP patwoNS_      // INTEGER 
// 	      );
// }


SEXP Beta24(SEXP Y_,          // REAL
	    SEXP X,           // REAL 
	    SEXP ID,          // INTEGER
	    SEXP B_,          // INTEGER
	    SEXP maxni,       // INTEGER
	    SEXP nn,          // INTEGER
	    SEXP M_,          // INTEGER
	    SEXP labelnp_,    // REAL
	    SEXP a_qG_,       // REAL
	    SEXP r_qG_,       // REAL
	    //SEXP max_aG_,   // REAL
	    SEXP mu_aG_,      // REAL
	    SEXP sd_aG_,      // REAL
	    SEXP mu_rG_,      // REAL
	    SEXP sd_rG_,      // REAL
	    SEXP mu_beta_,    // REAL
	    SEXP evalue_sigma_beta_, // REAL 
	    SEXP Inv_sigma_beta_,
	    SEXP a_D_,        // REAL
	    SEXP ib_D_,       // REAL
	    SEXP burnin_,     // INTEGER
	    SEXP FreqEvery,   // INTEGER
	    SEXP patwoNS_      // INTEGER 
	    )
{
  // Last modified:: Jan 6, 2013: the prior of beta is MVN(mu, Sigma_beta)
  // === Model ===
  // Nonparametric model with r.e changing between pre and new via Ji
  // Y_ij | Gij = gij ~ NB(size=exp(X_{ij}^T beta),prob=gij)
  // gij = g1 if j is in pre-scan 
  //     = g_new if j is in old-scan 
  // g_new = Ji * g1 + (1-Ji) * g2 
  // Ji ~ ber(qG)
  // qG ~ beta(a_qG,r_qG)
  // g1, g2 ~ sum_h=1^M pi_h beta(aG_h,rG_h)
  // beta ~ mvrnorm(mu_beta,sigma_beta)
  // log(aG) ~ norm(mu_aG,sd_aG) log(rG) ~ norm(mu_rG,sd_rG)
  const double *Y = REAL(Y_); 
  const double *labelnp = REAL(labelnp_);
  const double a_qG = REAL(a_qG_)[0];
  const double r_qG = REAL(r_qG_)[0];
  const int Ntot = length(ID); // # patients
  const int p = length(X)/Ntot; 
  const int N = INTEGER(nn)[0];
  const int M = INTEGER(M_)[0];
  //const double max_aG = REAL(max_aG_)[0];
  const double mu_aG = REAL(mu_aG_)[0];
  const double sd_aG= REAL(sd_aG_)[0];
  const double mu_rG= REAL(mu_rG_)[0];
  const double sd_rG= REAL(sd_rG_)[0];
  const int B = INTEGER(B_)[0];
  const double *mu_betas = REAL(mu_beta_); 
  const double *evalue_sigma_betas = REAL(evalue_sigma_beta_); // a vector of length p, containing eigen values of sigma (all positive)
  const double *Inv_sigma_beta = REAL(Inv_sigma_beta_); //  a vector of length p*p containing c(solve(sigma)[,1],solve(sigma)[,2],...)
  const int burnin = INTEGER(burnin_)[0];
  const double smallest = 1.0e-5;
  const double a_D = REAL(a_D_)[0];
  const double ib_D = REAL(ib_D_)[0];
  const int *patwoNS = INTEGER(patwoNS_);// If all patients have new scans, then patwoNS[0]=-1000.
  int NpatwoNS = length(patwoNS_);  // The number of patients without new scans
  if (patwoNS[0]==-1000) NpatwoNS = 0; 
  const int Nacc = M*2 + 1; // all covariates share the same proposal variance as beta ~ MVN
  // acc[0] = aGs[0],...,acc[M-1] = aGs[M-1],
  // acc[M] = rGs[0],...,acc[2M-1] = rGs[M-1],
  // acc[2M] = beta[0],
  int wID[INTEGER(maxni)[0]],wIDsize = 0, ih = 0, ip = 0, ipat = 0, ivec=0, idh=0, ibeta=0,iB=0;
  // wIDs[ipat][ivec] contains the position of the ivec^th repeated measure of the ipat^th patient at labelnp and Y
  //                  If ivec >= # repeated measures, wIDs[ipat][ivec] = -1000
  // wIDsizes[ipat] contains the number of repeated measure of ipat^th patient
  // sumYs[ipat][0] contains the sum of the CEL counts of ipat^th patients in the pre-scan period
  // sumYs[ipat][1] contains the sum of the CEL counts of ipat^th patients in the new-scan period
  int wIDs[N][INTEGER(maxni)[0]], wIDsizes[N];
  double sumYs[N][2];

  for (ipat = 0 ; ipat < N ; ipat++)
    {
      getwID( &wIDsize,wID,ID,ipat);
      wIDsizes[ipat] = wIDsize;
      sumYs[ipat][0] = 0;
      sumYs[ipat][1] = 0;
      for (ivec = 0 ; ivec < wIDsize; ivec++)
	{
	  wIDs[ipat][ivec]=wID[ivec];
	  if (labelnp[wIDs[ipat][ivec]]==0) sumYs[ipat][0] += Y[wIDs[ipat][ivec]];
	  if (labelnp[wIDs[ipat][ivec]]==1) sumYs[ipat][1] += Y[wIDs[ipat][ivec]];
	}
      for (ivec = wIDsize ; ivec < INTEGER(maxni)[0]; ivec++)
	wIDs[ipat][ivec]=-1000;
    }

  //h1s[N],h2s[N]
  double beta[p],canbeta[p],aGs[M],rGs[M],qG,g1s[N],g2s[N],js[N];
  double vs[M],weightH1[M],postprob[M];
  int h1s[N],h2s[N];
  double sizeij=0.0,cansizeij=0.0,w=0.0,canpij=0.0,pij=0.0,
    candi=0.0,MHrate=0.0,yij=0.0,sumgs=0.0;
  double  weipro=0.0,qG2=0.0,D=0.0;
  double att[Nacc], acc[Nacc],can[Nacc];
  for (ih = 0 ; ih < Nacc; ih ++)
    {
      att[ih] = 0.0;
      acc[ih] = 0.0;
      can[ih] = 5.0;
    }

  // http://stackoverflow.com/questions/8720550/how-to-return-array-of-structs-from-call-to-c-shared-library-in-r
  SEXP res = PROTECT(allocVector(VECSXP, 15)); // result is stored in res

  SEXP postqG = allocVector(REALSXP, B); 
  SET_VECTOR_ELT(res, 0, postqG);
  SEXP postD = allocVector(REALSXP, B); 
  SET_VECTOR_ELT(res, 1, postD);

  SEXP postgPre = allocVector(REALSXP, N*B); 
  SET_VECTOR_ELT(res, 2, postgPre); 
  SEXP postg2s = allocVector(REALSXP, N*B); 
  SET_VECTOR_ELT(res, 3, postg2s); 
  SEXP postgNew = allocVector(REALSXP, N*B); 
  SET_VECTOR_ELT(res, 4, postgNew); 
  SEXP postjs = allocVector(REALSXP, N*B); 
  SET_VECTOR_ELT(res, 5, postjs); 
  SEXP posth1s = allocVector(INTSXP, N*B);
  SET_VECTOR_ELT(res, 6, posth1s);
  SEXP posth2s = allocVector(INTSXP, N*B);
  SET_VECTOR_ELT(res, 7, posth2s);

  SEXP postweightH1 = allocVector(REALSXP,M*B);
  SET_VECTOR_ELT(res, 8, postweightH1);
  SEXP postvs = allocVector(REALSXP,M*B);
  SET_VECTOR_ELT(res, 9, postvs);
  SEXP postaGs = allocVector(REALSXP, M*B); 
  SET_VECTOR_ELT(res, 10, postaGs); 
  SEXP postrGs = allocVector(REALSXP, M*B); 
  SET_VECTOR_ELT(res, 11, postrGs); 

  SEXP postbetas = allocVector(REALSXP, p*B); 
  SET_VECTOR_ELT(res, 12, postbetas); 
  SEXP AR = allocVector(REALSXP,Nacc);
  SET_VECTOR_ELT(res, 13, AR);
  SEXP prp = allocVector(REALSXP,Nacc);
  SET_VECTOR_ELT(res, 14, prp);

  GetRNGstate();
  
  // ============ initialization ============ 
  for (ip=0;ip < p ; ip++) beta[ip]=0.0;
  for (ih=0;ih < M ; ih++)
    {
      aGs[ih] = exp(rnorm(mu_aG,sd_aG));
      rGs[ih] = exp(rnorm(mu_rG,sd_rG));
    }
  qG = 0.5;

  D = rgamma(a_D,1/ib_D);
  
  for (ih = 0 ; ih < (M-1) ; ih ++) vs[ih] = rbeta(1,D);
  vs[M-1]=1.0;

  weightH1[0] = vs[0];
  weipro = 1.0 ;
  for (ih=1; ih < M; ih++)
    {
      weipro = weipro*(1.0-vs[ih-1]);
      weightH1[ih] = vs[ih]*weipro;
    }

  for (ipat = 0 ; ipat < N ; ipat++)
    {
      h1s[ipat] = sample(weightH1,M);
      h2s[ipat] = sample(weightH1,M);
      g1s[ipat] = rbeta(aGs[h1s[ipat]],rGs[h1s[ipat]]);
      g2s[ipat]= rbeta(aGs[h2s[ipat]],rGs[h2s[ipat]]);
      if (g1s[ipat] > 1.0-smallest) g1s[ipat]= 1.0-smallest;
      if (g2s[ipat] > 1.0-smallest) g2s[ipat]= 1.0-smallest;
      js[ipat]=1.0;

      for (ih = 0 ; ih < NpatwoNS; ih++) 
	// If the patient does not have any new scans, h2s[ipat]=g2s[ipat]=js[ipat]=-1000 during the whole algorithm.
	{
	  if (ipat==patwoNS[ih])
	    {
	      js[ipat] = -1000;
	      g2s[ipat] = -1000;
	      h2s[ipat] = -1000;
	    }
	}
    }
  
  

  // ================ MCMC =================
  for (iB = 0 ; iB < B; iB++)
    {
      // STEP 1 and 2
      // Update g1s and g2s,  vectors of length N. 
      for (ipat = 0 ; ipat < N; ipat++)
	{
	  w = 0.0;
	  sumgs = 0.0;
	  for (ivec=0; ivec < wIDsizes[ipat]; ivec ++)
	    {
	      sizeij = getsizeij(ivec, X, wIDs[ipat], beta, Ntot, p);
	      w += sizeij*(1.0-labelnp[wIDs[ipat][ivec]]); // sum_{j in pre} exp(Xij%*%beta)
	      sumgs += sizeij*labelnp[wIDs[ipat][ivec]]; // sum_{j in new} exp(Xij%*%beta)
	    }
	  if (js[ipat]==0)  // gpre!=gnew
	    {
	      // update g1s
	      g1s[ipat] = rbeta(aGs[h1s[ipat]]+w,rGs[h1s[ipat]]+sumYs[ipat][0]);
	      // update g2s
	      g2s[ipat] = rbeta(aGs[h2s[ipat]]+sumgs,rGs[h2s[ipat]]+sumYs[ipat][1]);
	    }
	  if (js[ipat]==1) // gpre=gnew
	    {
	      // update g1s
	      g1s[ipat] = rbeta(aGs[h1s[ipat]]+w+sumgs,rGs[h1s[ipat]]+sumYs[ipat][0]+sumYs[ipat][1]);
	      // update g2s
	      g2s[ipat] = rbeta(aGs[h2s[ipat]],rGs[h2s[ipat]]);//prior

	    }
	  if (js[ipat]==-1000) // no new scans
	    {
	      g1s[ipat] = rbeta(aGs[h1s[ipat]]+w,rGs[h1s[ipat]]+sumYs[ipat][0]);
	      //g2s[ipat] = -1000; g2s[ipat] must already be -1000
	    }
	  if (g1s[ipat] > 1.0-smallest) g1s[ipat]= 1.0-smallest;
	  if (g2s[ipat] > 1.0-smallest) g2s[ipat]= 1.0-smallest;
	}
      
      // Step 3: update js, a vector of length N
      for (ipat = 0; ipat < N; ipat++)
        {
	  if (js[ipat]!=-1000)
	    {
	      w = qG;
	      sumgs = (1.0-qG) ;
	      for (ivec = 0 ; ivec < wIDsizes[ipat]; ivec++)
		{
		  if (labelnp[wIDs[ipat][ivec]]==1.0) // going through only the new scans
		    {
		      yij = Y[ wIDs[ipat][ivec] ];
		      sizeij = getsizeij(ivec, X, wIDs[ipat], beta, Ntot, p);
                  
		      // corresponds to Ji=1 gNew = gPre =g1s
		      canpij = getpij2(0.0,g1s[ipat],1.0);
		      w *= dnbinom(yij, sizeij, canpij, 0);
		      // corresponds to Ji=0 gNew = g2s != gPre =g1s
		      pij = getpij2(0.0,g2s[ipat],1.0);
		      sumgs *= dnbinom(yij, sizeij,  pij, 0);
		    }
		}
	      js[ipat] = rbinom(1.0,w/(sumgs+w));
	      //if (iB == 5000) Rprintf(" \v ipat %d js[ipat] %f g1s[ipat] %f g2s[ipat] %f w %f sumgs %f",
	      //ipat, js[ipat], g1s[ipat], g2s[ipat], w, sumgs);
	    }
        }  



      // STEP 4: update beta, a vector of length p, containing the coefficients b0, b1,...
      // beta ~ MVR(mu_beta,Sigma_beta); 

      att[2*M]++;
      // proposal distribution
      for (ibeta = 0; ibeta <p ; ibeta ++)
	canbeta[ibeta] = rnorm(beta[ibeta], can[2*M]); // the mean and sd must be double (not NumericVector)
      
      // prior
      MHrate = dmvnorm(canbeta, mu_betas,p, evalue_sigma_betas,Inv_sigma_beta,1)
	- dmvnorm(beta, mu_betas,p, evalue_sigma_betas,Inv_sigma_beta,1);


      for (ipat=0; ipat < N; ipat++)
	{
	  // if (printing)  Rprintf(" \v\v ipat: %d",ipat);
	  for ( ivec = 0; ivec < wIDsizes[ipat]; ivec++ )
	    {
	      yij = Y[ wIDs[ipat][ivec] ];
	      // yij is a vector b/c it is the first arg of dnbinom
	      pij = getpij2(g1s[ipat],
			    (g1s[ipat]*js[ipat]+g2s[ipat]*(1.0-js[ipat])),
			    labelnp[wIDs[ipat][ivec]]); 
	      //Rprintf("\v\v ipat %d js[ipat] %f ivec %d yij %f g1s[ipat] %f g2s[ipat] %f pij %f", 
	      //ipat,  js[ipat] , ivec,   yij,   g1s[ipat], g2s[ipat],  pij);
	      
	      cansizeij = getsizeij(ivec, X, wIDs[ipat], canbeta, Ntot, p); // the value of cansizeij is changed
	      
	      sizeij = getsizeij(ivec, X, wIDs[ipat], beta, Ntot, p); // the value of cansizeij is changed
	      
	      //Rprintf("canbeta %f cansizeij %f beta %f sizeij %f", canbeta[0],   cansizeij,   beta[0],   sizeij);
	      
	      MHrate += dnbinom(yij, cansizeij,  pij, 1)-dnbinom(yij, sizeij,  pij, 1); 
	      //Rprintf("\v dnbinom(yij, cansizeij,  pij, 1) %f",dnbinom(yij, cansizeij,  pij, 1));
	      //Rprintf("\v dnbinom(yij, sizeij,  pij, 1) %f",dnbinom(yij, sizeij,  pij, 1));
	      //Rprintf("\v summed MHrate %f", MHrate);
	    }
	  //Rprintf("\v");
	}    
      
      if(runif(0.0,1.0) < exp(MHrate))
	{
	  for (ibeta = 0; ibeta < p;ibeta++) beta[ibeta] = canbeta[ibeta];
	  acc[2*M]++;
	}
  
      // Step 5: update qG
      w=0.0;
      for (ipat = 0 ; ipat < N ; ipat ++) 
	{
	  if (js[ipat]==1)
	    w += js[ipat];
	}
      qG = rbeta(a_qG+w,r_qG+N-NpatwoNS-w);
      
      // Step 6 and 7: update aGs[ih] and rGs[ih], scalars for the base distribution of the random effect gi
      // gi ~ Beta(aG,rG)
      // log(aG[ih]) ~ N(mu_aG,sd_aG)
      // log(rG[ih]) ~ N(mu_rG,sd_rG)
      for (ih =0;ih < M;ih ++)
	{
	  candi=rnorm(aGs[ih],can[ih]);
	  if (candi > 1.0e-2) //&& candi < max_aG )
	    {
	      att[ih]++;
	      //MHrate = 0;
	      MHrate = dnorm(log(candi),mu_aG,sd_aG,1) - dnorm(log(aGs[ih]),mu_aG,sd_aG,1); // prior
	      for (ipat=0 ; ipat < N; ipat ++ ) 
		{
		  if (h1s[ipat]==ih) MHrate += dbeta(g1s[ipat],candi,rGs[ih],1)-dbeta(g1s[ipat],aGs[ih],rGs[ih],1); 
		  if (h2s[ipat]==ih&&js[ipat]!=-1000) MHrate += dbeta(g2s[ipat],candi,rGs[ih],1)-dbeta(g2s[ipat],aGs[ih],rGs[ih],1); 
		}
	      if (runif(0.0,1.0) < exp(MHrate) ) 
		{
		  aGs[ih] = candi;
		  acc[ih]++;
		}
	    }
      
	  candi=rnorm(rGs[ih],can[M+ih]);
	  if (candi > 1.0e-2 ) //&& candi < max_aG )
	    {
	      att[M+ih]++;
	      //MHrate = 0;
	      MHrate = dnorm(log(candi),mu_rG,sd_rG,1) - dnorm(log(rGs[ih]),mu_rG,sd_rG,1); // prior
	      for (ipat=0 ; ipat < N; ipat ++ ) 
		{
		  if (h1s[ipat]==ih)  MHrate += dbeta(g1s[ipat],aGs[ih],candi,1)-dbeta(g1s[ipat],aGs[ih],rGs[ih],1);   
		  if (h2s[ipat]==ih&&js[ipat]!=-1000) MHrate += dbeta(g2s[ipat],aGs[ih],candi,1)-dbeta(g2s[ipat],aGs[ih],rGs[ih],1); 
		}
	      if (runif(0.0,1.0) < exp(MHrate) ) 
		{
		  rGs[ih] = candi;
		  acc[M+ih]++;
		}
	    }
	}

      //  Steps 8 and 9: update h1s and h2s, vectors of length N, containing the cluster index of random effect 1 of each patient.
      //  the cluster index is between 1 and M.
      //  conjugate: posterior distribution of h1s has a categorical distribution with M categories. 
      for (ipat = 0 ; ipat < N ; ipat++)
	{
	  
	  // Given patient ipat, compute P(H_1i=ih;-) for each ih 
	  // P(H_1i=ih;-) = prod_{i=1}^{# repeated measure} dnbinom(yij; size=exp(b0+b1*x1i+..), prob=1/(gij+1))
	  for (ih=0; ih < M; ih++ ) postprob[ih] = dbeta(g1s[ipat],aGs[ih],rGs[ih],0)*weightH1[ih];
	  h1s[ipat] = sample(postprob,M); 
	  
	  if (js[ipat]!=-1000){
	    for (ih=0; ih < M; ih++ ) postprob[ih] = dbeta(g2s[ipat],aGs[ih],rGs[ih],0)*weightH1[ih];
	    h2s[ipat] = sample(postprob,M); 
	  }
	}

      // Step 10: update vs, a vector of length M, containing the latent variable v, used to construct pis. 
      // Once vs[ih] is updated, update the probabilities of the categorical distribution of H1 based on the formula:
      // weightH1[0] = vs[0]
      // weightH1[ih]= vs[ih]*(1-vh[ih-1])**(1-vh[ih-2])*...**(1-vh[0]) for ih =1,2,...M-1
      weipro = 1.0;
      for (ih = 0 ; ih < (M-1) ; ih ++) // vs[M-1] = 1 always!
	{
	  ip = 0;    // # patients with hs == ih
	  idh = 0;   // # patients with h2s == ih
	  ibeta = 0; // # patients with hs > ih
	  ivec = 0;  // # patients with h2s > ih
	  for (ipat = 0 ; ipat < N; ipat++ ) 
	    {
	      if ( h1s[ipat]==ih ) ip++ ;
	      if ( (h2s[ipat]==ih)&&js[ipat]!=-1000)  idh++; // h2s[ipat]=M+1 is not counted because they correspond to js[ipat]=-1000
	      if ( h1s[ipat] > ih ) ibeta++ ;
	      if ( (h2s[ipat] > ih)&&js[ipat]!=-1000)  ivec++;
	    }
	  
	  vs[ih] = rbeta(
			 (double) 1+ip+idh, 
			 D + (double)ibeta + (double)ivec
			 );
	  //Rprintf("ih %d sp1 %f sp2 %f",ih, (double) 1+ip+idh ,  D + (double)ibeta + (double)ivec);
	  if (ih > 0) weipro *= (1.0-vs[ih-1]);
	  weightH1[ih] = vs[ih]*weipro;
	}
      weightH1[M-1] = weipro*(1-vs[M-2]);


      // Step 11: update D
      // D ~ gamma(scale=a.D,shape=b.D) = gamma(scale=a.D,shape=1/ib.D)  
      w = 0.0;
      for (ih=0;ih<(M-1);ih++) w += log(1.0-vs[ih]);
      D = rgamma( (a_D+ (double) M-1.0), (1.0/(-w+ib_D)) ); // N+w is the number of random var following DP given js
      if (ISNAN(D)){ 
	D = 1.0; 
	Rprintf("\v ISNAN(D): a_D %f ib_D %f sum log(1-vs[ih]) %f ",a_D,ib_D,w);
      }
      if( D <0.01) D = 0.01; // if D is too small, rbeta(1,a,0) returns NaN for any a


      // http://stackoverflow.com/questions/8720550/how-to-return-array-of-structs-from-call-to-c-shared-library-in-r
      // SAVE the result from this iteration
      REAL(postqG)[iB] = qG;
      REAL(postD)[iB] = D;
      for (ih =0; ih < M ;ih++)
	{
	  idh = ih + iB*M;
	  REAL(postaGs)[idh] = aGs[ih];
	  REAL(postrGs)[idh] = rGs[ih];
	  REAL(postvs)[idh] = vs[ih];
	  REAL(postweightH1)[idh] =weightH1[ih];
	}
      for (ibeta = 0; ibeta < p; ibeta++) 
	{
	  idh = ibeta + iB*p;
	  REAL(postbetas)[idh] = beta[ibeta];
	}
      
      for (ipat = 0; ipat < N; ipat++ ) 
	{
	  idh = ipat + iB*N;
	  REAL(postgPre)[idh] = g1s[ipat];
	  REAL(postg2s)[idh] = g2s[ipat];
	  if (js[ipat]==-1000)
	    REAL(postgNew)[idh] = -1000;
	  else
	    REAL(postgNew)[idh] = g1s[ipat]*js[ipat]+g2s[ipat]*(1.0-js[ipat]);
	  REAL(postjs)[idh] = js[ipat];
	  INTEGER(posth1s)[idh] = h1s[ipat];
	  INTEGER(posth2s)[idh] = h2s[ipat];
	}

      if (iB < burnin ) //2*iB
	{
	  for(ih=0; ih < Nacc ; ih++)
	    {
	      // Keep the acceptance to range between 0.3 and 0.5
	      // The proposal variance is adjusted only during the burn-in period
	      if( ih < M ) can[ih]=FisherBeta(aGs[ih],rGs[ih],N)*3;  // update can of aGs[ih]
	      else if( M <= ih && ih < 2*M ) can[ih]=FisherBeta(rGs[ih-M],aGs[ih-M],N)*3;  // update can of rGs[ih]
	      else{ 
		if(att[ih]>50)  can[ih]=adjustcan(acc[ih]/att[ih],can[ih],0.05,10);
	      }
	      if (!R_FINITE(can[ih]))
		Rprintf("!");
	    }
	}else if (iB==burnin){
	for( ih=0; ih < Nacc; ih++){
	  att[ih] = 0.0;
	  acc[ih] = 0.0;
	}
      }  
      
      if(iB > 0 && iB % INTEGER(FreqEvery)[0] == 0) Rprintf("\v %d iterations are done...   ", iB);
    }
  for (ih=0;ih < Nacc; ih ++) 
    {
      REAL(AR)[ih] = acc[ih]/att[ih];
      REAL(prp)[ih] = can[ih];
    }
  
  PutRNGstate();
  UNPROTECT(1);
  return res;
}

// extern "C" {
//   SEXP Beta24Unif(SEXP Y_,          // REAL
// 		  SEXP X,           // REAL 
// 		  SEXP ID,          // INTEGER
// 		  SEXP B_,          // INTEGER
// 		  SEXP maxni,       // INTEGER
// 		  SEXP nn,          // INTEGER
// 		  SEXP M_,          // INTEGER
// 		  SEXP labelnp_,    // REAL
// 		  SEXP a_qG_,       // REAL
// 		  SEXP r_qG_,       // REAL
// 		  SEXP max_aG_,   // REAL
// 		  SEXP mu_beta_,    // REAL
// 		  SEXP evalue_sigma_beta_, // REAL 
// 		  SEXP Inv_sigma_beta_,
// 		  SEXP D_,        // REAL
// 		  SEXP burnin_,     // INTEGER
// 		  SEXP FreqEvery,   // INTEGER
// 		  SEXP patwoNS_,      // INTEGER 
// 		  SEXP initBeta_ // REAL
// 		  );
// }

SEXP Beta24Unif(SEXP Y_,          // REAL
		SEXP X,           // REAL 
		SEXP ID,          // INTEGER
		SEXP B_,          // INTEGER
		SEXP maxni,       // INTEGER
		SEXP nn,          // INTEGER
		SEXP M_,          // INTEGER
		SEXP labelnp_,    // REAL
		SEXP a_qG_,       // REAL
		SEXP r_qG_,       // REAL
		SEXP max_aG_,   // REAL
		//SEXP mu_aG_,      // REAL
		//SEXP sd_aG_,      // REAL
		//SEXP mu_rG_,      // REAL
		//SEXP sd_rG_,      // REAL
		SEXP mu_beta_,    // REAL
		SEXP evalue_sigma_beta_, // REAL 
		SEXP Inv_sigma_beta_,
		SEXP D_,        // REAL
		//SEXP ib_D_,       // REAL
		SEXP burnin_,     // INTEGER
		SEXP FreqEvery,   // INTEGER
		SEXP patwoNS_,      // INTEGER 
		SEXP initBeta_
		)
{
  // Last modified:: Jan 6, 2013: the prior of beta is MVN(mu, Sigma_beta)
  // === Model ===
  // Nonparametric model with r.e changing between pre and new via Ji
  // Y_ij | Gij = gij ~ NB(size=exp(X_{ij}^T beta),prob=gij)
  // gij = g1 if j is in pre-scan 
  //     = g_new if j is in old-scan 
  // g_new = Ji * g1 + (1-Ji) * g2 
  // Ji ~ ber(qG)
  // qG ~ beta(a_qG,r_qG)
  // g1, g2 ~ sum_h=1^M pi_h beta(aG_h,rG_h)
  // beta ~ mvrnorm(mu_beta,sigma_beta)
  // log(aG) ~ norm(mu_aG,sd_aG) log(rG) ~ norm(mu_rG,sd_rG)
  const double *Y = REAL(Y_); 
  const double *labelnp = REAL(labelnp_);
  const double a_qG = REAL(a_qG_)[0];
  const double r_qG = REAL(r_qG_)[0];
  const int Ntot = length(ID); // # patients
  const int p = length(X)/Ntot; 
  const int N = INTEGER(nn)[0];
  const int M = INTEGER(M_)[0];
  const double max_aG = REAL(max_aG_)[0];
  //const double mu_aG = REAL(mu_aG_)[0];
  //const double sd_aG= REAL(sd_aG_)[0];
  //const double mu_rG= REAL(mu_rG_)[0];
  //const double sd_rG= REAL(sd_rG_)[0];
  const int B = INTEGER(B_)[0];
  const double *mu_betas = REAL(mu_beta_); 
  const double *evalue_sigma_betas = REAL(evalue_sigma_beta_); // a vector of length p, containing eigen values of sigma (all positive)
  const double *Inv_sigma_beta = REAL(Inv_sigma_beta_); //  a vector of length p*p containing c(solve(sigma)[,1],solve(sigma)[,2],...)
  const int burnin = INTEGER(burnin_)[0];
  const double smallest = 1.0e-4,smallest2 = 0.01;
  const double max_D = REAL(D_)[0];
  //const double ib_D = REAL(ib_D_)[0];
  const int *patwoNS = INTEGER(patwoNS_);// If all patients have new scans, then patwoNS[0]=-1000.
  int NpatwoNS = length(patwoNS_);  // The number of patients without new scans
  if (patwoNS[0]==-1000) NpatwoNS = 0; 
  const int Nacc = M*2+2; // all covariates share the same proposal variance as beta ~ MVN
  // acc[0] = aGs[0],...,acc[M-1] = aGs[M-1],
  // acc[M] = rGs[0],...,acc[2M-1] = rGs[M-1],
  // acc[2M] = beta[0],
  // acc[2M+1] = D
  const double *initBeta = REAL(initBeta_);
  int wID[INTEGER(maxni)[0]],wIDsize = 0, ih = 0, ip = 0, ipat = 0, ivec=0, idh=0, ibeta=0,iB=0;
  // wIDs[ipat][ivec] contains the position of the ivec^th repeated measure of the ipat^th patient at labelnp and Y
  //                  If ivec >= # repeated measures, wIDs[ipat][ivec] = -1000
  // wIDsizes[ipat] contains the number of repeated measure of ipat^th patient
  // sumYs[ipat][0] contains the sum of the CEL counts of ipat^th patients in the pre-scan period
  // sumYs[ipat][1] contains the sum of the CEL counts of ipat^th patients in the new-scan period
  int wIDs[N][INTEGER(maxni)[0]], wIDsizes[N];
  double sumYs[N][2];

  for (ipat = 0 ; ipat < N ; ipat++)
    {
      getwID( &wIDsize,wID,ID,ipat);
      wIDsizes[ipat] = wIDsize;
      sumYs[ipat][0] = 0;
      sumYs[ipat][1] = 0;
      for (ivec = 0 ; ivec < wIDsize; ivec++)
	{
	  wIDs[ipat][ivec]=wID[ivec];
	  if (labelnp[wIDs[ipat][ivec]]==0) sumYs[ipat][0] += Y[wIDs[ipat][ivec]];
	  if (labelnp[wIDs[ipat][ivec]]==1) sumYs[ipat][1] += Y[wIDs[ipat][ivec]];
	}
      for (ivec = wIDsize ; ivec < INTEGER(maxni)[0]; ivec++)
	wIDs[ipat][ivec]=-1000;
    }

  //h1s[N],h2s[N]
  double beta[p],canbeta[p],aGs[M],rGs[M],qG,g1s[N],g2s[N],js[N];
  double vs[M],weightH1[M],postprob[M];
  int h1s[N],h2s[N];
  double sizeij=0.0,cansizeij=0.0,w=0.0,canpij=0.0,pij=0.0,
    candi=0.0,MHrate=0.0,yij=0.0,sumgs=0.0;
  double  weipro=0.0,qG2=0.0,D=0.0;
  double att[Nacc], acc[Nacc],can[Nacc];
  double logLcan, logL; // focused likelihood values
  for (ih = 0 ; ih < Nacc; ih ++)
    {
      att[ih] = 0.0;
      acc[ih] = 0.0;
      can[ih] = 5.0;
    }
  const double min_can[2] = {0.05, 0.35}; // minimum of the proposal variance for beta and D

  // http://stackoverflow.com/questions/8720550/how-to-return-array-of-structs-from-call-to-c-shared-library-in-r
  SEXP res = PROTECT(allocVector(VECSXP, 16)); // result is stored in res

  SEXP postqG = allocVector(REALSXP, B); 
  SET_VECTOR_ELT(res, 0, postqG);
  SEXP postD = allocVector(REALSXP, B); 
  SET_VECTOR_ELT(res, 1, postD);

  SEXP loglikL = allocVector(REALSXP, B); 
  SET_VECTOR_ELT(res, 2, loglikL);

  SEXP postgPre = allocVector(REALSXP, N*B); 
  SET_VECTOR_ELT(res, 3, postgPre); 
  SEXP postg2s = allocVector(REALSXP, N*B); 
  SET_VECTOR_ELT(res, 4, postg2s); 
  SEXP postgNew = allocVector(REALSXP, N*B); 
  SET_VECTOR_ELT(res, 5, postgNew); 
  SEXP postjs = allocVector(REALSXP, N*B); 
  SET_VECTOR_ELT(res, 6, postjs); 
  SEXP posth1s = allocVector(INTSXP, N*B);
  SET_VECTOR_ELT(res, 7, posth1s);
  SEXP posth2s = allocVector(INTSXP, N*B);
  SET_VECTOR_ELT(res, 8, posth2s);

  SEXP postweightH1 = allocVector(REALSXP,M*B);
  SET_VECTOR_ELT(res, 9, postweightH1);
  SEXP postvs = allocVector(REALSXP,M*B);
  SET_VECTOR_ELT(res, 10, postvs);
  SEXP postaGs = allocVector(REALSXP, M*B); 
  SET_VECTOR_ELT(res, 11, postaGs); 
  SEXP postrGs = allocVector(REALSXP, M*B); 
  SET_VECTOR_ELT(res, 12, postrGs); 

  SEXP postbetas = allocVector(REALSXP, p*B); 
  SET_VECTOR_ELT(res, 13, postbetas); 
  SEXP AR = allocVector(REALSXP,Nacc);
  SET_VECTOR_ELT(res, 14, AR);
  SEXP prp = allocVector(REALSXP,Nacc);
  SET_VECTOR_ELT(res, 15, prp);

  GetRNGstate();
  
  // ============ initialization ============ 
  // Initial values of the coefficient beta is given
  for (ip=0;ip < p ; ip++) beta[ip]=initBeta[ip];
  for (ih=0;ih < M ; ih++)
    {
      aGs[ih] = runif(smallest2,max_aG);
      rGs[ih] = runif(smallest2,max_aG);
    }
  qG = 0.5;
 
  D = runif(0.01,max_D);
  
  for (ih = 0 ; ih < (M-1) ; ih ++) vs[ih] = rbeta(1,D);
  vs[M-1]=1.0;

  weightH1[0] = vs[0];
  weipro = 1.0 ;
  for (ih=1; ih < M; ih++)
    {
      weipro = weipro*(1.0-vs[ih-1]);
      weightH1[ih] = vs[ih]*weipro;
    }

  for (ipat = 0 ; ipat < N ; ipat++)
    {
      h1s[ipat] = sample(weightH1,M);
      h2s[ipat] = sample(weightH1,M);
      g1s[ipat] = rbeta(aGs[h1s[ipat]],rGs[h1s[ipat]]);
      g2s[ipat]= rbeta(aGs[h2s[ipat]],rGs[h2s[ipat]]);
      if (g1s[ipat] > 1.0-smallest) g1s[ipat]= 1.0-smallest;
      if (g2s[ipat] > 1.0-smallest) g2s[ipat]= 1.0-smallest;
      js[ipat]=1.0;

      for (ih = 0 ; ih < NpatwoNS; ih++) 
	// If the patient does not have any new scans, h2s[ipat]=g2s[ipat]=js[ipat]=-1000 during the whole algorithm.
	{
	  if (ipat==patwoNS[ih])
	    {
	      js[ipat] = -1000;
	      g2s[ipat] = -1000;
	      h2s[ipat] = -1000;
	    }
	}
    }

  // ================ MCMC =================
  for (iB = 0 ; iB < B; iB++)
    {
      R_CheckUserInterrupt(); 
      //R_ProcessEvents();
      // STEP 1 and 2
      // Update g1s and g2s,  vectors of length N. 
      for (ipat = 0 ; ipat < N; ipat++)
	{
	  w = 0.0;
	  sumgs = 0.0;
	  for (ivec=0; ivec < wIDsizes[ipat]; ivec ++)
	    {
	      sizeij = getsizeij(ivec, X, wIDs[ipat], beta, Ntot, p);
	      w += sizeij*(1.0-labelnp[wIDs[ipat][ivec]]); // sum_{j in pre} exp(Xij%*%beta)
	      sumgs += sizeij*labelnp[wIDs[ipat][ivec]]; // sum_{j in new} exp(Xij%*%beta)
	    }
	  if (js[ipat]==0)  // gpre!=gnew
	    {
	      // update g1s
	      g1s[ipat] = rbeta(aGs[h1s[ipat]]+w,rGs[h1s[ipat]]+sumYs[ipat][0]);
	      // update g2s
	      g2s[ipat] = rbeta(aGs[h2s[ipat]]+sumgs,rGs[h2s[ipat]]+sumYs[ipat][1]);
	    }
	  if (js[ipat]==1) // gpre=gnew
	    {
	      // update g1s
	      g1s[ipat] = rbeta(aGs[h1s[ipat]]+w+sumgs,rGs[h1s[ipat]]+sumYs[ipat][0]+sumYs[ipat][1]);
	      // update g2s
	      g2s[ipat] = rbeta(aGs[h2s[ipat]],rGs[h2s[ipat]]);//prior

	    }
	  if (js[ipat]==-1000) // no new scans
	    {
	      g1s[ipat] = rbeta(aGs[h1s[ipat]]+w,rGs[h1s[ipat]]+sumYs[ipat][0]);
	      //g2s[ipat] = -1000; g2s[ipat] must already be -1000
	    }
	  if (g1s[ipat] > 1.0-smallest) g1s[ipat]= 1.0-smallest;
	  if (g2s[ipat] > 1.0-smallest) g2s[ipat]= 1.0-smallest;
	}
      
      // Step 3: update js, a vector of length N
      for (ipat = 0; ipat < N; ipat++)
        {
	  if (js[ipat]!=-1000)
	    {
	      w = qG;
	      sumgs = (1.0-qG) ;
	      for (ivec = 0 ; ivec < wIDsizes[ipat]; ivec++)
		{
		  if (labelnp[wIDs[ipat][ivec]]==1.0) // going through only the new scans
		    {
		      yij = Y[ wIDs[ipat][ivec] ];
		      sizeij = getsizeij(ivec, X, wIDs[ipat], beta, Ntot, p);
                  
		      // corresponds to Ji=1 gNew = gPre =g1s
		      canpij = getpij2(0.0,g1s[ipat],1.0);
		      w *= dnbinom(yij, sizeij, canpij, 0);
		      // corresponds to Ji=0 gNew = g2s != gPre =g1s
		      pij = getpij2(0.0,g2s[ipat],1.0);
		      sumgs *= dnbinom(yij, sizeij,  pij, 0);
		    }
		}
	      js[ipat] = rbinom(1.0,w/(sumgs+w));
	      //if (iB == 5000) Rprintf(" \v ipat %d js[ipat] %f g1s[ipat] %f g2s[ipat] %f w %f sumgs %f",
	      //ipat, js[ipat], g1s[ipat], g2s[ipat], w, sumgs);
	    }
        }  

      // STEP 4: update beta, a vector of length p, containing the coefficients b0, b1,...
      // beta ~ MVR(mu_beta,Sigma_beta); 

      att[2*M]++;
      // proposal distribution
      for (ibeta = 0; ibeta <p ; ibeta ++)
	canbeta[ibeta] = rnorm(beta[ibeta], can[2*M]); // the mean and sd must be double (not NumericVector)
      
      // prior
      MHrate = dmvnorm(canbeta, mu_betas,p, evalue_sigma_betas,Inv_sigma_beta,1) 
	- dmvnorm(beta, mu_betas,p, evalue_sigma_betas,Inv_sigma_beta,1);
      logL = 0.0; // the value of "focused" logLikelihood with current beta and the proposed beta are stored to evaluated the value of DIC
      logLcan = 0.0;
      for (ipat=0; ipat < N; ipat++)
	{
	  // if (printing)  Rprintf(" \v\v ipat: %d",ipat);
	  for ( ivec = 0; ivec < wIDsizes[ipat]; ivec++ )
	    {
	      yij = Y[ wIDs[ipat][ivec] ];
	      // yij is a vector b/c it is the first arg of dnbinom
	      pij = getpij2(g1s[ipat],
			    (g1s[ipat]*js[ipat]+g2s[ipat]*(1.0-js[ipat])),
			    labelnp[wIDs[ipat][ivec]]); 
	      //Rprintf("\v\v ipat %d js[ipat] %f ivec %d yij %f g1s[ipat] %f g2s[ipat] %f pij %f", 
	      //ipat,  js[ipat] , ivec,   yij,   g1s[ipat], g2s[ipat],  pij);
	      
	      cansizeij = getsizeij(ivec, X, wIDs[ipat], canbeta, Ntot, p); // the value of cansizeij is changed
	      
	      sizeij = getsizeij(ivec, X, wIDs[ipat], beta, Ntot, p); // the value of cansizeij is changed

	      //Rprintf("canbeta %f cansizeij %f beta %f sizeij %f", canbeta[0],   cansizeij,   beta[0],   sizeij);
	      logLcan += dnbinom(yij, cansizeij,  pij, 1);
	      logL += dnbinom(yij, sizeij,  pij, 1);

	      //Rprintf("\v dnbinom(yij, cansizeij,  pij, 1) %f",dnbinom(yij, cansizeij,  pij, 1));
	      //Rprintf("\v dnbinom(yij, sizeij,  pij, 1) %f",dnbinom(yij, sizeij,  pij, 1));
	      //Rprintf("\v summed MHrate %f", MHrate);
	    }
	  //Rprintf("\v");
	}    
      MHrate += logLcan - logL;
      if(runif(0.0,1.0) < exp(MHrate))
	{
	  for (ibeta = 0; ibeta < p;ibeta++) beta[ibeta] = canbeta[ibeta];
	  acc[2*M]++;
	  logL = logLcan;
	}
  
      // Step 5: update qG
      w=0.0;
      for (ipat = 0 ; ipat < N ; ipat ++) 
	{
	  if (js[ipat]==1)
	    w += js[ipat];
	}
      qG = rbeta(a_qG+w,r_qG+N-NpatwoNS-w);
      

      //Update aGs and rGs, scalars for the base distribution of the random effect gi
      //gi ~ Beta(aG,rG)
      //aG ~ unif(0,max_aG)
      //rG ~ unif(0,max_aG)
      for (ih = 0 ; ih < M ; ih ++)
      	{
      	  candi=rnorm(aGs[ih],can[ih]);
      	  att[ih]++;
      	  if (candi > smallest2 && candi < max_aG )
      	    {
      	      MHrate = 0;
      	      for (ipat=0 ; ipat < N; ipat ++ )  
      		{
      		  if (h1s[ipat]==ih)
      		    MHrate += dbeta(g1s[ipat],candi,rGs[ih],1)-dbeta(g1s[ipat],aGs[ih],rGs[ih],1); 
      		}
      	      if (runif(0.0,1.0) < exp(MHrate) ) 
      		{
      		  aGs[ih] = candi;
      		  acc[ih]++;
      		}
      	    }	  
      	  candi=rnorm(rGs[ih],can[ih+M]);
      	  att[M+ih]++;
      	  if (candi > smallest2 && candi < max_aG )
      	    {
      	      MHrate = 0;
      	      for (ipat=0 ; ipat < N; ipat ++ ) 
      		{
      		  if (h1s[ipat]==ih) 
      		    MHrate += dbeta(g1s[ipat],aGs[ih],candi,1)-dbeta(g1s[ipat],aGs[ih],rGs[ih],1);  
      		}
      	      if (runif(0.0,1.0) < exp(MHrate) ) 
      		{
      		  rGs[ih] = candi;
      		  acc[M+ih]++;
      		}
      	    }
      	}
     

      //  Steps 8 and 9: update h1s and h2s, vectors of length N, containing the cluster index of random effect 1 of each patient.
      //  the cluster index is between 1 and M.
      //  conjugate: posterior distribution of h1s has a categorical distribution with M categories. 
      for (ipat = 0 ; ipat < N ; ipat++)
	{
	  
	  // Given patient ipat, compute P(H_1i=ih;-) for each ih 
	  // P(H_1i=ih;-) = prod_{i=1}^{# repeated measure} dnbinom(yij; size=exp(b0+b1*x1i+..), prob=1/(gij+1))
	  for (ih=0; ih < M; ih++ ) postprob[ih] = dbeta(g1s[ipat],aGs[ih],rGs[ih],0)*weightH1[ih];
	  h1s[ipat] = sample(postprob,M); 
	  
	  if (js[ipat]!=-1000){
	    for (ih=0; ih < M; ih++ ) postprob[ih] = dbeta(g2s[ipat],aGs[ih],rGs[ih],0)*weightH1[ih];
	    h2s[ipat] = sample(postprob,M); 
	  }
	}

      // Step 10: update vs, a vector of length M, containing the latent variable v, used to construct pis. 
      // Once vs[ih] is updated, update the probabilities of the categorical distribution of H1 based on the formula:
      // weightH1[0] = vs[0]
      // weightH1[ih]= vs[ih]*(1-vh[ih-1])**(1-vh[ih-2])*...**(1-vh[0]) for ih =1,2,...M-1
      weipro = 1.0;
      for (ih = 0 ; ih < (M-1) ; ih ++) // vs[M-1] = 1 always!
	{
	  ip = 0;    // # patients with hs == ih
	  idh = 0;   // # patients with h2s == ih
	  ibeta = 0; // # patients with hs > ih
	  ivec = 0;  // # patients with h2s > ih
	  for (ipat = 0 ; ipat < N; ipat++ ) 
	    {
	      if ( h1s[ipat]==ih ) ip++ ;
	      if ( (h2s[ipat]==ih)&&js[ipat]!=-1000)  idh++; // h2s[ipat]=M+1 is not counted because they correspond to js[ipat]=-1000
	      if ( h1s[ipat] > ih ) ibeta++ ;
	      if ( (h2s[ipat] > ih)&&js[ipat]!=-1000)  ivec++;
	    }
	  
	  vs[ih] = rbeta(
			 (double) 1+ip+idh, 
			 D + (double)ibeta + (double)ivec
			 );
	  //Rprintf("ih %d sp1 %f sp2 %f",ih, (double) 1+ip+idh ,  D + (double)ibeta + (double)ivec);
	  if (ih > 0) weipro *= (1.0-vs[ih-1]);
	  weightH1[ih] = vs[ih]*weipro;
	}
      weightH1[M-1] = weipro*(1-vs[M-2]);


      // Step 11: update D
      // D ~ runif(0.01,maxD) //gamma(scale=a.D,shape=b.D) = gamma(scale=a.D,shape=1/ib.D)  
      candi = rnorm(D,can[2*M+1]);

      att[2*M+1]++;
      if (0.01 < candi && candi < max_D)
	{
	  MHrate = 0.0;
	  for (ih = 0 ; ih < (M-1);ih++)
	    {
	      MHrate += dbeta(vs[ih],1.0,candi,1)-dbeta(vs[ih],1.0,D,1);
	    }
	  if (runif(0.0,1.0) < exp(MHrate) ) 
	    {
	      D = candi;
	      acc[2*M+1]++;
	    }
	}

      // http://stackoverflow.com/questions/8720550/how-to-return-array-of-structs-from-call-to-c-shared-library-in-r
      // SAVE the result from this iteration
      REAL(postqG)[iB] = qG;
      REAL(postD)[iB] = D;
      REAL(loglikL)[iB] = logL;
      for (ih =0; ih < M ;ih++)
	{
	  idh = ih + iB*M;
	  REAL(postaGs)[idh] = aGs[ih];
	  REAL(postrGs)[idh] = rGs[ih];
	  REAL(postvs)[idh] = vs[ih];
	  REAL(postweightH1)[idh] = weightH1[ih];
	}
      for (ibeta = 0; ibeta < p; ibeta++) 
	{
	  idh = ibeta + iB*p;
	  REAL(postbetas)[idh] = beta[ibeta];
	}
      
      for (ipat = 0; ipat < N; ipat++ ) 
	{
	  idh = ipat + iB*N;
	  REAL(postgPre)[idh] = g1s[ipat];
	  REAL(postg2s)[idh] = g2s[ipat];
	  if (js[ipat]==-1000)
	    REAL(postgNew)[idh] = -1000;
	  else
	    REAL(postgNew)[idh] = g1s[ipat]*js[ipat]+g2s[ipat]*(1.0-js[ipat]);
	  REAL(postjs)[idh] = js[ipat];
	  INTEGER(posth1s)[idh] = h1s[ipat];
	  INTEGER(posth2s)[idh] = h2s[ipat];
	}
      
      if (iB < burnin ) //2*iB
	{
	  for(ih=0; ih < Nacc ; ih++)
	    {
	      // Keep the acceptance to range between 0.3 and 0.5
	      // The proposal variance is adjusted only during the burn-in period
	      if( ih < M ) can[ih]=FisherBeta(aGs[ih],rGs[ih],N)*3;  // update can of aGs[ih]
	      else if( M <= ih && ih < 2*M ) can[ih]=FisherBeta(rGs[ih-M],aGs[ih-M],N)*3;  // update can of rGs[ih]
	      else{ 
		// beta and D 
		if(att[ih]>50)  can[ih]=adjustcan(acc[ih]/att[ih],can[ih],min_can[ih-2*M],10);
	      }
	      if (!R_FINITE(can[ih]))
		Rprintf("!");
	    }
	}else if (iB==burnin){
	for( ih=0; ih < Nacc; ih++){
	  att[ih] = 0.0;
	  acc[ih] = 0.0;
	}
      }  
      
      if(iB > 0 && iB % INTEGER(FreqEvery)[0] == 0)
	{
	  Rprintf("\v %d iterations are done...   ", iB);
	  R_FlushConsole(); 
	  //R_ProcessEvents();
	}
    }
  for (ih=0;ih < Nacc; ih ++) 
    {
      REAL(AR)[ih] = acc[ih]/att[ih];
      REAL(prp)[ih] = can[ih];
    }
  
  PutRNGstate();
  UNPROTECT(1);
  return res;
}

// // declarations (necessary for C++ codes)
// extern "C" {
//   SEXP map_c(SEXP ma_, SEXP Npat_, SEXP B_);
// }




SEXP map_c(SEXP ma_, SEXP Npat_, SEXP B_)
{
  const int B = INTEGER(B_)[0];
  const int Npat = INTEGER(Npat_)[0];
  const int *ma = INTEGER(ma_); // B by Npat
  SEXP res;
  PROTECT(res = allocMatrix(INTSXP,Npat,Npat)); // Npat by Npat
  int count = 0;
  for (int iobs1=0; iobs1 < Npat; iobs1++)
    {
      for (int iobs2=iobs1+1; iobs2 < Npat; iobs2++)
	{
	  count = 0;
	  for (int iB=0; iB < B; iB++ )
	    {
	      if ( ma[iB+B*iobs1] == ma[iB+B*iobs2] ) 
		count++;
	    }
	  INTEGER(res)[iobs1+Npat*iobs2]=count;
	  INTEGER(res)[iobs2+Npat*iobs1]=count;
	}
    }
  UNPROTECT(1);
  return (res);
}

// // declarations (necessary for C++ codes)
// extern "C" {
//   SEXP pG1LeG2_c( SEXP inpg1s,SEXP inpg2s);
// }



SEXP pG1LeG2_c( SEXP inpg1s,SEXP inpg2s)
{

  const double *g1s = REAL(inpg1s);
  const double *g2s = REAL(inpg2s);
  // g1s: A vector of length B, containing B samples of g1 in the decreasing order 
  // g2s: A vector of length B, containing B samples of g2 in the decreasing order 
  // The goal of this function is to compute sum_{ib1=1}^B sum_{ib2=1}^B I{ g1_(ib1) <= g2_(ib2) }/(B*B)
  int ig2, prev_ig2, breaked = 0;
  int B = length(inpg1s);
  double temp=0.0;
  SEXP sumig2s;
  PROTECT(sumig2s = allocVector(REALSXP,1)); 

  prev_ig2=0;
  for (int ig1=0; ig1 < B; ig1++)
    {
      breaked = 0;
      for (ig2=prev_ig2; ig2 < B; ig2++)
	{
	  // given g1s[ig1],
	  // search for the smallest ig2 such that g1s[ig1] <= g2s[ig2]  
	  // then g1s[ig1] <= g2s[ig2], g1s[ig1] <= g2s[ig2+1],...,g1s[ig1] <= g2s[B-1]
	  // because g2s is sorted in the increasing way.
	  // That means (B-1)-ig2+1 elements of vector g2s are greater than g1s  
	  // i.e., (B-ig2) = sum_{ib2=1}^B I{ g1_(ib1) <= g2_(ib2) } 
	  if (g1s[ig1] <= g2s[ig2]){
	    breaked = 1;
	    break;  
	  }
	}
      // breaked = 0 means that for the given g1s[ig1], there is no element of g2s
      // that is greater than g1s[ig1].
      // This implies that g1s[ig1+1], g1s[ig2+2],...,g1s[B-1] are all greater than any 
      // values of g2s, because g1s is sorted in the increasing way.
      if ((ig2== B-1) && (breaked==0)) break;
      temp += B - ig2;
      prev_ig2 = ig2;
    }
  REAL(sumig2s)[0] = temp/(B*B);

  UNPROTECT(1);
  return(sumig2s);
  
}



/* SEXP test() */
/* { */


/*   double db = dnbinom(100, // response */
/* 		      5, // size */
/* 		      0.01, */
/* 		      1 // log */
/* 		      ); */
/*   // INF!! */
/*   Rprintf(" \v %f",db); */
/*   return R_NilValue; */
/* } */
