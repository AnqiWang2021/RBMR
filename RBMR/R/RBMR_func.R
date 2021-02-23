#' RBMR function for only one LD block
#'
#' When the number of LD block is 1, this function estimates the causal effect.
#'
#' @param bh1 A vector of SNP effects on the exposure variable, usually obtained from a GWAS.
#' @param bh2 A vector of SNP effects on the outcome varible, usually obtained from a GWAS.
#' @param se1 A vector of standard errors of \code{bh1}.
#' @param se2 A vector of standard errors of \code{bh2}.
#' @param R A matrix, LD matrix.
#' @param gamma A vector to intialize the true value of the mean of the exposure.
#' @param alpha A vector to intialize the true value of the mean of the direct effect.
#' @param sgga2 A constant to initialize the true value of the standard error of exposure \code{bh1}.
#' @param sgal2 A constant to initialize the true value of the standard error of direct effect \code{bh2}.
#' @param beta0 A constant, initialize the true value of causal effect.
#' @param constr 0 or 1, when constr is equal to 0, the function calculates the ELBO under the alternative hypothesis, when constr is equal to 1, the function calculates the ELBO under the null hypothesis.
#' @param alphag A constant to initialize the parameter of multivariate generalize t distribution
#' @param betag A constant to initialize the parameter of multivariate generalize t distribution
#' @param maxIter Maximum number of interations to solve the estimating equations.
#' @param epsStopLogLik Numerical precision.
#'
#' @param
#'
#'
#' @return A list
#' \describe{
#' \item{beta0}{Estimated causal effect.}
#' \item{sgal2}{Standard error of the direct effec.}
#' \item{sgga2}{Standard error of the exposure.}
#' \item{Iteration}{Number of iterations.}
#' \item{loglik}{The result of loglikelihood at convergence.}
#' \item{diff}{The difference between the value of loglikelihood at convergence and the value of the loglikelihood of previous itertaion.}
#' \item{tstat}{The value of the ELBO.}
#' \item{alpha_w}{The parameter of multivariate generalized t distribution.}
#' \item{beta_w}{The parameter of multivariate generalized t distribution.}
#' }
#' @export
#'

RBMR_func<-function(bh1,bh2,se1,se2,gamma,alpha,sgga2,sgal2,beta0,R,alphag,betag,constr,epsStopLogLik,maxIter)
{
  p<-length(bh1)

  sg2<-se1^2
  sG2<-se2^2

  GinvsG2<-bh2/sG2
  ginvsg2<-bh1/sg2


  insGRinsG<-diag(1/se2) %*% R %*% diag(1/se2)
  insgRinsg<-diag(1/se1) %*% R %*% diag(1/se1)

  diaginsGRinsG<-diag(insGRinsG)
  diaginsgRinsg<-diag(insgRinsg)

  Rins<-R%*%diag(1/se1)
  Rins2<-R%*%diag(1/se2)
  RinsGmu<-as.vector(R%*%diag(1/se2)%*%mu)
  Rinsgmu<-as.vector(R%*%diag(1/se1)%*%mu)
  RinsGmuA<-as.vector(R%*%diag(1/se2)%*%muA)

  v2<-vector()
  v2A<-vector()

  loglik<-rep(0,maxIter)

  Iteration<-1
  mu<-gamma
  muA<-alpha
  alpha_w<-alphag
  beta_w<-betag

  xi<-1

  for(iter in 1:maxIter){


    v2<-1/((beta0^2)*diaginsGRinsG+(xi^2)*diaginsgRinsg+rep(1/sgga2,p));

    for(j in 1:p){
      tmp1<-Rinsgmu-Rins[,j]*mu[j]
      tmp2<-RinsGmu-Rins2[,j]*mu[j]

      RinSmujj<-Rinsgmu[j]-R[j,j]*mu[j]/se1[j]
      RinSmujj2<-RinsGmu[j]-R[j,j]*mu[j]/se2[j]
      RinSmuAjj<-RinsGmuA[j]

      mu[j] =(beta0*GinvsG2[j] - beta0^2/se2[j]*RinSmujj2 - beta0/se2[j]*RinSmuAjj + xi*ginvsg2[j] - xi^2/ se1[j]*RinSmujj)*v2[j]
      Rinsgmu = tmp1 + Rins[,j]*mu[j]
      RinsGmu = tmp2 + Rins2[,j]*mu[j]
    }


    v2A<-1/(diaginsGRinsG + (alpha_w/beta_w)*rep(1/sgal2,p))
    for(k in 1:p){
      tmp3<-RinsGmuA - Rins2[,k]*muA[k]
      RinSmuAkk<-RinsGmuA[k] - R[k, k]*muA[k]/se2[k]
      muA[k]<-(GinvsG2[k] - beta0*(1/se2[k])*RinsGmu[k] - 1/se2[k]*RinSmuAkk)*v2A[k]
      RinsGmuA<-tmp3+Rins2[,k]*muA[k]
    }


    zong<-0
    for(g in 1:p){
      zong<-zong+muA[g]^2+v2A[g]
    }
    alpha_w<-alphag+p/2
    beta_w<-betag+zong/(2*sgal2)

    #M step
    #update beta0
    if(constr==1){
      beta0<-0
    }
    else{
      sig2b<-1/as.vector(t(mu)%*%insGRinsG%*%mu+t(v2)%*%diaginsGRinsG)
      beta0<-as.vector(as.vector(t(GinvsG2)%*%mu-t(muA)%*%insGRinsG%*%mu)*sig2b)
    }
    #update sgga2
    sgga2<-as.vector(t(mu)%*%mu+sum(v2))/p
    #update sgal2
    sgal2<-as.vector(t(muA)%*%(muA)+sum(v2A))/p

    #update xi
    xi<-as.vector((t(ginvsg2)%*%mu)/(t(mu)%*%insgRinsg%*%mu+t(v2)%*%diaginsgRinsg))

    #xi<-1
    #Reduction step
    xi2<-xi^2
    mu<-as.vector(xi*mu)
    beta0<-beta0/xi
    beta02<-beta0^2
    sgga2<-sgga2*xi2
    v2<-as.vector(xi2*v2)
    RinsGmu<-as.vector(xi*RinsGmu)
    Rinsgmu<-as.vector(xi*Rinsgmu)
    xi<-1
    xi2<-1

    #lower bound to check convergence
    term1<-as.vector(t(beta0*GinvsG2 + xi*ginvsg2)%*%mu) -
      0.5*as.vector(t(mu)%*%(beta02*insGRinsG + xi2*insgRinsg + 1/sgga2*diag(rep(1,p)))%*%mu) -
      0.5*as.vector(t(v2)%*%(beta02*diaginsGRinsG+ xi2*diaginsgRinsg + 1/sgga2))- 0.5*p*log(sgga2) + 0.5*sum(log(v2))

    term2<-as.vector(-beta0*t(muA)%*%insGRinsG%*%mu + t(GinvsG2)%*%muA - 0.5*t(muA)%*%(insGRinsG + 1/sgal2 *(alpha_w/beta_w)* diag(rep(1,p)))%*%muA)-
      0.5*as.vector(t(v2A)%*%(diaginsGRinsG +1/sgal2)) - 0.5*p*log(sgal2) + 0.5*sum(log(v2A))

    term3<-0.5*p*(digamma(alpha_w)-log(beta_w))

    low <- term1 + term2+ term3
    loglik[iter]<-low
    Iteration<-iter


    if(iter>1){
      if(abs(loglik[iter]-loglik[iter-1])<epsStopLogLik){
        break;
      }
    }
  }

  loglik_out<-vector()
  loglik_out<-loglik[1:Iteration]
  diff<-loglik_out[Iteration]-loglik_out[Iteration-1]

  #for loglikelihood ratio test
  MG <- beta02*insGRinsG + xi2*insgRinsg + 1./sgga2*diag(rep(1, p))
  MA <- insGRinsG + 1/sgal2*(alpha_w/beta_w)*diag(rep(1, p))
  SigG <- solve(MG)
  SigA <- solve(MA)
  t1 <- as.vector(t(beta0* GinvsG2 + xi*ginvsg2)%*%mu - 0.5*(t(mu)%*%MG%*%mu)) -
    0.5*p*(1 + log(sgga2))  + sum(log(diag(chol(SigG))))

  t2 <- as.vector(-beta0*t(mu)%*%insGRinsG%*%muA  + t(GinvsG2)%*%muA -
                    0.5 * t(muA)%*%MA%*%muA)- 0.5*p*(1 + log(sgal2))  + sum(log(diag(chol(SigA)))) + 0.5*p;

  t3<-0.5*p*(digamma(alpha_w)-log(beta_w))
  tstat = t1 + t2+ t3;

  output<-list(Iteration = Iteration,loglik=loglik_out,diff=diff,tstat=tstat,beta0=beta0,sgal2=sgal2,sgga2=sgga2,alpha_w=alpha_w,beta_w=beta_w)
  return (output)
}
