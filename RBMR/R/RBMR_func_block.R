#' RBMR function for multiple LD blocks
#'
#' When the number of LD block is greater than 1, this function gets the result of causal effect.
#'
#' @param F4Rblock A list contains the LD matrix of each block.
#' @param block_inf A matrix contains the information of the SNPS of the LD block matrix, the number of rows is the number of blocks, the number of columns is 2. Each row represents the start and end of the SNPs. That is, the first column is the starting position and the second column is the ending position.
#' @param nblocks Number of blocks of LD matrix.
#' @param bh1 A vector of SNP effects on the exposure variable, usually obtained from a GWAS.
#' @param bh2 A vector of SNP effects on the outcome varible, usually obtained from a GWAS.
#' @param se1 A vector of standard errors of \code{bh1}.
#' @param se2 A vector of standard errors of \code{bh2}.
#' @param gamma A vector to intialize the true value of the mean of the exposure.
#' @param alpha A vector to intialize the true value of the mean of the direct effect.
#' @param sgga2 A constant to initialize the true value of the standard error of exposure \code{bh1}.
#' @param sgal2 A constant to initialize the true value of the standard error of direct effect \code{bh2}.
#' @param beta0 A constant, initialize the true value of causal effect.
#' @param constr 0 or 1, when constr is equal to 0, the function calculates the ELBO under the alternative hypothesis, when constr is equal to 1, the function calculates ELBO under the null hypothesis.
#' @param epsStopLogLik Numerical precision
#' @param alphag A constant to initialize the parameter of multivariate generalize t distribution
#' @param IterMax Maximum number of interations to solve the estimating equations.
#' @param betag A constant to initialize the parameter of multivariate generalize t distribution
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
#' @examples
RBMR_func_block<-function(F4Rblock, block_inf, nblocks, bh1, bh2, se1, se2,
                          gamma, alpha, sgga2, sgal2, beta0, constr, epsStopLogLik, IterMax,alphag,betag){

  #define
  F4mu<-list();F4muA<-list();F4se1<-list();F4se2<-list()
  F4sg2<-list();F4sG2<-list();F4GinvsG2<-list();F4ginvsg2<-list()
  F4insGRinsG<-list();F4insgRinsg<-list();F4diaginsGRinsG<-list();F4diaginsgRinsg<-list()
  F4Rins<-list();F4Rins2<-list();F4RinsGmu<-list();F4Rinsgmu<-list();F4RinsGmuA<-list()
  NB<-vector();F4v2<-list();F4v2A<-list();betade<-vector();betanu<-vector();
  alpha_w_b<-vector();beta_w_b<-vector();low_b<-vector()

  mu<-gamma
  muA<-alpha

  p<-length(bh1)
  for(nn in 1:nblocks)
  {
    se1_block<-se1[(block_inf[nn,1]+1):(block_inf[nn,2]+1)]
    se2_block<-se2[(block_inf[nn,1]+1):(block_inf[nn,2]+1)]
    sg2_block<-se1_block^2
    sG2_block<-se2_block^2
    mu_block<-mu[(block_inf[nn,1]+1):(block_inf[nn,2]+1)]
    muA_block<-muA[(block_inf[nn,1]+1):(block_inf[nn,2]+1)]
    bh1_block<-bh1[(block_inf[nn,1]+1):(block_inf[nn,2]+1)]
    bh2_block<-bh2[(block_inf[nn,1]+1):(block_inf[nn,2]+1)]

    R_block<-F4Rblock[nn,1]
    NB[nn] = length(se1_block)
    F4mu[[nn]] = mu_block
    F4muA[[nn]] = muA_block
    F4se1[[nn]] = se1_block
    F4se2[[nn]] = se2_block

    F4sg2[[nn]] = sg2_block;
    F4sG2[[nn]] = sG2_block;

    F4GinvsG2[[nn]] = bh2_block / sG2_block;
    F4ginvsg2[[nn]] = bh1_block / sg2_block;

    F4insGRinsG[[nn]] = diag(1. / se2_block,NB[nn],NB[nn])%*%R_block[[1]]%*%diag(1. / se2_block,NB[nn],NB[nn]);
    F4insgRinsg[[nn]] = diag(1. / se1_block,NB[nn],NB[nn])%*%R_block[[1]]%*%diag(1. / se1_block,NB[nn],NB[nn]);

    F4diaginsGRinsG[[nn]] = diag(F4insGRinsG[[nn]]);
    F4diaginsgRinsg[[nn]] = diag(F4insgRinsg[[nn]]);

    F4Rins[[nn]] = R_block[[1]]%*%diag(1 / se1_block,NB[nn],NB[nn]);
    F4Rins2[[nn]] = R_block[[1]]%*%diag(1 / se2_block,NB[nn],NB[nn]);

    F4RinsGmu[[nn]] = R_block[[1]]%*%diag(1 / se2_block,NB[nn],NB[nn])%*%mu_block;
    F4Rinsgmu[[nn]] = R_block[[1]]%*%diag(1 / se1_block,NB[nn],NB[nn])%*%mu_block;
    F4RinsGmuA[[nn]] = R_block[[1]]%*%diag(1. / se2_block,NB[nn],NB[nn])%*%muA_block;
  }

  xi<-1
  loglik<-numeric(0)
  Iteration<-1
  alpha_w<-alphag
  beta_w<-betag

  for(iter in 1:IterMax){
    for(nb1 in 1:nblocks){


      v2 <- 1/((beta0^2)*F4diaginsGRinsG[[nb1]] + (xi^2)*F4diaginsgRinsg[[nb1]] + 1/ sgga2)
      F4v2[[nb1]]<-v2
      p_block<-NB[nb1]

      for(j in 1:(p_block)){
        tmp1<-F4Rinsgmu[[nb1]] - F4Rins[[nb1]][,j]*F4mu[[nb1]][j];
        tmp2 <- F4RinsGmu[[nb1]] - F4Rins2[[nb1]][,j]*F4mu[[nb1]][j];

        RinSmujj <- F4Rinsgmu[[nb1]][j] - F4Rblock[nb1,1][[1]][j, j]*F4mu[[nb1]][j] / F4se1[[nb1]][j]
        RinSmujj2 <- F4RinsGmu[[nb1]][j] - F4Rblock[nb1,1][[1]][j, j]*F4mu[[nb1]][j] / F4se2[[nb1]][j]
        RinSmuAjj = F4RinsGmuA[[nb1]][j]

        F4mu[[nb1]][j] <- (beta0*F4GinvsG2[[nb1]][j] - beta0^2/F4se2[[nb1]][j]*RinSmujj2  - beta0/F4se2[[nb1]][j]*RinSmuAjj+
                             xi*F4ginvsg2[[nb1]][j] - xi^2/ F4se1[[nb1]][j]*RinSmujj)*v2[j];
        F4Rinsgmu[[nb1]] <- tmp1 + F4Rins[[nb1]][,j]*F4mu[[nb1]][j];
        F4RinsGmu[[nb1]] <- tmp2 + F4Rins2[[nb1]][,j]*F4mu[[nb1]][j];

      }


      v2A <- 1/ (1/ (F4sG2[[nb1]])+(alpha_w/beta_w)*1/ sgal2)
      F4v2A[[nb1]]<-v2A
      pa_block<-NB[nb1]
      for(k in 1:pa_block){

        tmp3 <- F4RinsGmuA[[nb1]] - F4Rins2[[nb1]][,k]*F4muA[[nb1]][k];
        RinSmuAkk = F4RinsGmuA[[nb1]][k] - F4Rblock[nb1,1][[1]][k, k]*F4muA[[nb1]][k] / F4se2[[nb1]][k];
        F4muA[[nb1]][k] = (F4GinvsG2[[nb1]][k] - beta0*( 1/ F4se2[[nb1]][k])*F4RinsGmu[[nb1]][k] - 1 / F4se2[[nb1]][k]*RinSmuAkk)*v2A[k];
        F4RinsGmuA[[nb1]] = tmp3 + F4Rins2[[nb1]][,k]*F4muA[[nb1]][k];
      }

      #estimate the alpha_w and beta_w
      zong<-0
      pp_block<-NB[nb1]
      for(g in 1:pp_block){
        zong<-zong+F4muA[[nb1]][g]^2+v2A[g]
      }
      alpha_w_b[nb1]<-alphag+pp_block/2
      beta_w_b[nb1]<-betag+zong/(2*sgal2)
    }

    #M step
    #update beta0
    if (constr == 1){
      beta0 = 0;
    }
    else {
      for (jj in 1: nblocks)
      {
        betade[jj] = as.vector(t(F4mu[[jj]])%*%F4insGRinsG[[jj]]%*%F4mu[[jj]] + t(F4v2[[jj]])%*%F4diaginsGRinsG[[jj]]);
        betanu[jj] = as.vector(t(F4GinvsG2[[jj]])%*%F4mu[[jj]]-t(F4muA[[jj]])%*%F4insGRinsG[[jj]]%*%F4mu[[jj]]);
      }
      beta0 = sum(betanu) / sum(betade);
    }

    sgga2_b<-vector()
    sgal2_b<-vector()
    xi_nu<-vector();xi_de<-vector()

    for(kk in 1:nblocks){
      sgga2_b[kk] = sum(t(F4mu[[kk]])%*%F4mu[[kk]]) + sum(F4v2[[kk]]);
      sgal2_b[kk] = sum(t(F4muA[[kk]])%*%F4muA[[kk]]) + sum(F4v2A[[kk]]);
      xi_nu[kk] = as.vector(t(F4ginvsg2[[kk]])%*%F4mu[[kk]]);
      xi_de[kk] = as.vector(t(F4mu[[kk]])%*%F4insgRinsg[[kk]]%*%F4mu[[kk]] + t(F4v2[[kk]])%*%F4diaginsgRinsg[[kk]]);
    }

    #update sgga2
    sgga2<-sum(sgga2_b)/p
    #update sgal2
    sgal2<-sum(sgal2_b)/p
    #update alpha_w
    alpha_w<-sum(alpha_w_b)/p
    #update beta_w
    beta_w<-sum(beta_w_b)/p
    #update xi
    xi<-sum(xi_nu)/sum(xi_de)
    #xi<-1

    #Reduction step
    xi2<-xi*xi
    beta0<-beta0/xi
    beta02<-beta0*beta0
    sgga2<-sgga2*xi2

    for(ll in 1:nblocks){
      F4v2[[ll]] = F4v2[[ll]]*xi2;
      F4mu[[ll]] = F4mu[[ll]]*xi;
      F4RinsGmu[[ll]] = F4RinsGmu[[ll]]*xi;
      F4Rinsgmu[[ll]] = F4Rinsgmu[[ll]]*xi;
    }

    xi<-1
    xi2<-1

    #lower bound to check convergence
    low_b<-vector()
    for(ll in 1:nblocks)
    {
      p_b <- NB[ll];
      term1 = as.vector(t(beta0*F4GinvsG2[[ll]] + xi*F4ginvsg2[[ll]])%*%F4mu[[ll]]) -
        0.5*as.vector(t(F4mu[[ll]])%*%(beta02*F4insGRinsG[[ll]] + xi2*F4insgRinsg[[ll]] + 1./sgga2*diag(rep(1,p_b),nrow=p_b,ncol=p_b))%*%F4mu[[ll]]) -
        0.5*as.vector(t(F4v2[[ll]])%*%(beta02*F4diaginsGRinsG[[ll]]  + xi2*F4diaginsgRinsg[[ll]] + 1/sgga2)) + 0.5*sum(log(F4v2[[ll]]))-0.5*p_b*log(sgga2);

      term2 = as.vector(- beta0*t(F4muA[[ll]])%*%F4insGRinsG[[ll]]%*%F4mu[[ll]] + t(F4GinvsG2[[ll]])%*%F4muA[[ll]] -
                          0.5*t(F4muA[[ll]])%*%(F4insGRinsG[[ll]] + 1/sgal2 * (alpha_w/beta_w)*diag(rep(1,p_b),nrow=p_b,ncol=p_b))%*%F4muA[[ll]]) -
        0.5* as.vector(t(F4v2A[[ll]])%*%(F4diaginsGRinsG[[ll]] + 1/sgal2)) + 0.5*sum(log(F4v2A[[ll]]))-0.5*p_b*log(sgal2);
      term3<-0.5*p_b*(digamma(alpha_w)-log(beta_w))

      low_b[ll] = term1 + term2 + term3;
    }

    low<-sum(low_b)
    loglik<-c(loglik,low)
    Iteration<-iter

    if(iter>1){
      if(abs(loglik[iter]-loglik[iter-1])<epsStopLogLik){
        break
      }
    }

  }
  loglik_out<-vector()
  loglik_out<-loglik[1:Iteration]
  diff<-loglik_out[Iteration]-loglik_out[Iteration-1]

  tstat_b<-vector()

  for(bb in 1:nblocks){
    pb <- NB[bb];
    MG <- beta02*F4insGRinsG[[bb]] + xi2*F4insgRinsg[[bb]] + 1./sgga2*diag(rep(1,pb),nrow=pb,ncol=pb);
    MA <- F4insGRinsG[[bb]] + 1./sgal2*diag(rep(1,pb),nrow=pb,ncol=pb);
    SigG <- solve(MG);
    SigA <- solve(MA);

    t1 = as.vector(t(beta0* F4GinvsG2[[bb]] + xi*F4ginvsg2[[bb]])%*%F4mu[[bb]]) -
      0.5*(t(F4mu[[bb]])%*%MG%*%F4mu[[bb]])  + sum(log(diag(chol(SigG))))-0.5*pb*(1+log(sgga2));

    t2 = as.vector(-beta0*t(F4mu[[bb]])%*%F4insGRinsG[[bb]]%*%F4muA[[bb]]  + t(F4GinvsG2[[bb]])%*%F4muA[[bb]] -
                     0.5 * t(F4muA[[bb]])%*%MA%*%F4muA[[bb]])  + sum(log(diag(chol(SigA))))-0.5*pb*(1+log(sgal2));

    t3 = 0.5*pb*(digamma(alpha_w)-log(beta_w))

    tstat_b[bb] = t1 + t2+t3;
  }
  tstat<-sum(tstat_b)
  output<-list(Iteration = Iteration,loglik=loglik_out,diff=diff,tstat=tstat,beta0=beta0,sgal2=sgal2,sgga2=sgga2,alpha_w=alpha_w,beta_w=beta_w)
  return (output)
}
