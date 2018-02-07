
library(MASS)
library(MCMCpack)
library(mvtnorm)
library(cluster)

rm(list=ls())
list.files()

###########################################################################################################
##### sally.multouts() is a function which takes the (1) the data vector, Y, stacked by outcomes, 
#####         (2) covariate matrices for each of the ones specified in Luo et al. 2014: Sf, SDf, SD, SO  
#####         (3/4) the outcome labels y.jlabs, and the observation labels, y.ilabs, for the Y vector 
#####         (5/6)  priors as a list, priors, and the initial values as a list, inits
#####         (7) the number of iterations to perform the mcmc sampler (8) the FIXED domain assignments
###########################################################################################################

sally.multouts <- function(Y, Sf, SDf, SD, SO,  y.jlabs, y.ilabs, priors, inits, niter, printIter, dom.assns, ycorrs ){
  #### y.jlabs = distinguish outcomes
  #### y.ilabs = distinguish subjects
  #### printIter = TRUE to show the iteration for the overall MCMC
  #### printSumms = TRUE to show the clustering process for each outcome at each iteration
  #### sim.seed = random seed used in simulating everything but the covariate matrices in the function multOutData() 
  
  ##########################################
  ########## Define the following:
  nobs.vec <- table(y.jlabs)  # number of observations for each outcome (same if all have complete obs on all outcomes)
  pF <- ncol(Sf) 
  pDF<- ncol(SDf) 
  pD <- ncol(SD) 
  pO <- ncol(SO)
  J <- nlevels(factor(y.jlabs)) # number of outcomes
  n <- nlevels(factor(y.ilabs))# number of subjects/outcome
  y.mislabs <- apply(matrix(Y,ncol=1), 1, is.na) 
  nmis <- sum(y.mislabs)
  
  
  ##########################################
  ########## Matrices to Store Parameter Draws: 
  ## Domain assignments:
  d.draws <- dom.assns
  
  ## Non-domain specific parameters:
  K <- nlevels(factor(d.draws))
  betaf.draws <- matrix(NA, nrow=niter, ncol=pF)
  bO.draws <- matrix(NA, nrow=niter, ncol=pO*J)
  sig2bO.draws <- matrix(NA, nrow=niter, ncol=pO)
  r.draws <- matrix(NA, nrow=niter, ncol=n)
  sig2r.draws <- rep(NA, niter)
  sig2eps.draws <- matrix(NA, nrow=niter, ncol=J)
  sig2rD.draws <- matrix(NA, nrow=niter, ncol=K)
  sig2bD.draws <- matrix(NA, nrow=niter, ncol=pD)
  
  ## Domain-specific parameters: 
  betaDf.draws <- list()
  bD.draws <- list()
  rD.draws <- list()
  
  ##########################################
  ########## Matrix to Store Posterior Predictive Draws:
  if(!is.null(ycorrs)) Ypred <- matrix(NA, nrow=n*J, ncol=niter/ycorrs) 

    
  ##########################################
  ########## Initial Values (entered from inits=list()) :
  betaf.draws[1,] <- inits$betaF
  bO.draws[1,] <- inits$bO
  colnames(bO.draws) <- rep(1:J,each=pO)
  sig2bO.draws[1,] <- inits$sig2bO
  r.draws[1,] <-  inits$r
  sig2r.draws[1] <- inits$sig2r
  sig2eps.draws[1,] <- inits$sig2eps
  sig2rD.draws[1,] <- inits$sig2rD
  sig2bD.draws[1,] <-inits$sig2bD
  
  betaDf.draws[[1]] <- inits$betaDf
  bD.draws[[1]] <- inits$bD
  rD.draws[[1]] <- inits$rD
  
  
  ##########################################
  ########## Prior Values (entered from inits=list()) :
  pr.beta0f <- priors$beta0f
  pr.Sigma0Invf <- priors$SigInv.betaf
  pr.beta0Df <- priors$beta0Df
  pr.Sigma0InvDf <- priors$SigInv.betaDf 
  pr.A0.bO <- priors$sig2.bO[1]
  pr.B0.bO <- priors$sig2.bO[2]
  pr.A0.bD <- priors$sig2.bD[1]
  pr.B0.bD <- priors$sig2.bD[2] 
  pr.A0.eps <- priors$sig2.eps[1]
  pr.B0.eps <- priors$sig2.eps[2]
  pr.A0.r <- priors$sig2.r[1]
  pr.B0.r <- priors$sig2.r[2]
  pr.A0.rD <- priors$sig2.rD[1]
  pr.B0.rD <- priors$sig2.rD[2] 
  pr.alpha <- priors$alpha
  
  ##########################################
  ##########  Membership Matrix - j rows represent the outcomes, d columns represent 
  #  the membership of the jth outcome in the dth domain
  #   ---> this matrix will NO LONGER need to be updated at each iteration
  Lambda <- cbind(1*(d.draws==1),model.matrix( ~ as.factor(d.draws)  ,contrasts='contr.sum')[,-1])
  RDvec <- kronecker(Lambda,diag(1,n))%*%unlist(rD.draws[[1]])
  Rvec <- t(t(kronecker(c(rep(1,J)),r.draws[1,])))

  ##########################################
  ##########  Covariate matrices NO LONGER change with changing domains:   
  Fmat.F <- kronecker(rep(1, J),Sf)  
  Zmat.O <- kronecker(diag(1,J), SO) 
  Fmat.DF <- kronecker(Lambda,SDf) 
  Zmat.D <- kronecker(Lambda, SD)
  
  ## label outcomes according to domain assignment:
  d.labs <- rep(d.draws,nobs.vec)  
  ## number of outcomes in each domain:
  njD <- table(d.draws)  

  for(t in 2:niter){
    
    if(printIter==TRUE) cat('\n t=',t,' K=',K,';\n')
    ess.test <- 0
     
   
    
    ########################
    ### ---> Draw missing values
    means.Y <- Fmat.F%*%(betaf.draws[t-1,]) + Fmat.DF%*%unlist(betaDf.draws[[t-1]]) + Zmat.D%*%unlist(bD.draws[[t-1]]) + Zmat.O%*%(bO.draws[t-1,]) + RDvec + Rvec
    sds.Y <- rep(sqrt(sig2eps.draws[t-1,]),each=n)
    Ymis <- rnorm(nmis, mean=means.Y[y.mislabs], sd=sds.Y[y.mislabs])
    Y[y.mislabs] <- Ymis    
    
    ################################################
    #### ---> Gibbs Steps for Domain-Specific Parameters:
    ## EDIT: added in D_epsilon,d for means.Df..
    means.Df <- lapply(seq_along(1:K), FUN=function(k) pr.Sigma0InvDf%*%pr.beta0Df + t(kronecker(rep(1,njD[k]),SDf))%*%diag(rep(1/sig2eps.draws[(t-1),unique(y.jlabs[d.labs==k])],nobs.vec[unique(y.jlabs[d.labs==k])]))%*%(as.matrix(Y[d.labs==k],ncol=1) 
                                                        - kronecker(rep(1,njD[k]),(Sf%*%betaf.draws[(t-1),])) - kronecker(rep(1,njD[k]),(SD%*%(bD.draws[[t-1]][[k]]))) - kronecker(diag(1,njD[k]),SO)%*%bO.draws[(t-1),(colnames(bO.draws) %in% y.jlabs[d.labs==k])] 
                                                        - as.matrix(kronecker(rep(1,njD[k]),r.draws[(t-1),]) - kronecker(rep(1,njD[k]),(rD.draws[[t-1]][[k]])),ncol=1)) )
    var.betaDf <- lapply(seq_along(1:K), FUN=function(k) solve(pr.Sigma0InvDf + sum(1/1/sig2eps.draws[(t-1),unique(y.jlabs[d.labs==k])])*t(SDf)%*%SDf) )#t(kronecker(rep(1,njD[k]),SDf))%*%diag(rep(1/sig2eps.draws[(t-1),unique(y.jlabs[d.labs==k])],nobs.vec[unique(y.jlabs[d.labs==k])]))%*%kronecker(rep(1,njD[k]),SDf)) )
    betaDf.draws[[t]] <- lapply(X=seq_along(1:K), FUN = function(k) mvrnorm(1, mu=unlist(var.betaDf[[k]])%*%unlist(means.Df[[k]]) , Sigma =unlist(var.betaDf[[k]])))  

    
    means.bD <- lapply(seq_along(1:K), FUN=function(k) t(kronecker(rep(1,njD[k]),SD))%*%diag(rep(1/sig2eps.draws[(t-1),unique(y.jlabs[d.labs==k])],nobs.vec[unique(y.jlabs[d.labs==k])]))%*%
                         (as.matrix(Y[d.labs==k],ncol=1) - kronecker(rep(1,njD[k]),(Sf%*%betaf.draws[(t-1),])) - kronecker(rep(1,njD[k]),(SDf%*%betaDf.draws[[t]][[k]])) 
                          - kronecker(diag(1,njD[k]),SO)%*%bO.draws[(t-1),(colnames(bO.draws) %in% y.jlabs[d.labs==k])] - as.matrix(kronecker(rep(1,njD[k]),r.draws[(t-1),]) 
                           - kronecker(rep(1,njD[k]),(rD.draws[[t-1]][[k]])),ncol=1)) )
    var.bD <- lapply(seq_along(1:K), FUN=function(k) solve(diag(1/sig2bD.draws[(t-1),],pD) + t(kronecker(rep(1,njD[k]),SD))%*%diag(rep(1/sig2eps.draws[(t-1),unique(y.jlabs[d.labs==k])],nobs.vec[unique(y.jlabs[d.labs==k])]))%*%kronecker(rep(1,njD[k]),SD)) )
    bD.draws[[t]] <- lapply(X=seq_along(1:K), FUN = function(k) mvrnorm(1, mu=unlist(var.bD[[k]])%*%unlist(means.bD[[k]]) , Sigma =unlist(var.bD[[k]])))  
    
    var.rD <- lapply(seq_along(1:K), FUN=function(k) solve((1/sig2rD.draws[t-1,k]) + sum(1/sig2eps.draws[t-1,d.draws==k])) )
    separ.rD <- by((Y - Fmat.F%*%(betaf.draws[t-1,]) - Fmat.DF%*%unlist(betaDf.draws[[t]]) - Zmat.D%*%unlist(bD.draws[[t]])  - Zmat.O%*%(bO.draws[t-1,]) - Rvec),
                   INDICES=factor(y.jlabs),FUN=function(x) 1*x)
    means.rDbyout <- lapply(seq(1:J), FUN=function(j) (1/sig2eps.draws[t-1,j])*unlist(separ.rD[[j]]) )
    means.rDmat <- matrix(unlist(means.rDbyout),ncol=J,byrow=F)
    means.rD <- lapply(seq_along(1:K), FUN=function(k) 
      if( sum(d.draws==k)==1){
        means.rDmat[,k] 
      }else{
        rowSums(means.rDmat[,d.draws==k])
      }        )
    rD.draws[[t]] <- lapply(seq_along(1:K), FUN=function(k) rnorm(n, mean=unlist(var.rD[[k]])%*%unlist(means.rD[[k]]) , sd=sqrt(unlist(var.rD[[k]])))) 
    RDvec <- kronecker(Lambda,diag(1,n))%*%unlist(rD.draws[[t]])

    
    ###################
    #### ---> Gibbs Steps for non-domain specific parameters:
    ## Changed to sum over the n-subjects
    VbetafInv <- solve(pr.Sigma0Invf + sum(1/sig2eps.draws[t-1,])*t(Sf)%*%Sf)
    means.betaf1 <- VbetafInv%*%(pr.Sigma0Invf%*%pr.beta0f +(rep(1/sig2eps.draws[t-1,],each=n)*rep(Sf,J))%*%(Y - Fmat.DF%*%unlist(betaDf.draws[[t]]) - Zmat.D%*%unlist(bD.draws[[t]]) - Zmat.O%*%(bO.draws[(t-1),]) - Rvec - RDvec))
    betaf.draws[t,] <- rnorm(1, mean=means.betaf1,  sd=sqrt(VbetafInv)) 
    
    VbOInv <- by(sig2eps.draws[t-1,], INDICES=factor(1:J), FUN=function(x) solve(diag(1/sig2bO.draws[t-1,],pO)+(1/x)*t(SO)%*%SO)) 
    means.bO <- lapply(seq_along(1:J), FUN=function(j) unlist(VbOInv[[j]])%*%(1/sig2eps.draws[t-1,j])*t(SO)%*%
                         (Y[y.jlabs==j]- Sf%*%(betaf.draws[t,]) - SDf%*%betaDf.draws[[t]][[d.draws[j]]] - SD%*%bD.draws[[t]][[d.draws[j]]] - r.draws[t-1,] - rD.draws[[t]][[d.draws[j]]]))                     
    bO.draws[t,] <- unlist(lapply(seq_along(1:J), FUN=function(j) mvrnorm(1, mu=unlist(means.bO[[j]]), Sigma=unlist(VbOInv[[j]])) ))

    
    ## EDIT: r.draws <- rnorm ( ... aggregate(mean.diff, by=list(outcomes=y.ilabs), FUN=function(x) sum(x/c(1:J)))[,'V1'] )
    VrInv <- solve(1/sig2r.draws[t-1] + sum(1/sig2eps.draws[t-1,]))
    mean.diffr <- (Y - Fmat.F%*%(betaf.draws[t,]) - Fmat.DF%*%unlist(betaDf.draws[[t]]) - Zmat.D%*%unlist(bD.draws[[t]])  - Zmat.O%*%(bO.draws[t,]) - RDvec)
    mean.diffrlist <- by(mean.diffr, INDICES=factor(y.jlabs), FUN=function(x) 1*x)
    means.rbyout <-  lapply(seq(1:J), FUN=function(j) (1/sig2eps.draws[t-1,j])*unlist(mean.diffrlist[[j]]))
    means.r <- apply(matrix(unlist(means.rbyout),nrow=J,byrow=T),2,function(x) VrInv*sum(x))
    r.draws[t,] <- rnorm(n, mean=means.r, sd=sqrt(VrInv))

    
    Rvec <- t(t(kronecker(c(rep(1,J)),r.draws[t,]))) 
    
    
    ### ---> Variance parameters:
    bOj.sums <- apply(matrix(bO.draws[t,], nrow=pO),1,function(x) sum(x^2))
    sig2bO.draws[t,] <- rinvgamma(pO, shape=pr.A0.bO+J/2, scale=pr.B0.bO+.5*bOj.sums)
    
    bDd.sums <- apply(matrix(unlist(bD.draws[[t]]), nrow=pD),1,function(x) sum(x^2))
    sig2bD.draws[t,] <- rinvgamma(pD, shape=pr.A0.bD + K/2, scale=pr.B0.bD + .5*bDd.sums)
    
    sig2r.draws[t] <- rinvgamma(1, shape=pr.A0.r + n/2, scale=pr.B0.r + .5*t(r.draws[t,])%*%r.draws[t,])
    
    sig2rDsums <- matrix(apply(matrix(unlist(rD.draws[[t]]), nrow=K, byrow=T),1,function(x) sum(x^2)),nrow=1)
    sig2rD.draws[t,] <- apply(sig2rDsums, 2, function(x) rinvgamma(1, shape=pr.A0.rD + n/2, scale=pr.B0.rD + .5*x))
    
    mean.diff2 <- (Y - Fmat.F%*%(betaf.draws[t,]) - Fmat.DF%*%unlist(betaDf.draws[[t]]) - Zmat.D%*%unlist(bD.draws[[t]])  - Zmat.O%*%(bO.draws[t,]) - RDvec - Rvec)
    sig2eps.draws[t,] <- rinvgamma(J, shape=pr.A0.eps + n/2, scale=pr.B0.eps + .5*aggregate(mean.diff2, by=list(outcomes=y.jlabs), FUN=function(x) sum(x^2))[,'V1'])
    
    ##### Draws of Y from the posterior predictive distribution for model checking:
    if(!is.null(ycorrs)){
      ## predictive distribution of y
      if((t%%ycorrs)==0){
 	  means.Y <- Fmat.F%*%(betaf.draws[t-1,]) + Fmat.DF%*%unlist(betaDf.draws[[t-1]]) + Zmat.D%*%unlist(bD.draws[[t-1]]) + Zmat.O%*%(bO.draws[t-1,]) + RDvec + Rvec
        sds.Y <- rep(sqrt(sig2eps.draws[t-1,]),each=n)
        Ymis <- rnorm(nmis, mean=means.Y[y.mislabs], sd=sds.Y[y.mislabs])
        Ypred[,t/ycorrs] <- rnorm(n*J, mean=means.Y, sd=sds.Y)
      }
    }
  }    
  if(is.null(ycorrs)) Ypred <- NULL
  return( list( betaf.draws=betaf.draws,bO.draws=bO.draws, sig2r.draws=sig2r.draws, sig2eps.draws=sig2eps.draws, sig2bO.draws=sig2bO.draws,
          sig2bD.draws=sig2bD.draws, sig2rD.draws=sig2rD.draws, betaDf.draws=betaDf.draws,  bD.draws=bD.draws, Ypred=Ypred ))
}


  