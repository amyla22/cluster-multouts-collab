
##############################################################################################################
#####     FUNCTION TO DRAW THE DOMAIN INDICATOR VARIABLE IN NONCONJUGATE MODEL USING SPLIT-MERGE PROCEDURE 
##############################################################################################################  


draw.JN.d <- function(Y, tune,  d.olddraws, y.jlabs, fixed.mean, alpha, K, mean.diffnoDs, SD, SDf,
                      sig2eps=sig2eps.draws[t,],sig2rD=sig2rD.draws[[t]],sig2bD=sig2bD.draws[t,], 
                      r.D=rD.draws[[t]], b.D=bD.draws[[t]], beta.Df = betaDf.draws[[t]],
                      pr.Sigma0InvDf, pr.beta0Df, pr.A0.rD, pr.B0.rD, printSumms ){
  
  ## Y = vector of all outcomes (Y1, ... YJ); Yj = (y1, ..., yn)
  ## tune = (tuning param) # of intermediate RGSS to modify the launch_split and launch_merge states
  ## d.olddraws = vector of domain assignments for all outcomes prior to updating
  ## y.jlabs = outcome labels for Y
  ## fixed.mean =  current value of Xbeta + bO + r (all non-domain specific parameters)
  ## alpha = concentration parameter for the DP 
  ## K = current number of groups
  ## mean.diffnoDs = current value of Y - fixed.mean
  ## SD, SDf = domain-specific covariate matrices
  ## sig2eps, sig2bD = current value of non-domain specific variances 
  ## sig2rD, r.D, b.D, beta.Df = current value of all domain specific
  ## pr.Sigma0InvDf, pr.beta0Df, pr.A0.rD, pr.B0.rD = hyperparameter values for domain-specific parameters
  ## printSumms = if true, print summaries pertaining to clustering algorithm
  ## accepted = 1-accepted a split, 2-rejected a split, 3-accepted a merge, 4-rejected a merge
  
  
  qniter <- tune[1]
  rniter <- tune[2]
  
  pDF <- ncol(SDf); pD <- ncol(SD)
  n <- nrow(SDf)
  J <- length(unique(y.jlabs))
  
  outs <- sample(1:J, size=2, replace=F)
  # outs=c(1, 13); d_1=2; d_13=3
  d.Lsplit <- d.Lmerge <-  d.olddraws[which(d.olddraws==d.olddraws[outs[1]]|d.olddraws==d.olddraws[outs[2]])]
  ## S = set out outcomes in same group as outcome 1 and outcome 2 sampled in 'outs'
  S <- which(d.olddraws==d.olddraws[outs[1]]|d.olddraws==d.olddraws[outs[2]])[!which(d.olddraws==d.olddraws[outs[1]]|d.olddraws==d.olddraws[outs[2]])%in%(outs)] 
  ## S.ij = S + outcome i and j from 'outs'
  S.ij <- which(d.olddraws==d.olddraws[outs[1]]|d.olddraws==d.olddraws[outs[2]])
  ## Kcur = total current number of domains
  Kcur <- length(unique(d.olddraws))
  
  
  ################ Step 3 in Paper ###################
  ##### CONSTRUCT THE LAUNCH STATES:
  if(d.olddraws[outs[1]]==d.olddraws[outs[2]]){ 
    ### If i and j are in the same domain, 
    #    - split them into two by changing outcome 1 and keeping 2 the same for launch_split
    #    - make no changes for the launch_merge
    d.Lsplit[S.ij==outs[1]] <-  di_new <-  Kcur + 1 # need for launch_split
    dj <- d.olddraws[outs[2]]
    cat('Split: outs=',outs,' dij=',c(di_new,dj))
    
  }else{
    
    ### If i and j are in separate domains, 
    #    - let them keep their original labels for launch_split
    #    - assign all outcomes to the j^th domain (outcome 2) for launch_merge
    di_new <- d.olddraws[outs[1]] # need for launch_split
    dj <- d.olddraws[outs[2]]
    d.Lmerge <- rep(dj, length(S.ij))
    cat('Merge: outs=',outs,' dij=',c(di_new,dj))
    
  }
  ## Allocate the remaining outcomes in S to either of the two components w.p. 1/2
  d.Lsplit[S.ij%in%S] <- sample(c(di_new,dj), size=sum(S.ij%in%S),replace=T)
  
  # set.seed(97)
  ## Draw new parameters from the PRIOR for each component of luanch_split and only 1 for the launch_merge:
  # (1) sig2rD.new ~ IG(a, b)Ind(0,10) where a and b are hyperparameters and Ind(0,10) is the truncation
  sig2rD.Lsplit <- list(rinvgamma(1, shape=pr.A0.rD, scale=pr.B0.rD), rinvgamma(1, shape=pr.A0.rD, scale=pr.B0.rD))
  sig2rD.Lmerge <- rinvgamma(1, shape=pr.A0.rD, scale=pr.B0.rD)
  
  # (2) rD.Lsplit|sig2rD.Lsplit ~ N(0, sqrt(sig2rD.Lsplit)) 
  rD.Lsplit <- list(rnorm(n, mean=0, sd=sqrt(sig2rD.Lsplit[[1]])), rnorm(n, mean=0, sd=sqrt(sig2rD.Lsplit[[1]])))
  rD.Lmerge <- rnorm(n, mean=0, sd=sqrt(sig2rD.Lmerge))
  
  # (3) bD.Lsplit ~ N(0, sqrt(sig2bD.curr)) 
  bD.Lsplit <- list(rnorm(pD, mean=0, sd=sqrt(sig2bD)), rnorm(pD, mean=0, sd=sqrt(sig2bD)))
  bD.Lmerge <- rnorm(pD, mean=0, sd=sqrt(sig2bD) )
  
  # (4) betaDF.Lsplit ~ N(pr.beta0Df, pr.Sigma0InvDf) -- independent in the prior
  betaDf.Lsplit <- list(matrix(rmvnorm(1, mean=pr.beta0Df, sigma = solve(pr.Sigma0InvDf)),ncol=1),
                       matrix(rmvnorm(1, mean=pr.beta0Df, sigma = solve(pr.Sigma0InvDf)),ncol=1))
  betaDf.Lmerge <- matrix(rmvnorm(1, mean=pr.beta0Df, sigma = solve(pr.Sigma0InvDf)),ncol=1)

   
  #### UPDATE THE LAUNCH_SPLIT STATES WITH q AND r RESTRICTED GIBBS SAMPLING SCANS (RGSS):
  meandiffs <- Y - fixed.mean
  for(i in 1:qniter){
    # cat('\ni=',i,'dsplit=',d.Lsplit)
    # This will force outcome i's domain to always be in the first component of the parameter list 
    #  (recall that in launch split if d_i=d_j originally then d_i = K+1)
    d.ij <- c(di_new, dj)
    
    sig2rD.Lsplit <- lapply(seq(1:2), FUN=function(k)  
      rinvgamma(1, shape=pr.A0.rD + n/2, scale=pr.B0.rD + .5*t(unlist(rD.Lsplit[[k]]))%*%unlist(rD.Lsplit[[k]])))
    
    vrDLsplit <- lapply(seq_along(1:2), FUN=function(k) solve((1/sig2rD.Lsplit[[k]]) + sum(1/sig2eps[S.ij[d.Lsplit==d.ij[k]]])) )
    meansrDLsplit <- lapply(seq_along(1:2), FUN=function(k) diag(rep(1/sig2eps[S.ij[d.Lsplit==d.ij[k]]],each=n))%*%
                              matrix(meandiffs[y.jlabs%in%S.ij[d.Lsplit==d.ij[k]]],ncol=1) - kronecker(rep(1,sum(d.Lsplit==d.ij[k])),(SDf%*%(betaDf.Lsplit[[k]])))
                            - kronecker(rep(1,sum(d.Lsplit==d.ij[k])),(SD%*%(bD.Lsplit[[k]]))))         
    rD.Lsplit <- lapply(seq_along(1:2), FUN=function(k) rnorm(n, mean=unlist(vrDLsplit[[k]])%*%t(unlist(meansrDLsplit[[k]]) ), sd=sqrt(unlist(vrDLsplit[[k]])))) 
    
    meansDfLsplit <- lapply(seq_along(1:2), FUN=function(k) pr.Sigma0InvDf%*%pr.beta0Df + t(kronecker(rep(1,sum(d.Lsplit==d.ij[k])),SDf))%*%
                              diag(rep(1/sig2eps[S.ij[d.Lsplit==d.ij[k]]],each=n))%*%(matrix(meandiffs[y.jlabs%in%S.ij[d.Lsplit==d.ij[k]]],ncol=1)
                                                                                      - kronecker(rep(1,sum(d.Lsplit==d.ij[k])),(SD%*%(bD.Lsplit[[k]]))) - matrix(kronecker(rep(1,sum(d.Lsplit==d.ij[k])),(rD.Lsplit[[k]])),ncol=1 ))  )                                                       
    vDfLsplit <- lapply(seq_along(1:2), FUN=function(k) solve(pr.Sigma0InvDf + sum(1/sig2eps[S.ij[d.Lsplit==d.ij[k]]])*t(SDf)%*%SDf)) ## optimize other sampler!!!
    betaDf.Lsplit <- lapply(X=seq_along(1:2), FUN = function(k) mvrnorm(1, mu=unlist(vDfLsplit[[k]])%*%unlist(meansDfLsplit[[k]]) , Sigma =unlist(vDfLsplit[[k]])))  
    
    meansbDLsplit <- lapply(seq_along(1:2), FUN=function(k) t(kronecker(rep(1,sum(d.Lsplit==d.ij[k])),SD))%*%diag(rep(1/sig2eps[S.ij[d.Lsplit==d.ij[k]]],each=n))%*%
                              (matrix(meandiffs[y.jlabs%in%S.ij[d.Lsplit==d.ij[k]]],ncol=1)
                               - kronecker(rep(1,sum(d.Lsplit==d.ij[k])),(SDf%*%(betaDf.Lsplit[[k]]))) - matrix(kronecker(rep(1,sum(d.Lsplit==d.ij[k])),(rD.Lsplit[[k]])),ncol=1 ))  )                                                       
    vbDLsplit<- lapply(seq_along(1:2), FUN=function(k) solve(diag(1/sig2bD,pD) + sum(1/sig2eps[S.ij[d.Lsplit==d.ij[k]]])*t(SD)%*%SD) )
    bD.Lsplit <- lapply(X=seq_along(1:2), FUN = function(k) rnorm(pD, mean=unlist(vbDLsplit[[k]])%*%unlist(meansbDLsplit[[k]]) , sd = sqrt(vbDLsplit[[k]])))  
    
    if(length(d.Lsplit)>2){
      lmeans <- lapply(seq_along(1:2),FUN=function(k) 
        matrix(fixed.mean[y.jlabs%in%S],ncol=1 ) + kronecker(rep(1,length(S)),SDf)%*%betaDf.Lsplit[[k]] + 
          kronecker(rep(1,length(S)),SD)%*%bD.Lsplit[[k]] + matrix(kronecker(rep(1,length(S)),rD.Lsplit[[k]]),ncol=1) )
      lsds <- rep(sqrt(sig2eps[S]),each=n)
      ldens <- lapply(seq_along(1:2), FUN=function(k) by(dnorm(Y[y.jlabs%in%S],mean=lmeans[[k]],sd=lsds,log=T),  
                                                         INDICES = factor(y.jlabs[y.jlabs%in%S]), FUN = sum))
      # Question: Can we use log-densities here?
      phi <- -max(ldens[[1]],ldens[[2]])
      dens <- lapply(seq_along(1:2), FUN = function(k) exp(ldens[[k]]+phi))
      # cat('\n1-dens:',c(dens[[1]]),'\n2-dens',c(dens[[2]]))
      if(sum(dens[[1]]==0 & dens[[2]]==0)>0){
        ## problem is 1 (1) 1 1 1 1 1 1 1 1 3 2 2 2 2 2 2 2 2 (3) we propose to merge these 
        #   but when looking at whether the other (3) fits in any of the proposed splits, we get an NaN
        ## we set the more negative l-likelihood to 0 so that it always groups into the domain
        #   with the less negative log likelihood
        fix <- which(dens[[1]]==0 & dens[[2]]==0)
        probi <- (dens[[1]]*sum(d.Lsplit==d.ij[1]))/(dens[[1]]*sum(d.Lsplit==d.ij[1]) +dens[[2]]*sum(d.Lsplit==d.ij[2]))
        fix1 <- ldens[[1]][fix]<ldens[[2]][fix]
        cat('fix1=',fix1)
        if(sum(fix1)>0){
          dens[[1]][fix] <- dens[[2]][fix] <- 1
          probi <- (dens[[1]]*sum(d.Lsplit==d.ij[1]))/(dens[[1]]*sum(d.Lsplit==d.ij[1]) +dens[[2]]*sum(d.Lsplit==d.ij[2]))
          probi[fix[fix1==T]] <- 0
          probi[fix[fix1!=T]] <- 1
        }else{
          probi[fix[fix1!=T]] <- 1
        }
        # cat('\nfix[fix1==T]',fix[fix1==T],'\nfix[fix1!=T]',fix[fix1!=T],'prob.i=',probi)
      }else{
        probi <- (dens[[1]]*sum(d.Lsplit==d.ij[1]))/(dens[[1]]*sum(d.Lsplit==d.ij[1]) +dens[[2]]*sum(d.Lsplit==d.ij[2]))
      }
      # cat('\nfix[fix1!=T]',fix[fix1!=T],'prob.i=',probi)
      d.Lsplit[S.ij%in%S] <- unlist(lapply(seq_along(S), FUN=function(x) sample(d.ij, size=1, prob=c(probi[[x]],1-probi[[x]]))))
      
    }
  }
  
  for(i in 1:rniter){
    sig2rD.Lmerge <- rinvgamma(1, shape=pr.A0.rD + n/2, scale=pr.B0.rD + .5*t(rD.Lmerge)%*%rD.Lmerge)
    # while(sig2rD.Lmerge>500) sig2rD.Lmerge <- rinvgamma(1, shape=pr.A0.rD + n/2, scale=pr.B0.rD + .5*t(rD.Lmerge)%*%rD.Lmerge)
    
    vrD.Lmerge <- solve((1/sig2rD.Lmerge) + sum(1/sig2eps[S.ij]))
    meanrD.Lmerge <-  lapply( 1:length(S.ij), function(k) (1/sig2eps[S.ij[k]])*
                                by(Y[y.jlabs%in%S.ij]- fixed.mean[y.jlabs%in%S.ij] - kronecker(rep(1,length(S.ij)),SDf)%*%betaDf.Lmerge - t(t(rep(SD,length(S.ij))))%*%bD.Lmerge ,
                                   INDICES = factor(y.jlabs[y.jlabs%in%S.ij]), FUN= function(x) 1*x)[[k]])
    rD.Lmerge <- rnorm(n, mean=vrD.Lmerge%*%rowSums(matrix(unlist(meanrD.Lmerge),nrow=n,byrow=F)), sd=sqrt(vrD.Lmerge))
    
    vbD.Lmerge <- solve(1/sig2bD + sum(apply(matrix((1/sig2eps[S.ij]),ncol=length(S.ij)),2,function(x) x*t(SD)%*%SD)))
    meanbD.Lmerge <- t(SD)%*%rowSums(matrix(unlist(lapply( 1:length(S.ij), function(k) (1/sig2eps[S.ij[k]])*
                                                             by(Y[y.jlabs%in%S.ij]- fixed.mean[y.jlabs%in%S.ij] - kronecker(rep(1,length(S.ij)),SDf)%*%betaDf.Lmerge - matrix(kronecker(rep(1,length(S.ij)),rD.Lmerge),ncol=1) , INDICES = factor(y.jlabs[y.jlabs%in%S.ij]), FUN= function(x) 1*x)[[k]])),nrow=n,byrow=F))    
    bD.Lmerge <- rnorm(pD, mean=vbD.Lmerge%*%meanbD.Lmerge, sd = sqrt(vbD.Lmerge))
    
    vDfLmerge <- solve(pr.Sigma0InvDf + Reduce('+',lapply(1:length(S.ij),FUN=function(k) (1/sig2eps[S.ij[k]])*t(SDf)%*%SDf)))
    meanbetaDf.Lmerge <- pr.Sigma0InvDf%*%pr.beta0Df + t(SDf)%*%
      Reduce('+',lapply( 1:length(S.ij), function(k) (1/sig2eps[S.ij[k]])*
                           by(Y[y.jlabs%in%S.ij]- fixed.mean[y.jlabs%in%S.ij] - kronecker(rep(1,length(S.ij)),SD)%*%bD.Lmerge - matrix(kronecker(rep(1,length(S.ij)),rD.Lmerge),ncol=1) , INDICES = factor(y.jlabs[y.jlabs%in%S.ij]), FUN= function(x) 1*x)[[k]]))  
    betaDf.Lmerge <- matrix(rmvnorm(1, mean=vDfLmerge%*%meanbetaDf.Lmerge, sigma=vDfLmerge),ncol=1)
  }
  
  
  
  if(d.olddraws[outs[1]]==d.olddraws[outs[2]]){
    ### PROPOSE A SPLIT:
    d.split <- d.Lsplit
    ## Assign all k in S to either K+1 or based on launch_split 
    ##  We can re-use the final probability calculated in the launch_split procedure:   lmeans.split <- lapply(seq_along(1:2),FUN=function(k) matrix(fixed.mean[y.jlabs%in%S],ncol=1 ) + kronecker(rep(1,length(S)),SDf)%*%betaDf.Lsplit[[k]] + kronecker(rep(1,length(S)),SD)%*%bD.Lsplit[[k]] + matrix(kronecker(rep(1,length(S)),rD.Lsplit[[k]]),ncol=1) )
    #    probi <- (dens[[1]]*sum(d.Lsplit==d.ij[1]))/(dens[[1]]*sum(d.Lsplit==d.ij[1]) +dens[[2]]*sum(d.Lsplit==d.ij[2]))
    
    if(length(d.Lsplit)>2){
      d.split[S.ij%in%S] <- unlist(lapply(seq_along(S), FUN=function(x) sample(d.ij, size=1, prob=c(probi[[x]],1-probi[[x]]))))
    }
    sig2rD.split <- lapply(seq(1:2), FUN=function(k)  
      rinvgamma(1, shape=pr.A0.rD + n/2, scale=pr.B0.rD + .5*t(unlist(rD.Lsplit[[k]]))%*%unlist(rD.Lsplit[[k]])))
    
    
    vrDsplit <- lapply(seq_along(1:2), FUN=function(k) solve((1/sig2rD.split[[k]]) + sum(1/sig2eps[S.ij[d.split==d.ij[k]]])) )
    meansrDsplit <- lapply(seq_along(1:2), FUN=function(k) diag(rep(1/sig2eps[S.ij[d.split==d.ij[k]]],each=n))%*%
                             matrix(meandiffs[y.jlabs%in%S.ij[d.split==d.ij[k]]],ncol=1) - kronecker(rep(1,sum(d.split==d.ij[k])),(SDf%*%(betaDf.Lsplit[[k]])))
                           - kronecker(rep(1,sum(d.split==d.ij[k])),(SD%*%(bD.Lsplit[[k]]))))         
    rD.split <- lapply(seq_along(1:2), FUN=function(k) rnorm(n, mean=unlist(vrDsplit[[k]])%*%t(unlist(meansrDsplit[[k]]) ), sd=sqrt(unlist(vrDsplit[[k]])))) 
    
    meansDfsplit <- lapply(seq_along(1:2), FUN=function(k) pr.Sigma0InvDf%*%pr.beta0Df + t(kronecker(rep(1,sum(d.split==d.ij[k])),SDf))%*%
                             diag(rep(1/sig2eps[S.ij[d.split==d.ij[k]]],each=n))%*%(matrix(meandiffs[y.jlabs%in%S.ij[d.split==d.ij[k]]],ncol=1)
                                                                                    - kronecker(rep(1,sum(d.split==d.ij[k])),(SD%*%(bD.Lsplit[[k]]))) 
                                                                                    - matrix(kronecker(rep(1,sum(d.split==d.ij[k])),(rD.split[[k]])),ncol=1 ))  )                                                       
    vDfsplit <- lapply(seq_along(1:2), FUN=function(k) solve(pr.Sigma0InvDf + sum(1/sig2eps[S.ij[d.split==d.ij[k]]])*t(SDf)%*%SDf)) ## optimize other sampler!!!
    betaDf.split <- lapply(X=seq_along(1:2), FUN = function(k) mvrnorm(1, mu=unlist(vDfsplit[[k]])%*%unlist(meansDfsplit[[k]]) , Sigma =unlist(vDfsplit[[k]])))  
    
    meansbDsplit <- lapply(seq_along(1:2), FUN=function(k) t(kronecker(rep(1,sum(d.split==d.ij[k])),SD))%*%diag(rep(1/sig2eps[S.ij[d.split==d.ij[k]]],each=n))%*%
                             (matrix(meandiffs[y.jlabs%in%S.ij[d.split==d.ij[k]]],ncol=1)
                              - kronecker(rep(1,sum(d.split==d.ij[k])),(SDf%*%(betaDf.split[[k]]))) - matrix(kronecker(rep(1,sum(d.split==d.ij[k])),(rD.split[[k]])),ncol=1 ))  )                                                       
    vbDsplit<- lapply(seq_along(1:2), FUN=function(k) solve(diag(1/sig2bD,pD) + sum(1/sig2eps[S.ij[d.split==d.ij[k]]])*t(SD)%*%SD) )
    bD.split <- lapply(X=seq_along(1:2), FUN = function(k) rnorm(pD, mean=unlist(vbDsplit[[k]])%*%unlist(meansbDsplit[[k]]) , sd = sqrt(vbDsplit[[k]])))  
    
    
    
    lqcurr.split <- (dmvnorm(beta.Df[[dj]], mean=betaDf.Lmerge, sigma=vDfLmerge, log=T) + dnorm(b.D[[dj]], mean=bD.Lmerge, sd=sqrt(vbD.Lmerge), log=T)
                     + sum(dnorm(r.D[[dj]], mean=rD.Lmerge, sd=sqrt(vrD.Lmerge),log=T)) + actuar::dinvgamma(sig2rD[[dj]], shape=pr.A0.rD + n/2, scale=pr.B0.rD + .5*t(rD.Lmerge)%*%rD.Lmerge,log=T ) )
    cat('betadf=',dmvnorm(beta.Df[[dj]], mean=betaDf.Lmerge, sigma=vDfLmerge, log=T),'\nbd=', dnorm(b.D[[dj]], mean=bD.Lmerge, sd=sqrt(vbD.Lmerge), log=T)
        ,'\nrd=' ,sum(dnorm(r.D[[dj]], mean=rD.Lmerge, sd=sqrt(vrD.Lmerge),log=T)) ,'\ns2rd', actuar::dinvgamma(sig2rD[[dj]], shape=pr.A0.rD + n/2, scale=pr.B0.rD + .5*t(rD.Lmerge)%*%rD.Lmerge,log=T ))
    
    lpr.curr <- (log(factorial(length(S.ij) -1) )
                 + dmvnorm(beta.Df[[dj]],  mean=pr.beta0Df, sigma=solve(pr.Sigma0InvDf), log=T) + dnorm(b.D[[dj]], mean=0, sd=sqrt(sig2bD), log=T)
                 + sum(dnorm(r.D[[dj]],mean=rep(0,n), sd=sqrt(sig2rD[[dj]]), log=T) ) + actuar::dinvgamma(sig2rD[[dj]], shape=pr.A0.rD , scale=pr.B0.rD, log=T ) ) 
    
    lpr.split <- (log(factorial(sum(d.split==d.ij[1]) -1) + factorial(sum(d.split==d.ij[2]) -1) ) 
                  + Reduce(sum, lapply(seq_along(1:2), function(k) dmvnorm(betaDf.split[[k]], mean=pr.beta0Df, sigma=solve(pr.Sigma0InvDf), log=T) 
                                       + dnorm(bD.split[[k]], mean=0, sd=sqrt(sig2bD), log=T)
                                       + sum(dnorm(rD.split[[k]], mean=rep(0,n), sd=sqrt(sig2rD.split[[k]]), log=T) )
                                       + actuar::dinvgamma(sig2rD.split[[k]], shape=pr.A0.rD, scale=pr.B0.rD , log=T) )) )
    
    if(length(d.Lsplit)>2){
      splmeans <- lapply(seq_along(1:2),FUN=function(k) matrix(fixed.mean[y.jlabs%in%S],ncol=1 ) + kronecker(rep(1,length(S)),SDf)%*%betaDf.split[[k]] + kronecker(rep(1,length(S)),SD)%*%bD.split[[k]] + matrix(kronecker(rep(1,length(S)),rD.split[[k]]),ncol=1) )
      splsds <- rep(sqrt(sig2eps[S]),each=n)
      spldens <- lapply(seq_along(1:2), FUN=function(k) by(dnorm(Y[y.jlabs%in%S],mean=splmeans[[k]],sd=splsds,log=T),  INDICES = factor(y.jlabs[y.jlabs%in%S]), FUN = sum))
      
      phi <- -max(spldens[[1]],spldens[[2]])
      dens <- lapply(seq_along(1:2), FUN = function(k) exp(spldens[[k]]+phi))
      # probi <- (dens[[1]]*sum(d.Lsplit==d.ij[1]))/(dens[[1]]*sum(d.Lsplit==d.ij[1]) + dens[[2]]*sum(d.Lsplit==d.ij[2]))
      if(sum(dens[[1]]==0 & dens[[2]]==0)>0){
        
        fix <- which(dens[[1]]==0 & dens[[2]]==0)
        probi <- (dens[[1]]*sum(d.Lsplit==d.ij[1]))/(dens[[1]]*sum(d.Lsplit==d.ij[1]) +dens[[2]]*sum(d.Lsplit==d.ij[2]))
        fix1 <- spldens[[1]][fix]<spldens[[2]][fix]
        if(sum(fix1)>0){
          dens[[1]][fix] <- dens[[2]][fix] <- 1
          probi <- (dens[[1]]*sum(d.Lsplit==d.ij[1]))/(dens[[1]]*sum(d.Lsplit==d.ij[1]) +dens[[2]]*sum(d.Lsplit==d.ij[2]))
          probi[fix[fix1==T]] <- 0
          probi[fix[fix1!=T]] <- 1        
        }else{
          probi[fix[fix1!=T]] <- 1
        }
        
      }else{
        probi <- (dens[[1]]*sum(d.Lsplit==d.ij[1]))/(dens[[1]]*sum(d.Lsplit==d.ij[1]) +dens[[2]]*sum(d.Lsplit==d.ij[2]))
      }
      probs <- list(probi, (1-probi))
      
      lqsplit.curr <- (Reduce(sum, lapply(seq_along(1:2), function(k) dmvnorm(betaDf.split[[k]], mean=betaDf.Lsplit[[k]], sigma=vDfLsplit[[k]], log=T) 
                                          + dnorm(bD.split[[k]], mean=bD.Lsplit[[k]], sd=sqrt(vbDLsplit[[k]]), log=T)
                                          + sum(dnorm(rD.split[[k]], mean=rD.Lsplit[[k]], sd=sqrt(vrDLsplit[[k]]), log=T) )
                                          + actuar::dinvgamma(sig2rD.split[[k]], shape=pr.A0.rD + n/2, scale=pr.B0.rD + .5*t(rD.Lsplit[[k]])%*%rD.Lsplit[[k]], log=T ) ))
                       + Reduce(sum, lapply(seq(1,length(S)), FUN = function(k) log(probs[[apply(matrix(d.split[S.ij%in%S],nrow=1),2 , function(x) which(d.ij==x))[k]]][[k]])
                       ) ) )
    } else {
      lqsplit.curr <- Reduce(sum, lapply(seq_along(1:2), function(k) dmvnorm(betaDf.split[[k]], mean=betaDf.Lsplit[[k]], sigma=vDfLsplit[[k]], log=T) 
                                         + dnorm(bD.split[[k]], mean=bD.Lsplit[[k]], sd=sqrt(vbDLsplit[[k]]), log=T)
                                         + sum(dnorm(rD.split[[k]], mean=rD.Lsplit[[k]], sd=sqrt(vrDLsplit[[k]]), log=T) )
                                         + actuar::dinvgamma(sig2rD.split[[k]], shape=pr.A0.rD + n/2, scale=pr.B0.rD + .5*t(rD.Lsplit[[k]])%*%rD.Lsplit[[k]], log=T ) ))
      
    }               
    
    llik.curr <- sum(dnorm(Y[y.jlabs%in%S.ij], mean=matrix(fixed.mean[y.jlabs%in%S.ij],ncol=1 ) + kronecker(rep(1,length(S.ij)),SDf)%*%beta.Df[[dj]] 
                           + kronecker(rep(1,length(S.ij)),SD)%*%b.D[[dj]] + matrix(kronecker(rep(1,length(S.ij)),r.D[[dj]]),ncol=1), 
                           sd = rep(sqrt(sig2eps[S.ij]),each=n) ,log=T))
    
    newlabsdij <- as.numeric(factor(d.split,levels=d.ij))
    
    llik.split <-  Reduce(sum, lapply(seq(1,2), FUN = function(k)  dnorm(Y[y.jlabs[y.jlabs%in%(S.ij[newlabsdij==k])]], mean=matrix(fixed.mean[y.jlabs[y.jlabs%in%(S.ij[newlabsdij==k])]],ncol=1 ) + 
                                                                           kronecker(rep(1,sum(newlabsdij==k)),SDf)%*%betaDf.split[[k]]
                                                                         + kronecker(rep(1,sum(newlabsdij==k)),SD)%*%bD.split[[k]] + matrix(kronecker(rep(1,sum(newlabsdij==k)),rD.split[[k]]),ncol=1), 
                                                                         sd = rep(sqrt(sig2eps[unique(y.jlabs[y.jlabs%in%(S.ij[newlabsdij==k])])]),each=n) ,log=T)))
    cat(' lqcur.spl=',lqcurr.split,' lqspl.cur=',lqsplit.curr,' lpr.spl=',lpr.split,' lpr.cur=',lpr.curr,' llik.spl=',llik.split, ' llik.cur=',llik.curr)
    
    acc.split <- min( 1,  exp((lqcurr.split - lqsplit.curr) + (lpr.split - lpr.curr) + (llik.split - llik.curr)) )
    
    u <- runif(1)
    if(u < acc.split){
      accepted <- 1
      d.olddraws[S.ij] <- d.split
      # Assign the 2nd component (dj) to the previous dj component and place new di at the end:
      beta.Df[[dj]] <- betaDf.split[[2]]; beta.Df[[length(beta.Df)+1]] <- betaDf.split[[1]]
      b.D[[dj]] <- bD.split[[2]]; b.D[[length(b.D)+1]] <- bD.split[[1]]
      r.D[[dj]] <- rD.split[[2]]; r.D[[length(r.D)+1]] <- rD.split[[1]]
      sig2rD[[dj]] <- sig2rD.split[[2]]; sig2rD <- c(sig2rD, sig2rD.split[[1]])
    } else{
      accepted <- 2
      d.olddraws <- d.olddraws
    }
    
  } else {
    
    ### PROPOSE A MERGE:
    ## Assign all k in S.ij to dj as in launch_merge 
    d.merge <- d.Lmerge
    
    ## Perform one final update for phi_merge
    sig2rD.merge <- rinvgamma(1, shape=pr.A0.rD + n/2, scale=pr.B0.rD + .5*t(rD.Lmerge)%*%rD.Lmerge)
    # while(sig2rD.merge>500) sig2rD.merge <- rinvgamma(1, shape=pr.A0.rD + n/2, scale=pr.B0.rD + .5*t(rD.Lmerge)%*%rD.Lmerge)
    vrDmerge <- solve((1/sig2rD.merge) + sum(1/sig2eps[S.ij]))
    meanrDmerge <-  lapply( 1:length(S.ij), function(k) (1/sig2eps[S.ij[k]])*
                              by(Y[y.jlabs%in%S.ij]- fixed.mean[y.jlabs%in%S.ij] - kronecker(rep(1,length(S.ij)),SDf)%*%betaDf.Lmerge - t(t(rep(SD,length(S.ij))))%*%bD.Lmerge , INDICES = factor(y.jlabs[y.jlabs%in%S.ij]), FUN= function(x) 1*x)[[k]])
    rD.merge <- rnorm(n, mean=vrDmerge%*%rowSums(matrix(unlist(meanrDmerge),nrow=n,byrow=F)), sd=sqrt(vrDmerge))
    
    vbDmerge <- solve(1/sig2bD + sum(apply(matrix((1/sig2eps[S.ij]),ncol=length(S.ij)),2,function(x) x*t(SD)%*%SD)))
    meanbDmerge <- t(SD)%*%rowSums(matrix(unlist(lapply( 1:length(S.ij), function(k) (1/sig2eps[S.ij[k]])*
                                                           by(Y[y.jlabs%in%S.ij]- fixed.mean[y.jlabs%in%S.ij] - kronecker(rep(1,length(S.ij)),SDf)%*%betaDf.Lmerge - matrix(kronecker(rep(1,length(S.ij)),rD.merge),ncol=1) , INDICES = factor(y.jlabs[y.jlabs%in%S.ij]), FUN= function(x) 1*x)[[k]])),nrow=n,byrow=F))    
    bD.merge <- rnorm(pD, mean=vbDmerge%*%meanbDmerge, sd = sqrt(vbDmerge))
    
    vbetaDfmerge <- solve(pr.Sigma0InvDf + Reduce('+',lapply(1:length(S.ij),FUN=function(k) (1/sig2eps[S.ij[k]])*t(SDf)%*%SDf)))
    meanbetaDfmerge <- pr.Sigma0InvDf%*%pr.beta0Df + t(SDf)%*%
      Reduce('+',lapply( 1:length(S.ij), function(k) (1/sig2eps[S.ij[k]])*
                           by(Y[y.jlabs%in%S.ij]- fixed.mean[y.jlabs%in%S.ij] - kronecker(rep(1,length(S.ij)),SD)%*%bD.merge - matrix(kronecker(rep(1,length(S.ij)),rD.merge),ncol=1) , INDICES = factor(y.jlabs[y.jlabs%in%S.ij]), FUN= function(x) 1*x)[[k]]))  
    betaDf.merge <- matrix(rmvnorm(1, mean=vbetaDfmerge%*%meanbetaDfmerge, sigma=vbetaDfmerge),ncol=1)
    
    if(length(d.Lsplit)>2){
      currmeans <- lapply(seq_along(1:2),FUN=function(k) matrix(fixed.mean[y.jlabs%in%S],ncol=1 ) + kronecker(rep(1,length(S)),SDf)%*%beta.Df[[d.ij[k]]] + kronecker(rep(1,length(S)),SD)%*%b.D[[d.ij[k]]] + matrix(kronecker(rep(1,length(S)),r.D[[d.ij[k]]]),ncol=1) )
      currsds <- rep(sqrt(sig2eps[S]),each=n)
      lcurrdens <- lapply(seq_along(1:2), FUN=function(k) by(dnorm(Y[y.jlabs%in%S],mean=currmeans[[k]],sd=currsds,log=T),  INDICES = factor(y.jlabs[y.jlabs%in%S]), FUN = sum))
      
      phi <- -max(lcurrdens[[1]],lcurrdens[[2]])
      dens <- lapply(seq_along(1:2), FUN = function(k) exp(lcurrdens[[k]]+phi))
      # probi <- (dens[[1]]*sum(d.Lsplit==d.ij[1]))/(dens[[1]]*sum(d.Lsplit==d.ij[1]) +dens[[2]]*sum(d.Lsplit==d.ij[2]))
      
      if(sum(dens[[1]]==0 & dens[[2]]==0)>0){
        
        fix <- which(dens[[1]]==0 & dens[[2]]==0)
        probi <- (dens[[1]]*sum(d.Lsplit==d.ij[1]))/(dens[[1]]*sum(d.Lsplit==d.ij[1]) +dens[[2]]*sum(d.Lsplit==d.ij[2]))
        fix1 <- lcurrdens[[1]][fix]<lcurrdens[[2]][fix]
        if(sum(fix1)>0){
          dens[[1]][fix] <- dens[[2]][fix] <- 1
          probi <- (dens[[1]]*sum(d.Lsplit==d.ij[1]))/(dens[[1]]*sum(d.Lsplit==d.ij[1]) +dens[[2]]*sum(d.Lsplit==d.ij[2]))
          probi[fix[fix1==T]] <- 0
          probi[fix[fix1!=T]] <- 1        
        }else{
          probi[fix[fix1!=T]] <- 1
        }
      }else{
        probi <- (dens[[1]]*sum(d.Lsplit==d.ij[1]))/(dens[[1]]*sum(d.Lsplit==d.ij[1]) +dens[[2]]*sum(d.Lsplit==d.ij[2]))
      }
      probs <- list(probi, (1-probi))     
      cat('probi=', probi)
      lqcurr.merge <- (Reduce(sum, lapply(seq_along(1:2), function(k) dmvnorm(beta.Df[[d.ij[k]]], mean=betaDf.Lsplit[[k]], sigma=vDfLsplit[[k]], log=T) 
                                          + dnorm(b.D[[d.ij[k]]], mean=bD.Lsplit[[k]], sd=sqrt(vbDLsplit[[k]]), log=T)
                                          + sum(dnorm(r.D[[d.ij[k]]], mean=rD.Lsplit[[k]], sd=sqrt(vrDLsplit[[k]]), log=T) )
                                          + actuar::dinvgamma(sig2rD[[d.ij[k]]], shape=pr.A0.rD + n/2, scale=pr.B0.rD + .5*t(rD.Lsplit[[k]])%*%rD.Lsplit[[k]],log=T ) ))
                       + Reduce(sum, lapply(seq(1,length(S)), FUN = function(k) log(probs[[apply(matrix(d.olddraws[S],nrow=1),2 , function(x) which(d.ij==x))[k]]][[k]])
                       ) ))
    } else {
      lqcurr.merge <- (Reduce(sum, lapply(seq_along(1:2), function(k) dmvnorm(beta.Df[[d.ij[k]]], mean=betaDf.Lsplit[[k]], sigma=vDfLsplit[[k]], log=T) 
                                          + dnorm(b.D[[d.ij[k]]], mean=bD.Lsplit[[k]], sd=sqrt(vbDLsplit[[k]]), log=T)
                                          + sum(dnorm(r.D[[d.ij[k]]], mean=rD.Lsplit[[k]], sd=sqrt(vrDLsplit[[k]]), log=T) )
                                          + actuar::dinvgamma(sig2rD[[d.ij[k]]], shape=pr.A0.rD + n/2, scale=pr.B0.rD + .5*t(rD.Lsplit[[k]])%*%rD.Lsplit[[k]],log=T ) )))
    }
    
    # cat('betadfi=',dmvnorm(beta.Df[[d.ij[1]]], mean=betaDf.Lsplit[[1]], sigma=vDfLsplit[[1]], log=T),'\nbetadfj=',dmvnorm(beta.Df[[d.ij[2]]], mean=betaDf.Lsplit[[1]], sigma=vDfLsplit[[1]], log=T),
    #     '\nbdi=', dnorm(b.D[[d.ij[1]]], mean=bD.Lsplit[[1]], sd=sqrt(vbDLsplit[[1]]), log=T),'\nbdj=', dnorm(b.D[[d.ij[2]]], mean=bD.Lsplit[[2]], sd=sqrt(vbDLsplit[[2]]), log=T),
    #     '\nrdi=',sum(dnorm(r.D[[d.ij[1]]], mean=rD.Lsplit[[2]], sd=sqrt(vrDLsplit[[1]]), log=T)) ,'\nrdj=' ,sum(dnorm(r.D[[d.ij[2]]], mean=rD.Lsplit[[2]], sd=sqrt(vrDLsplit[[2]]), log=T)) ,
    #     '\ns2rdi=', actuar::dinvgamma(sig2rD[[d.ij[1]]], shape=pr.A0.rD + n/2, scale=pr.B0.rD + .5*t(rD.Lsplit[[1]])%*%rD.Lsplit[[1]],log=T ),'\ns2rdi=', actuar::dinvgamma(sig2rD[[d.ij[2]]], shape=pr.A0.rD + n/2, scale=pr.B0.rD + .5*t(rD.Lsplit[[2]])%*%rD.Lsplit[[2]],log=T ))
    #     
    lqmerge.curr <- (dmvnorm(t(betaDf.merge), mean=t(betaDf.Lmerge), sigma=vDfLmerge, log=T) + dnorm(bD.merge, mean=bD.Lmerge, sd=sqrt(vbD.Lmerge), log=T)
                     + sum(dnorm(rD.merge, mean=rD.Lmerge, sd=sqrt(vrD.Lmerge))) + actuar::dinvgamma(sig2rD.merge, shape=pr.A0.rD + n/2, scale=pr.B0.rD + .5*t(rD.Lmerge)%*%rD.Lmerge,log=T ) )
    
    lpr.merg <- log(factorial(length(S.ij) -1) ) + (dmvnorm(t(betaDf.merge), mean=pr.beta0Df, sigma=solve(pr.Sigma0InvDf), log=T) + dnorm(bD.merge, mean=0, sd=sqrt(sig2bD), log=T)
                                                    + sum(dnorm(rD.merge, mean=rep(0,n), sd=sqrt(sig2rD.merge), log=T)) + actuar::dinvgamma(sig2rD.merge, shape=pr.A0.rD, scale=pr.B0.rD , log=T) )
    lpr.curr <- (log(factorial(sum(d.olddraws==d.ij[1])-1) + factorial(sum(d.olddraws==d.ij[2])-1) )
                 + Reduce(sum, lapply(seq_along(1:2), function(k) dmvnorm(beta.Df[[d.ij[k]]],  mean=pr.beta0Df, sigma=solve(pr.Sigma0InvDf), log=T)
                                      + dnorm(b.D[[d.ij[k]]], mean=0, sd=sqrt(sig2bD), log=T)
                                      + sum(dnorm(r.D[[d.ij[k]]],mean=rep(0,n), sd=sqrt(sig2rD[[d.ij[k]]]), log=T) )
                                      + actuar::dinvgamma(sig2rD[[d.ij[k]]], shape=pr.A0.rD , scale=pr.B0.rD, log=T ) )) )
    llik.merg <- sum(dnorm(Y[y.jlabs%in%S.ij], mean=matrix(fixed.mean[y.jlabs%in%S.ij],ncol=1 ) + kronecker(rep(1,length(S.ij)),SDf)%*%betaDf.merge 
                           + kronecker(rep(1,length(S.ij)),SD)%*%bD.merge + matrix(kronecker(rep(1,length(S.ij)),rD.merge),ncol=1), 
                           sd = rep(sqrt(sig2eps[S.ij]),each=n) ,log=T))
     
    llik.curr <-  Reduce(sum, lapply(seq(1,length(S.ij)), FUN = function(k)  dnorm(Y[y.jlabs==S.ij[k]], mean=matrix(fixed.mean[y.jlabs==S.ij[k]],ncol=1 ) + SDf%*%matrix(unlist(beta.Df[d.olddraws[S.ij[k]]]),ncol=1)
                                                                       + SD%*%unlist(b.D[d.olddraws[S.ij[k]]]) + unlist(r.D[d.olddraws[S.ij[k]]]), 
                                                                       sd = rep(sqrt(sig2eps[S.ij[k]]),each=n) ,log=T)))
    
    cat(' lqcur.mer=',lqcurr.merge,' lqmer.cur=',lqmerge.curr,' lpr.mer=',lpr.merg,' lpr.cur=',lpr.curr,' llik.mer=',llik.merg, ' llik.cur=',llik.curr)
    
    acc.merge <- min( 1,  exp((lqcurr.merge - lqmerge.curr) + (lpr.merg - lpr.curr) + (llik.merg - llik.curr)) )
    
    u <- runif(1)
    if(u < acc.merge){
      accepted <- 3
      d.olddraws[S.ij] <- d.merge
      d.olddraws <- as.numeric(factor(d.olddraws))
      beta.Df <- beta.Df[-d.ij[1]]
      b.D <- b.D[-d.ij[1]]
      r.D <- r.D[-d.ij[1]]
      sig2rD <- sig2rD[-d.ij[1]]
    } else{
      accepted <- 4
      d.olddraws <- d.olddraws
    }
    
  }
  #print(accepted)
  return( list(d.new=d.olddraws, betaDf.new=beta.Df, bD.new=b.D, rD.new=r.D, sig2rD.new=sig2rD, accepted=accepted)) 
}


