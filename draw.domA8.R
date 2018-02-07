


###########################################################################################################
#####     FUNCTION TO DRAW THE DOMAIN INDICATOR VARIABLE IN NONCONJUGATE MODEL USING GIBBS SAMPLING 
#####       Based on Neal (2000) algorithm 8
###########################################################################################################  


draw.domA8 <- function( Y, n, J, d.olddraws, y.jlabs, fixed.mean, alpha, SD, SDf, sig2eps, sig2rD, sig2bD, 
                        r.D, b.D, beta.Df, pr.Sigma0InvDf, pr.beta0Df,pr.A0.rD, pr.B0.rD, m, printSumms ){
  
  ## Y = vector of all outcomes (Y1, ... YJ); Yj = (y1, ..., yn)
  ## n  = number of subjects
  ## d.olddraws = vector of domain assignments for all outcomes prior to updating
  ## y.jlabs = outcome labels for Y
  ## fixed.mean =  current value of Xbeta + bO + r (all non-domain specific parameters)
  ## alpha = concentration parameter for the DP 
  ## SD, SDf = domain-specific covariate matrices
  ## sig2eps, sig2bD = current value of non-domain specific variances 
  ## sig2rD, r.D, b.D, beta.Df = current value of all domain specific
  ## pr.Sigma0InvDf, pr.beta0Df, pr.A0.rD, pr.B0.rD = hyperparameter values for domain-specific parameters
  ## printSumms = if true, print summaries pertaining to clustering algorithm
  ## m = tuning parameter in Neal's algorithm (number of draws in Monte Carlo simulation)
  
  for(j in 1:J){
    
    K=max(d.olddraws)
    if(printSumms==TRUE) cat('j=',j,': K.i=',K,' d.j=',d.olddraws[j]) 
    ## Tabulate the group sizes
    # j=1
    n.minus.j <- table(d.olddraws[-j])
    Kminus <- length(n.minus.j)
    h <- Kminus + m
    d.new.minus.j <- d.olddraws[-j]
    
    if( Kminus < K ){  
      ## The case where d_j != d_i for any other i:
      d.new.minus.j <- as.numeric(factor(d.olddraws[-j],labels=1:nlevels(factor(d.olddraws[-j]))))
      dj <- Kminus+1
      
      ## Need to draw new phi_c for m-1 possible new states:
      sig2rDnew <- as.list(rinvgamma(m-1, shape=pr.A0.rD , scale=pr.B0.rD ))
      sig2rDall <- c(sig2rD[-d.olddraws[j]], sig2rD[d.olddraws[j]], sig2rDnew )
      r.Dall <- c(r.D[-d.olddraws[j]], r.D[d.olddraws[j]], lapply(seq(1,(m-1)), FUN=function(x) rnorm(n, mean= 0, sd= sqrt(sig2rDnew[[x]]))))
      b.Dall <- c(b.D[-d.olddraws[j]], b.D[d.olddraws[j]], as.list(rnorm(m-1, mean= rep(0,pD), sd = sqrt(sig2bD))))
      beta.Dfall <- c(beta.Df[-d.olddraws[j]],beta.Df[d.olddraws[j]], lapply(seq(1,(m-1)), FUN=function(x) rmvnorm(1, mean=pr.beta0Df, sigma=solve(pr.Sigma0InvDf))))
    }  else {
      ## The case where d_j == d_i for some other i:
      
      ## Need to draw new phi_c for m possible new states:
      sig2rDnew <- as.list(rinvgamma(m, shape=pr.A0.rD , scale=pr.B0.rD ))
      sig2rDall <- c(sig2rD, sig2rDnew )
      r.Dall <- c(r.D, lapply(seq(1,m), FUN=function(x) rnorm(n, mean= 0, sd= sqrt(sig2rDnew[[x]]))))
      b.Dall <- c(b.D, as.list(rmvnorm(m, mean= rep(0,pD), sigma= diag(sig2bD,pD))))
      beta.Dfall <- c(beta.Df, lapply(seq(1,m), FUN=function(x) rmvnorm(1, mean=pr.beta0Df, sigma=solve(pr.Sigma0InvDf))))
    } 
    
    
    ## Log-likelihoods, and normalizing constant b_j for determining probability of membership for
    ##    d_j=d | all other d_i
    
    currmeans <- lapply(seq_along(1:h),FUN=function(d) matrix(fixed.mean[y.jlabs==j] + SDf%*%matrix(unlist(beta.Dfall[[d]]),pDF,1) + 
                                                                SD%*%matrix(unlist(b.Dall[[d]]),pD,1) + matrix(unlist(r.Dall[[d]],n,1))) )
    
    log.liks <- lapply(seq_along(1:h), FUN=function(d) sum(dnorm(Y[y.jlabs==j],mean=currmeans[[d]],sd=rep(sqrt(sig2eps[j]),n),log=T)))
    #     
    phi= -max(unlist(log.liks))
    e2phi <- lapply(X=1:h, FUN=function(d) exp(phi+log.liks[[d]]))
    weights <- c(as.vector((n.minus.j)),rep(alpha/m,m))/(J - 1 + alpha)
    
    
    ## Normalize using e^(c+phi) / sum_c e^(c+phi) :
    probs.j <- lapply(X=seq(1:h), FUN=function(d) weights[d]*e2phi[[d]]/sum(weights%*%unlist(e2phi)))
    if(printSumms==TRUE) print(unlist(probs.j))
    
    
    d.new.j <- sample(1:h,1,prob=probs.j)
    if(printSumms==TRUE) cat(' d.new.j=', d.new.j,'\n')
    
    if( Kminus<K ){ 
      # i.e. we were looking at a singleton to begin with
      if(d.new.j > (Kminus+1)){
        #  i.e. stay a singleton but at a different state
        r.D[[d.olddraws[j]]] <- r.Dall[[d.new.j]]
        b.D[[d.olddraws[j]]] <- b.Dall[[d.new.j]]
        beta.Df[[d.olddraws[j]]] <- beta.Dfall[[d.new.j]]
        sig2rD[[d.olddraws[j]]] <- sig2rDall[[d.new.j]]
        # I don't think we need to reassign the domain a new number 
      }else{
        if(d.new.j <= Kminus){
          # i.e. join an existing domain
          r.D <- r.D[-d.olddraws[j]]
          b.D <- b.D[-d.olddraws[j]]
          beta.Df <- beta.Df[-d.olddraws[j]]
          sig2rD <- sig2rD[-d.olddraws[j]]
          d.olddraws[j] <- unique(d.olddraws[-j][which(d.new.minus.j==d.new.j)])
          d.olddraws <- as.numeric(factor(d.olddraws))
          # I don't think we need to re-order the 
        }
        # else, do nothing because d_j stayed a singleton with it's original state
      } 
    }
    
    if(Kminus==K){
      if(d.new.j > Kminus){
        # i.e. an outcome breaks out to it's own domain
        beta.Df[[length(beta.Df)+1]] <- beta.Dfall[[d.new.j]]
        r.D[[length(r.D)+1]] <- r.Dall[[d.new.j]]
        b.D <- c(b.D,b.Dall[[d.new.j]])
        sig2rD <- c(sig2rD,sig2rDall[[d.new.j]])
        d.olddraws[j] <- Kminus+1
      }else{
        # i.e. an outcome joins a new domain or stays in its original domain
        d.olddraws[j] <- d.new.j
        
      }
    }
    
    
    if(printSumms==TRUE) cat('d.draws=',d.olddraws,'\n\n')
    ### We can now update the d.olddraws and K
    K <- nlevels(as.factor(d.olddraws))
    
    log.liks <- probs.j <- currmeans <-  NULL
    beta.Dfall <- b.Dall <- r.Dall <- sig2rDall <- NULL
  }
  return(list(d.new=d.olddraws, betaDf.new=beta.Df, bD.new=b.D, rD.new=r.D, sig2rD.new=sig2rD))
}

