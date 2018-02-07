###########################################################################################################
##### prelabels() is a function which creates the parameter names for a given set of parameters which 
#####           may be indexed by j or require a squared symbol
###########################################################################################################

prelabels <- function(symbol,j,sq=TRUE) {
  if(sq==TRUE){
    paste(paste(paste(symbol,j,sep='['),']',sep=''),'^2',sep='')
  } else{
    paste(paste(symbol,j,sep='['),']',sep='')
  }
}


###########################################################################################################
##### plot.params() is a function which creates the plots for multiple (or single parameters) all at once
###########################################################################################################

plot.params <- function(draws, par.name, lim.y, plotdims=NULL, plotcolors=NULL){
  if(is.matrix(draws)){
    par(mar=c(2, 2.5, 1.5, 1.5))
    if(is.null(plotdims)){
      dimRows <- ncol(draws)/2
      par(mfrow=c(dimRows,2))
    }else{
      par(mfrow=plotdims)
    }
    niter = nrow(draws)
    if(sum(lim.y)==0){
      if(sum(grepl('sig',par.name))==0){
        lim.y=c(min(draws),max(draws))
      }else{
        lim.y=c(0,max(draws))
      }
    }
    if(is.null(plotcolors)){
      pcolor = 1
    } else {
      pcolor = plotcolors
    }
    for(j in 1:ncol(draws)){
      max.y = max(draws[,j])
      plot(1:niter, draws[,j], type='l',xlab='',ylab='',ylim=lim.y,cex.axis=.9, col=pcolor)
      text(niter*.9,max(lim.y)*.95 ,parse(text=paste(par.name[j])))
    }
  }else{
    if(sum(lim.y)==0){
      if(sum(grepl('sig',par.name))==0){
        lim.y=c(min(draws),max(draws))
      }else{
        lim.y=c(0,max(draws))
      }
    }
    niter = length(draws)
    plot(1:niter, draws, type='l',xlab='',ylab='',ylim=lim.y,cex.axis=.9)
    text(niter*.1,max(lim.y)*.95 ,parse(text=paste(par.name)))
  }
}

###########################################################################################################
#####  covsbydomain.plot() is a function to summarize the domain-specific fixed covariate effects visually 
##        with Means and 95% posterior intervals for each domain within each covariate
###########################################################################################################
covsbydomain.plot <- function( param.draws, ndom, covnames, fullcovnames, xrange){
  ## param.draws = matrix of all covariate effects of interest for all domains
  ## ndom = number of domains as discovered in the MCMC clustering
  ## covnames = predictor variable names as found in param.draws
  ## pcolor = colors associated with domains (if NULL these will be assigned per rainbow())
  ## fullcovnames = names to appear in title
  ## xrange = range of xvalues to be plotted over across all covariates for consistency
  p <- length(covnames)
  par(mfrow=c(ceiling(p/2), 2))
  par(mar=c(2, 4, 1.5, 1.5))
  
  
  pcolor <- 1:ndom
    #c('slateblue','seagreen','deeppink','darkorchid4','yellow','seagreen','navyblue','lavender','magenta','yellow')
  
  for(i in 1:p){
    ypts <- ndom:1
    xpts <- apply(param.draws[,which(colnames(param.draws)==covnames[i])], 2, mean )
    lower <- apply(param.draws[,which(colnames(param.draws)==covnames[i])],2,function(x) quantile(x,probs=0.025)) 
    upper <- apply(param.draws[,which(colnames(param.draws)==covnames[i])],2,function(x) quantile(x,probs=0.975))
    plot(x=xpts, y = ypts, xlim = xrange, ylim=c(0,ndom+1),  xlab='', ylab='', type='n',yaxt='n')
    title(ylab="Domain", line=2, cex.lab=1.2)
    axis(2,at=1:ndom,labels=ypts,las=2,cex.axis=1)
    
    for(d in 1:ndom){
      points(xpts[d],ypts[d], col=pcolor[d], pch=16)
      segments(lower[d], ypts[d] ,upper[d], ypts[d] ,col=pcolor[d],lwd=2)
    }
    abline(v=0); abline(v=.1); abline(v=-.1)
    if(is.null(fullcovnames)){
      title(main=covnames[i],cex=.8)
    }else{
      title(main=fullcovnames[i],cex=.8)
    }
  }
  
}


###########################################################################################################
#####  outcomesindomain.plot() is a function to summarize the results visually with Means/Medians and 95% posterior intervals for the 
##        Bayesian multiple outcome results with the option to compare side-by-side with the frequentist results.
###########################################################################################################

outcomesindomain.plot <- function(param.draws,nexpos,ndoms=K,nouts=J, legendloc="topleft",longnames=NULL, outcomenames=NULL,
                                  domcolors, legendsize=1, plottitle, compareBoth=NULL, xrange=NULL, dom.assn, orig.assn, fixedEffect=FALSE){
  ##### This function draws a plot to visualize/compare the posterior intervals
  ##### and estimates of a single parameter across different simulation scenarios
  ## param.draws = matrix of the draws of the slope parameters (overall fixed slope, domain-specific RE and outcome-specific RE)
  ## nexpos = number of exposures considered (e.g. mercury in the SCDS)
  ## longnames = names to appear in the legend if domains are assigned names
  ## outcomenames = names to appear on the lhs plot y-axis if not set to NULL
  ## ndoms = number of domains
  ## dom.assn = domain assignment  <------- think about/ask Tanzy
  ## nouts = number of outcomes
  ## range = the amount of 'space' to add to each end of the most extreme 
  #   posterior intervals -- this is mainly for appearance
  ## intercepts could plot the true values
  ## xaxis = a vector that can paste a word to numbers to label the x-axis
  #   if not sequentially labeled, need to use label word + at least 3 numbers:
  #    e.g. ("Prior",c(2,4,8))
  ## xrange = range of x axis
  ## orig.assn = determines colors for original domain assignments
  ## fixedEffect = FALSE if plotting the overall mercury effects for each outcome
  #   if TRUE the function can plot overall domain-specific fixed effects (e.g. sex, age, etc.)
  #   if TRUE, must specify the xrange

  if(!is.null(compareBoth)){
    par(mfrow=c(1,2))
    par(mar=c(2, 4, 1.5, 0))
  } 
  outcome.vec <- seq(1,nouts)
   
  if(fixedEffect==FALSE){
    betaf <- param.draws[,grep("betaf",colnames(param.draws))]
    bDmat <- param.draws[,grep("bD",colnames(param.draws))]
    bOmat <- param.draws[,grep("bO",colnames(param.draws))]
    expcoef.outs <- matrix(NA, nrow=nrow(param.draws), ncol=J)
    for( j in 1:J){
      expcoef.outs[,j] <- betaf + bOmat[,j] + bDmat[,dom.assn[j]]
    }
  }else{
    expcoef.outs <- matrix(NA, nrow=nrow(param.draws), ncol=nouts)
    for( j in 1:nouts){
      expcoef.outs[,j] <- param.draws[,dom.assn[j]]
    }
  }
  lower <- apply(expcoef.outs,2,function(x) quantile(x,probs=0.025)) 
  upper <- apply(expcoef.outs,2,function(x) quantile(x,probs=0.975))
  means.all <- apply(expcoef.outs,2,function(x) quantile(x,probs=0.5))
  
  if(is.null(outcomenames)){
    ylabels <- c(paste('y',seq(1:nouts),sep='_'))
  }else{
    ylabels <- outcomenames  
  }
  if(is.null(xrange)){
    xseq <- seq(min(lower,quantile(betaf,probs=.025))-.1,max(upper,quantile(betaf,probs=.975))+.1,length.out=(nouts+1))
  } else {
    xseq <- seq(xrange[1],xrange[2], length.out=(nouts+1))
  }
  
  plot(x=xseq,y=seq(1,nouts+3,length.out=(nouts+1))-1,type="n",xlab="",ylab="",yaxt="n",cex.axis=.8)
  axis(2,at=outcome.vec,labels=ylabels,las=2,cex.axis=.48)
  segments(0,outcome.vec[1],0,outcome.vec[length(outcome.vec)],col='lightgray')
  
  ### longnames = NULL if just assigning Domain 1:K, otherwise it is set to the domain names
  if(is.null(longnames))   longnames <- 1:ndoms
  if(is.null(domcolors)){
    pcolor = rainbow(ndoms)
  } else { 
    pcolor = domcolors
  }
  legend(legendloc[1], inset=.005, ncol=ndoms, legend=longnames , fill=pcolor,cex=legendsize,bty='n', title="Domain")
  title(main=plottitle[1],cex.main=as.numeric(plottitle[3]))
#   points(median(betaf),outcome.vec[1],cex=.8,pch=16)
#   segments(quantile(betaf,probs=.025), outcome.vec[1],quantile(betaf,probs=.975),outcome.vec[1]) 
  
  for(j in 1:(nouts)){
    segments(xrange[1], outcome.vec[j], xrange[2],outcome.vec[j],col='lightgray') 
    points(means.all[j], outcome.vec[j],col=pcolor[dom.assn[j]],pch=12+dom.assn[j],cex=1)
    segments(lower[j], outcome.vec[j],upper[j],outcome.vec[j],col=pcolor[dom.assn[j]],lwd=2) 
  }
  pcolor2 =c('green','cyan2','purple','red','orange')
  
  if(!is.null(compareBoth)){
    par(mar=c(2, 2.5, 1.5, 1.5))
    SEs <- compareBoth[,2]
    betaEsts <- compareBoth[,1]
    plot(x=xseq,y=seq(1,nouts+3,length.out=(nouts+1))-1,type="n",xlab="",ylab="",yaxt="n",cex.axis=.8)
    segments(0,outcome.vec[1],0,outcome.vec[length(outcome.vec)],col='lightgray')
    #axis(2,at=outcome.vec,labels=ylabels,las=2,cex.axis=.68)
    title(main=plottitle[2],cex.main=as.numeric(plottitle[3]))
    for(j in 1:(nouts)){
      segments(xrange[1], outcome.vec[j],xrange[2],outcome.vec[j],col='lightgray',cex=1) 
      points(betaEsts[j], outcome.vec[j],col=pcolor2[orig.assn[j]],pch=14+orig.assn[j], cex=1)
      segments(betaEsts[j]-SEs[j]*qnorm(.975), outcome.vec[j],betaEsts[j]+SEs[j]*qnorm(.975),outcome.vec[j],col=pcolor2[orig.assn[j]],lwd=2) 
    }
    if(max(orig.assn)==4){
      orig.doms <- c('cognition','memory','motor','behavior')
    }else{
      orig.doms <- c('cognition','memory','motor','behavior', 'attention')
    }
    legend(legendloc[1], inset=.005, ncol=ceiling(max(orig.assn)/2), legend=orig.doms, fill=pcolor2,cex=legendsize*1,bty='n', title="Original Domain")
  }
}



# --------------------------------------------------------------------------------------------------------------- #
############### Function to summarize the results into tables with Means/Medians, SEs and 95% posterior intervals.
summarize.sims <- function(results1, results2, dom.assn, outvarnames, sig2rDeq){
  
  nouts <- length(dom.assn)
  Kpost <- max(dom.assn)
  
  if(sig2rDeq==FALSE){
    param.names1 <- c("$\\beta_f$", paste("$b_{\\mathcal{D},",1:Kpost,"}$",sep=''), "$\\sigma_{b,\\mathcal{D}}$", 
                      "$\\sigma_{b,\\mathcal{O}}$","$\\sigma_r$",paste("$\\sigma_{r,\\mathcal{D},",1:Kpost,"}$",sep=''),
                      paste("$\\sigma_\\epsilon,",1:nouts,"$",sep=''))
  }else{
    param.names1 <- c("$\\beta_f$", paste("$b_{\\mathcal{D},",1:Kpost,"}$",sep=''), "$\\sigma_{b,\\mathcal{D}}$", 
                    "$\\sigma_{b,\\mathcal{O}}$","$\\sigma_r$","$\\sigma_{r,\\mathcal{D}}$",
                    paste("$\\sigma_\\epsilon,",1:nouts,"$",sep=''))
  }
  nparams1 <- length(param.names1)
  
  table1 <- matrix(NA,nrow=nparams1,ncol=4)
  colnames(table1) <- c("Parameter", "Estimate", "SE","95\\% Interval")
  
  
  ci.lower1 <- apply(results1, 2, function(x) quantile(x, probs=.025))
  ci.upper1 <- apply(results1, 2, function(x) quantile(x, probs=.975))
  
  sdparams <- grep('sigma', param.names1)
  nonsds <- grep('sigma', param.names1, invert=TRUE)
  
  table1[,"Parameter"] = param.names1
  table1[nonsds,"Estimate"] = round(apply(X=results1[,nonsds],MARGIN=2,FUN=mean),3)
  table1[sdparams,"Estimate"] = round(apply(X=results1[,sdparams],MARGIN=2,FUN=median),3)
  table1[,"95\\% Interval"] = paste('(',round( ci.lower1,3),', ',round(ci.upper1,3),')', sep='')
  table1[,"SE"] = round(apply(X=results1,MARGIN=2,FUN=sd),3)
  
  
  table2 <- matrix(NA,nrow=nouts,ncol=6)
  colnames(table2) <- c("Outcome","Original Domain","New Domain", "MeHg Effect Estimate", "SE","95\\% Interval")
  ci.lower2 <- apply(results2, 2, function(x) quantile(x, probs=.025))
  ci.upper2 <- apply(results2, 2, function(x) quantile(x, probs=.975))
  
  # Parse domain names from the labels from the outcome variables:
  domnames <- apply(matrix(names(outvarnames),ncol=1),1,function(x) substring(x, 1,nchar(x)-1) )
  betaf <- param.draws[,grep("betaf",colnames(param.draws))]
  bDmat <- param.draws[,grep("bD",colnames(param.draws))]
  bOmat <- param.draws[,grep("bO",colnames(param.draws))]
  coef.outs <- matrix(NA, nrow=nrow(param.draws), ncol=J)
  for(j in 1:J){
    coef.outs[,j] <- betaf + bOmat[,j] + bDmat[,dom.assn[j]]
  }
  ci.lower2 <- apply(coef.outs, 2, function(x) quantile(x, probs=.025))
  ci.upper2 <- apply(coef.outs, 2, function(x) quantile(x, probs=.975))
  
  table2[,"Outcome"] = outvarnames
  table2[,'Original Domain'] = domnames
  table2[,'New Domain'] = dom.assn
  table2[,"MeHg Effect Estimate"] = round(apply(X=coef.outs,MARGIN=2,FUN=median),3)
  table2[,"95\\% Interval"] = paste('(',round( ci.lower2,3),', ',round(ci.upper2,3),')', sep='')
  table2[,"SE"] = round(apply(X=coef.outs,MARGIN=2,FUN=sd),3)
  
  
  return(list(table1, table2))
}

###########################################################################################################
##### is.same() is a function which takes the final posterior clustering vector and a row of d.draws and 
#####    returns true if their levels "match". This replaces the identical() function, which will not work
#####    for the labelings since they may not be equal in value but are in fact equal in terms of domains.
###########################################################################################################

is.same <- function(rowd, dom.assns){
  nsame <- 0
  
  for(i in 1:length(rowd)){
    if(as.numeric(reorder(factor(rowd),dom.assns))[i]==dom.assns[i]){
      nsame <- nsame + 1
    }
  }
  return(nsame)
}

###########################################################################################################
##### summAll() is a function which takes the final posterior draws from the MCMC (after the posterior 
#####     clustering, which may have reduced our number of draws removed after the burnin) and produces key
#####     summary output.
###########################################################################################################


summAll <- function(mcmc.results,dom.assns){
  varcols <- grep('sig2',colnames(mcmc.results))
  vardraws <- mcmc.results[,varcols]
  coefdraws <- mcmc.results[,-varcols]
  
  resultsMat <- matrix(c(mean(coefdraws[,1]),sd(coefdraws[,1]),quantile(coefdraws[,1],probs=c(0.025,0.975))),nrow=1)
  for(j in 2:ncol(coefdraws)){
    resultsMat <- cbind(resultsMat, cbind(mean(coefdraws[,j]),sd(coefdraws[,j]),quantile(coefdraws[,j],probs=c(0.025)),quantile(coefdraws[,j],probs=c(0.975))))
  }
  for(j in 1:ncol(vardraws)){
    resultsMat <- cbind(resultsMat, cbind(median(vardraws[,j]),sd(vardraws[,j]),quantile(vardraws[,j],probs=c(0.025)),quantile(vardraws[,j],probs=c(0.975))))
  }
  resultsMat <- cbind(resultsMat, matrix(dom.assns,nrow=1))
  colnames(resultsMat) <- c(paste(rep(colnames(mcmc.results),each=4),rep(c('est','se','low','up'),ncol(mcmc.results)),sep='.'),
                            paste('d',1:length(dom.assns),sep='-')) 
  return(resultsMat)
}

###########################################################################################################
##### dissMatrix() is a function which takes the matrix of group assignments at each iteration and 
#####         calculates a dissimilarity matrix, which is a summary of the proportion of times the 
#####         observations are NOT grouped together (i.e. high probabilities = very dissimilar)
###########################################################################################################

library(cluster)  # need for function pam() for the posterior clustering of group assignments

dissMatrix <- function(foodat,printF=FALSE){
  nvars <- ncol(foodat)
  niter <- nrow(foodat)
  propMat <- matrix(0, nrow=nvars, ncol=nvars)
  for(i in 1:niter){
    if(printF==TRUE) cat('i=',i,'; ')
    distMat <- matrix(0, nrow=nvars, ncol=nvars)
    for(j in 1:nvars){
      distMat[j,which(foodat[i,]==foodat[i,j])] <- 1  
    }
    propMat <- propMat + distMat
  }
  propMat <- propMat/niter
  dissMat <- 1 - propMat
  return(dissMat)
}

###########################################################################################################
##### dahlLSclust() is a function which takes the matrix of group assignments at each iteration and 
#####         computes least squared error between the association matrix at each iteration, t, and
#####         the estimated similarity matrix, which is a summary of the proportion of times the 
#####         observations are grouped together (i.e. 1 - dissimilarity matrix)
###########################################################################################################

dahlLSclust <- function(d.draws, S.hat ){
  nouts <- ncol(d.draws)
  niter <- nrow(d.draws)
  sums.t <- rep(NA, niter)
  for(t in 1:niter){
    ## Create association matrix at each iteration:
    assocMat <- matrix(0, nrow=nouts, ncol=nouts)
    for(j in 1:nouts){
      assocMat[j,which(d.draws[t,]==d.draws[t,j])] <- 1  
    }
    sumij <- 0
    for(i in 1:nouts){
      sumij <- sumij + sum((assocMat[i,] - S.hat[i,])^2  )
    }
    sums.t[t] <- sumij
  } 
  return(sums.t)
}
###########################################################################################################
##### renumber() is a function which takes the matrix of group assignments at each iteration and 
#####         relabels them in increasing order to make the groups more comparable across iterations
###########################################################################################################

renumber <- function(domdraws, Kpost=NULL){
  assn1st <- apply(domdraws, 1, function(x) seq_along(x)[!duplicated(x)] ) # index of 1st unique outcomes
  domorder <- apply(domdraws, 1, function(x) x[seq_along(x)[!duplicated(x)]] ) # value of unique outcomes in order
  relabds <- matrix(NA, nrow=nrow(domdraws), ncol=ncol(domdraws))
  J <- ncol(domdraws)
  
  for(i in 1:nrow(domdraws)){
    # i = 1
    if(is.null(Kpost)){
      Kpost <- length(domorder[[i]])
      relabds[i, ] <- factor(domdraws[i,], levels = domorder[[i]], labels = 1:Kpost)
      Kpost <- NULL
    } else {
     # Kpost <- length(domorder[,i])
      relabds[i, ] <- factor(domdraws[i,], levels = domorder[,i], labels = 1:Kpost)
    }
  }
  
  p=1
  relabs2gps <- relabds
  patterns <- list()
  groups <- list()
  k=1
  while(nrow(relabs2gps)>0){
    p1 <- relabs2gps[1,] #assn1st[,p]
    patterns[[k]] <- p1
    relabs2gps <- matrix(relabs2gps[-1,],ncol=J)
    gp <- which(apply(relabs2gps,1,function(x) sum(x==p1)==J))
    groups[[k]] <- c(p,gp+p)
    if(length(gp)>0)  relabs2gps <- matrix(relabs2gps[-gp,],ncol=J)
    p <- p + length(gp) + 1
    k <- k + 1
  }
  
  return(list(newds=relabds,gpindices=groups,gppatterns=patterns))
}

###########################################################################################################
##### minmax.corrs() is a function which takes the correlation matrix observed/sampled draws of the Y 
#####         outcomes and identifies highly (+/-) correlated outcomes returning the values and names
###########################################################################################################


minmax.corrs <- function(corrMat, roundto){
  #### corrMat = correlation matrix to evaluate
  #### roundto = number of digits to round the correlation values to
  
  minmaxMat <- matrix(NA, ncol=9, nrow=nrow(corrMat))
  colnames(minmaxMat) <- c('Outcome', 'Min. (-) Corr.', 'Min. (-) Corr. Outcome','Max. (+) Corr.', 'Max. (+) Corr. Outcome',
                           'Min. Corr.', 'Min. Corr. Outcome', 'Max. Corr.', 'Max. Corr. Outcome')
  minmaxMat[,'Outcome'] <- colnames(corrMat)
  p <- ncol(corrMat)
  for(i in 1:(p)){
    minmaxMat[i,'Min. (-) Corr.'] <- round(min(corrMat[i,]),roundto)
    minmaxMat[i,'Min. (-) Corr. Outcome'] <- colnames(corrMat)[which(corrMat[i,]==min(corrMat[i,]))]
    minmaxMat[i,'Max. (+) Corr.'] <- round(max(corrMat[i,-i]),roundto)
    minmaxMat[i,'Max. (+) Corr. Outcome'] <- colnames(corrMat)[which(corrMat[i,]==max(corrMat[i,-i]))]
    minmaxMat[i,'Min. Corr.'] <- round(corrMat[i,which(min(abs(corrMat[i,]))==abs(corrMat[i,]))],roundto)
    minmaxMat[i,'Min. Corr. Outcome'] <- colnames(corrMat)[which(abs(corrMat[i,])==min(abs(corrMat[i,])))]
    minmaxMat[i,'Max. Corr.'] <- round(corrMat[i,which(abs(corrMat[i,])==max(abs(corrMat[i,-i])))],roundto)
    minmaxMat[i,'Max. Corr. Outcome'] <- colnames(corrMat)[which(abs(corrMat[i,])==max(abs(corrMat[i,-i])))]
    
  }
  return(minmaxMat)
}
# 
# postysMat <- NULL
# postysvec <- NULL
# for(t in 1:5){
#   for(j in 1:J){
#     postysvec <- c(postysvec, rnorm(n, mean=Sf%*%betaf.hat + SDf%*%betaDf.hat[domain.assns[j],] 
#                                     + SD%*%bD.hat[domain.assns[j],] + SO%*%bO.hat[j],sd=rep(sqrt(eps.hat[j]),n)))
#   }
#   postysMat <- cbind(postysMat,postysvec) 
#   postysvec <- NULL
# }


###############################################################################################################
##### postcorrs() is a function which takes the matrix of the draws from the posterior predictive distribution 
#####         of Y and identifies each outcomes most highly (+/-) correlated outcomes returning the values and 
#####         names of these outcomes in 2 separate list of length J. The real data is also evaluated in the 
#####         same way with the results returned as a matrix. This function uses minmax.corrs()
###############################################################################################################

postcorrs <- function(postysMat, outnames, realYs ){
  #### postysMat = posterior predictive draws
  #### outnames = names of the J outcomes
  #### realYs = observed values of Y
  
  ysMat <- matrix(postysMat[,1],nrow=n,ncol=J,byrow=F)
  corrMat <- cor(ysMat)
  colnames(corrMat) <- outnames
  minmaxMat <- minmax.corrs(corrMat,roundto=10)
  outslist <- lapply(1:nrow(minmaxMat), FUN = function(j) c(as.numeric(minmaxMat[j,c(2,4,6,8)])))
  outnameslist <- lapply(1:nrow(minmaxMat), FUN = function(j) c(minmaxMat[j,c(3,5,7,9)]))
  
  realYsMat <- matrix(realYs,nrow=n,ncol=J,byrow=F)
  realcorrMat <- cor(realYsMat,use='pairwise.complete.obs')
  colnames(realcorrMat) <- outnames
  realminmaxMat <- minmax.corrs(realcorrMat, roundto = 10)
  actualsMat <- realminmaxMat[,c(1,2,4,6,8,3,5,7,9)]
  
  
  for(c in 2:ncol(postysMat)){
    ysMat <- matrix(postysMat[,c],nrow=n,ncol=J,byrow=F)
    corrMat <- cor(ysMat)
    colnames(corrMat) <- outnames
    minmaxMat <- minmax.corrs(corrMat,roundto=10)
    outslist <- lapply(1:nrow(minmaxMat), FUN = function(j) rbind(outslist[[j]],as.numeric(minmaxMat[j,c(2,4,6,8)])))
    outnameslist <- lapply(1:nrow(minmaxMat), FUN = function(j) rbind(outslist[[j]],minmaxMat[j,c(3,5,7,9)]))
  }
  
  return(list(minmaxCorrs=outslist, minmaxnames=outnameslist, actualminmax=actualsMat)) 
}


