

source('~/clustMultOutwImpSplitMerge.R')
## location of 9-yr SCDS dat
dat1 <- read.table('/projects/32885/sey2/Main/9Years/Amy/sasdata/missx.txt')


predictor.varnames <- c('MERC_S','CSEX','MAGE','HOME','KBIT_STD','M8HOLLSC','child_age')

## Subset further to ensure the ELIGIBLE subjects with complete covariate data have at least two outcomes in each domain,
##  where domains are defined as in Thurston (2009):
outsin09domains <- list(cognition=c('M8VERBIQ','M8PERFIQ','M8NOCUE'), memory=c('nnlog_m8rhitss','M8M5SS','M8SHRTSS','M8LONGSS','M8DSSC','M8WRMLSS'),motor=c('nlog_m8domsec','nlog_m8nonsec','M8DOMAVG','M8NONAVG','nlog_m8ttimea','nlog_m8ttimeb','TOT_ST','M8VMISS'),behavior=c('M8INT_T','M8EX_T','nlog_m8trt4') )
cat('The number of subjects left is',nrow(dat1),'which is equal to 533 from the 2014 paper.')
## Center and Scale outcomes:
predictorvars <- scale(dat1[,predictor.varnames])
head(predictorvars)

## Define outcome variables
outcome.varnames <- unlist(outsin09domains)
outcome.varnames

## Center and Scale outcomes:
outcomevars <- scale(dat1[,outcome.varnames])


J <- length(unlist(outsin09domains)) # outcomes
n <- nrow(dat1) # subjects

pF <- 1 # no. of overall fixed effects covariates (MeHg exposure)
pDF <- 6 # no. of covariates w/ domain specific fixed effects 
# sex, mother's age, HOME score (stimulation in home environment), KBIT (mother's IQ), 
# Hollingshead SES (paternal education and employment), and child's age at testing
pD <- 1 # no. of covariates w/ domain-specific random effects
pO <- 1 # no. of covariates w/ outcome-specific random effects

Y <- NULL
for(j in 1:J){
  Y <- c(Y, outcomevars[,unlist(outsin09domains)[j]])
}
Sf <- SD <- SO <- matrix(predictorvars[,'MERC_S'] , ncol=1)
SDf <- predictorvars[,-which(colnames(predictorvars)=='MERC_S')]
y.ilabs <- rep(1:n, J)
y.jlabs <- rep(1:J, each=n)


#### Priors from Sally's paper:
prior.NI <- list(beta0f=rep(0,pF), SigInv.betaf=diag(.00001,pF),beta0Df=rep(0,pDF), SigInv.betaDf=diag(.00001,pDF), 
                 sig2.bO=c(0.5,0.00005), sig2.bD=c(0.5,0.00005), sig2.eps=c(0.5,0.0005),sig2.r=c(0.5,0.00005), 
                 sig2.rD=c(0.5,0.00005), alpha=1)

domAssn <- rep(1:10,length.out=J)
Kstart=length(unique(domAssn))
set.seed(956)
inits <- list(D.init=Kstart, dvec=domAssn, betaF = rep(0,pF), bO = rep(0,pO*J), 
              r = rep(0,n), sig2bO=rep(.05,pO), sig2bD=rep(.05,pD), sig2eps=rep(.05,J), sig2r=.05,
              betaDf = lapply(1:Kstart, function(i) as.vector(matrix(0,nrow=Kstart,ncol=pDF)[i,,drop=FALSE])), 
              bD = lapply(1:Kstart, function(i) as.vector(matrix(0,nrow=Kstart,ncol=pD)[i,,drop=FALSE])), 
              sig2rD = lapply(1:Kstart, function(i) as.vector(matrix(0.05,nrow=Kstart,ncol=1)[i,,drop=FALSE])), 
              rD=lapply(1:Kstart, function(i) as.vector(matrix(0,nrow=Kstart,ncol=n)[i,,drop=FALSE])))

clust1 <- cluster.multouts( Y=Y, Sf=Sf, SDf=SDf, SD=SD, SO=SO, y.jlabs=y.jlabs, y.ilabs=y.ilabs, 
                            priors=prior.NI, inits=inits, niter=30000, printIter=TRUE, tune=c(30,30),
                            m=6,  printSumms=FALSE, ycorrs=NULL)

save(clust1, file='clusteringDraws.Rda')

## Remove burnin draws:
nburn = 1e3
d.drawsMinusBurnin <- clust1$d.draws[(nburn + 1):niter,] 
## Calculate number of groups at each iteration:
K.draws <- apply(d.drawsMinusBurnin,1,max)

cat('The sampler visited states which had ',paste('d=',unique(K.draws),' ',sep=''),'domains \nafter the ',nburn,'-draw burnin period.',sep='')
print(table(K.draws))

#### Clustering via least-squares with all draws after the burnin.
d.drawsNewAllK <- renumber(d.drawsMinusBurnin, Kpost = NULL)
grpinds <- d.drawsNewAllK$gpindices
grppatts <- d.drawsNewAllK$gppatterns

S.hat <- 1 - dissMatrix(d.drawsMinusBurnin)
LSclustresults <- dahlLSclust(d.drawsMinusBurnin, S.hat)
LSminAll.1 <- which(LSclustresults==min(LSclustresults))[1]
dom.assnsLS <- d.drawsNewAllK$newds[LSminAll.1,]
write.table(dom.assnsLS, file = 'final_clust_all.txt')


###### ------ Alternatively, find posterior grouping using only draws of Kpost groups ----- #######

Kpost <- as.numeric(names(which.max(table(K.draws))))
Krows7 <- which(K.draws==Kpost)

#### Clustering via least-squares does not require an estimate of K-posterior 
###   and we can incorporate all draws after the burnin.
d.drawsNewAllK7 <- renumber(d.drawsMinusBurnin[Krows7,], Kpost = Kpost)
grpinds <- d.drawsNewAllK7$gpindices
grppatts <- d.drawsNewAllK7$gppatterns

S.hat <- 1 - dissMatrix(d.drawsMinusBurnin[Krows7,])
LSclustresults <- dahlLSclust(d.drawsMinusBurnin[Krows7,], S.hat)
LSmin7 <- which(LSclustresults==min(LSclustresults))[1]
dom.assnsLS.7 <- d.drawsNewAllK7$newds[LSmin7,]
write.table(dom.assnsLS.7, file = 'final_clust_Kpost.txt')
