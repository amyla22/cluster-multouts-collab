
source('~/clustMultOutwImpSplitMerge.R')

dat1 <- read.table('/projects/32885/sey2/Main/9Years/Amy/sasdata/missx.txt')


predictor.varnames <- c('MERC_S','CSEX','MAGE','HOME','KBIT_STD','M8HOLLSC','child_age')

## Subset further to ensure the ELIGIBLE subjects with complete covariate data have at least two outcomes in each domain,
##  where domains are defined as in Thurston (2009):
outsin09domains <- list(cognition=c('M8VERBIQ','M8PERFIQ','M8NOCUE'), memory=c('nnlog_m8rhitss','M8M5SS','M8SHRTSS','M8LONGSS','M8DSSC','M8WRMLSS'),motor=c('nlog_m8domsec','nlog_m8nonsec','M8DOMAVG','M8NONAVG','nlog_m8ttimea','nlog_m8ttimeb','TOT_ST','M8VMISS'),behavior=c('M8INT_T','M8EX_T','nlog_m8trt4') )
cat('The number of subjects left is',nrow(dat1),'which is equal to 533 from the 2014 paper.')
## Center and Scale outcomes:
num.predictorvars <- scale(dat1[,predictor.varnames[-2]])
predictorvars <- cbind(MERC_S=num.predictorvars[,1],CSEX=dat1[,'CSEX'],num.predictorvars[,2:6])
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


#### In Sally's paper:
priors <- list(beta0f=rep(0,pF), SigInv.betaf=diag(.00001,pF),beta0Df=rep(0,pDF), SigInv.betaDf=diag(.00001,pDF), 
                 sig2.bO=c(0.5,0.00005), sig2.bD=c(0.5,0.00005), sig2.eps=c(0.5,0.0005),sig2.r=c(0.5,0.00005), 
                 sig2.rD=c(0.5,0.00005), alpha=1)

### Choose the overall a posteriori clustering OR 
domAssn <- read.table('fixed_clust_all.txt')
#    the estimate from those draws with Kpost (posterior estimate of # of groups) groups
# domAssn <- read.table('fixed_clust_Kpost.txt')

Kstart=length(unique(domAssn))
set.seed(956)
inits <- list(betaF = rep(0,pF), bO = rep(0,pO*J), 
              r = rep(0,n), sig2bO=rep(.05,pO), sig2bD=rep(.05,pD), sig2eps=rep(.05,J), sig2rD=rep(.05,K), sig2r=.05,
              betaDf = lapply(1:K, function(i) as.vector(matrix(0,nrow=K,ncol=pDF)[i,,drop=FALSE])), 
              bD = lapply(1:K, function(i) as.vector(matrix(0,nrow=K,ncol=pD)[i,,drop=FALSE])), 
              rD=lapply(1:K, function(i) as.vector(matrix(0,nrow=K,ncol=n)[i,,drop=FALSE])))

fixedrun1 <- sally.multouts( Y=Y, Sf=Sf, SDf=SDf, SD=SD, SO=SO, y.jlabs=y.jlabs, y.ilabs=y.ilabs, 
                          priors=priors, inits=inits, niter=10e3, 
                          printIter=TRUE, dom.assns=domAssn, ycorrs=25, thin=1)

save(fixedrun1, file='fixedModelDraws.Rda')