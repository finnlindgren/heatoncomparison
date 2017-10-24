rm(list=ls())
## Libraries
library(mcmc)
library(MASS)
library(KernSmooth)
library(fields)
library(pscl)
library(mvtnorm)
library(spBayes)
library(MCMCpack)
library(Mposterior)
library(parallel)
library(doParallel)
library(foreach)

## Satellite Data
load("SatelliteTemps.RData")
## object is called "sat.temps"
n <- dim(sat.temps)[1]
na.obs <- which(is.na(sat.temps$Temp))
length(na.obs)
## [1] 44431
## Nothing else missing in Lat-Lon
##range(dat.nomiss$Temp)
##[1] 24.37 55.41
## Training Data
dat.miss <- sat.temps[na.obs,] 
## Test Data     
dat.nomiss <- sat.temps[-na.obs,] 


## Useful Quantities
## Total number of training samples
n.sample <- nrow(dat.nomiss) 
## Total number of test samples
n.test <- nrow(dat.miss)  
## Number of cores or subsets
n.core <- 30  
per.core <- floor(nrow(dat.nomiss)/n.core)  
## Number of observations in different subsets
n.part <- c(rep(per.core,n.core-1),n.sample-per.core*(n.core-1))  
## This is same as n.core
n.split <- length(n.part)
## Number of MCMC iterations
mcmc.sample <- 2000  
## Burn in
n.burn <- mcmc.sample*0.5  
## Divide predicted data and predict them
## independently in each subset
p.sub <- c(0,seq(1000,44000,by=1000),n.test)


## GP regression on subsamples

a = 1:n.sample 
sample.loc <- dat.nomiss[,1:2] 
y <- dat.nomiss[,3]
## Training response
Y.train <- y 
## Training coordinates
coords.train <- sample.loc
## Test coordinates
coords.pred <- dat.miss[,1:2] 
## Test sample size
n.test <- nrow(dat.miss) 
index.part <- list()
X.part <- list()
Y.part <- list()
coords.part <- list()

for(i in 1:n.split){
  beg<-Sys.time()
  index.part[[i]] <- sample(a,n.part[i],replace=FALSE)  
  ## Response in i th subset
  Y.part[[i]] <- Y.train[index.part[[i]]]   
  ## Coordinates in i th subset
  coords.part[[i]] <- coords.train[index.part[[i]],] 
  ## Predictor in i th subset
  X.part[[i]] <- cbind(rep(1,n.part[i]),coords.part[[i]][,2],coords.part[[i]][,1]) 
  a <- setdiff(a,index.part[[i]])
}

## Predictors in the test set
X.pred <- cbind(rep(1,n.test),coords.pred[,2],coords.pred[,1])  

##########################################
##Partitioned GP function
##Works with subset i, for i=1,...,n.core
##########################################

partitioned_GP <- function(i){
  ## Model fitting

  library(spBayes)
  starting <- list("phi"=3, "sigma.sq"=5, "tau.sq"=1) 
  tuning <- list("phi"=0.01, "sigma.sq"=0.01, "tau.sq"=0.01)
  priors.1 <- list("beta.Norm"=list(rep(0,3), diag(1000,3)),
                   "phi.Unif"=c(3/10, 3/0.1), "sigma.sq.IG"=c(2, 2),
                   "tau.sq.IG"=c(2, 0.1))
  cov.model <- "exponential"
  ## Response in subset i
  ZZ <- Y.part[[i]]
  ## Predictor in subset i 
  XX <- X.part[[i]]
  ## Coordinates in subset i 
  CC <- coords.part[[i]]
  CC <- cbind(CC$Lon,CC$Lat)
  ## GP computation in each subset
  m.1 <- spLM(ZZ~XX-1, coords=CC, starting=starting,
              tuning=tuning, priors=priors.1, cov.model=cov.model,
              n.samples=mcmc.sample,verbose=FALSE)
  ## Recover all MCMC samples in each subset
  m.1.samp <- spRecover(m.1, start=n.burn+1, verbose=FALSE)
  ## List of MCMC iterates from subset i
  subAtom <- cbind(m.1.samp$p.beta.recover.samples,m.1.samp$p.theta.recover.samples) 
  ## Delete this quantity from the memory
  rm(m.1.samp)
  ## Garbage cleaning 
  gc()

  ## Prediction
  Y.pred <- t(spPredict(m.1, pred.covars=X.pred[(p.sub[1]+1):p.sub[2],],
                           pred.coords=coords.pred[(p.sub[1]+1):p.sub[2],],
                           start=0.5*mcmc.sample)$p.y.predictive.samples)

  for(l in 2:(length(p.sub)-1)){
     m.1.pred <- spPredict(m.1, pred.covars=X.pred[(p.sub[l]+1):p.sub[l+1],],
                           pred.coords=coords.pred[(p.sub[l]+1):p.sub[l+1],],
                           start=0.5*mcmc.sample)$p.y.predictive.samples
     Y.pred <- cbind(Y.pred,t(m.1.pred))
  }
  ## MCMC samples for both parameter estimates and prediction after burn-in
  hh <- list(subAtom,Y.pred)
  names(hh) <- c("atoms","predictions") 
  return(hh)
}

####### Parallelization ########
## Number of clusters for parallel implementation
cl<-makeCluster(n.core)  
registerDoParallel(cl)

## Start time
strt<-Sys.time()
## Parallelized subset computation of GP in different cores
obj <- foreach(i=1:n.core) %dopar% partitioned_GP(i)  
## Total time for parallelized inference
final.time <- Sys.time()-strt  
stopCluster(cl)

######## combine ###############

subAtomList <- list()
for(i in 1:length(n.part)){
  ## MCMC samples to run Weiszfeld algorithm
  subAtomList[[i]] <- obj[[i]]$atoms[seq(1,mcmc.sample-n.burn,10),]
}
beg.rec <- Sys.time()
## Combination using Weiszfeld algorithm
medPosterior <- findWeiszfeldMedian(subAtomList, sigma = 0.1, maxit = 100, tol = 1e-5) 
## Posterior wts. using Weiszfeld's algorithm
wts <- medPosterior$weiszfeldWts
## Time for combining subset posteriors
recomb.time <- Sys.time()-beg.rec 

######## prediction############## 
Y.med <- numeric()
Y.lower_qnt <- numeric()
Y.upper_qnt <- numeric()
beg <- Sys.time()
## Number of posterior predictive samples
n.pred.sample <- nrow(obj[[1]]$predictions)

for(j in 1:ncol(obj[[1]]$predictions)){
   ## Posterior predictive samples for j th predicted point from first subsample
  a.sub.i<- obj[[1]]$predictions[,j] 
  a.wts <- wts[1]*rep(1/n.pred.sample,n.pred.sample)
  for(i in 2:n.split){
     ## Posterior predictive samples for j th predicted point from i th subsample
     a.sub.i <- c(a.sub.i,obj[[i]]$predictions[,j])  
     ## Corresponding empirical weights
     a.wts <- c(a.wts, wts[i]*rep(1/n.pred.sample,n.pred.sample))
  }
  new.atom <- sort(a.sub.i,decreasing=F,index.return=T)$ix
  ## MCMC samples for the meta posterior
  atoms <- a.sub.i[new.atom]
  ## Corresponding weights
  prob <- a.wts[new.atom]
  id11 <- min(which(cumsum(prob)>=0.025))
  id12 <- min(which(cumsum(prob)>=0.975))
  id13 <- min(which(cumsum(prob)>=0.5))
  ## Posterior 2.5% quantile for j th predicted point
  Y.lower_qnt[j] <- atoms[id11]  
  ## Posterior 97.5% quantile for j th predicted point         
  Y.upper_qnt[j] <- atoms[id12]  
  ## Posterior 50% quantile for j th predicted point         
  Y.med[j] <- atoms[id13]               
}

## Time for calculating quantiles
quantile.calculation.time <- Sys.time()-beg  

## Total time for subset running, combining subset posteriors and calculating quantiles
Time <- final.time+recomb.time+quantile.calculation.time   

## Save Results
res.df <- data.frame(Y.med=Y.med,Y.lower_qnt=Y.lower_qnt,Y.upper_qnt=Y.upper_qnt)
write.table(x=res.df,file="./MetaKrigingSatResults.txt",quote=FALSE,row.names=FALSE)
save(file="MetaSatTime.RData",list=c("Time"))
