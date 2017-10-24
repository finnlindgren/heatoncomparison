#######################################################################
## Author: Matt Heaton
## Date: 10 May 2017
## Description: Code to fit partitioned spatial model
## Instructions: Please set folder to current directory and source it
## Requirements: You need to have the following packages installed
##               -LatticeKrig
##               -plyr
##               -parallel
#######################################################################


############################################### 
## Load Libraries and Source Clust Functions ##
###############################################
rm(list=ls())
library(LatticeKrig)
library(plyr)
library(parallel)
n.cores <- 55

######################
## Read in the data ##
######################
data.set <- "sim"
if(data.set=="sat"){
  read.path <- '../../Data/SatelliteTemps.RData'
  write.path <- './SatelliteResults.RData'
  load(read.path)
  spatdat <- sat.temps
  nrows <- 500
  rm(list="sat.temps")
} else if(data.set=="sim"){
  read.path <- '../../Data/SimulatedTemps.RData'
  write.path <- './SimulationResults.RData'
  load(read.path)
  spatdat <- sim.data
  nrows <- 500
  rm(list="sim.data")
} else if(data.set=="test"){
  read.path <- '../../Data/SmallTestData.RData'
  write.path <- './SmallTestResults.RData'
  load(read.path)
  spatdat <- code.test
  names(spatdat) <- c("Lon","Lat","Temp","TrueObs")
  nrows <- 100
} else {
  stop('data.set must be either sat, sim or test')
}

##############################################
## Define a few variables for ease up front ##
##############################################
lon <- matrix(spatdat$Lon,nrow=nrows)
lat <- matrix(spatdat$Lat,nrow=nrows)
u.lon <- unique(spatdat$Lon)
u.lat <- unique(spatdat$Lat)
minDist <- min(abs(u.lon[1]-u.lon[2]),abs(u.lat[1]-u.lat[2]))
maxDist <- sqrt((u.lon[1]-u.lon[length(u.lon)])^2+(u.lat[1]-u.lat[length(u.lat)])^2)
obs <- which(!is.na(spatdat$Temp))
N <- length(obs)

###################################
## Define Global Basis Functions ##
###################################
n.knots <- 4 ## n.knots^2 actually
xlim <- range(spatdat$Lon)
ylim <- range(spatdat$Lat)
xlocs <- seq(xlim[1],xlim[2],length=n.knots+1)[-1]
xlocs <- xlocs-0.5*mean(diff(xlocs))
ylocs <- seq(ylim[1],ylim[2],length=n.knots+1)[-1]
ylocs <- ylocs-0.5*mean(diff(ylocs))
knot.locs <- expand.grid(xlocs,ylocs)
B <- rdist(cbind(spatdat$Lon,spatdat$Lat),knot.locs)
rl <- 2*min(rdist(knot.locs)[lower.tri(rdist(knot.locs))])
B <- ((1-(B/rl)^2)^2)*(B<=rl)
# plot(spatdat$Lon,spatdat$Lat,pch=19,cex=.2)
# points(knot.locs[,1],knot.locs[,2],pch=19,col="red")
# image.plot(lon,lat,matrix(B[,2],nrow=nrows))

##############
## X-matrix ##
##############
if(data.set=="sat"){
  X <- cbind(1,spatdat$Lon,spatdat$Lat,B)
} else if(data.set=="sim") {
  X <- cbind(1,B)
} else {
  X <- cbind(1,B)
}
P <- ncol(X)

##########################################
## Cluster Residuals via Agg clustering ##
##########################################
n.grps <- 5
n.clust <- n.grps^2
xlocs <- seq(xlim[1],xlim[2],length=n.grps+1)[-1]
xlocs <- xlocs-0.5*mean(diff(xlocs))
ylocs <- seq(ylim[1],ylim[2],length=n.grps+1)[-1]
ylocs <- ylocs-0.5*mean(diff(ylocs))
part.centroids <- expand.grid(xlocs,ylocs)
D <- rdist(cbind(spatdat$Lon,spatdat$Lat),part.centroids)
clust <- apply(D,1,which.min)
spatdat$cluster <- clust
# plot(spatdat$Lon,spatdat$Lat,pch=19,cex=.2)
# points(part.centroids[,1],part.centroids[,2],pch=19,col="red")
# image.plot(lon,lat,matrix(spatdat$cluster,nrow=nrows))
range(table(spatdat$cluster[obs]))
summary(lm(Temp~X+as.factor(cluster),data=spatdat))$r.sq

######################################################
## Grid search for omega (nugget) and alpha (decay) ##
######################################################
n.om <- 20
n.alpha <- 10
om.grid <- seq(0,0.99,length=n.om)
alpha.grid <- 3/seq(minDist,maxDist,length=n.alpha)
oa.grid <- expand.grid(om.grid,alpha.grid)

############################################
## Function to Evaluate Likelihood pieces ##
############################################
setup.MLE <- function(xlist){
  
  ## Set up variables
  nc <- nrow(yc)
  D <- rdist(sc)
  Rc <- xlist$omega*exp(-xlist$alpha*D)+(1-xlist$omega)*diag(nc)
  if(xlist$omega==0){ #if omega=0 then no spatial correlation and Rc=I
    Rc.chol <- Rc
    Rc.inv <- Rc
  } else {
    Rc.chol <- chol(Rc)
    Rc.inv <- chol2inv(Rc.chol)
  }
  
  ## Calculate critical quantities
  log.det <- 2*sum(diag(Rc.chol))
  Xp.Rcinv.y <- t(Xc)%*%Rc.inv%*%yc
  Xp.Rcinv.X <- t(Xc)%*%Rc.inv%*%Xc
  yp.Rcinv.y <- t(yc)%*%Rc.inv%*%yc
  
  ## Add to Sum
  xlist$Xp.Rcinv.y.sum <- xlist$Xp.Rcinv.y.sum+Xp.Rcinv.y
  xlist$Xp.Rcinv.X.sum <- xlist$Xp.Rcinv.X.sum+Xp.Rcinv.X
  xlist$yp.Rcinv.y.sum <- xlist$yp.Rcinv.y.sum+yp.Rcinv.y
  xlist$logdet.sum <- xlist$logdet.sum + log.det
  
  ## Return new xlist
  return(xlist)
  
}

############################################
## Function to Evaluate MLEs and log-like ##
############################################
get.MLEs <- function(xlist){
  
  ## Calculate beta.hat
  beta.hat <- chol2inv(chol(xlist$Xp.Rcinv.X.sum))%*%xlist$Xp.Rcinv.y.sum
  
  ## Calculate sigma2.hat
  ss <- (xlist$yp.Rcinv.y.sum-2*t(beta.hat)%*%xlist$Xp.Rcinv.y.sum+
           t(beta.hat)%*%xlist$Xp.Rcinv.X.sum%*%beta.hat)
  sig2.hat <- ss/(N-P)
  
  ## Calculate log-like
  loglike <- -(N/2)*log(sig2.hat)-0.5*xlist$logdet.sum-0.5*(ss/sig2.hat)
  
  ## Return info
  return(list(beta.hat=beta.hat,sig2.hat=sig2.hat,omega.hat=xlist$omega,alpha.hat=xlist$alpha,loglike=loglike))
  
}

####################################
## Create lists w/unique om/alpha ##
####################################
oa.list <- vector("list",length=nrow(oa.grid))
for(i in 1:length(oa.list)){
  oa.list[[i]]$omega <- oa.grid[i,1]
  oa.list[[i]]$alpha <- oa.grid[i,2]
  oa.list[[i]]$Xp.Rcinv.y.sum <- matrix(0,P,1)
  oa.list[[i]]$Xp.Rcinv.X.sum <- matrix(0,P,P)
  oa.list[[i]]$yp.Rcinv.y.sum <- 0
  oa.list[[i]]$logdet.sum <- 0
}

###########################################
## Find the MLEs via parallel processing ##
###########################################

## Calculate critical quantities in each cluster
#chk.list <- oa.list[[34]]
comp.time <- system.time({
for(i in 1:n.clust){
  clust.obs <- which(spatdat$cluster==i & !is.na(spatdat$Temp))
  yc <- matrix(spatdat$Temp[clust.obs],ncol=1)
  Xc <- matrix(X[clust.obs,],ncol=P)
  sc <- cbind(spatdat$Lon[clust.obs],spatdat$Lat[clust.obs])
  oa.list <- mclapply(oa.list,setup.MLE,mc.cores=n.cores)
  #chk.list <- setup.MLE(chk.list)
  print(paste0(100*round(i/n.clust,2)," Percent Complete"))
}
})
tot.time <- comp.time[3]
comp.time <- comp.time[3]/(nrow(oa.grid))
#chk.MLEs <- get.MLEs(chk.list)
#chk.lm <- summary(lm(Temp~X-1,data=spatdat))

## Find MLEs and log.like
all.mles <- lapply(oa.list,get.MLEs)

## Select the model that is the MLE
the.mle <- which.max(sapply(all.mles,function(xlist){return(xlist$loglike)}))
fitted.model <- all.mles[[the.mle]]
# image.plot(lon,lat,matrix(X%*%fitted.model$beta.hat,nrow=nrows))
# image.plot(lon,lat,matrix(spatdat$Temp,nrow=nrows))

###################################################
## Function to calculate predictions in parallel ##
###################################################
get.preds <- function(xlist){
  
  ## Setup Variables
  sc.obs <- which(!is.na(xlist$yc))
  D <- rdist(xlist$sc)
  R <- fitted.model$omega.hat*exp(-fitted.model$alpha.hat*D)+
    (1-fitted.model$omega.hat)*diag(nrow(D))
  mn <- xlist$Xc%*%fitted.model$beta.hat
  
  ## Calculate Predictive Distribution
  R12.R22inv <- R[-sc.obs,sc.obs]%*%chol2inv(chol(R[sc.obs,sc.obs]))
  pred.mn <- mn[-sc.obs,]+R12.R22inv%*%(xlist$yc[sc.obs,]-mn[sc.obs,])
  pred.var <- diag(R[-sc.obs,-sc.obs]-R12.R22inv%*%R[sc.obs,-sc.obs])
  pred.var <- fitted.model$sig2.hat*pred.var
  
  ## Return predmn and variances
  all.pred.mns <- xlist$yc
  all.pred.mns[-sc.obs] <- pred.mn
  all.pred.vars <- matrix(0,nrow=nrow(xlist$yc),ncol=1)
  all.pred.vars[-sc.obs] <- pred.var
  return(data.frame(pred.mn=all.pred.mns,pred.sd=sqrt(all.pred.vars)))
  
}
  
#############################
## Predictions in Parallel ##
#############################
clustdata <- vector("list",length=n.clust)
for(i in 1:n.clust){
  clust.obs <- which(spatdat$cluster==i)
  clustdata[[i]]$yc <- matrix(spatdat$Temp[clust.obs],ncol=1)
  clustdata[[i]]$Xc <- X[clust.obs,]
  clustdata[[i]]$sc <- cbind(spatdat$Lon[clust.obs],spatdat$Lat[clust.obs])
}
pred.time <- system.time({
  preds <- mclapply(clustdata,get.preds,mc.cores=n.cores)
})
tot.time <- tot.time+pred.time[3]

#######################################
## Divide out Predictions to Dataset ##
#######################################
spatdat$pred <- spatdat$Temp
spatdat$pred.se <- 0
for(i in 1:n.clust){
  spatdat$pred[spatdat$cluster==i] <- preds[[i]]$pred.mn
  spatdat$pred.se[spatdat$cluster==i] <- preds[[i]]$pred.sd
}
spatdat$lowlim <- spatdat$pred-1.96*spatdat$pred.se
spatdat$uplim <- spatdat$pred+1.96*spatdat$pred.se
# image.plot(lon,lat,matrix(spatdat$pred,nrow=nrows))
# image.plot(lon,lat,matrix(spatdat$Temp,nrow=nrows))

###############################
## Check Predictive Accuracy ##
###############################
if(data.set=="sat"){
  load('../../Data/AllSatelliteTemps.RData')
  tst.set <- which(is.na(all.sat.temps$MaskTemp) & !is.na(all.sat.temps$TrueTemp))
  tst.obs <- all.sat.temps$TrueTemp[tst.set]
} else if(data.set=='sim'){
  load('../../Data/AllSimulatedTemps.RData')
  tst.set <- which(is.na(all.sim.data$MaskTemp))
  tst.obs <- all.sim.data$TrueTemp[tst.set]
} else {
  tst.set <- which(is.na(spatdat$Temp))
  tst.obs <- spatdat$TrueObs[tst.set]
}
ll <- spatdat$lowlim[tst.set]
pred <- spatdat$pred[tst.set]
ul <- spatdat$uplim[tst.set]
print(mean(ll<=tst.obs & ul>=tst.obs)) #CVG
print(mean(pred-tst.obs)) #Bias
print(sqrt(mean((pred-tst.obs)^2))) #RMSE

###########################
## Save your predictions ##
###########################
if(data.set=="sat"){
  sat.results <- spatdat
  save(file=write.path,list=c("sat.results","tot.time","comp.time","pred.time"))
} else if(data.set=="sim"){
  sim.results <- spatdat
  save(file=write.path,list=c("sim.results","tot.time","comp.time","pred.time"))
} else {
  test.results <- spatdat
  save(file=write.path,list=c("test.results","tot.time","comp.time","pred.time"))
}




  