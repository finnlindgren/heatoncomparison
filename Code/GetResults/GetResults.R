####################
## Load Libraries ##
####################
library(LatticeKrig)
library(R.matlab)
library(xtable)
library(gridExtra)
library(ggplot2)

#######################################
## Load Simulated and Satellite Data ##
#######################################
load('../../Data/AllSatelliteTemps.RData')
load('../../Data/AllSimulatedTemps.RData')
Lon <- matrix(all.sat.temps$Lon,nrow=500)
Lat <- matrix(all.sat.temps$Lat,nrow=500)
Sat.true <- matrix(all.sat.temps$TrueTemp,nrow=500)
Sim.true <- matrix(all.sim.data$TrueTemp,nrow=500)
sim.miss <- which(is.na(all.sim.data$MaskTemp))
sat.miss <- which(is.na(all.sat.temps$MaskTemp))
sim.preds <- 1:length(sim.miss)
sat.preds <- 1:length(sat.miss)

####################################################
## Should I include the long-distance prediction? ##
####################################################
kp.long <- TRUE
if(kp.long==FALSE){
  ld.locs <- which(Lat>36.49353 & Lon > -94.0070 & Lon< -91.83842)
  sat.preds <- sat.preds[!sat.miss%in%ld.locs]
  sim.preds <- sim.preds[!sim.miss%in%ld.locs]
  sim.miss <- sim.miss[!sim.miss%in%ld.locs]
  sat.miss <- sat.miss[!sat.miss%in%ld.locs]
}

######################
## Get NNGP Results ##
######################
nngp.sim <- read.csv('../NNGP/nngp-resp-sim-pred.csv')
nngp.sat <- read.csv('../NNGP/nngp-resp-real-pred.csv')
names(nngp.sim) <- c("Pred","Lower","Upper")
names(nngp.sat) <- c("Pred","Lower","Upper")
nngp.sat$sd <- (nngp.sat$Upper-nngp.sat$Lower)/(2*qnorm(0.975))
nngp.sim$sd <- (nngp.sim$Upper-nngp.sim$Lower)/(2*qnorm(0.975))
nngp.sim.conj <- read.csv('../NNGP/new-nngp-results/updated-nngp-conj-sim-pred.csv')
nngp.sat.conj <- read.csv('../NNGP/new-nngp-results/updated-nngp-conj-real-pred.csv')
names(nngp.sim.conj) <- c("Pred","Lower","Upper")
names(nngp.sat.conj) <- c("Pred","Lower","Upper")
nngp.sat.conj$sd <- (nngp.sat.conj$Upper-nngp.sat.conj$Lower)/(2*qnorm(0.975))
nngp.sim.conj$sd <- (nngp.sim.conj$Upper-nngp.sim.conj$Lower)/(2*qnorm(0.975))


##########################
## Get Tapering Results ##
##########################
load('../Tapering/SatelliteTemps_tapering_0.2.RData')
taper.sat <- list(pred=Ypred,unc=preduncertainty,timing=timing)
load('../Tapering/SimulatedTemps_tapering_0.25.RData')
taper.sim <- list(pred=Ypred,unc=preduncertainty,timing=timing)
rm(list=c("delta","predCIlower","scholCall","res","scholCobs",
          "shobs","shpred","shpredobs","verbose","Yconditional",
          "Ypred","preduncertainty"))

#########################
## Get Gapfill Results ##
#########################
load('../GapFill/SatelliteTemps_gapfill.RData')
gapfill.sat <- data.frame(pred=c(prediction),low=c(ciLo),up=c(ciUp))
gapfill.sat$pse <- (gapfill.sat$up-gapfill.sat$low)/(2*qnorm(0.975))
load('../GapFill/SimulatedTemps_gapfill.RData')
gapfill.sim <- data.frame(pred=c(prediction),low=c(ciLo),up=c(ciUp))
gapfill.sim$pse <- (gapfill.sim$up-gapfill.sim$low)/(2*qnorm(0.975))

###################################
## Get Spatial Partition Results ##
###################################
load("../Partition/SatelliteResults.RData")
partition.sat <- sat.results
partition.sat.tottime <- tot.time
partition.sat.liketime <- comp.time
partition.sat.predtime <- pred.time
load("../Partition/SimulationResults.RData")
partition.sim <- sim.results
partition.sim.tottime <- tot.time
partition.sim.liketime <- comp.time
partition.sim.predtime <- pred.time
rm(list=c("sat.results","sim.results","tot.time","comp.time","pred.time"))

#####################
## Get FRK Results ##
#####################
load("../FRK/FRKSatResults.RData")
load("../FRK/FRKSimResults.RData")

#############################
## Get LatticeKrig Results ##
#############################
load('../LatticeKrig/SatelliteTempsNC40nlevel4finalResults.rda')
LK.sat <- finalResults
load('../LatticeKrig/simulatedNC30nlevel4finalResults.rda')
LK.sim <- finalResults

#######################
## Get SPDEs Results ##
#######################
load('../SPDEs/SatTempsFitted.RData')
SPDE.sat <- reconstruction3
load('../SPDEs/SimTempsFitted.RData')
SPDE.sim <- reconstruction3
rm(list=c("reconstruction3","reconstruction1","reconstruction2","indata"))

##############################
## Get Meta Kriging Results ##
##############################
# meta.sat <- read.table("../MetaKriging/Meta_Kriging_RealData_30_subsets.txt") #UCSC Results
meta.sat <- read.table("../MetaKriging/MetaKrigingSatResults.txt",header=TRUE)
meta.sat$SD <- (meta.sat$Y.upper_qnt-meta.sat$Y.lower_qnt)/(2*qnorm(0.975))
#meta.sim <- read.table("../MetaKriging/Meta_Kriging_Simulation_30_subsets.txt") #UCSC Results
meta.sim <- read.table("../MetaKriging/MetaKrigingSimResults.txt",header=TRUE)
meta.sim$SD <- (meta.sim$Y.upper_qnt-meta.sim$Y.lower_qnt)/(2*qnorm(0.975))

#####################
## Get MRA Results ##
#####################
mra.sat <- readMat("../MRA/ResultsMrasatellite.mat")
mra.sat$predMat <- matrix(0,nrow=nrow(Lon),ncol=ncol(Lon))
mra.sat$sdMat <- mra.sat$predMat
for(rw in 1:nrow(Lon)){ 
  ##prediction locations from MRA are in different order so need to reorder
  the.lon <- unique(Lon[rw,])
  mra.vals <- mra.sat$predMean[mra.sat$predloc[,1]==the.lon]
  mra.sd <- sqrt(mra.sat$predVariance[mra.sat$predloc[,1]==the.lon])
  mra.latvals <- mra.sat$predloc[mra.sat$predloc[,1]==the.lon,2]
  mra.vals <- mra.vals[order(mra.latvals,decreasing=TRUE)]
  mra.sd <- mra.sd[order(mra.latvals,decreasing=TRUE)]
  mra.sat$predMat[rw,] <- mra.vals
  mra.sat$sdMat[rw,] <- mra.sd
}
mra.sim <- readMat("../MRA/ResultsMrasimulated.mat")
mra.sim$predMat <- matrix(0,nrow=nrow(Lon),ncol=ncol(Lon))
mra.sim$sdMat <- mra.sim$predMat
for(rw in 1:nrow(Lon)){ 
  ##prediction locations from MRA are in different order so need to reorder
  the.lon <- unique(Lon[rw,])
  mra.vals <- mra.sim$predMean[mra.sim$predloc[,1]==the.lon]
  mra.sd <- sqrt(mra.sim$predVariance[mra.sim$predloc[,1]==the.lon])
  mra.latvals <- mra.sim$predloc[mra.sim$predloc[,1]==the.lon,2]
  mra.vals <- mra.vals[order(mra.latvals,decreasing=TRUE)]
  mra.sd <- mra.sd[order(mra.latvals,decreasing=TRUE)]
  mra.sim$predMat[rw,] <- mra.vals
  mra.sim$sdMat[rw,] <- mra.sd
}

############################
## Get Pred Proc. Results ##
############################
sat.pp <- read.csv("../PredProc/pp-real-pred.csv")
sat.pp$sd <- (sat.pp[,3]-sat.pp[,2])/(2*qnorm(0.975))
sim.pp <- read.csv("../PredProc/pp-sim-pred.csv")
sim.pp$sd <- (sim.pp[,3]-sim.pp[,2])/(2*qnorm(0.975))

######################
## Get LAGP Results ##
######################
load("../LocalApproxGP/SatLAGP.RData")
sat.lagp <- list(mean=pmean,lower=plower,upper=pupper)
sat.lagp$sd <- (sat.lagp$upper-sat.lagp$lower)/(2*qnorm(0.975))
load("../LocalApproxGP/SimLAGP.RData")
sim.lagp <- list(mean=pmean,lower=plower,upper=pupper)
sim.lagp$sd <- (sim.lagp$upper-sim.lagp$lower)/(2*qnorm(0.975))


##############################
## Plot Each For Comparison ##
##############################

## Setup to Plot Satellite Predictions
preds.NNGP <- Sat.true
preds.NNGP[sat.miss] <- nngp.sat$Pred
preds.NNGP.conj <- Sat.true
preds.NNGP.conj[sat.miss] <- nngp.sat.conj$Pred
preds.part <- matrix(partition.sat$pred,nrow=500,ncol=300)
preds.FRK <- Sat.true
preds.FRK[sat.miss] <- FRKreal$results$mu
preds.LK <- Sat.true
preds.LK[sat.miss] <- LK.sat$yHat
preds.gapfill <- Sat.true
preds.gapfill[sat.miss] <- gapfill.sat$pred[sat.miss]
preds.taper <- Sat.true
preds.taper[sat.miss] <- taper.sat$pred
preds.SPDE <- Sat.true
preds.SPDE[sat.miss] <- SPDE.sat$Temp[sat.miss]
preds.meta <- Sat.true
preds.meta[sat.miss]  <- meta.sat$Y.med
preds.mra <- Sat.true
preds.mra[sat.miss] <- mra.sat$predMat[sat.miss]
preds.pp <- Sat.true
preds.pp[sat.miss] <- sat.pp[,1]
preds.lagp <- Sat.true
preds.lagp[sat.miss] <- sat.lagp$mean
zlimits <- range(c(preds.NNGP,preds.FRK,preds.part,preds.LK,
                   preds.gapfill,preds.taper,preds.SPDE,preds.meta,
                   preds.mra,preds.pp,preds.lagp,Sat.true),na.rm=TRUE)

## Plot em
preds.df <- data.frame(Lon=c(Lon),Lat=c(Lat),Pred=c(preds.FRK))
FRK.plot <- ggplot(preds.df) + geom_raster(aes(Lon,Lat,fill=Pred)) +
  theme_bw() + scale_fill_distiller(palette="Spectral",name="",limits=zlimits)+ ggtitle(label="(a) FRK")
preds.df <- data.frame(Lon=c(Lon),Lat=c(Lat),Pred=c(preds.gapfill))
gapfill.plot <- ggplot(preds.df) + geom_raster(aes(Lon,Lat,fill=Pred)) +
  theme_bw() + scale_fill_distiller(palette="Spectral",name="",limits=zlimits)+ ggtitle(label="(b) Gapfill")
preds.df <- data.frame(Lon=c(Lon),Lat=c(Lat),Pred=c(preds.LK))
LK.plot <- ggplot(preds.df) + geom_raster(aes(Lon,Lat,fill=Pred)) +
  theme_bw() + scale_fill_distiller(palette="Spectral",name="",limits=zlimits)+ ggtitle(label="(c) Lattice Krig")
preds.df <- data.frame(Lon=c(Lon),Lat=c(Lat),Pred=c(preds.lagp))
lagp.plot <- ggplot(preds.df) + geom_raster(aes(Lon,Lat,fill=Pred)) +
  theme_bw() + scale_fill_distiller(palette="Spectral",name="",limits=zlimits)+ ggtitle(label="(d) LAGP")
preds.df <- data.frame(Lon=c(Lon),Lat=c(Lat),Pred=c(preds.meta))
meta.plot <- ggplot(preds.df) + geom_raster(aes(Lon,Lat,fill=Pred)) +
  theme_bw() + scale_fill_distiller(palette="Spectral",name="",limits=zlimits)+ ggtitle(label="(e) Metakriging")
preds.df <- data.frame(Lon=c(Lon),Lat=c(Lat),Pred=c(preds.mra))
mra.plot <- ggplot(preds.df) + geom_raster(aes(Lon,Lat,fill=Pred)) +
  theme_bw() + scale_fill_distiller(palette="Spectral",name="",limits=zlimits)+ ggtitle(label="(f) MRA")
preds.df <- data.frame(Lon=c(Lon),Lat=c(Lat),Pred=c(preds.NNGP.conj))
NNGP.conj.plot <- ggplot(preds.df) + geom_raster(aes(Lon,Lat,fill=Pred)) +
  theme_bw() + scale_fill_distiller(palette="Spectral",name="",limits=zlimits)+ ggtitle(label="(g) NNGP-Conjugate")
preds.df <- data.frame(Lon=c(Lon),Lat=c(Lat),Pred=c(preds.NNGP))
NNGP.plot <- ggplot(preds.df) + geom_raster(aes(Lon,Lat,fill=Pred)) +
  theme_bw() + scale_fill_distiller(palette="Spectral",name="",limits=zlimits)+ ggtitle(label="(h) NNGP-Response")
preds.df <- data.frame(Lon=c(Lon),Lat=c(Lat),Pred=c(preds.part))
Partition.plot <- ggplot(preds.df) + geom_raster(aes(Lon,Lat,fill=Pred)) +
  theme_bw() + scale_fill_distiller(palette="Spectral",name="",limits=zlimits)+ ggtitle(label="(i) Partitioning")
preds.df <- data.frame(Lon=c(Lon),Lat=c(Lat),Pred=c(preds.pp))
pp.plot <- ggplot(preds.df) + geom_raster(aes(Lon,Lat,fill=Pred)) +
  theme_bw() + scale_fill_distiller(palette="Spectral",name="",limits=zlimits) + ggtitle(label="(j) Pred. Proc.")
preds.df <- data.frame(Lon=c(Lon),Lat=c(Lat),Pred=c(preds.SPDE))
SPDE.plot <- ggplot(preds.df) + geom_raster(aes(Lon,Lat,fill=Pred)) +
  theme_bw() + scale_fill_distiller(palette="Spectral",name="",limits=zlimits)+ ggtitle(label="(k) SPDEs")
preds.df <- data.frame(Lon=c(Lon),Lat=c(Lat),Pred=c(preds.taper))
taper.plot <- ggplot(preds.df) + geom_raster(aes(Lon,Lat,fill=Pred)) +
  theme_bw() + scale_fill_distiller(palette="Spectral",name="",limits=zlimits)+ ggtitle(label="(l) Tapering")
preds.df <- data.frame(Lon=c(Lon),Lat=c(Lat),Pred=c(Sat.true))
true.plot <- ggplot(preds.df) + geom_raster(aes(Lon,Lat,fill=Pred)) +
  theme_bw() + scale_fill_distiller(palette="Spectral",name="",limits=zlimits)+ ggtitle(label="(m) Truth")
gridExtra::grid.arrange(FRK.plot,gapfill.plot,LK.plot,lagp.plot,meta.plot,
                        mra.plot,NNGP.conj.plot,NNGP.plot,Partition.plot,pp.plot,
                        SPDE.plot,taper.plot,ncol=2)

## Plot Simulated Predictions
preds.NNGP <- Sim.true
preds.NNGP[sim.miss] <- nngp.sim$Pred
preds.NNGP.conj <- Sim.true
preds.NNGP.conj[sim.miss] <- nngp.sim.conj$Pred
preds.part <- matrix(partition.sim$pred,nrow=500,ncol=300)
preds.FRK <- Sim.true
preds.FRK[sim.miss] <- FRKsim$results$mu
preds.LK <- Sim.true
preds.LK[sim.miss] <- LK.sim$yHat
preds.gapfill <- Sim.true
preds.gapfill[sim.miss] <- gapfill.sim$pred[sim.miss]
preds.taper <- Sim.true
preds.taper[sim.miss] <- taper.sim$pred
preds.SPDE <- Sim.true
preds.SPDE[sim.miss] <- SPDE.sim$Temp[sim.miss]
preds.meta <- Sim.true
preds.meta[sim.miss] <- meta.sim$Y.med
preds.mra <- Sim.true
preds.mra[sim.miss] <- mra.sim$predMat[sim.miss]
preds.pp <- Sim.true
preds.pp[sim.miss] <- sim.pp[,1]
preds.lagp <- Sim.true
preds.lagp[sim.miss] <- sim.lagp$mean
zlimits <- range(c(preds.NNGP,preds.FRK,preds.part,preds.LK,preds.NNGP.conj,
                   preds.gapfill,preds.taper,preds.SPDE,preds.meta,
                   preds.mra,preds.pp,preds.lagp,Sim.true))

## Plot em
preds.df <- data.frame(Lon=c(Lon),Lat=c(Lat),Pred=c(preds.FRK))
FRK.plot <- ggplot(preds.df) + geom_raster(aes(Lon,Lat,fill=Pred)) +
  theme_bw() + scale_fill_distiller(palette="Spectral",name="",limits=zlimits)+ ggtitle(label="(a) FRK")
preds.df <- data.frame(Lon=c(Lon),Lat=c(Lat),Pred=c(preds.gapfill))
gapfill.plot <- ggplot(preds.df) + geom_raster(aes(Lon,Lat,fill=Pred)) +
  theme_bw() + scale_fill_distiller(palette="Spectral",name="",limits=zlimits)+ ggtitle(label="(b) Gapfill")
preds.df <- data.frame(Lon=c(Lon),Lat=c(Lat),Pred=c(preds.LK))
LK.plot <- ggplot(preds.df) + geom_raster(aes(Lon,Lat,fill=Pred)) +
  theme_bw() + scale_fill_distiller(palette="Spectral",name="",limits=zlimits)+ ggtitle(label="(c) Lattice Krig")
preds.df <- data.frame(Lon=c(Lon),Lat=c(Lat),Pred=c(preds.lagp))
lagp.plot <- ggplot(preds.df) + geom_raster(aes(Lon,Lat,fill=Pred)) +
  theme_bw() + scale_fill_distiller(palette="Spectral",name="",limits=zlimits)+ ggtitle(label="(d) LAGP")
preds.df <- data.frame(Lon=c(Lon),Lat=c(Lat),Pred=c(preds.meta))
meta.plot <- ggplot(preds.df) + geom_raster(aes(Lon,Lat,fill=Pred)) +
  theme_bw() + scale_fill_distiller(palette="Spectral",name="",limits=zlimits)+ ggtitle(label="(e) Metakriging")
preds.df <- data.frame(Lon=c(Lon),Lat=c(Lat),Pred=c(preds.mra))
mra.plot <- ggplot(preds.df) + geom_raster(aes(Lon,Lat,fill=Pred)) +
  theme_bw() + scale_fill_distiller(palette="Spectral",name="",limits=zlimits)+ ggtitle(label="(f) MRA")
preds.df <- data.frame(Lon=c(Lon),Lat=c(Lat),Pred=c(preds.NNGP.conj))
NNGP.conj.plot <- ggplot(preds.df) + geom_raster(aes(Lon,Lat,fill=Pred)) +
  theme_bw() + scale_fill_distiller(palette="Spectral",name="",limits=zlimits)+ ggtitle(label="(g) NNGP-Conjugate")
preds.df <- data.frame(Lon=c(Lon),Lat=c(Lat),Pred=c(preds.NNGP))
NNGP.plot <- ggplot(preds.df) + geom_raster(aes(Lon,Lat,fill=Pred)) +
  theme_bw() + scale_fill_distiller(palette="Spectral",name="",limits=zlimits)+ ggtitle(label="(h) NNGP-Response")
preds.df <- data.frame(Lon=c(Lon),Lat=c(Lat),Pred=c(preds.part))
Partition.plot <- ggplot(preds.df) + geom_raster(aes(Lon,Lat,fill=Pred)) +
  theme_bw() + scale_fill_distiller(palette="Spectral",name="",limits=zlimits)+ ggtitle(label="(i) Partitioning")
preds.df <- data.frame(Lon=c(Lon),Lat=c(Lat),Pred=c(preds.pp))
pp.plot <- ggplot(preds.df) + geom_raster(aes(Lon,Lat,fill=Pred)) +
  theme_bw() + scale_fill_distiller(palette="Spectral",name="",limits=zlimits) + ggtitle(label="(j) Pred. Proc.")
preds.df <- data.frame(Lon=c(Lon),Lat=c(Lat),Pred=c(preds.SPDE))
SPDE.plot <- ggplot(preds.df) + geom_raster(aes(Lon,Lat,fill=Pred)) +
  theme_bw() + scale_fill_distiller(palette="Spectral",name="",limits=zlimits)+ ggtitle(label="(k) SPDEs")
preds.df <- data.frame(Lon=c(Lon),Lat=c(Lat),Pred=c(preds.taper))
taper.plot <- ggplot(preds.df) + geom_raster(aes(Lon,Lat,fill=Pred)) +
  theme_bw() + scale_fill_distiller(palette="Spectral",name="",limits=zlimits)+ ggtitle(label="(l) Tapering")
preds.df <- data.frame(Lon=c(Lon),Lat=c(Lat),Pred=c(Sim.true))
true.plot <- ggplot(preds.df) + geom_raster(aes(Lon,Lat,fill=Pred)) +
  theme_bw() + scale_fill_distiller(palette="Spectral",name="",limits=zlimits)+ ggtitle(label="(m) Truth")
gridExtra::grid.arrange(FRK.plot,gapfill.plot,LK.plot,lagp.plot,meta.plot,
                        mra.plot,NNGP.conj.plot,NNGP.plot,Partition.plot,pp.plot,
                        SPDE.plot,taper.plot,ncol=2)

##########################################
## Numerical Results for Satellite Data ##
##########################################
satres.df <- data.frame(Method=c("NNGP","NNGPc","Taper","GapFill","Partition","FRK","LK","SPDE","Meta","MRA","PP","LAGP"))

## MAE
satres.df$MAE <- c(mean(abs(nngp.sat$Pred[sat.preds]-Sat.true[sat.miss]),na.rm=TRUE),
                   mean(abs(nngp.sat.conj$Pred[sat.preds]-Sat.true[sat.miss]),na.rm=TRUE),
                    mean(abs(taper.sat$pred[sat.preds]-Sat.true[sat.miss]),na.rm=TRUE),
                    mean(abs(gapfill.sat$pred[sat.miss]-Sat.true[sat.miss]),na.rm=TRUE),
                    mean(abs(partition.sat$pred[sat.miss]-Sat.true[sat.miss]),na.rm=TRUE),
                    mean(abs(FRKreal$results$mu[sat.preds]-Sat.true[sat.miss]),na.rm=TRUE),
                    mean(abs(LK.sat$yHat[sat.preds]-Sat.true[sat.miss]),na.rm=TRUE),
                    mean(abs(SPDE.sat$Temp[sat.miss]-Sat.true[sat.miss]),na.rm=TRUE),
                   mean(abs(meta.sat$Y.med[sat.preds]-Sat.true[sat.miss]),na.rm=TRUE),
                   mean(abs(mra.sat$predMat[sat.miss]-Sat.true[sat.miss]),na.rm=TRUE),
                   mean(abs(sat.pp[sat.preds,1]-Sat.true[sat.miss]),na.rm=TRUE),
                   mean(abs(sat.lagp$mean[sat.preds]-Sat.true[sat.miss]),na.rm=TRUE)
                   )

## RMSE
satres.df$RMSE <- c(sqrt(mean((nngp.sat$Pred[sat.preds]-Sat.true[sat.miss])^2,na.rm=TRUE)),
                    sqrt(mean((nngp.sat.conj$Pred[sat.preds]-Sat.true[sat.miss])^2,na.rm=TRUE)),
                    sqrt(mean((taper.sat$pred[sat.preds]-Sat.true[sat.miss])^2,na.rm=TRUE)),
                    sqrt(mean((gapfill.sat$pred[sat.miss]-Sat.true[sat.miss])^2,na.rm=TRUE)),
                    sqrt(mean((partition.sat$pred[sat.miss]-Sat.true[sat.miss])^2,na.rm=TRUE)),
                    sqrt(mean((FRKreal$results$mu[sat.preds]-Sat.true[sat.miss])^2,na.rm=TRUE)),
                    sqrt(mean((LK.sat$yHat[sat.preds]-Sat.true[sat.miss])^2,na.rm=TRUE)),
                    sqrt(mean((SPDE.sat$Temp[sat.miss]-Sat.true[sat.miss])^2,na.rm=TRUE)),
                    sqrt(mean((meta.sat$Y.med[sat.preds]-Sat.true[sat.miss])^2,na.rm=TRUE)),
                    sqrt(mean((mra.sat$predMat[sat.miss]-Sat.true[sat.miss])^2,na.rm=TRUE)),
                    sqrt(mean((sat.pp[sat.preds,1]-Sat.true[sat.miss])^2,na.rm=TRUE)),
                    sqrt(mean((sat.lagp$mean[sat.preds]-Sat.true[sat.miss])^2,na.rm=TRUE))
                    )

## CRPS
crps <- function(predlist,trueobs) {
  z <- as.numeric((trueobs - predlist$mean) / predlist$sd)
  scores <- predlist$sd * (z *(2 * pnorm(z, 0, 1) - 1) +
                             2 * dnorm(z, 0, 1) - 1/sqrt(pi))
  return(scores)
}
satres.df$CRPS <- c(mean(crps(list(mean=nngp.sat$Pred[sat.preds],sd=nngp.sat$sd[sat.preds]),Sat.true[sat.miss]),na.rm=TRUE),
                    mean(crps(list(mean=nngp.sat.conj$Pred[sat.preds],sd=nngp.sat.conj$sd[sat.preds]),Sat.true[sat.miss]),na.rm=TRUE),
                    mean(crps(list(mean=taper.sat$pred[sat.preds],sd=taper.sat$unc[sat.preds]),Sat.true[sat.miss]),na.rm=TRUE),
                    mean(crps(list(mean=gapfill.sat$pred[sat.miss],sd=gapfill.sat$pse[sat.miss]),Sat.true[sat.miss]),na.rm=TRUE),
                    mean(crps(list(mean=partition.sat$pred[sat.miss],sd=partition.sat$pred.se[sat.miss]),Sat.true[sat.miss]),na.rm=TRUE),
                    mean(crps(list(mean=FRKreal$results$mu[sat.preds],sd=FRKreal$results$sd[sat.preds]),Sat.true[sat.miss]),na.rm=TRUE),
                    mean(crps(list(mean=LK.sat$yHat[sat.preds],sd=LK.sat$standError[sat.preds]),Sat.true[sat.miss]),na.rm=TRUE),
                    mean(crps(list(mean=SPDE.sat$Temp[sat.miss],sd=SPDE.sat$SD[sat.miss]),Sat.true[sat.miss]),na.rm=TRUE),
                    mean(crps(list(mean=meta.sat$Y.med[sat.preds],sd=meta.sat$SD[sat.preds]),Sat.true[sat.miss]),na.rm=TRUE),
                    mean(crps(list(mean=mra.sat$predMat[sat.miss],sd=mra.sat$sdMat[sat.miss]),Sat.true[sat.miss]),na.rm=TRUE),
                    mean(crps(list(mean=sat.pp[sat.preds,1],sd=sat.pp$sd[sat.preds]),Sat.true[sat.miss]),na.rm=TRUE),
                    mean(crps(list(mean=sat.lagp$mean[sat.preds],sd=sat.lagp$sd[sat.preds]),Sat.true[sat.miss]),na.rm=TRUE)
)

## Interval Score
intscore <- function(x, y, alpha=0.05) {
  hw <- -qnorm(alpha/2) * x$sd
  scores <- 2 * hw + (2/alpha) * (((x$mean - hw) - y) * (y < x$mean - hw) +
                                    (y - (x$mean + hw)) * (y > x$mean + hw))
  return(scores)
}
satres.df$intscore <- c(mean(intscore(list(mean=nngp.sat$Pred[sat.preds],sd=nngp.sat$sd[sat.preds]),Sat.true[sat.miss]),na.rm=TRUE),
                        mean(intscore(list(mean=nngp.sat.conj$Pred[sat.preds],sd=nngp.sat.conj$sd[sat.preds]),Sat.true[sat.miss]),na.rm=TRUE),
                    mean(intscore(list(mean=taper.sat$pred[sat.preds],sd=taper.sat$unc[sat.preds]),Sat.true[sat.miss]),na.rm=TRUE),
                    mean(intscore(list(mean=gapfill.sat$pred[sat.miss],sd=gapfill.sat$pse[sat.miss]),Sat.true[sat.miss]),na.rm=TRUE),
                    mean(intscore(list(mean=partition.sat$pred[sat.miss],sd=partition.sat$pred.se[sat.miss]),Sat.true[sat.miss]),na.rm=TRUE),
                    mean(intscore(list(mean=FRKreal$results$mu[sat.preds],sd=FRKreal$results$sd[sat.preds]),Sat.true[sat.miss]),na.rm=TRUE),
                    mean(intscore(list(mean=LK.sat$yHat[sat.preds],sd=LK.sat$standError[sat.preds]),Sat.true[sat.miss]),na.rm=TRUE),
                    mean(intscore(list(mean=SPDE.sat$Temp[sat.miss],sd=SPDE.sat$SD[sat.miss]),Sat.true[sat.miss]),na.rm=TRUE),
                    mean(intscore(list(mean=meta.sat$Y.med[sat.preds],sd=meta.sat$SD[sat.preds]),Sat.true[sat.miss]),na.rm=TRUE),
                    mean(intscore(list(mean=mra.sat$predMat[sat.miss],sd=mra.sat$sdMat[sat.miss]),Sat.true[sat.miss]),na.rm=TRUE),
                    mean(intscore(list(mean=sat.pp[sat.preds,1],sd=sat.pp$sd[sat.preds]),Sat.true[sat.miss]),na.rm=TRUE),
                    mean(intscore(list(mean=sat.lagp$mean[sat.preds],sd=sat.lagp$sd[sat.preds]),Sat.true[sat.miss]),na.rm=TRUE)
)                   

## Coverage
cvg <- function(x, y, alpha=0.05) {
  hw <- -qnorm(alpha/2) * x$sd
  scores <- y >= (x$mean - hw) & y <= (x$mean + hw)
  return(scores)
}
satres.df$CVG <- c(mean(cvg(list(mean=nngp.sat$Pred[sat.preds],sd=nngp.sat$sd[sat.preds]),Sat.true[sat.miss]),na.rm=TRUE),
                   mean(cvg(list(mean=nngp.sat.conj$Pred[sat.preds],sd=nngp.sat.conj$sd[sat.preds]),Sat.true[sat.miss]),na.rm=TRUE),
                        mean(cvg(list(mean=taper.sat$pred[sat.preds],sd=taper.sat$unc[sat.preds]),Sat.true[sat.miss]),na.rm=TRUE),
                        mean(cvg(list(mean=gapfill.sat$pred[sat.miss],sd=gapfill.sat$pse[sat.miss]),Sat.true[sat.miss]),na.rm=TRUE),
                        mean(cvg(list(mean=partition.sat$pred[sat.miss],sd=1.1*partition.sat$pred.se[sat.miss]),Sat.true[sat.miss]),na.rm=TRUE),
                        mean(cvg(list(mean=FRKreal$results$mu[sat.preds],sd=FRKreal$results$sd[sat.preds]),Sat.true[sat.miss]),na.rm=TRUE),
                        mean(cvg(list(mean=LK.sat$yHat[sat.preds],sd=LK.sat$standError[sat.preds]),Sat.true[sat.miss]),na.rm=TRUE),
                        mean(cvg(list(mean=SPDE.sat$Temp[sat.miss],sd=SPDE.sat$SD[sat.miss]),Sat.true[sat.miss]),na.rm=TRUE),
                        mean(cvg(list(mean=meta.sat$Y.med[sat.preds],sd=meta.sat$SD[sat.preds]),Sat.true[sat.miss]),na.rm=TRUE),
                        mean(cvg(list(mean=mra.sat$predMat[sat.miss],sd=mra.sat$sdMat[sat.miss]),Sat.true[sat.miss]),na.rm=TRUE),
                        mean(cvg(list(mean=sat.pp[sat.preds,1],sd=sat.pp$sd[sat.preds]),Sat.true[sat.miss]),na.rm=TRUE),
                        mean(cvg(list(mean=sat.lagp$mean[sat.preds],sd=sat.lagp$sd[sat.preds]),Sat.true[sat.miss]),na.rm=TRUE)
)                   

## Rank Methods
method.rank <- satres.df[,-1]
method.rank$CVG <- abs(method.rank$CVG-0.95)
for(metric in 1:ncol(method.rank)){
  method.order <- order(method.rank[,metric])
  method.rank[method.order,metric] <- 1:nrow(method.rank)
}
method.rank <- rowMeans(method.rank)
satres.df$RunTime <- c(2571,2.06*60,7995.7,83.527,4798.6,134.146,1675.154,7219.541,173311,936.78,38428.51,136)/60
satres.df$MeanRank <- method.rank
satres.df <- satres.df[order(method.rank),]
satres.df <- satres.df[order(substr(satres.df$Method,1,2)),]
rownames(satres.df) <- as.character(1:nrow(satres.df))
print(satres.df)
xtable(satres.df,digits=2)

##########################################
## Numerical Results for Simulated Data ##
##########################################
simres.df <- data.frame(Method=c("NNGP","NNGPc","Taper","GapFill","Partition","FRK","LK","SPDE","Meta","MRA","PP","LAGP"))

## MAE
simres.df$MAE <- c(mean(abs(nngp.sim$Pred[sim.preds]-Sim.true[sim.miss]),na.rm=TRUE),
                   mean(abs(nngp.sim.conj$Pred[sim.preds]-Sim.true[sim.miss]),na.rm=TRUE),
                   mean(abs(taper.sim$pred[sim.preds]-Sim.true[sim.miss]),na.rm=TRUE),
                   mean(abs(gapfill.sim$pred[sim.miss]-Sim.true[sim.miss]),na.rm=TRUE),
                   mean(abs(partition.sim$pred[sim.miss]-Sim.true[sim.miss]),na.rm=TRUE),
                   mean(abs(FRKsim$results$mu[sim.preds]-Sim.true[sim.miss]),na.rm=TRUE),
                   mean(abs(LK.sim$yHat[sim.preds]-Sim.true[sim.miss]),na.rm=TRUE),
                   mean(abs(SPDE.sim$Temp[sat.miss]-Sim.true[sim.miss]),na.rm=TRUE),
                   mean(abs(meta.sim$Y.med[sim.preds]-Sim.true[sim.miss]),na.rm=TRUE),
                   mean(abs(mra.sim$predMat[sim.miss]-Sim.true[sim.miss]),na.rm=TRUE),
                   mean(abs(sim.pp[sim.preds,1]-Sim.true[sim.miss]),na.rm=TRUE),
                   mean(abs(sim.lagp$mean[sim.preds]-Sim.true[sim.miss]),na.rm=TRUE)
)

## RMSE
simres.df$RMSE <- c(sqrt(mean((nngp.sim$Pred[sim.preds]-Sim.true[sim.miss])^2,na.rm=TRUE)),
                    sqrt(mean((nngp.sim.conj$Pred[sim.preds]-Sim.true[sim.miss])^2,na.rm=TRUE)),
                    sqrt(mean((taper.sim$pred[sim.preds]-Sim.true[sim.miss])^2,na.rm=TRUE)),
                    sqrt(mean((gapfill.sim$pred[sim.miss]-Sim.true[sim.miss])^2,na.rm=TRUE)),
                    sqrt(mean((partition.sim$pred[sim.miss]-Sim.true[sim.miss])^2,na.rm=TRUE)),
                    sqrt(mean((FRKsim$results$mu[sim.preds]-Sim.true[sim.miss])^2,na.rm=TRUE)),
                    sqrt(mean((LK.sim$yHat[sim.preds]-Sim.true[sim.miss])^2,na.rm=TRUE)),
                    sqrt(mean((SPDE.sim$Temp[sim.miss]-Sim.true[sim.miss])^2,na.rm=TRUE)),
                    sqrt(mean((meta.sim$Y.med[sim.preds]-Sim.true[sim.miss])^2,na.rm=TRUE)),
                    sqrt(mean((mra.sim$predMat[sim.miss]-Sim.true[sim.miss])^2,na.rm=TRUE)),
                    sqrt(mean((sim.pp[sim.preds,1]-Sim.true[sim.miss])^2,na.rm=TRUE)),
                    sqrt(mean((sim.lagp$mean[sim.preds]-Sim.true[sim.miss])^2,na.rm=TRUE))
)

## CRPS
simres.df$CRPS <- c(mean(crps(list(mean=nngp.sim$Pred[sim.preds],sd=nngp.sim$sd[sim.preds]),Sim.true[sim.miss]),na.rm=TRUE),
                    mean(crps(list(mean=nngp.sim.conj$Pred[sim.preds],sd=nngp.sim.conj$sd[sim.preds]),Sim.true[sim.miss]),na.rm=TRUE),
                    mean(crps(list(mean=taper.sim$pred[sim.preds],sd=taper.sim$unc[sim.preds]),Sim.true[sim.miss]),na.rm=TRUE),
                    mean(crps(list(mean=gapfill.sim$pred[sim.miss],sd=gapfill.sim$pse[sim.miss]),Sim.true[sim.miss]),na.rm=TRUE),
                    mean(crps(list(mean=partition.sim$pred[sim.miss],sd=partition.sim$pred.se[sim.miss]),Sim.true[sim.miss]),na.rm=TRUE),
                    mean(crps(list(mean=FRKsim$results$mu[sim.preds],sd=FRKsim$results$sd[sim.preds]),Sim.true[sim.miss]),na.rm=TRUE),
                    mean(crps(list(mean=LK.sim$yHat[sim.preds],sd=LK.sim$standError[sim.preds]),Sim.true[sim.miss]),na.rm=TRUE),
                    mean(crps(list(mean=SPDE.sim$Temp[sim.miss],sd=SPDE.sat$SD[sim.miss]),Sim.true[sim.miss]),na.rm=TRUE),
                    mean(crps(list(mean=meta.sim$Y.med[sim.preds],sd=meta.sim$SD[sim.preds]),Sim.true[sim.miss]),na.rm=TRUE),
                    mean(crps(list(mean=mra.sim$predMat[sim.miss],sd=mra.sim$sdMat[sim.miss]),Sim.true[sim.miss]),na.rm=TRUE),
                    mean(crps(list(mean=sim.pp[sim.preds,1],sd=sim.pp$sd[sim.preds]),Sim.true[sim.miss]),na.rm=TRUE),
                    mean(crps(list(mean=sim.lagp$mean[sim.preds],sd=sim.lagp$sd[sim.preds]),Sim.true[sim.miss]),na.rm=TRUE)
)

## Interval Score
simres.df$intscore <- c(mean(intscore(list(mean=nngp.sim$Pred[sim.preds],sd=nngp.sim$sd[sim.preds]),Sim.true[sim.miss]),na.rm=TRUE),
                        mean(intscore(list(mean=nngp.sim.conj$Pred[sim.preds],sd=nngp.sim.conj$sd[sim.preds]),Sim.true[sim.miss]),na.rm=TRUE),
                        mean(intscore(list(mean=taper.sim$pred[sim.preds],sd=taper.sim$unc[sim.preds]),Sim.true[sim.miss]),na.rm=TRUE),
                        mean(intscore(list(mean=gapfill.sim$pred[sim.miss],sd=gapfill.sim$pse[sim.miss]),Sim.true[sim.miss]),na.rm=TRUE),
                        mean(intscore(list(mean=partition.sim$pred[sim.miss],sd=partition.sim$pred.se[sim.miss]),Sim.true[sim.miss]),na.rm=TRUE),
                        mean(intscore(list(mean=FRKsim$results$mu[sim.preds],sd=FRKsim$results$sd[sim.preds]),Sim.true[sim.miss]),na.rm=TRUE),
                        mean(intscore(list(mean=LK.sim$yHat[sim.preds],sd=LK.sim$standError[sim.preds]),Sim.true[sim.miss]),na.rm=TRUE),
                        mean(intscore(list(mean=SPDE.sim$Temp[sim.miss],sd=SPDE.sat$SD[sim.miss]),Sim.true[sim.miss]),na.rm=TRUE),
                        mean(intscore(list(mean=meta.sim$Y.med[sim.preds],sd=meta.sim$SD[sim.preds]),Sim.true[sim.miss]),na.rm=TRUE),
                        mean(intscore(list(mean=mra.sim$predMat[sim.miss],sd=mra.sim$sdMat[sim.miss]),Sim.true[sim.miss]),na.rm=TRUE),
                        mean(intscore(list(mean=sim.pp[sim.preds,1],sd=sim.pp$sd[sim.preds]),Sim.true[sim.miss]),na.rm=TRUE),
                        mean(intscore(list(mean=sim.lagp$mean[sim.preds],sd=sim.lagp$sd[sim.preds]),Sim.true[sim.miss]),na.rm=TRUE)
)

## Coverage
simres.df$CVG <- c(mean(cvg(list(mean=nngp.sim$Pred[sim.preds],sd=nngp.sim$sd[sim.preds]),Sim.true[sim.miss]),na.rm=TRUE),
                   mean(cvg(list(mean=nngp.sim.conj$Pred[sim.preds],sd=nngp.sim.conj$sd[sim.preds]),Sim.true[sim.miss]),na.rm=TRUE),
                   mean(cvg(list(mean=taper.sim$pred[sim.preds],sd=taper.sim$unc[sim.preds]),Sim.true[sim.miss]),na.rm=TRUE),
                   mean(cvg(list(mean=gapfill.sim$pred[sim.miss],sd=gapfill.sim$pse[sim.miss]),Sim.true[sim.miss]),na.rm=TRUE),
                   mean(cvg(list(mean=partition.sim$pred[sim.miss],sd=partition.sim$pred.se[sim.miss]),Sim.true[sim.miss]),na.rm=TRUE),
                   mean(cvg(list(mean=FRKsim$results$mu[sim.preds],sd=FRKsim$results$sd[sim.preds]),Sim.true[sim.miss]),na.rm=TRUE),
                   mean(cvg(list(mean=LK.sim$yHat[sim.preds],sd=LK.sim$standError[sim.preds]),Sim.true[sim.miss]),na.rm=TRUE),
                   mean(cvg(list(mean=SPDE.sim$Temp[sim.miss],sd=SPDE.sat$SD[sim.miss]),Sim.true[sim.miss]),na.rm=TRUE),
                   mean(cvg(list(mean=meta.sim$Y.med[sim.preds],sd=meta.sim$SD[sim.preds]),Sim.true[sim.miss]),na.rm=TRUE),
                   mean(cvg(list(mean=mra.sim$predMat[sim.miss],sd=mra.sim$sdMat[sim.miss]),Sim.true[sim.miss]),na.rm=TRUE),
                   mean(cvg(list(mean=sim.pp[sim.preds,1],sd=sim.pp$sd[sim.preds]),Sim.true[sim.miss]),na.rm=TRUE),
                   mean(cvg(list(mean=sim.lagp$mean[sim.preds],sd=sim.lagp$sd[sim.preds]),Sim.true[sim.miss]),na.rm=TRUE)
)

## Rank Methods
method.rank <- simres.df[,-1]
method.rank$CVG <- abs(method.rank$CVG-0.95)
for(metric in 1:ncol(method.rank)){
  method.order <- order(method.rank[,metric])
  method.rank[method.order,metric] <- 1:nrow(method.rank)
}
method.rank <- rowMeans(method.rank)
simres.df$MeanRank <- method.rank
simres.df$RunTime <- c(130.80,1.99*60,37.56,1534.743,173333.2,814.20,2703.34,4653.71,38353.51,8300.606,11301.38,137)/60
simres.df <- simres.df[order(method.rank),]
simres.df <- simres.df[order(substr(simres.df$Method,1,2)),]
rownames(simres.df) <- as.character(1:nrow(simres.df))
print(simres.df)
xtable(simres.df,digits=2)




