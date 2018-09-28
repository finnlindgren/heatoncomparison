####################
## Load Libraries ##
####################
library(LatticeKrig)
library(R.matlab)
library(xtable)
library(gridExtra)
library(ggplot2)
library(ggpubr)
library(dplyr)

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
## Find how close nearest obs is to pred location ##
####################################################
# nmiss <- length(sat.miss)
# grp <- ceiling(seq(1,nmiss,length=100))
# satdists <- rep(NA,length(sat.miss))
# simdists <- rep(NA,length(sim.miss))
# for(g in 2:length(grp)){
#   if(g==length(grp)){
#     gseq <- grp[g-1]:(grp[g])
#   } else {
#     gseq <- grp[g-1]:(grp[g]-1)
#   }
#   mindist <- rdist(all.sat.temps[sat.miss[gseq],c("Lon","Lat")],
#                    all.sat.temps[-sat.miss,c("Lon","Lat")]) %>% apply(MARGIN=1,FUN=min)
#   satdists[gseq] <- mindist
#   mindist <- rdist(all.sim.data[sim.miss[gseq],c("Lon","Lat")],
#                    all.sim.data[-sim.miss,c("Lon","Lat")]) %>% apply(MARGIN=1,FUN=min)
#   simdists[gseq] <- mindist
# }
# save(file="./MinDists.RData",list=c("simdists","satdists"))
load("MinDists.RData")

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
sat.pp <- read.csv("../PredProc/pp2-real-pred.csv")
sat.pp$sd <- (sat.pp[,3]-sat.pp[,2])/(2*qnorm(0.975))
sim.pp <- read.csv("../PredProc/pp2-sim-pred.csv")
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
## Get Circ. Embed. Results ##
##############################
ce.sim <- read.table("../CirculantEmbedding/CESimResults.txt",header=TRUE)
ce.sat <- read.table("../CirculantEmbedding/CESatResults.txt",header=TRUE)

#######################################################
## Create a Data frame with each method result in it ##
#######################################################
nbreaks <- 5
satpred.df <- data.frame(Obs=Sat.true[sat.miss],
                         NNGP=nngp.sat$Pred,
                         NNGP.se=nngp.sat$sd,
                         NNGPc=nngp.sat.conj$Pred,
                         NNGPc.se=nngp.sat.conj$sd,
                         Taper=taper.sat$pred,
                         Taper.se=taper.sat$unc,
                         Gap=gapfill.sat$pred[sat.miss],
                         Gap.se=gapfill.sat$pse[sat.miss],
                         Part=partition.sat$pred[sat.miss],
                         Part.se=partition.sat$pred.se[sat.miss],
                         FRK=FRKreal$results$mu,
                         FRK.se=FRKreal$results$sd,
                         LK=LK.sat$yHat,
                         LK.se=LK.sat$standError,
                         SPDE=SPDE.sat$Temp[sat.miss],
                         SPDE.se=SPDE.sat$SD[sat.miss],
                         Meta=meta.sat$Y.med,
                         Meta.se=meta.sat$SD,
                         MRA=mra.sat$predMat[sat.miss],
                         MRA.se=mra.sat$sdMat[sat.miss],
                         PP=sat.pp[sat.preds,1],
                         PP.se=sat.pp$sd,
                         LAGP=sat.lagp$mean,
                         LAGP.se=sat.lagp$sd,
                         CircEmbed=ce.sat$pred[sat.miss],
                         CircEmbed.se=ce.sat$se[sat.miss],
                         Dists=cut(satdists,breaks=nbreaks))
satpred.df <- satpred.df %>% filter(!is.na(Obs))

simpred.df <- data.frame(Obs=Sim.true[sim.miss],
                         NNGP=nngp.sim$Pred,
                         NNGP.se=nngp.sim$sd,
                         NNGPc=nngp.sim.conj$Pred,
                         NNGPc.se=nngp.sim.conj$sd,
                         Taper=taper.sim$pred,
                         Taper.se=taper.sim$unc,
                         Gap=gapfill.sim$pred[sim.miss],
                         Gap.se=gapfill.sim$pse[sim.miss],
                         Part=partition.sim$pred[sim.miss],
                         Part.se=partition.sim$pred.se[sim.miss],
                         FRK=FRKsim$results$mu,
                         FRK.se=FRKsim$results$sd,
                         LK=LK.sim$yHat,
                         LK.se=LK.sim$standError,
                         SPDE=SPDE.sim$Temp[sim.miss],
                         SPDE.se=SPDE.sim$SD[sim.miss],
                         Meta=meta.sim$Y.med,
                         Meta.se=meta.sim$SD,
                         MRA=mra.sim$predMat[sim.miss],
                         MRA.se=mra.sim$sdMat[sim.miss],
                         PP=sim.pp[sim.preds,1],
                         PP.se=sim.pp$sd,
                         LAGP=sim.lagp$mean,
                         LAGP.se=sim.lagp$sd,
                         CircEmbed=ce.sim$pred[sim.miss],
                         CircEmbed.se=ce.sat$se[sim.miss],
                         Dists=cut(simdists,breaks=nbreaks))
simpred.df <- simpred.df %>% filter(!is.na(Obs))

#######################################
## Calculate RMSE & CRPS by Distance ##
#######################################
crps <- function(pred,pred.se,trueobs) {
  z <- as.numeric((trueobs - pred) / pred.se)
  scores <- pred.se * (z *(2 * pnorm(z, 0, 1) - 1) +
                             2 * dnorm(z, 0, 1) - 1/sqrt(pi))
  return(scores)
}

sat.res <- satpred.df %>% group_by(Dists) %>%
  summarize(NNGP.RMSE = sqrt(mean((NNGP-Obs)^2)),
            NNGPc.RMSE = sqrt(mean((NNGPc-Obs)^2)),
            Taper.RMSE = sqrt(mean((Taper-Obs)^2)),
            Gap.RMSE = sqrt(mean((Gap-Obs)^2)),
            Part.RMSE = sqrt(mean((Part-Obs)^2)),
            FRK.RMSE = sqrt(mean((FRK-Obs)^2)),
            LK.RMSE = sqrt(mean((LK-Obs)^2)),
            SPDE.RMSE = sqrt(mean((SPDE-Obs)^2)),
            Meta.RMSE = sqrt(mean((Meta-Obs)^2)),
            MRA.RMSE = sqrt(mean((MRA-Obs)^2)),
            PP.RMSE = sqrt(mean((PP-Obs)^2)),
            LAGP.RMSE = sqrt(mean((LAGP-Obs)^2)),
            CE.RMSE = sqrt(mean((CircEmbed-Obs)^2)),
            NNGP.CRPS = mean(crps(NNGP,NNGP.se,Obs)),
            NNGPc.CRPS = mean(crps(NNGPc,NNGPc.se,Obs)),
            Taper.CRPS = mean(crps(Taper,Taper.se,Obs)),
            Gap.CRPS = mean(crps(Gap,Gap.se,Obs),na.rm=TRUE),
            Part.CRPS = mean(crps(Part,Part.se,Obs)),
            FRK.CRPS = mean(crps(FRK,FRK.se,Obs)),
            LK.CRPS = mean(crps(LK,LK.se,Obs)),
            SPDE.CRPS = mean(crps(SPDE,SPDE.se,Obs)),
            Meta.CRPS = mean(crps(Meta,Meta.se,Obs)),
            MRA.CRPS = mean(crps(MRA,MRA.se,Obs)),
            PP.CRPS = mean(crps(PP,PP.se,Obs)),
            LAGP.CRPS = mean(crps(LAGP,LAGP.se,Obs)),
            CE.CRPS = mean(crps(CircEmbed,CircEmbed.se,Obs))
            )
Methods <- names(sat.res)[substr(names(sat.res),nchar(names(sat.res))-3,nchar(sat.res))=="RMSE"]
Methods <- sub("[.].*","",Methods)
sat.res <- data.frame(Dists=rep(sat.res$Dists,
                                 length(Methods)),
                       Method=rep(Methods,each=nlevels(sat.res$Dists)),
                       RMSE=unlist(c(sat.res[,substr(names(sat.res),nchar(names(sat.res))-3,nchar(sat.res))=="RMSE"])),
                      CRPS=unlist(c(sat.res[,substr(names(sat.res),nchar(names(sat.res))-3,nchar(sat.res))=="CRPS"]))
                      )
rownames(sat.res) <- 1:nrow(sat.res)

ggplot(sat.res,aes(x=Dists,y=RMSE,group=Method,color=Method))+
  geom_line()+geom_point()+xlab("Distance to Nearest Observation")+
  ggtitle("(a) RMSE")

ggplot(sat.res,aes(x=Dists,y=CRPS,group=Method,color=Method))+
  geom_line()+geom_point()+xlab("Distance to Nearest Observation")

p1 <- ggplot(filter(sat.res,Method%in%c("NNGPc","MRA","SPDE","LK","Part","FRK")),aes(x=Dists,y=RMSE,group=Method,color=Method))+
  geom_line()+geom_point()+xlab("Distance to Nearest Observation")+
  ggtitle("(a) RMSE")

p2 <- ggplot(filter(sat.res,Method%in%c("NNGPc","MRA","SPDE","LK","Part","FRK")),aes(x=Dists,y=CRPS,group=Method,color=Method))+
  geom_line()+geom_point()+xlab("Distance to Nearest Observation")+
  ggtitle("(b) CRPS")

ggarrange(p1,p2)

##########################################################
## Calculate RMSE & CRPS by Distance for Simulated Data ##
##########################################################
sim.res <- simpred.df %>% group_by(Dists) %>%
  summarize(NNGP.RMSE = sqrt(mean((NNGP-Obs)^2)),
            NNGPc.RMSE = sqrt(mean((NNGPc-Obs)^2)),
            Taper.RMSE = sqrt(mean((Taper-Obs)^2)),
            Gap.RMSE = sqrt(mean((Gap-Obs)^2)),
            Part.RMSE = sqrt(mean((Part-Obs)^2)),
            FRK.RMSE = sqrt(mean((FRK-Obs)^2)),
            LK.RMSE = sqrt(mean((LK-Obs)^2)),
            SPDE.RMSE = sqrt(mean((SPDE-Obs)^2)),
            Meta.RMSE = sqrt(mean((Meta-Obs)^2)),
            MRA.RMSE = sqrt(mean((MRA-Obs)^2)),
            PP.RMSE = sqrt(mean((PP-Obs)^2)),
            LAGP.RMSE = sqrt(mean((LAGP-Obs)^2)),
            CE.RMSE = sqrt(mean((CircEmbed-Obs)^2)),
            NNGP.CRPS = mean(crps(NNGP,NNGP.se,Obs)),
            NNGPc.CRPS = mean(crps(NNGPc,NNGPc.se,Obs)),
            Taper.CRPS = mean(crps(Taper,Taper.se,Obs)),
            Gap.CRPS = mean(crps(Gap,Gap.se,Obs),na.rm=TRUE),
            Part.CRPS = mean(crps(Part,Part.se,Obs)),
            FRK.CRPS = mean(crps(FRK,FRK.se,Obs)),
            LK.CRPS = mean(crps(LK,LK.se,Obs)),
            SPDE.CRPS = mean(crps(SPDE,SPDE.se,Obs)),
            Meta.CRPS = mean(crps(Meta,Meta.se,Obs)),
            MRA.CRPS = mean(crps(MRA,MRA.se,Obs)),
            PP.CRPS = mean(crps(PP,PP.se,Obs)),
            LAGP.CRPS = mean(crps(LAGP,LAGP.se,Obs)),
            CE.CRPS = mean(crps(CircEmbed,CircEmbed.se,Obs))
            )
sim.res <- data.frame(Dists=rep(sim.res$Dists,
                                length(Methods)),
                      Method=rep(Methods,each=nlevels(sim.res$Dists)),
                      RMSE=unlist(c(sim.res[,substr(names(sim.res),nchar(names(sim.res))-3,nchar(sim.res))=="RMSE"])),
                      CRPS=unlist(c(sim.res[,substr(names(sim.res),nchar(names(sim.res))-3,nchar(sim.res))=="CRPS"]))
                      )
rownames(sim.res) <- 1:nrow(sim.res)

p1 <- ggplot(filter(sim.res,Method%in%c("NNGPc","MRA","SPDE","LK","Part")),aes(x=Dists,y=RMSE,group=Method,color=Method))+
  geom_line()+geom_point()+xlab("Distance to Nearest Observation")+
  ggtitle("(a) RMSE")

ggplot(sim.res,aes(x=Dists,y=RMSE,group=Method,color=Method))+
  geom_line()+geom_point()+xlab("Distance to Nearest Observation")

ggplot(sim.res,aes(x=Dists,y=CRPS,group=Method,color=Method))+
  geom_line()+geom_point()+xlab("Distance to Nearest Observation")

p2 <- ggplot(filter(sim.res,Method%in%c("NNGPc","MRA","SPDE","LK","Part")),aes(x=Dists,y=CRPS,group=Method,color=Method))+
  geom_line()+geom_point()+xlab("Distance to Nearest Observation")+
  ggtitle("(b) CRPS")

ggarrange(p1,p2)

