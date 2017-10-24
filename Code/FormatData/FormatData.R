##
## Test to see about downloaded data
##

library(raster)
library(rgdal)
library(LatticeKrig)
library(maps)
library(gridExtra)
library(ggplot2)
source("./fitMaternGP.R")
library(RandomFields)

########################################
## Format the MODIS Data for Analysis ##
########################################
temp <- raster('../../Data/MOD11A1.MRTWEB.A2016217.041.LST_Day_1km.tif')
tod <- raster('../../Data/MOD11A1.MRTWEB.A2016217.041.Day_view_time.tif')
#plot(temp)
t2 <- t(as.matrix(temp))
tod <- t(as.matrix(tod))
#t2 <- (9/5)*t2-459.67
#image.plot(1:1200,1:1200,t(t2))
t2[t2==0] <- NA
t2 <- t2*0.02-273.15
x.grid <- (seq(-10007555,-8895604,length=1200))/100000 #LON coordinates
y.grid <- rev(seq(3335852,4447802,length=1200))/100000 #LAT coordinates
xy.grid <- expand.grid(x.grid,y.grid)
xgmat <- matrix(xy.grid[,1],nrow=length(x.grid))
ygmat <- matrix(xy.grid[,2],nrow=length(x.grid))
# map('state')
# image.plot(xgmat,ygmat,t2,add=TRUE)
# map('state',add=TRUE)

## Plot it
par(mfrow=c(1,2))
sub.len.lat <- 300-1
sub.len.lon <- 500-1
lon.cut <- 450
lat.cut <- 800
sub.t2 <- t2[lon.cut:(lon.cut+sub.len.lon),lat.cut:(lat.cut+sub.len.lat)]
sub.tod <- tod[lon.cut:(lon.cut+sub.len.lon),lat.cut:(lat.cut+sub.len.lat)]
subx <- x.grid[lon.cut:(lon.cut+sub.len.lon)]
suby <- y.grid[lat.cut:(lat.cut+sub.len.lat)]
subxy <- expand.grid(subx,suby)
subxmat <- matrix(subxy[,1],nrow=length(subx))
subymat <- matrix(subxy[,2],nrow=length(subx))
#map('state')
image.plot(subxmat,subymat,sub.t2)
map('state',add=TRUE)

## Read in the mask to create the training dataset
mask <- raster('../../Data/MOD11A1.MRTWEB.A2016218.041.Clear_day_cov.tif')
mask <- t(as.matrix(mask))
mask[mask==0] <- NA
mask <- mask[lon.cut:(lon.cut+sub.len.lon),lat.cut:(lat.cut+sub.len.lat)]
masked <- sub.t2
masked[is.na(mask)] <- NA
image.plot(subxmat,subymat,masked)
map('state',add=TRUE)
sum(is.na(t2[lon.cut:(lon.cut+sub.len.lon),lat.cut:(lat.cut+sub.len.lat)]))/prod(dim(sub.t2))
sum(!is.na(masked))
sum(is.na(masked))/(prod(dim(sub.t2)))

## Create satellite dataframes
all.sat.temps <- data.frame(Lon=c(subxmat),Lat=c(subymat),MaskTemp=c(masked),TrueTemp=c(sub.t2))
sat.temps <- all.sat.temps[,-4]
names(sat.temps) <- c("Lon","Lat","Temp")

###############################################
## Explore the Data to Create Simulated Data ##
###############################################

## Fit a Matern 1/2 model to sample of true data to get parameters
sub.data <- sample(1:nrow(sat.temps),2500)
sub.data <- all.sat.temps[sub.data,]
sub.data <- sub.data[!is.na(sub.data[,4]),]
quilt.plot(sub.data[,1],sub.data[,2],sub.data[,4],nx=nrow(subxmat),ny=ncol(subxmat))
#my.GP <- fit.Matern.MC(matrix(sub.data[,4],ncol=1),matrix(1,nrow=nrow(sub.data)),as.matrix(sub.data[,1:2]),nu=1/2)
## alpha=1.264009, sigma2=16.40771, tau2=0.8635636, mean=44.49105

## Simulate
xlocs <- unique(sat.temps$Lon)
ylocs <- unique(sat.temps$Lat)
minloc <- c(xlocs[1],ylocs[1])
maxloc <- c(xlocs[length(xlocs)],ylocs[length(ylocs)])
maxD <- sqrt(sum((minloc-maxloc)^2))
myRFmod <- RMmatern(nu=0.5,scale=1/0.75,var=16.40771)+
  RMnugget(var=0.05)+RMtrend(mean=44.49105)
tst <- RFsimulate(myRFmod,x=xlocs,y=ylocs)
tst <- as.matrix(tst)
par(mfrow=c(1,2))
image.plot(subxmat,subymat,tst)
tst.mask <- tst
tst.mask[is.na(sat.temps$Temp)] <- NA
image.plot(subxmat,subymat,tst.mask)

## Create Simulated Dataframes
sim.data <- as.data.frame(cbind(sat.temps$Lon,sat.temps$Lat,c(tst.mask)))
names(sim.data) <- c("Lon","Lat","Temp")
all.sim.data <- as.data.frame(cbind(sat.temps$Lon,sat.temps$Lat,c(tst.mask),c(tst)))
names(all.sim.data) <- c("Lon","Lat","MaskTemp","TrueTemp")

#######################################################
## Simulate a Smaller Test Dataset for code checking ##
#######################################################
xlocs <- seq(0,1,length=100)
ylocs <- seq(0,1,length=100)
xmat <- as.matrix(expand.grid(xlocs,ylocs))
myRFmod <- RMmatern(nu=0.5,scale=1/1.5,var=16.40771)+
  RMnugget(var=0.05)+RMtrend(mean=44.49105)
tst <- RFsimulate(myRFmod,x=xlocs,y=ylocs)
tst <- as.matrix(tst)
ymat <- matrix(xmat[,2],nrow=100)
xmat <- matrix(xmat[,1],nrow=100)
par(mfrow=c(1,2))
image.plot(xmat,ymat,tst)
holdout <- sample(100^2,round(0.29*prod(dim(xmat))))
masked <- c(tst)
masked[holdout] <- NA
code.test <- data.frame(Lon=c(xmat),Lat=c(ymat),MaskedData=masked,FullData=c(tst))
image.plot(xmat,ymat,matrix(masked,nrow=100))

###################################################
## Write Files (WARNING: WILL OVERRIDE DATASETS) ##
###################################################
#save(file='~/Google Drive/Research/LargeSpatial/Datasets/SmallTestData.RData',list="code.test")
#save(file='../../Data/SimulatedTemps.RData',list="sim.data")
#save(file='../../Data/AllSimulatedTemps.RData',list="all.sim.data")
#save(file='../../Data/AllSatelliteTemps.RData',list="all.sat.temps")
#save(file='../../Data/SatelliteTemps.RData',list="sat.temps")

#########################
## Plot Final Datasets ##
#########################
rm(list=ls())
load("../../Data/AllSatelliteTemps.RData")
sat.gg <- ggplot(all.sat.temps) + geom_raster(aes(Lon,Lat,fill=TrueTemp)) +
  theme_bw() + scale_fill_distiller(palette="Spectral",name="") +
  ggtitle("(a)")
satmiss.gg <- ggplot(all.sat.temps) + geom_raster(aes(Lon,Lat,fill=MaskTemp)) +
  theme_bw() + scale_fill_distiller(palette="Spectral",name="") +
  ggtitle("(b)")
load("../../Data/AllSimulatedTemps.RData")
sim.gg <- ggplot(all.sat.temps) + geom_raster(aes(Lon,Lat,fill=TrueTemp)) +
  theme_bw() + scale_fill_distiller(palette="Spectral",name="") +
  ggtitle("(c)")
simmiss.gg <- ggplot(all.sim.data) + geom_raster(aes(Lon,Lat,fill=MaskTemp)) +
  theme_bw() + scale_fill_distiller(palette="Spectral",name="") +
  ggtitle("(d)")
gridExtra::grid.arrange(sat.gg,satmiss.gg,sim.gg,simmiss.gg,ncol=2)












