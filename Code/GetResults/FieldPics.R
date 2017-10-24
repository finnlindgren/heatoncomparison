####################
## Load Libraries ##
####################
library(LatticeKrig)
library(ggplot2)
library(gridExtra)
fig.path <- "~/Documents/Professional/Research/Active Research/Spatial/LargeNComparison/LaTeXFiles/Figures/"

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
names(all.sat.temps)
zlims <- range(c(all.sat.temps$TrueTemp,all.sim.data$TrueTemp),na.rm=TRUE)


######################
## Plot of Datasets ##
######################
satplot <- ggplot(all.sat.temps) + geom_raster(aes(Lon,Lat,fill=TrueTemp)) +
  theme_bw() + scale_fill_distiller(palette="Spectral",name="",limits=zlims) + 
  ggtitle("(a)")
satmaskplot <- ggplot(all.sat.temps) + geom_raster(aes(Lon,Lat,fill=MaskTemp)) +
  theme_bw() + scale_fill_distiller(palette="Spectral",name="",limits=zlims) +
  ggtitle("(b)")
simplot <- ggplot(all.sim.data) + geom_raster(aes(Lon,Lat,fill=TrueTemp)) +
  theme_bw() + scale_fill_distiller(palette="Spectral",name="",limits=zlims) + 
  ggtitle("(c)")
simmaskplot <- ggplot(all.sim.data) + geom_raster(aes(Lon,Lat,fill=MaskTemp)) +
  theme_bw() + scale_fill_distiller(palette="Spectral",name="",limits=zlims) + 
  ggtitle("(d)")
postscript(file=paste0(fig.path,"DataFigure.eps"),width=11,heigh=8.5)
gridExtra::grid.arrange(satplot,satmaskplot,simplot,simmaskplot,ncol=2)
dev.off()






