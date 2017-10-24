# LatticeKrig model for big N challenge
# fits LatticeKrig model for different data sets and
# predicts at missing locations
#
########################################
######## Specify settings  #############
########################################
dataCase<-"SatelliteTemps"  # Choices are "test","simulated", and "SatelliteTemps"
runMode <-"parallel"  # Choices are "sequential" and "parallel" [Use "parallel" for cluster.]
M <- 100  # Number of conditional draws (used for standard errors) [Final model: 100]
nCores <- 55 # Number of cores (only relevant for parallel run mode) [Depends on processors/memory on cluster]
plotting <- "no" # Choices are "yes" or "no" [Use "yes" for testing on local machine, "no" for cluster.]


########################################
######## Load files and libraries  #####
########################################
library( LatticeKrig)
source("functions.R")

if( runMode=="parallel"){
  library(doParallel)
  library(foreach)
  source("LKrig.sim.conditional.foreach.R")
}
 
  
################################################
#### Load data and specify model parameters ####
################################################
  
# switch among the three cases
if( dataCase=="test"){
  load("../../Data/SmallTestData.RData")
  lon <- code.test$Lon
  lat <- code.test$Lat
  temps <- code.test$MaskedData
  NC<- 16
  nlevel<- 2
  a.wght<- 4.4
  nu<- .5
}

if( dataCase=="simulated"){
  load("../../Data/SimulatedTemps.RData")
  lon <- sim.data$Lon
  lat <- sim.data$Lat
  temps <- sim.data$Temp
  NC<- 30
  nlevel<-4
  a.wght<-  4.4
  nu<- .1
}

if( dataCase=="SatelliteTemps"){
  load("../../Data/SatelliteTemps.RData")
  lon <- sat.temps$Lon
  lat <- sat.temps$Lat
  temps <- sat.temps$Temp
  NC<- 40
  nlevel<- 4
  a.wght<-  10.25
  nu<- .1
 }

########################
#### Fit the model ####
########################

# create locations and responses (y)
dataObject<- makeData(lon,lat, temps)

# setup LKrig object
LKinfoFinal<- LKrigSetup( dataObject$x,
                         NC = NC,
                         nlevel = nlevel,
                         a.wght = a.wght,
                             nu = nu)
# determine fit
  comp.time <- system.time(fitFinal<- LatticeKrig( dataObject$x, dataObject$y, LKinfo= LKinfoFinal))
  if( plotting=="yes"){
  print( fitFinal) }# Display of fit

###########################################################
#### Conditional simulation to determine prediction SE ####
###########################################################
# M is the number of independent draws from conditional
  set.seed(234)
  
  if( runMode=="parallel"){
    comp.time <- comp.time+system.time(outputSim<- LKrig.sim.conditional.foreach(fitFinal,
                                              x.grid = dataObject$xMissing,
                                              M = M, nCores=nCores))}
  if( runMode=="sequential"){
    comp.time <- comp.time + system.time(outputSim<- LKrig.sim.conditional(fitFinal,
                                    x.grid = dataObject$xMissing,
                                    M = M))}
  

# NOTE: prediction SE needs to include uncertainty due to nugget
# random component ( estimated as fitFinal$sigma.MLE )
  
  comp.time <- comp.time + system.time(standardError<- sqrt( 
      apply( outputSim$g.draw, 1, "var") +
      fitFinal$sigma.MLE^2))

  
  yHat<- outputSim$ghat
  CI95Lower<- yHat - 1.96* standardError
  CI95Upper<- yHat + 1.96* standardError


########################################################
####  Save everything important to an R binary file ####
########################################################
  
  finalResults<- list(x=dataObject$xMissing,
                      yHat=yHat, 
                      standError=standardError)
  
  save( finalResults, fitFinal, comp.time, file=paste0(dataCase,"NC",NC,"nlevel",nlevel,"finalResults.rda") )

  
###################################
#### Summary figure with fit  ####
##################################
  if( plotting=="yes"){
  pdf(paste0(dataCase,"NC",NC,"nlevel",nlevel,"SummaryFit.pdf"), width=6, height=8)
    makeFig() # pdf panel of fit and prediction
  dev.off()

  cat("All Done", fill=TRUE)}







