#######################################################################
## Author: Andrew Zammit-Mangion
## Date: 8 May 2017
## Description: Reproducible code for FRK
## Instructions: Please set folder to current directory and source it
## Requirements: The folder needs to have the files
##               -SmallTestData.RData
##               -SimulatedTemps.RData
##               -SatelliteTemps.RData
##
##               You need to have the following packages installed
##               -FRK
##               -sp
##               -ggplot2
##               -gridExtra
##               -gstat
##               -INLA
##               -splancs
#######################################################################

## Load packages
library(FRK)
library(sp)
library(ggplot2)
library(gridExtra)
library(INLA)
library(splancs)

## runFRK: Main program function
## Takes a data frame containing Lon, Lat, Temp containing NAs for missing data
## and returns a list with (i) a results data frame containing prediction, prediction
## errors and 95% coverage intervals at unobserved locations, (ii) the time taken
## to run FRK, and (iii) spatial plots of the data, predictions and prediction errors.
runFRK <- function(df, f = Temp ~ Lon + Lat, tol = 0.01) {

    t1 <- proc.time()                     # start timer

    ## Make BAUs as SpatialPixels
    BAUs <- df                            # assign BAUs
    BAUs$Missing <- is.na(BAUs$Temp)      # mark which BAUs contain missing data
    BAUs$Temp <- NULL                     # remove data from BAUs
    BAUs$fs <- 1                          # set fs variation to unity
    coordinates(BAUs)  <- ~Lon+Lat        # convert to SpatialPointsDataFrame
    gridded(BAUs) <- TRUE                 # convert to SpatialPixelsDataFrame

    ## Make Data as SpatialPoints
    dat <- subset(df,!is.na(Temp))        # no missing data in data frame
    coordinates(dat)  <- ~Lon+Lat         # Convert to SpatialPointsDataFrame

    basis <- auto_basis(plane(),          # we are on the plane
                        data = dat,       # data around which to make basis
                        regular = 0,      # irregular basis
                        nres = 3,         # 3 resolutions
                        scale_aperture = 1)   # aperture scaling of basis functions 
    
    ## Remove basis functions in problematic region
    if(nrow(df) == 150000) {
        basis_df <- data.frame(basis)
        rmidx <- which(basis_df$loc2 > 36.5 &
                           basis_df$loc1 > -94.5 &
                           basis_df$res == 3)
        basis <- remove_basis(basis,rmidx)
    }

    ## Estimate using ML
    S <- FRK(f = f,                       # formula for SRE model
             data = dat,                  # data
             basis = basis,               # Basis
             BAUs = BAUs,                 # BAUs
             tol = tol)                   # EM iterations
             

    ## Predict
    BAUs_pred <- SRE.predict(S)           # predict over all BAUs
    BAUs_pred_df <- data.frame(BAUs_pred) # convert to data frame

    ## Compute variance of predicted observations and conf. intervals
    BAUs_pred_df$sd_obs <- sqrt(BAUs_pred_df$sd^2 + S@Ve[1,1])
    BAUs_pred_df$conflo <- BAUs_pred_df$mu - 1.96*BAUs_pred_df$sd_obs
    BAUs_pred_df$confhi <- BAUs_pred_df$mu + 1.96*BAUs_pred_df$sd_obs

    ## Only return missing data and relevant variables
    Res <- subset(BAUs_pred_df,Missing==TRUE)
    Res$Missing <- Res$fs <- Res$var <- NULL

    t2 <- proc.time()                     # stop timer

    ## Plots for illustration purposes
    ## Data
    gdata <- ggplot(df) + geom_raster(aes(Lon,Lat,fill=Temp)) +
        theme_bw() + scale_fill_distiller(palette="Spectral") +
        ggtitle("Data")
    ## Predictions
    g1 <- ggplot(BAUs_pred_df) + geom_raster(aes(Lon,Lat,fill=mu)) +
        theme_bw() + scale_fill_distiller(palette="Spectral",name="Pred.") + ggtitle("Pred.")
    ## Prediction errors
    g2 <- ggplot(BAUs_pred_df) + geom_raster(aes(Lon,Lat,fill=sqrt(sd^2 + S@Ve[1,1]))) +
        theme_bw() + scale_fill_distiller(palette="Spectral",name="Error") + ggtitle("Pred. Error")

    ## Return in a list
    list(results = Res,
         SREmodel = S,
         tot_time = t2 - t1,
         plotData = gdata,
         plotPred = g1,
         plotSE = g2)
}


## Small Test Data experiment
load("../../Data/SmallTestData.RData")
X <- code.test                           # save code.test in X
code.test$Temp <- code.test$MaskedData   # change data frame to mimic other data frames
code.test$MaskedData <-
    code.test$FullData <- NULL

set.seed(25)
FRKtest <- runFRK(code.test)             # run FRK with this data frame
X <- merge(FRKtest$results,
               subset(X,is.na(MaskedData)))  # merge FRK results with original data
mean(X$mu - X$FullData)                               # Compute bias
sqrt(mean((X$mu - X$FullData)^2))                     # Compute RMSPE
mean(X$FullData > X$conflo & X$FullData < X$confhi)   # Compute empiriral 95% coverage

## Simulated Data experiment
load("../../Data/SimulatedTemps.RData")
set.seed(25)
FRKsim <- runFRK(sim.data,f <- Temp ~ 1, tol=0.1)
save(file='./FRKSimResults.RData',list="FRKsim")

## Satellite Data experiment
load("../../Data/SatelliteTemps.RData")
set.seed(25)
FRKreal <- runFRK(sat.temps,tol=0.02)
save(file='./FRKSatResults.RData',list="FRKreal")

## Plot results
# gridExtra::grid.arrange(FRKtest$plotData,FRKtest$plotPred,FRKtest$plotSE,ncol=3)
# gridExtra::grid.arrange(FRKsim$plotData,FRKsim$plotPred,FRKsim$plotSE,ncol=3)
# gridExtra::grid.arrange(FRKreal$plotData,FRKreal$plotPred,FRKreal$plotSE,ncol=3)
