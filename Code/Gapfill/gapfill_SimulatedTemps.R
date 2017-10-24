rm(list = ls())
## all packages are available on CRAN
library("gapfill") 
library("doParallel");
registerDoParallel(40) # run 40 tasks in parallel
library("abind")

## load data
load("../../Data/SimulatedTemps.RData")

## rearrange data as matrix
data <- with(sim.data,
             array(Temp, c(length(unique(Lon)), length(unique(Lat)))))
dim <- dim(data)
nx <- dim[1]; ny <- dim[2]

## display data
## Image(data)

## augment data: since the gapfill method is designed for
## spatio-temproal data, we artificially create 9 additional
## and similar images by shifting the given image
data_augmented <- abind(data,
                        data[c(1,1:(nx-1)),],
                        data[,c(1,1:(ny-1))],
                        data[c(2:nx,nx),],
                        data[,c(2:ny,ny)],

                        data[c(2:nx,nx),c(2:ny,ny)],
                        data[c(2:nx,nx),c(1,1:(ny-1))],
                        data[c(1,1:(nx-1)),c(2:ny,ny)],
                        data[c(1,1:(nx-1)),c(1,1:(ny-1))],
                        
                        data[c(3:nx,nx,nx),],
                        data[,c(3:ny,ny,ny)],
                        data[c(1,1,1:(nx-2)),],
                        data[,c(1,1,1:(ny-2))],
                        along = 3)
dim(data_augmented) <- c(nx, ny, dim(data_augmented)[3], 1)


## predict missing values
out <- Gapfill(data_augmented,
               ## only predict missing values in first (original) image
               subset = which(is.na(data_augmented[,,1,1])),
               ## use parallel processing via R package foreach
               dopar = TRUE,
               ## tuning parameters of the algorithm
               initialSize = c(2L, 2L, 100L, 100L),
               nTargetImage = 2,
               nQuant = 3,
               ## restrict values to the following range
               clipRange = range(data_augmented[,,1,1], na.rm=TRUE), 
               ## return prediction interval
               nPredict = 3, predictionInterval = TRUE  
               )

## extract prediction and prediction interval
prediction <- out$fill[,,1,1,1]
ciLo <- out$fill[,,1,1,2]
ciUp <- out$fill[,,1,1,3]

## save results to file
save(prediction, ciLo, ciUp, file = "./data/SimulatedTemps_gapfill.RData")


## summary images
jpeg("./figs/SimulatedTemps.jpeg", 1000, 700)
Image(abind(data, prediction, along = 3))
dev.off()

jpeg("./figs/SimulatedTemps_ci.jpeg", 1000, 700)
Image(abind(ciLo, ciUp, along = 3))
dev.off()

jpeg("./figs/SimulatedTemps_ci.jpeg", 500, 700)
Image(ciUp-ciLo)
dev.off()


sessionInfo()
