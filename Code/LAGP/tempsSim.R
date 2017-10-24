## Submission to Matthew Heaton's "Large N Competition"
## by Robert B. Gramacy and Furong Sun, Virginia Tech

## change this to 2x the number of cores on the node
## (assumes hyper-threading; only one node is needed)
nth <- 40

## Here is what we need
library(LatticeKrig)
library(tgp)
library(laGP)
library(plgp)

## simulated problem
load("../../Data/SimulatedTemps.RData")
pdf("temps_sim.pdf")
data <- sim.data
nr <- 500; nc <- 300
Lon <- matrix(data$Lon, nrow=nr, ncol=nc)
Lat <- matrix(data$Lat, nrow=nr, ncol=nc)
Temp <- matrix(data$Temp, nrow=nr, ncol=nc)

## extract the training/testing problem
nas <- which(is.na(data$Temp))
X <- cbind(data$Lon, data$Lat)[-nas,] # training input
N <- nrow(X)
y <- data$Temp[-nas] # training resposne
XX <- cbind(data$Lon, data$Lat)[nas,] # test input
NN <- nrow(XX)

## code inputs and predictive grid to unit cube
maxX <- apply(rbind(X, XX), 2, max)
minX <- apply(rbind(X, XX), 2, min)
for (j in 1:ncol(X)){
     X[,j] <- X[,j] - minX[j]	
     X[,j] <- X[,j]/(maxX[j]-minX[j])  
     XX[,j] <- XX[,j] - minX[j]	
     XX[,j] <- XX[,j]/(maxX[j]-minX[j])
}

##
## time for fitting
##

## macro-scale analysis on a maximum entropy sub-design of size n=100
n <- 100
sub <- dopt.gp(100, Xcand=X)

## priors for the global (subset) GP
da <- darg(list(mle=TRUE, max=10), X)
ga <- garg(list(mle=TRUE, max=10), y)

## fit the global GP
gpsepi <- newGPsep(sub$XX, y[sub$fi], d=da$start, g=ga$start, dK=TRUE)
that <- mleGPsep(gpsepi, param="both", tmin=c(da$min, ga$min), 
  tmax=c(da$max, ga$max), ab=c(da$ab, ga$ab), maxit=200)

## predictions from the global GP on the test set
psub <- predGPsep(gpsepi, XX, lite=TRUE)

## calculation of residuals at the full set of input locations
pX <- predGPsep(gpsepi, X, lite=TRUE)
yresid <- y - pX$mean

## done with the psub GP, clean up
deleteGPsep(gpsepi)

## Now for laGP on residuals from the (subset) global GP

## scale the inputs according to the macro-analysis lengthscales
scale <- sqrt(that$theta[1:2])
Xs <- X; XXs <- XX
for(j in 1:ncol(Xs)){
 	  Xs[,j] <- Xs[,j] / scale[j]
	  XXs[,j] <- XXs[,j] / scale[j]
}

## local analysis on residuals (and scaled inputs) from local analysis
out <- aGPsep(Xs, yresid, XXs, d=list(start=1, max=20), g=that$theta[3], 
  omp.threads=nth, verb=0)

## 
## Calulations are basically done, just have to process some of the
## outputs for intervals and visualization
##

## 2x2 grid for visualization
par(mfrow=c(2,2))

## original data
image.plot(Lon, Lat, Temp, zlim=range(data$Temp, na.rm=TRUE),
           xlab="Longitude", ylab="Latitude",
           main="Original Data with Subset")
## with space-filling sub-design for global GP
points(data$Lon[-nas][sub$fi], data$Lat[-nas][sub$fi], pch=19)

## visualization of global GP predictive surface
p <- data$Temp
p[nas] <- psub$mean
PTemp <- matrix(p, nrow=nr, ncol=nc)
image.plot(Lon, Lat, PTemp, zlim=range(data$Temp, na.rm=TRUE),
	         xlab="Longitude", ylab="Latitude",
	         main="Global GP on Subset")

## visualiation of the global-local hybrid
p2 <- data$Temp
p2[nas] <- out$mean + psub$mean
PTemp2 <- matrix(p2, nrow=nr, ncol=nc)
image.plot(Lon, Lat, PTemp2, zlim=range(data$Temp, na.rm=TRUE),
	         xlab="Longitude", ylab="Latitude",
	         main="Global-Local laGP")

##
## Assessing predictive uncertainty for the global-local hybrid
## requires combining uncertainty from two models: uncertainty
## in the global gp mean used to define residuals, and the 
## full predictive uncertainty from the laGP on the residuals.
##

## Deriving the uncertainty in the global GP mean requires 
## removing the nugget from psub variance.
library(plgp)
K <- covar.sep(sub$XX, d=that$theta[1:2], g=that$theta[3])
psi <- drop(t(y[sub$fi]) %*% solve(K) %*% y[sub$fi])
s2.f <- psub$s2 * nrow(sub$XX) / psi
eps <- sqrt(.Machine$double.eps)
s2.f <- s2.f - that$theta[3] + eps
s2.f <- s2.f * psi / nrow(sub$XX)

## combine variances from psub (from residuals) and local GP
stdev <- sqrt(s2.f + out$var)
p3 <- data$Temp
p3[nas] <- stdev
p3[-nas] <- NA
PTemp3 <- matrix(p3, nrow=nr, ncol=nc)
image.plot(Lon, Lat, PTemp3, xlab="Longitude", ylab="Latitude",
 	         main="Predictive SD based on Global-Local laGP")

## outputs for missing values in terms of means and 95% interval
pmean <- out$mean + psub$mean
pupper <- pmean + 1.96*stdev
plower <- pmean - 1.96*stdev

save(file="./SimLAGP.RData",list=c("pmean","pupper","plower"))

##
## That's all folks
##
dev.off()
