library(INLA) ## See http://www.r-inla.org/download/ for installation instructions.
library(fields) ## Only used for plotting
library(pals) ## Palettes for plotting

## Input:
data.path <- "../../Data" ## Path to data files
dataset <- 2
## 1: SatelliteTemps.RData
## 2: SimulatedTemps.RData
## 3: SmallTestData.RData

## Model settings:
spde.alpha <- 1.5 ## 1.5 gives an approximation of an exponential covariance

## inla() options:
num.threads = 2 ## Limit memory usage by limiting the number of openmp threads.
openmp.strategy = "large" ## "large" turns off nested parallelism

## Output:
## 
## Model 1, "Trend": intercept + longitude + latitude
## Model 2, "TrendField": intercept + longitude + latitude + gmrf/spde
##
## Posterior reconstructions for all lattice points, in the same order as the input data:
##   reconstructionsX: data.frame(Lon, Lat, Temp, SD)
##   predX: data.frame(mean, sd), where mean=Temp and sd=SD:
##          special objects used when calling the scoring methods from assessment.R
##   tempsX: list(mean, sd), same information as predX, but in matrix format for plotting
##   temps, temps0, temps00: Data in matrix format, tempsX, but only list(mean)
##
## resultX: inla model outputs
## times: Overall running times for each method
## scores.observed: Assessment scores evaluated for the observed data;
##                  Only to be used for sanity checking.
## scores.unobserved: Assessment scores evaluated for the holdout(unbserved) data
##                    for dataset 3.
##                    Set Temp0 below to the full data for dataset 2 for comparisons.
##
## The posterior predictive distributions F_i=N(m_i, s_i^2) include the observation noise
## variance, so are appropriate to evaluate on "new" observations y.
##
## Mean scores: S({F_i}, {y_i}) := 1/N \sum_{i=1}^N S(F_i, y_i)
##
## All the scores are negatively oriented, i.e. small values are
## better, and all are proper, i.e. the expected score under y~G
## cannot be made smaller than when F=G (the predictive distribution
## is the same as the truth).
## For more scoring details, see Gneiting & Raftery, JASA, 2007.
##
## MSE: mean squared error, S(F,y) = (y-m)^2
## MAE: mean absolute error, S(F,y) = |y-m|
## IGN: ignorance score, S(F,y) = (y-m)^2 / s^2 / 2 + log(s)
## CRPS: contiunuous ranked probability score, S(F,y) = CRPS(N(m,s^2), y)
## Interval: interval score, S(F,y) = IntervalScore([m-q*s, m+q*s], y)
## Coverage: prediction interval coverage (not a score) = I{y in [m-q*s, m+q*s]}
## q is the 1-alpha/2 quantile of N(0,1), and alpha=0.05
##
## Remark: (y-m)^2 + s^2 is not a proper score, and should not be
## used. It's frighteningly common to see it used in the wild. Resist
## the temptation. (See G&R above for more details.)


## Unify data storage to allow the same code to run for all three data sets.
## Temp: properly observed observations
## Temp0: all observations, including holdout data
## Temp00: all observations not in Temp (Temp00 is only used when plotting results)
if (dataset == 1) {
    load(file.path(data.path, "SatelliteTemps.RData"))
    n.row <- 500
    n.col <- 300
    indata <- sat.temps
    indata$Temp0 <- indata$Temp
    extend <- -0.5
    max.edge <- 360/38820*10
} else if (dataset == 2) {
    load(file.path(data.path, "SimulatedTemps.RData"))
    n.row <- 500
    n.col <- 300
    indata <- sim.data
    indata$Temp0 <- indata$Temp ## Change to the full dataset for comparison
    extend <- -0.5
    max.edge <- 360/38820*10
} else {
    load(file.path(data.path, "SmallTestData.RData"))
    n.row <- 100
    n.col <- 100
    indata <- code.test
    indata$Temp <- indata$MaskedData
    indata$Temp0 <- indata$FullData
    extend <- -0.5
    max.edge <- 1/100*10
}

## Construct prior parameter distributions heuristically:
## Medium range, relative half-spread factor
prior.range <- c(max(c(diff(range(indata$Lon)),
                       diff(unique(indata$Lat)))) * 1/5, 5)
## Medium sd, relative half-spread factor
prior.sigma <- c(sd(indata$Temp, na.rm=TRUE) / 2, 4)

## Extract available but "not observed" values.
indata$Temp00 <- indata$Temp0
indata$Temp00[!is.na(indata$Temp)] <- NA

## Construct centred covariate versions
LonCentre <- mean(range(indata$Lon))
LatCentre <- mean(range(indata$Lat))
indata$LonC <- indata$Lon - LonCentre
indata$LatC <- indata$Lat - LatCentre

lon <- matrix(indata$Lon,nrow=n.row)
lat <- matrix(indata$Lat,nrow=n.row)
temps <- list(mean=matrix(indata$Temp,nrow=n.row))
temps0 <- list(mean=matrix(indata$Temp0,nrow=n.row))
temps00 <- list(mean=matrix(indata$Temp00,nrow=n.row))

## Run a linear trend model

time1 <- system.time({
    formula1 <- Temp ~ 1 + LonC + LatC
    result1 <- inla(formula1, data=indata, family="gaussian",
                    control.predictor=list(compute=TRUE),
                    control.inla=list(strategy="gaussian", int.strategy="eb"),
                    verbose=TRUE, num.threads=num.threads,
                    control.compute=list(openmp.strategy=openmp.strategy))
    
    reconstruction1 <- data.frame(
        Lon=indata$Lon,
        Lat=indata$Lat,
        Temp=result1$summary.linear.predictor[,"mean"],
        SD=sqrt(result1$summary.linear.predictor[,"sd"]^2
                + 1/result1$summary.hyperpar[1,"0.5quant"]))
    
    temps1 <- list(mean=matrix(reconstruction1$Temp,nrow=n.row),
                   sd=matrix(reconstruction1$SD,nrow=n.row))
})


## Run the full model.
time2 <- system.time({
time2.prep <- system.time({                                                                                            
    loc2 <- cbind(indata$Lon, indata$Lat)
    mesh2 <- inla.mesh.create(loc=loc2,
                              extend=list(offset=extend),
                              refine=list(min.angle=30))
    proj2 <- inla.mesh.projector(mesh2, loc2)
    col.idx <- apply(proj2$proj$bary, 1, which.max)
    idx2 <- mesh2$graph$tv[proj2$proj$t + (col.idx-1)*nrow(mesh2$graph$tv)]
    
    ## lognormal prior
    ## Centre the parameterisation at range=1, sigma=1
    nu <- spde.alpha - 1
    kappa.zero <- sqrt(8*nu) / 1
    tau.zero <- (gamma(nu) / (gamma(spde.alpha) * 4*pi * kappa.zero^(2*nu)) )^0.5 / 1
    ## sigma^2 = Gamma(0.5)/( Gamma(1.5) (4\pi)^(dim/2) * kappa^(2*0.5) * tau^2 )
    ## tau = [ Gamma(0.5)/( Gamma(1.5) (4\pi)^(dim/2) * kappa^(2*0.5) ) ]^0.5 / sigma
    ## kappa = sqrt(8*0.5) / range
    ## tau = [ Gamma(0.5)/( Gamma(1.5) (4\pi)^(dim/2) * (sqrt(8*0.5)/range)^(2*0.5) ) ]^0.5 / sigma
    B.tau   <- cbind(log(tau.zero), nu, -1)
    B.kappa <- cbind(log(kappa.zero), -1, 0)

    theta.prior.mean <- log(c(prior.range[1], prior.sigma[2]))
    ## Half-width of prior prediction interval, on log-scale: 2*sd = log(rel)
    theta.prior.prec <- 4 / log(c(prior.range[2], prior.sigma[2]))^2
    spde2 <- inla.spde2.matern(mesh2, alpha=spde.alpha,
                               B.tau = B.tau, B.kappa = B.kappa,
                               theta.prior.mean = theta.prior.mean,
                               theta.prior.prec = theta.prior.prec,
                               constr=FALSE)
    
    hyper.range.initial <- log(prior.range[1])
    hyper.sigma.initial <- log(prior.sigma[1])
    hyper.family.initial <- 2
      
    spde2$f$hyper.default$theta1$initial <- hyper.range.initial
    spde2$f$hyper.default$theta2$initial <- hyper.sigma.initial

    formula2 <- Temp ~ 1 + LonC + LatC + f(field, model=spde)
    
    ok <- !is.na(indata$Temp) | TRUE ## Include all; NA data are predicted
    data2 <- list(Temp=indata$Temp[ok], LonC=indata$LonC[ok], LatC=indata$LatC[ok],
                  field=idx2[ok], spde=spde2)
})
  time2.run <- system.time({  
    result2 <- inla(formula2, data=data2, family="gaussian",
                    control.family=list(hyper=list(theta=list(initial=hyper.family.initial))),
                    control.predictor=list(compute=TRUE),
                    control.inla=list(strategy="gaussian", int.strategy="eb",
                                      force.diagonal=TRUE, stupid.search=FALSE),
                    verbose=TRUE, num.threads=num.threads,
                    control.compute=list(openmp.strategy=openmp.strategy))
  })
    times2.post <- system.time({
    reconstruction2 <- data.frame(Lon=indata$Lon,
                                  Lat=indata$Lat,
                                  Temp=indata$Temp*NA,
                                  SD=indata$Temp*NA)
    reconstruction2$Temp[ok] <- result2$summary.linear.predictor[,"mean"]
    reconstruction2$SD[ok] <- 
        sqrt(result2$summary.linear.predictor[,"sd"]^2
             + 1/result2$summary.hyperpar[1,"0.5quant"])

    temps2 <- list(mean=matrix(reconstruction2$Temp,nrow=n.row),
                   sd=matrix(reconstruction2$SD,nrow=n.row))
    
    ## Convert model parameters to humanly readable form
    result2.field <- inla.spde.result(result2, "field", spde=spde2)
    param2 <- rbind("nugget.sd"=vapply(2:7, function(x) {
      if (x==1) {
        1
      } else if (x==8) {
        0
      } else if (x==2) {
          inla.emarginal(function(x) x^(-0.5),
                                      result2$marginals.hyperpar[[1]])
                     } else if (x==3) {
                       m <- inla.emarginal(function(x) x^(-0.5),
                                           result2$marginals.hyperpar[[1]])
                       inla.emarginal(function(x) (x^(-0.5) - m)^2,
                                      result2$marginals.hyperpar[[1]])^0.5
                     } else {
                       kk <- c(5,4,3,6)[x-3]
                       result2$summary.hyperpar[1,kk]^(-0.5)
                     }
                   }, 1.0))
    param2 <- rbind(param2, "field.range"=vapply(2:7, function(k) {
      if (k==1 | k==8) {
        result2.field$summary.log.range.nominal[1, k]
      } else if (k==2) {
        inla.emarginal(function(x) x,
                       result2.field$marginals.range.nominal[[1]])
      } else if (k==3) {
        m <- inla.emarginal(function(x) x,
                            result2.field$marginals.range.nominal[[1]])
        inla.emarginal(function(x) (x - m)^2,
                       result2.field$marginals.range.nominal[[1]])^0.5
      } else {
        exp(result2.field$summary.log.range.nominal[,k])
      }
    }, 1.0))
    param2 <- rbind(param2,
                   "field.sd"=vapply(2:7, function(k) {
                     if (k==1 | k==8) {
                       result2.field$summary.log.variance.nominal[1, k]
                     } else if (k==2) {
                       inla.emarginal(function(x) x^0.5, result2.field$marginals.variance.nominal[[1]])
                     } else if (k==3) {
                       m <- inla.emarginal(function(x) x^0.5,
                                           result2.field$marginals.variance.nominal[[1]])
                       inla.emarginal(function(x) (x^0.5 - m)^2,
                                      result2.field$marginals.variance.nominal[[1]])^0.5
                     } else {
                       exp(result2.field$summary.log.variance.nominal[,k] / 2)
                     }
                   }, 1.0))
    colnames(param2) <- colnames(result2.field$summary.log.range.nominal)[-c(1,8)]
    param2 <- as.data.frame(param2)
    })
})


## Get a rough kriging&variance time estimate, and sanity check that the posterior
## median paramters give a similar result to the Empirical Bayes estimate.
time3 <- system.time({
  time3.setup <- system.time({
    if (TRUE) {
      ### lognormal prior
      Qx <- inla.spde2.precision(spde2,
                               theta=(c(result2$summary.hyperpar[2,"0.5quant"],
                                        result2$summary.hyperpar[3,"0.5quant"])))
    } else {
      ### pcprior
      Qx <- inla.spde2.precision(spde2,
                                 theta=log(c(result2$summary.hyperpar[2,"0.5quant"],
                                             result2$summary.hyperpar[3,"0.5quant"])))
    }
    ok3 <- !is.na(indata$Temp) ## Include only observed values. Prediction is separate
    covar <- cbind(1,
                   mesh2$loc[idx2[ok3],1,drop=FALSE] - LonCentre,
                   mesh2$loc[idx2[ok3],2,drop=FALSE] - LatCentre)
    Qcovar <- Diagonal(3, c(0, 1e-8, 1e-8))
    Qprior <- rBind(cBind(Qx, Matrix(0, nrow(Qx), ncol(covar))),
                    cBind(Matrix(0, ncol(covar), nrow(Qx)), Qcovar))
    Afield <- sparseMatrix(i=1:sum(ok3), j=idx2[ok3], x=1, dims=c(sum(ok3), spde2$n.spde))
    Aobs <- cBind(Afield, covar)
    Qobs <- Diagonal(sum(ok3), result2$summary.hyperpar[1,"0.5quant"])
    Qpost <- Qprior + t(Aobs) %*% Qobs %*% Aobs
  
    Apred <- cBind(inla.spde.make.A(mesh2, loc2),
                   1,
                   loc2[,1] - LonCentre,
                   loc2[,2] - LatCentre)
  })
  
  time3.mean <- system.time({
    x3.mean <- Apred %*% inla.qsolve(Qpost, as.matrix(t(Aobs) %*% (Qobs %*% indata$Temp[ok3])))
  })
  time3.sd <- system.time({
    ## Split into spare (s) and dense (d) parts:
    ## A S A' = [As Ad] [Sss Ssd; Sds Sdd] [As'; Ad']
    ##        = As Sss As' + As Ssd Ad' + Ad Sds As' + Ad Sdd Ad'
    ## (A S A')_kk = \sum_{i=s,j=s} As_ki Sss_ij As_kj
    ##               +\sum_{i=s,j=d} As_ki Ssd_ij Ad_kj
    ##               +\sum_{i=d,j=s} Ad_ki Sds_ij As_kj
    ##               +\sum_{i=d,j=d} Ad_ki Sdd_ij Ad_kj
    idx.sparse <- 1:mesh2$n
    idx.dense <- mesh2$n + 1:ncol(covar)
    ## Make sure the non-zero pattern covers A' A:
    Spost <- inla.qinv(Qpost + (t(Apred) %*% Apred) * 0)
    ## Calculate (A inverse(Q) A')_kk : 
    As <- Apred[,idx.sparse,drop=FALSE]
    Ad <- Apred[,idx.dense,drop=FALSE]
    AsSs <- As %*% Spost[idx.sparse, idx.sparse, drop=FALSE] 
    AsSd <- As %*% Spost[idx.sparse, idx.dense, drop=FALSE] 
    AdSd <- Ad %*% Spost[idx.dense, idx.dense, drop=FALSE] 
    x3.sd <- (rowSums(AsSs * As) + 2*rowSums(AsSd * Ad) + rowSums(AdSd * Ad))^0.5
  })
  
  reconstruction3 <- data.frame(Lon=indata$Lon,
                                Lat=indata$Lat,
                                Temp=indata$Temp*NA,
                                SD=indata$Temp*NA)
  reconstruction3$Temp <- x3.mean
  reconstruction3$SD <- 
    sqrt(x3.sd^2 + 1/result2$summary.hyperpar[1,"0.5quant"])

  temps3 <- list(mean=matrix(reconstruction3$Temp,nrow=n.row),
                 sd=matrix(reconstruction3$SD,nrow=n.row))
})
  
## Save Output ##
# set.names <- c("SatTemps","SimTemps","TestData")
# write.name <- paste0("./",set.names[dataset],'Fitted.RData')
# save(file=write.name,list=c("reconstruction1","reconstruction2","reconstruction3","indata",
#                             "time1","time3","result1","result2"))



#### Assess predictions #####

# source("assessment.R")
# 
# calc.scores <- function(pred, y, name) {
#   s <- data.frame(c(
#     MSE=mean(sqerr(pred, y)),
#     MAE=mean(abserr(pred, y)),
#     IGN=mean(ign(pred, y)),
#     CRPS=mean(crps(pred, y)),
#     Interval=mean(inter(pred, y)),
#     Coverage=mean(cover(pred, y))))
#   colnames(s) <- name
#   s
# }
# 
# pred1 <- pred.obj(reconstruction1$Temp, reconstruction1$SD)
# pred2 <- pred.obj(reconstruction2$Temp, reconstruction2$SD)
# pred3 <- pred.obj(reconstruction3$Temp, reconstruction3$SD)
# 
# ok.observed <- !is.na(indata$Temp)
# scores.observed <-
#   cbind(calc.scores(pred1[ok.observed,], indata$Temp0[ok.observed], "Trend"),
#         calc.scores(pred2[ok.observed,], indata$Temp0[ok.observed], "TrendField"),
#         calc.scores(pred3[ok.observed,], indata$Temp0[ok.observed], "Postprocessed"))
# scores.unobserved <-
#   cbind(calc.scores(pred1[!ok.observed,], indata$Temp0[!ok.observed], "Trend"),
#         calc.scores(pred2[!ok.observed,], indata$Temp0[!ok.observed], "TrendField"),
#         calc.scores(pred3[!ok.observed,], indata$Temp0[!ok.observed], "Postprocessed"))
# 
# ## Summarise running time
# 
# times <- data.frame("Trend"=time1[1:5], "TrendField"=time2[1:5])
# times.krig <- data.frame("Total"=time3[1:5],
#                          "setup"=time3.setup[1:5],
#                          "mean"=time3.mean[1:5],
#                          "sd"=time3.sd[1:5])
# 
# ####################################################
# ## Model estimate summaries
# 
# print(summary(result1))
# print(summary(result2))
# print(param2)
# 
# ####################################################
# ## Print summary statistics
# 
# print(scores.observed)
# print(scores.unobserved)
# print(times)
# 
# cat(paste("Note: The kriging timings include overhead\n",
#           "      from writing and reading sparse matrices\n",
#           "      to disk, which would not be part of\n",
#           "      a normal inla() run, where the matrices are\n",
#           "      only transferred once, and the symbolic\n",
#           "      Cholesky factorisation is also constructed\n",
#           "      only once.", sep=""))
# ## Times to
# ##   setup the kriging problem matrices (setup)
# ##   compute the kriging point estimate (mean)
# ##   compute the kriging standard deviations (sd)
# ## The calculations needed for the point estimate also
# ##   would provide the log-likelihood for no extra cost
# print(times.krig)
# 
# ####################################################
# 
# 
# pal <- pals::ocean.haline(100)
# dpal <- pals::ocean.balance(100)
# 
# zlim <- range(c(range(temps$mean, na.rm = TRUE),
#                 range(temps1$mean, na.rm = TRUE),
#                 range(temps2$mean, na.rm = TRUE)))
# par(mfrow=c(2,3))
# image.plot(lon,lat,temps$mean,xlab="Longitude",ylab="Latitude",col=pal,asp=1,zlim=zlim)
# image.plot(lon,lat,temps1$mean,xlab="Longitude",ylab="Latitude",col=pal,asp=1,zlim=zlim)
# image.plot(lon,lat,temps2$mean,xlab="Longitude",ylab="Latitude",col=pal,asp=1,zlim=zlim)
# if (any(!is.na(temps00$mean))) {
#   image.plot(lon,lat,temps00$mean,xlab="Longitude",ylab="Latitude",col=pal,asp=1,zlim=zlim)
# } else {
#   plot(mesh2)
# }
# image.plot(lon,lat,temps1$sd,xlab="Longitude",ylab="Latitude",col=pal,asp=1)
# image.plot(lon,lat,temps2$sd,xlab="Longitude",ylab="Latitude",col=pal,asp=1)
# par(mfrow=c(1,1))
# 
# zlim <- c(-1, 1) * max(c(max(abs(temps1$mean - temps0$mean), na.rm = TRUE),
#                          max(abs(temps2$mean - temps0$mean), na.rm = TRUE)))
# par(mfrow=c(1+any(!is.na(temps00$mean)),2))
# if (any(!is.na(temps00$mean))) {
#   zlim <- c(-1, 1) * max(c(max(zlim),
#                            max(abs(temps1$mean - temps00$mean), na.rm = TRUE),
#                            max(abs(temps2$mean - temps00$mean), na.rm = TRUE)))
# }
# image.plot(lon,lat,temps1$mean-temps0$mean,xlab="Longitude",ylab="Latitude",col=dpal,asp=1,zlim=zlim)
# image.plot(lon,lat,temps2$mean-temps0$mean,xlab="Longitude",ylab="Latitude",col=dpal,asp=1,zlim=zlim)
# if (any(!is.na(temps00$mean))) {
#   image.plot(lon,lat,temps1$mean-temps00$mean,xlab="Longitude",ylab="Latitude",col=dpal,asp=1,zlim=zlim)
#   image.plot(lon,lat,temps2$mean-temps00$mean,xlab="Longitude",ylab="Latitude",col=dpal,asp=1,zlim=zlim)
# }
# par(mfrow=c(1,1))
