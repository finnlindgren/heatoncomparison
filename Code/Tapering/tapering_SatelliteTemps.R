args <- commandArgs(TRUE)

if (length(args)==0) { 
  
  rm(list = ls())
  delta <- 0.2     # taper range    
  verbose <- list(print=TRUE, plot=FALSE)
  
} else {
  
  verbose <- list(print=TRUE, plot=FALSE)
  for(i in 1:length(args)){
    eval(parse(text=args[[i]]))
  }
  
}
cat( '\n\n\nTaper range: ',delta,'\n\n')

timing <- list(`Setting up structure and load data` = system.time({
    
    source('tapering_functions.R')

    # load CRAN packages and possibly other elements 
    library("spam")
    library("fields")
    library("gapfill")      # for visualization:

    set.seed( 14)
    taperpar <- c(98121332,29962933,29131276,532192428,1226520061)  # taper range theta=0.2


    if(verbose$print) {
        cat( '\n\n\nTaper range: ',delta,'\n\n')
        print( args <- commandArgs())
        print(taperpar)
        cat( '\n\n\n')
    }

    load("../../Data/SatelliteTemps.RData")  # Load data and prepare variables.
    Lat <- unique(sat.temps$Lat)
    Lon <- unique(sat.temps$Lon)
    nx <- length( Lon)
    ny <- length( Lat)
    n <- nx*ny
    
    datamat <- array(sat.temps$Temp, c(nx, ny))  # Arrange in matrix form for better display
                                                #  dim(datamat)  is   500 300
    notNA <- !is.na(sat.temps[,3])
    isNA <- !notNA
    nnotNA <- sum(notNA)
    nisNA <- n-nnotNA
    
    xobs <- sat.temps[ notNA, 1:2]
    xpred <- sat.temps[ isNA, 1:2]
    Y <- sat.temps[ notNA, 3]
    
})[1]  )

timing$`Estimation of trend`= system.time({

  # Simple model for the mean. Easy to extend with `poly` and `lm` construction:
    X <- as.data.frame( poly( scale( sat.temps[,1:2]), degree=2, raw=TRUE) )
    drift <- lm( Y ~ `1.0`+`0.1`, data=cbind(Y, X[ notNA,]))
    mupred <- predict( drift, X)
    
    if(verbose$print) 
        summary( drift)

    if (verbose$plot) {
        plotImage( datamat)
        quilt.plot( sat.temps[ notNA,-4], nx=nx, ny=ny)
        plotImage( mupred)
    }
    
})[1]

######################################################################################
# estimation
timing$`Estimation of covariance structure`= system.time({
    
    residmat <- datamat
    residmat[notNA] <- drift$resid
    
    
    vg <- vgMatrix(residmat,xlag=1:20,ylag=1:21)
    
    if(verbose$plot)  image.plot(vg)
#

    dlat <- as.vector( nearest.dist( cbind(Lon[1], Lat), cbind(Lon[1], Lat[1]), miles=FALSE, delta=10))[-1]
# length 300
    dlon <- as.vector( nearest.dist( cbind(Lon, Lat[1]), cbind(Lon[1], Lat[1]), miles=FALSE, delta=10))[-1]
# length 500


# We assume different ranges (geometric anisotropy), determine the different scalings:
    empirical <- extractDirVg( vg, dlon, dlat)
    emp <- empirical$emp[c(1,3,5)]
    lag <- empirical$lag[c(1,3,5)]
    
    
    res <- optim( c(.05,8,.1), fit, method = "L-BFGS-B", 
                  lower=c(.025, 3, 0), 
                  upper=c(.3,  14, 1),  hessian=TRUE)
   # Fitting is suboptimal. Weighting might improve the situation a bit.
      
   if(verbose$print) {
        cat('\n\n\n')
        print(res)
        cat('\n\n\n')
    }
    
    if (verbose$plot) {
        plot( lag[[1]], emp[[1]], ylim=c(0,max(unlist(emp))), type='l')
        for (i in 2:3)   lines( lag[[i]], emp[[i]], col=i) 
        
        h <- seq(0,to=1.7,l=100)
        lines( h, vg.exp( h, res$par[1:3]), col=2, lwd=2)
    }

    
} )[1]

######################################################################################

timing$`Construction of Sigma` <- system.time({
    spam.options(nearestdistnnz=c(taperpar[1], 400))   
    hobs <- nearest.dist( xobs, miles=FALSE, delta=delta)
    
    
    if( verbose$print) {
        cat('\nDistance matrix (upper only):')
        shobs <- summary(hobs)
        print( c( log(shobs$nnz, 2)))
    }
    
    hobs <- hobs + t(hobs)
    Cobsobs <- cov.exp( hobs, theta=c( res$par[1:3])) * cov.wend1( hobs, theta=c( delta,1,0))
                       
})[1]

timing$`Cholesky decomposion` <- system.time({
    
    cholCobsobs <- chol( Cobsobs, memory=list(nnzR=taperpar[4]))
    if( verbose$print) {
        cat('\nCholesky factor:')
        scholCobs <- summary(cholCobsobs)
        print( c( log(scholCobs$nnzR,2)))
    }
    
})[1]


timing$`Prediction` <- system.time({ 
    spam.options(nearestdistnnz=c(taperpar[2], 400))   
    hpredobs <- nearest.dist( xpred, xobs, method='eucli', miles=FALSE, delta=delta)
    if( verbose$print) {
        cat('\nDistance matrix (pred-obs):')
        shpredobs <- summary(hpredobs)
    }
    
    Cpredobs <-  cov.exp( hpredobs, theta=c( res$par[1:3])) * cov.wend1( hpredobs, theta=c( delta,1,0))
    Ypred <- c( Cpredobs %*% solve.spam( cholCobsobs, drift$resid) +  mupred[ isNA])
    
})[1]

timing$`Uncertainty calculation` <- system.time({ 
                                  
    M <- 250  # some large number 
    
    spam.options(nearestdistnnz=c(taperpar[3], 400))
    hpred <- nearest.dist( xpred, method='eucli', miles=FALSE, delta=delta)
    if( verbose$print) {
        cat('\nDistance matrix (pred; upper only):')
        shpred <- summary(hpred)
    }
    hpred <- hpred + t(hpred)
    
    Cpredpred <- cov.exp( hpred, theta=c( res$par[1:3])) * cov.wend1( hpred, theta=c( delta,1,0))
    
    Call <- rbind( cbind(Cobsobs, t(Cpredobs)), cbind( Cpredobs, Cpredpred))
    
    cholCall <- chol( Call, memory=list(nnzR=taperpar[5]))
    
    if( verbose$print){
        cat('\nCholesky factor:')
        scholCall <- summary( cholCall)
        print( c( log(scholCall$nnzR,2)))
    }
    
    samples <- t( rmvnorm.spam( M, Sigma=Call, Rstruct=cholCall))

    Yconditional <- samples[(nnotNA+1):n,] - Cpredobs %*% solve.spam( cholCobsobs, samples[1:nnotNA,]) # + mupred[ isNA]   # last is not necessary here.
    
    preduncertainty <- apply( Yconditional, 1, sd)
    
    predCIlower <- Ypred + apply( Yconditional, 1, quantile, probs=0.025)
    predCIupper <- Ypred + apply( Yconditional, 1, quantile, probs=0.975)
    
})[1]

if (verbose$print){
    cat( '\n\n\n')
    print(timing)
    cat( '\n\n\nTotal time: ', sum( unlist(timing)),'s, ',
        sum( unlist(timing))/3600,'h.\n\n\n')
}

if (verbose$plot) {
    plotImage( Ypred, all=F)

    # quilt.plot( rbind(xobs,xpred),apply( samples,1,mean))
    # quilt.plot( rbind(xobs,xpred),apply( samples,1,sd))
    quilt.plot( xpred, preduncertainty)

    plotImage( apply( Yconditional, 1, mean) )
    plotImage(predCIlower, all=FALSE )
    plotImage(predCIupper, all=FALSE )
    plotImage(predCIupper-predCIlower, all=FALSE )
    plotImage((predCIupper-predCIlower)/preduncertainty, all=FALSE )

}



## summary images
jpeg(paste0("figs/SatelliteTemps_fit_",delta,".jpeg"), 1000, 700)
plotImage( Ypred)
dev.off()

jpeg(paste0("figs/SatelliteTemps_SE_",delta,".jpeg"), 1000, 700)
plotImage( preduncertainty, all=FALSE)
dev.off()

jpeg(paste0("figs/SatelliteTemps_summary_",delta,".jpeg"), 1000, 700)
hist( Ypred, n=120)
dev.off()



## save results to file
save(Ypred, preduncertainty, res, timing, delta, verbose,
     Yconditional, predCIlower, predCIlower,
     shobs, scholCobs, shpred, scholCobs, scholCall,
     file=paste0("./SatelliteTemps_tapering_",delta,".RData"))

sessionInfo()
