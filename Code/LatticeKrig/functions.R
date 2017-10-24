 
#### creates a list that has the different parts of the data
#### doTask operates on this. 
makeData<- function(lon,lat,temps){
  #full set of locations
  xFull<- cbind( c(lon), c(lat))
  xGrid<-  sort(unique( lon))
  yGrid<- sort( unique( lat))
  # find locations that do not have missing values  
  ind<- !is.na( temps)
  x<-xFull[c(ind),]
  y<- temps[ind]
  # these are locations to predict
  xMissing<- xFull[c(!ind),]
  #  divide up sample into  10/90 for cross-validation 
  set.seed( 122)
  N<- length( y)
  outSample<- sample( 1:N, round(.1*N), replace=FALSE)
  # fit to xIn, yIn  validate by predicting to xCV
  yIn<- y[-outSample]
  xIn<- x[-outSample,]
  xCV<- x[outSample,]
  yCV<- y[outSample]
  return( 
    list( x = x,
          y = y, 
        yIn = yIn,
        xIn = xIn, 
        xCV = xCV, 
        yCV = yCV, 
   xMissing = xMissing,
      xGrid = xGrid,
      yGrid = yGrid)
   )
}

doTask<- function(k,dataObject,batchTable){
  # this function has a call amenable to using Rmpi
  # (e.g. ID is k that points to pars in a table)
  # instead of a for loop.
  require(LatticeKrig)
  NCTemp<- batchTable$NC[k]
  nlevelTemp<- batchTable$nlevel[k]
  a.wghtTemp<- batchTable$a.wght[k]
  nuTemp<- batchTable$nu[k]
  tick<-  proc.time()[3]
  # this call to LatticeKrig will find the MLE over 
  # the sill and nugget variance parameters
  # for fixed a.wght (equivalent to the range parameter)
  out<- LatticeKrig( dataObject$xIn,dataObject$yIn,
                     NC=NCTemp,
                     nlevel=nlevelTemp,
                     a.wght=a.wghtTemp,
                     nu= nuTemp
  )
  tock<-  proc.time()[3]
  timeLK<- tock- tick  
  # predict at omitted points 
  predictCV<- predict( out, dataObject$xCV)
  MSECV<- mean( (predictCV - dataObject$yCV)^2)
  # output list with summaries of results  
  summaryList<- list( rhoMLE = out$rho.MLE.FULL,
                   lambdaMLE =  out$lambda,
                    sigmaMLE = out$sigma.MLE.FULL,
                  lnProfLike = out$lnProfileLike.FULL,
                       effDf = out$eff.df
               )
  temp<-  
    list(      ID = k, 
               pars = batchTable[k,],
               summary = summaryList,
               MSECV   = MSECV,
               timeLK = timeLK
    )
  cat( "MSECV: ", MSECV, "time: ", timeLK, fill=TRUE )
  
  #############################################################
  # May want to omit this output block.
  # This is for debugging and to recover if the for loop is interupted.
  # save this result as R binary file
  # when these each take an hour good idea to save them
  # if loop is interupted
  #  ID<- substring(format(k+1000),2,4) # IDs will alpha sort nicely
  #  outName<- paste0("ResultsTask",ID,".rda")
  #  cat( "Saving ID = " , k ," to ", outName, fill=TRUE) 
  #  save( temp, file = outName)  
  #############################################################
  return( temp)
}


#### completely optional checks
###  But as Reagan said of the Russians: trust but verify
checkData<-function(dataObject){
  
  # testing:
  indX<- cbind( match( dataObject$xIn[,1], dataObject$xGrid ),
                match( dataObject$xIn[,2], dataObject$yGrid )
  )
  indCV<- cbind( match( dataObject$xCV[,1], dataObject$xGrid ),
                 match( dataObject$xCV[,2], dataObject$yGrid )
  )
  indMissing<- cbind( match( dataObject$xMissing[,1], dataObject$xGrid ),
                      match( dataObject$xMissing[,2], dataObject$yGrid )
  )             
  hold<- matrix( 0, length( dataObject$xGrid), length( dataObject$yGrid)) 
  hold[indX]<- 1
  hold[indCV]<- 2 + hold[indCV]
  hold[indMissing]<- 4 + hold[indMissing]
  #image.plot( xGrid,yGrid, hold, col=rainbow(3))
  return(
    table( c( hold))
  )
}
# 1    2    4 
# 6390  710 2900 
makeFig<- function(){
    set.panel( 2,2)
    indX<- cbind( match( dataObject$x[,1], dataObject$xGrid ),
                  match( dataObject$x[,2], dataObject$yGrid )
    )
    m1<- length( dataObject$xGrid)
    m2<- length( dataObject$yGrid)
    Z<- matrix( NA, m1, m2)
    ##### data
    Z[indX]<- dataObject$y
    zRange<- range(dataObject$y)
    image.plot( Z, zlim = zRange)
    title("Raw data field")
    ##### residuals
    Z[indX]<- fitFinal$residuals
    image.plot( Z)
    title("Residual field")
    
    ##### predicted
    Z<- matrix( NA, m1, m2)
    indX<- cbind( match( dataObject$xMissing[,1], dataObject$xGrid ),
                  match( dataObject$xMissing[,2], dataObject$yGrid )
    )
    Z[indX]<- outputSim$ghat
    image.plot( Z, zlim = zRange)
    title("predicted")
    ##### SE
    Z[indX]<- standardError
    image.plot( Z)
    title("standard error")
}    


findCov<- function( LKinfo, NG=100){
  grid.info<- LKinfo$latticeInfo$grid.info
  NCx<- (grid.info$xmax - grid.info$xmin)/grid.info$delta + 1
  NCy<- (grid.info$ymax - grid.info$ymin)/grid.info$delta + 1
  NC<- max( NCx, NCy)
#
# if this is big model don't compute over whole domain to get a localized covariance function.
# this will be an approximation but seems adequate.   
  if( NC > 64){
  NC0 <- 64
  scaleNC<- (NC0 - 1)/(NC - 1)
  }
  else{
    scaleNC<- 1.0
  }
#  
  
  x0<- colMeans( grid.info$range)
  rX<- scaleNC*diff(grid.info$range[,1] )/2
  rY<- scaleNC*diff(grid.info$range[,2] )/2
  xLoc0<-    cbind( c( -rX, rX) + x0[1],
                    c( -rY, rY) + x0[2]
                    )

  LKinfo0<- LKrigSetup( xLoc0,
                        nlevel = LKinfo$nlevel,
                        NC = NC,
                        a.wght = LKinfo$a.wght,
                        nu = LKinfo$nu)
  

d0<- seq( 0, rX, length.out= NG)
xGrid<- cbind( x0[1] + d0, rep( x0[2],NG) )
  
  hold0<- c(
    LKrig.cov( rbind(x0), xGrid , LKinfo0)
  )
return( list( d=d0, cov=hold0, xGrid= xGrid))
}
