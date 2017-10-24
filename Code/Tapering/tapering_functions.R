

plotImage <- function(x, all=TRUE, xlab='Lon', ylab='Lat', ...) {
    # Plotting the data  
    if(!is.matrix(x))
        if (length( x)==n) {
            datamat[] <- x
        } else if (length( x)==nisNA) {
            datamat[isNA] <- x
        } else  if (length( x)==nnotNA) {
            datamat[notNA] <- x
        }  else   image.plot( x, xlab='Lon', ylab='Lat', ...)  # fall back
    if (!all)     datamat[notNA] <- NA
    gapfill::Image(t(datamat), xlab=xlab, ylab=ylab, ...)
}

vgMatrix <- function( x, xlag=1:50, ylag=1:50){
    # Rudimentary estimation of the variogram. Result is a matrix.
    dimx <- dim(x)[1]
    dimy <- dim(x)[2]
    nx <- length(xlag)
    ny <- length(ylag)
    vg <- matrix(0, nx, ny)
    for( i in 1:nx) {
      for( j in 1:ny) {
        xind <- 1:(dimx-xlag[i])
        yind <- 1:(dimy-ylag[j])
        # Here is the actual estimation (Materons's estimator), very rudimentary...
        vg[i,j] <- mean( (x[ xind, yind] - x[ xind+xlag[i], yind+ylag[j]])^2, na.rm=T)
      }
    }
    return(vg)
}

extractDirVg <- function(vgmat, dlon, dlat) {
    # From a 'varigram' matrix, we extract individual directional variograms.
    dvg <- 1:(2*floor( min( dim(vgmat)-1)/2))
    dvg2 <- (1:floor( min( dim(vgmat)-1)/2))
    h <- emp <- list()
    h[[1]] <- dlon[dvg]
    emp[[1]] <- vgmat[dvg,1]
    h[[2]] <- sqrt( dlon[2*dvg2+1]^2+dlat[dvg2+1]^2)
    emp[[2]] <- vgmat[cbind(2*dvg2+1,dvg2+1)]
    h[[3]] <- sqrt( dlon[dvg]^2+dlat[dvg]^2)
    emp[[3]] <- vgmat[cbind(dvg,dvg)]
    h[[4]] <- sqrt( dlon[dvg2+1]^2+dlat[2*dvg2+1]^2)
    emp[[4]] <- vgmat[cbind(dvg2+1,2*dvg2+1)]
    h[[5]] <- dlat[dvg]
    emp[[5]] <- vgmat[1,dvg]
    return( list( lag=h, emp=emp))
}


vg.exp <- function( h, theta)
    # Exponential variogram, analogue to cov.exp, simplified here...
    theta[3] + theta[2] * (1- exp(-h/theta[1]))


fit <- function(theta) {
    # Calculate (unweighted) sums of squares for 'optim'.
    # theta=(range, sill, nugget, anisotropy angle, anisotropy ratio)
    ss <- sum( (emp[[1]]-vg.exp(lag[[1]], theta[1:3]))^2)
    for( i in 2:length(emp))
        ss <- ss + sum( (emp[[i]]-vg.exp(lag[[i]], theta[1:3]) )^2)
    return(ss)
    }
