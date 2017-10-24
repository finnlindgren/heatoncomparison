library(LatticeKrig)
library(parallel)

fit.Matern.GP <- function(y,X,s,nu){
  
  ## Create a Sequence for Spatial Range
  D <- rdist(s)
  max.dist <- max(D)
  sr.seq <- seq(0.001,max.dist,length=10)
  for(i in 1:length(sr.seq)){
    sr.seq[i] <- 1/Matern.cor.to.range(sr.seq[i],nu=nu,cor.target=0.05)
  }
  
  ## Create a Sequence for %Spatial
  pct.spatial <- seq(0,.95,length=10)
  
  ## log-likelihoods
  ll <- matrix(0,nrow=length(sr.seq),ncol=length(pct.spatial))
  s2.hats <- ll
  beta.hats <- matrix(0,nrow=ncol(X),ncol=length(sr.seq)*length(pct.spatial))
  the.it <- 0
  for(pct in pct.spatial){
    for(sr in sr.seq){
      the.it <- the.it + 1
      R <- pct*Matern(D,nu=nu,alpha=sr)+(1-pct)*diag(nrow(s))
      R.chol <- t(chol(R))
      first.piece <- forwardsolve(R.chol,X)
      XpRinvX <- t(first.piece)%*%first.piece
      last.piece <- forwardsolve(R.chol,y)
      beta.hat <- solve(XpRinvX)%*%t(first.piece)%*%last.piece
      beta.hats[,the.it] <- beta.hat
      ss <- forwardsolve(R.chol,y-X%*%beta.hat)
      ss <- t(ss)%*%ss
      s2.hat <- ss/nrow(s)
      s2.hats[sr==sr.seq,pct==pct.spatial] <- s2.hat
      log.like <- -0.5*nrow(s)*log(s2.hat)-sum(log(diag(R.chol)))-0.5*ss/(s2.hat)
      ll[sr==sr.seq,pct==pct.spatial] <- log.like
    }
  }
  
  ## Find maximum likelihood estimate
  the.mle <- which(ll==max(ll),arr.ind=TRUE)
  sr <- sr.seq[the.mle[1,1]]
  pct <- pct.spatial[the.mle[1,2]]
  the.mle <- which(ll==max(ll))
  s2.hat <- s2.hats[the.mle]
  beta.hat <- beta.hats[,the.mle]
  
  ## Return List
  return(list(alpha=sr,sigma2=pct*s2.hat,tau2=(1-pct)*s2.hat,beta.hat=beta.hat,log.like=ll[the.mle]))
  
}

## Test to see if fit.Matern.GP works
# om <- 0.95
# s2 <- 1
# mu <- 0
# nu <- 1/2
# s <- runif(250,0,1)
# D <- rdist(s)
# V <- s2*(om*Matern(D,alpha=3/.5,nu=nu)+(1-om)*diag(length(s)))
# y <- mu+t(chol(V))%*%rnorm(length(s))
# plot(s,y,pch=19)
# tst <- fit.Matern.GP(y,matrix(1,nrow=length(s),ncol=1),matrix(s,ncol=1),nu)
# tst

fit.Matern.MC <- function(y,X,s,nu,num.cores=1){
  n <- nrow(X)
  D <- rdist(s)
  max.dist <- max(D)
  dist <- seq(max.dist/100000,max.dist,length=10)
  alpha <- sapply(1:length(dist), function(i) 1/Matern.cor.to.range(dist[i],nu=nu,cor.target=0.05))
  a2M <- function(a){Matern(D,alpha=a,nu=nu)}
  M.list <- mclapply(alpha,a2M,mc.cores=num.cores)
  omega <- seq(0,0.99,length=20)
  aMw.list <- vector('list',length(M.list)*length(omega))
  i <- 1
  for(a in 1:length(alpha)){
    for(w in 1:length(omega)){
      aMw.list[[i]] <- list(alpha[a],M.list[[a]],omega[w])
      i <- i+1
    }
  }
  aMw <- aMw.list[[1]]
  aMw2lik <- function(aMw){
    M <- aMw[[2]]
    w <- aMw[[3]]
    R <- w*M+(1-w)*diag(n)
    R.chol <- t(chol(R))
    first.piece <- forwardsolve(R.chol,X)
    XpRinvX <- t(first.piece)%*%first.piece
    last.piece <- forwardsolve(R.chol,y)
    B_hat <- solve(XpRinvX)%*%t(first.piece)%*%last.piece
    ss <- forwardsolve(R.chol,y-X%*%B_hat)
    sigma2_hat <- t(ss)%*%ss/n
    ll <- -n/2*log(sigma2_hat)-sum(log(diag(R.chol)))-0.5*sum(ss^2)/sigma2_hat
    return(c(B_hat,sigma2_hat,ll))
  }
  B_hat.sigma2.ll.list <- mclapply(aMw.list,aMw2lik,mc.cores=num.cores)
  B_hat <- sapply(1:length(B_hat.sigma2.ll.list), function(t) B_hat.sigma2.ll.list[[t]][1])
  sigma2_hat <- sapply(1:length(B_hat.sigma2.ll.list), function(t) B_hat.sigma2.ll.list[[t]][2])
  ll <- sapply(1:length(B_hat.sigma2.ll.list), function(t) B_hat.sigma2.ll.list[[t]][3])
  
  max.ll <- which(ll==max(ll))
  B_MLE <- B_hat[max.ll]
  sigma2_MLE <- sigma2_hat[max.ll]
  alpha_MLE <- aMw.list[[max.ll]][[1]]
  omega_MLE <- aMw.list[[max.ll]][[3]]
  list(alpha=alpha_MLE,tau2=(1-omega_MLE)*sigma2_MLE,beta.hat=B_MLE,sigma2=omega_MLE*sigma2_MLE,log.like=ll[max.ll])
}
	