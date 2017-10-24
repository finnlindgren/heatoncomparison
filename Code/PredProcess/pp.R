##export OMP_NUM_THREADS=6

rm(list=ls())
library(spBayes)

coords <- as.matrix(read.table("../../coords.mod"))
X <- as.matrix(read.table("../../x.mod"))
y <- as.matrix(read.table("../../y.mod"))[,1]

coords.ho <- as.matrix(read.table("../../coords.ho"))
X.ho <- as.matrix(read.table("../../x.ho"))
y.ho <- as.matrix(read.table("../../y.ho"))[,1]

n.samples <- 10000

starting <- list("phi"=1.5, "sigma.sq"=6.12, "tau.sq"=0.0002)

tuning <- list("phi"=0.001, "sigma.sq"=0.0001, "tau.sq"=0.0001)

priors <- list("beta.Flat",
               "phi.Unif"=c(0.6, 30), "sigma.sq.IG"=c(2, 5),
               "tau.sq.IG"=c(2, 0.0001))

cov.model <- "exponential"

n.report <- 500
verbose <- TRUE

m.1 <- spLM(y~X-1, coords=coords, starting=starting, knots=c(25, 25),
            tuning=tuning, priors=priors, cov.model=cov.model, 
            n.samples=n.samples, verbose=verbose, n.report=10, modified.pp=TRUE)

##recover beta and spatial random effects
burn.in <- 5000
m.1 <- spRecover(m.1, thin=10, start=burn.in)

round(summary(m.1$p.theta.recover.samples)$quantiles[,c(3,1,5)],2)
round(summary(m.1$p.theta.recover.samples)$quantiles[,c(3,1,5)],4)
round(summary(m.1$p.beta.recover.samples)$quantiles[,c(3,1,5)],2)

set.seed(1)

beta.s <- m.1$p.beta.recover.samples
theta.s <- m.1$p.theta.recover.samples
w.str.s <- m.1$p.wStr.recover.samples
n.sub <- nrow(beta.s)
n.pred <- nrow(coords.ho)

pred.samples <- matrix(0, n.pred, n.sub)

D.knots <- iDist(m.1$knot.coords)
D.knots.pred <- iDist(coords.ho, m.1$knot.coords)

for(s in 1:n.sub){
    pred.samples[,s] <- rnorm(n.pred, X.ho%*%beta.s[s,]+exp(-theta.s[s,"phi"]*D.knots.pred)%*%chol2inv(chol(exp(-theta.s[s,"phi"]*D.knots)))%*%w.str.s[,s], sqrt(theta.s[s,"tau.sq"]))
    print(s)                                                                                                             
}

save.image(file="mod")
##load("mod")

a <- read.csv("../../nngp-results/nngp-resp-sim-pred.csv")
plot(a[,1], apply(pred.samples, 1, mean))

quants <- function(x){quantile(x, prob=c(0.5, 0.025, 0.975))}

pred.q <- apply(pred.samples, 1, quants)

write.csv(file="pp-real-pred.csv", t(pred.q), row.names=F)
