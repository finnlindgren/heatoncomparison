rm(list=ls())
library(spNNGP)
library(MBA)
library(fields)

coords <- as.matrix(read.table("../data/coords.mod"))
X <- as.matrix(read.table("../data/x.mod"))
y <- as.matrix(read.table("../data/y.mod"))[,1]

coords.ho <- as.matrix(read.table("../data/coords.ho"))
X.ho <- as.matrix(read.table("../data/x.ho"))
y.ho <- as.matrix(read.table("../data/y.ho"))[,1]

cov.model <- "exponential"

sigma.sq <- 6.5

sigma.sq.IG <- c(2, sigma.sq)

g <- 5
theta.alpha <- as.matrix(expand.grid(seq(7, 9, length.out=g), seq(0.00001/sigma.sq, 0.001/sigma.sq, length.out=g)))

colnames(theta.alpha) <- c("phi", "alpha")

m.c <- spConjNNGP(y~X-1, coords=coords, n.neighbors = 15,
                  k.fold = 5, score.rule = "crps",
                  n.omp.threads = 3,
                  theta.alpha = theta.alpha, sigma.sq.IG = sigma.sq.IG,
                  cov.model = cov.model)

m.c$run.time[3]/60

m.c$sigma.sq.hat
tau.sq <- m.c$theta.alpha[2]*m.c$sigma.sq.hat
tau.sq

crps.surf <- mba.surf(m.c$k.fold.scores[,c("phi","alpha","crps")], no.X=100, no.Y=100)$xyz.est
image.plot(crps.surf, xlab="phi", ylab="alpha=tau^2/sigma^2", main="CRPS (lower is better)")
points(m.c$theta.alpha, col="white", pch=2)

##prediction
theta.alpha <- as.vector(m.c$theta.alpha)
names(theta.alpha) <- c("phi", "alpha")

m.p <- spConjNNGP(y~X-1, coords=coords, n.neighbors = 15,
                  X.0=X.ho, coords.0=coords.ho,
                  n.omp.threads = 10,
                  theta.alpha = theta.alpha, sigma.sq.IG = sigma.sq.IG,
                  cov.model = cov.model)

m.p$run.time[3]/60

##total runtime for x-validation and prediction
total.time <- m.c$run.time[3]/60+m.p$run.time[3]/60
total.time

##estimates
round(m.p$sigma.sq.hat,2)

m.p$theta.alpha
tau.sq <- m.p$theta.alpha[2]*m.p$sigma.sq.hat
tau.sq

round(m.p$beta.ha,2)

##confidence interval: (y(s)|Y - hat(y(s)) / sqrt((a-1)/a*Varhat(y(s)) ~ t-distributions with 2a degrees of freedom 

##So, the alpha^th quantile for y(s)|Y will be given by hat(y(s)) + sqrt((a-1)/a*Varhat(y(s))) * t_{2a,alpha}. Then 95% CI will be (2.5 qntl, 97.5 qntl)
a <- m.p$ab[2]
t <- qt(0.975, 2*a)

me <- sqrt((a-1)/a*m.p$y.0.hat.var)*t

pred <- cbind(m.p$y.0.hat, m.p$y.0.hat-me, m.p$y.0.hat+me)

colnames(pred) <- c("50%","2.5%","97.5%")

write.csv(pred, "updated-nngp-conj-real-pred.csv", row.names=F)

##read check
r.pred <- read.csv("../nngp-results/nngp-resp-real-pred.csv", header=TRUE)

c.pred <- read.csv("updated-nngp-conj-real-pred.csv", header=TRUE)

par(mfrow=c(1,3))
plot(c.pred[,1],r.pred[,1])
plot(c.pred[,2],r.pred[,2])
plot(c.pred[,3],r.pred[,3])
