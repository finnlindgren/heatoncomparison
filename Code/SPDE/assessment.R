pred.obj <- function(mean, sd) {
    obj <- data.frame(mean=as.numeric(mean), sd=as.numeric(sd))
    class(obj) <- c("pred.obj", "data.frame")
    obj
}


ign <- function(x, ...) {
    UseMethod("ign")
}
crps <- function(x, ...) {
    UseMethod("crps")
}
sqerr <- function(x, ...) {
  UseMethod("sqerr")
}
abserr <- function(x, ...) {
  UseMethod("abserr")
}
inter <- function(x, ...) {
  UseMethod("inter")
}
cover <- function(x, ...) {
  UseMethod("cover")
}
pit <- function(x, ...) {
    UseMethod("pit")
}
ign.pred.obj <- function(x, y, ...) {
    ## Compute normalised prediction residuals
    z <- as.numeric((y - x$mean) / x$sd)
    ## Compute the log-score:
    scores <- z^2 / 2 + log(x$sd)
    attr(scores, "score") <- mean(scores, na.rm=TRUE)
    scores
}
sqerr.pred.obj <- function(x, y, ...) {
  ## Compute prediction residuals
  z <- as.numeric((y - x$mean))
  ## Compute the sqerr-score:
  scores <- z^2
  attr(scores, "score") <- mean(scores, na.rm=TRUE)
  scores
}
abserr.pred.obj <- function(x, y, ...) {
  ## Compute prediction residuals
  z <- as.numeric((y - x$mean))
  ## Compute the absolute error score:
  scores <- abs(z)
  attr(scores, "score") <- mean(scores, na.rm=TRUE)
  scores
}
crps.pred.obj <- function(x, y, ...) {
    ## Compute normalised prediction residuals
    z <- as.numeric((y - x$mean) / x$sd)
    ## Compute the crps-score:
    scores <- x$sd * (z *(2 * pnorm(z, 0, 1) - 1) +
                      2 * dnorm(z, 0, 1) - 1/sqrt(pi))
    attr(scores, "score") <- mean(scores, na.rm=TRUE)
    scores
}
inter.pred.obj <- function(x, y, alpha=0.05, ...) {
  ## Compute the Interval score:
  hw <- -qnorm(alpha/2) * x$sd
  scores <- 2 * hw + (2/alpha) * (((x$mean - hw) - y) * (y < x$mean - hw) +
                                    (y - (x$mean + hw)) * (y > x$mean + hw))
  attr(scores, "score") <- mean(scores, na.rm=TRUE)
  scores
}
cover.pred.obj <- function(x, y, alpha=0.05, ...) {
  ## Compute the coverage:
  hw <- -qnorm(alpha/2) * x$sd
  scores <- (x$mean - hw <= y) & (y <= x$mean + hw)
  attr(scores, "score") <- mean(scores, na.rm=TRUE)
  scores
}
pit.pred.obj <- function(x, y, ...) {
    ## Compute normalised prediction residuals
    z <- as.numeric((y - x$mean) / x$sd)
    ## Compute the P(Y <= y) probabilities:
    scores <- pnorm(z, 0, 1)
    attr(scores, "score") <- NA
    scores
}
