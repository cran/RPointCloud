## higherD.r
## Copyright (C) 2022-4 Kevin R. Coombes and Jake Reed
## LICENSE: Perl Artistic License 2.0

expoSSE <- function(lambda, X0, pdf) {
  theo <- dexp(X0, lambda)
  weird <- theo - pdf
  sum(weird^2)
}

getDelta <- function(mu) {
  digits <- 0.15*mu
  val <- round(c(mu-digits, mu+digits))
  if(val[1] < 0) val[1]  <- 0
  if(val[2] < 1) val[2] <- 1
  val
}

ExpoFit  <- function(edata, resn = 200) {
  if (min(edata) < 0) {
    stop("Exponential duistributions cannot yield negative values.\n")
  }
  ## Use 'hist' to get the empirical density function
  bks <- seq(0, max(edata), length = resn)
  h0 <- hist(edata, breaks=bks, plot = FALSE)
  X0 <- h0$mids
  pdf <- h0$density

  ## Compute the crude recip parameter estimate
  mu  <- mean(edata)    # mean
  delta <- getDelta(1/mu)
  recip = 1/mu # rate is appromixately mean reciprocal
  oo <- optimize(expoSSE, delta, X0 = X0, pdf = pdf)
  lambda <- oo$minimum
  val <- new("ExpoFit",
             X0 = X0, pdf = pdf, mu = mu, h0 = h0,
             edata = edata, lambda = lambda)
  val
}

setOldClass("histogram")
setClass("ExpoFit", slots = c(X0 = "numeric",
                              pdf = "numeric",
                              mu = "numeric",
                              h0 = "histogram",
                              edata = "numeric",
                              lambda = "numeric"))

setMethod("plot", c("ExpoFit", "missing"), function(x, y, ...) {
  delta <- getDelta(recip <- 1/x@mu)
  ss <- seq(delta[1], delta[2], length=300)
  tt <- sapply(ss, expoSSE, X0 = x@X0, pdf = x@pdf)
  plot(ss, tt, type='l', xlab="Lambda", ylab="SS_Error", lwd=2)
  points(x@lambda, expoSSE(x@lambda, x@X0, x@pdf),
         pch=16, col="red", cex=1.5)
  points(recip, expoSSE(1/x@mu, x@X0, x@pdf),
         pch=16, col="blue", cex=1.5)
  legend("top", c("lambda", "1/mu"), col=c("red", "blue"), pch=16)
  invisible(x)
})

EBexpo <- function(edata, resn = 200) {
  if (inherits(edata, "ExpoFit")) {
    object <- edata
  } else {
    object <- ExpoFit(edata, resn)
  }
  ## fit the model
  X0 <- object@X0
  theo <- dexp(X0, object@lambda)
  weird <- theo - object@pdf
  Y0 <- log(object@pdf) - log(theo)
  click <- !is.infinite(Y0)
  Y <- Y0[click]
  X <- X0[click]
  Z <- lm(Y ~ bs(X, df=5))
  YP <- predict(Z, data.frame(X=X0)) # gives a warning that we are currently ignoring
  unravel <- exp(YP + log(theo))
  val  <- new("EBexpo",
              expo = object,
              theoretical.pdf = theo,
              unravel = unravel)
  val
}

setClass("EBexpo", slots = c(expo = "ExpoFit",
                             theoretical.pdf = "numeric",
                             unravel = "numeric"))

setMethod("hist", "EBexpo", function(x, xlab="", ylab="Prob(Interesting | X)", main="", ...) {
  top <- max(c(x@unravel, x@theoretical.pdf))
  xvals <- x@expo@X0
  hist(x@expo@edata, probability=TRUE, breaks=100, ylim=c(0, top),
       xlim=c(min(xvals), max(xvals)), xlab=xlab, main=main)
  lines(xvals, x@theoretical.pdf, col="red", lwd=2)
  lines(xvals, x@unravel, col="blue", lwd=2)
  legend("topright", c('Empirical', 'Theoretical'),
         col=c("blue", "red"), lwd=2)
  invisible(x)
})


### posterior probability of difference (p0 = prior)
probDiff <- function(p0, object) {1 - p0*object@theoretical.pdf/object@unravel}

### maximum value of p0 that keeps posterior non-negative
estPrior <- function(object, minp, maxp, resn = 1000) {
  sap <- sapply(sip <- seq(minp, maxp, length=resn),
                function(x) min(probDiff(x, object)))
  prior <- sip[max(which(sap >= 0))]
  prior
}

cutoff <- function(target, prior, object) {
  pd <- probDiff(prior, object)
  X0  <- object@expo@X0
  X0[min(which(pd > target))]
}
setMethod("plot", c("EBexpo", "missing"), function(x, prior, post = c(0.5, 0.8, 0.9), ...) {
  cuts <- sapply(post, cutoff, prior = prior, object = x)
  f1 <- function(i) {
    lines(c(0, cuts[i], cuts[i]),
          c(post[i], post[i], 0), col="gray")
  }
  pd <- probDiff(prior, x)
  plot(x@expo@X0, pd, type='l', lwd=2, xaxs='i',
       xlab="Duration",
       ylab="Probability(Interesting | Duration)")
  abline(h = c(0, 1))
  sapply(1:length(post), f1)
  invisible(x)
})

### Extract a sub-histogram from the tail of the distibution
tailHisto <- function(H, target) {
  L <- length(H$breaks)
  L1 <- L - 1
  start <- which(H$breaks > target)
  if (length(start) == 0) {
    stop("The 'target' is bigger than all breaks.\n")
  }
  start <- start[1]
  awkward <- structure(list(breaks = H$breaks[start:L],
                            counts = H$counts[start:L1],
                            density = H$density[start:L1],
                            mids = H$mids[start:L1],
                            xname = H$xname,
                            equidist = H$equidist),
                       class="histogram")
  awkward
}
