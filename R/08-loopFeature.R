## 08-loopFeature.R
## Copyright (C) 2023-4 Kevin R. Coombes, RB McGee, and Jake Reed
## LICENSE: Perl Artistic License 2.0

setClass("LoopFeature",
         slots = c(
           data = "matrix",
           fitted = "matrix",
           RSS = "numeric",
           V = "numeric",
           Kstat = "numeric"
         ))

LoopFeature <- function(circMeans) {
  steps <- 360/nrow(circMeans)
  theta <- seq(steps, 360, steps)/180*pi
  eta <- function(param, y) {
    a <- param[1]
    b <- param[2]
    h <- param[3]
    r <- h + a*cos(theta) + b*sin(theta)
    x <- (r-y)^2
    sum(x, na.rm = TRUE)
  }

  fitted <- sapply(colnames(circMeans), function(N) {
    opi <- optim(c(0, 1, 1), eta, y=circMeans[,N])
    fit <- opi$par[3] + opi$par[1]*cos(theta) + opi$par[2]*sin(theta)
    fit
  })

  RSS <- apply((fitted - circMeans)^2, 2, sum, na.rm = TRUE)
  spread <- apply(rbind(circMeans, fitted), 2, function(X) {
    diff(range(X, na.rm = TRUE))
  })
  V <- apply(circMeans, 2, sd, na.rm=TRUE)^2
  Kstat <- (RSS/nrow(fitted))/V
  new("LoopFeature",
      data = circMeans,
      fitted = fitted,
      RSS = RSS,
      V = V,
      Kstat = Kstat
      )
}

setMethod("plot", c("LoopFeature", "missing"), function(x, y, ...) {
  plot(x, colnames(x@fitted))
})


setMethod("plot", c("LoopFeature", "character"), function(x, y, ...) {
  if (missing(y) | (length(y) == 1 && y == "all")) y <- colnames(x@fitted)
  for (N in y) {
    rg <- range(c(x@data[,N], x@fitted[,N]), na.rm = TRUE)
    steps <- 360/nrow(x@data)
    degrees <- seq(0, 345, steps)
    plot(degrees, x@fitted[,N], type = "b", main = N,
         ylim = rg, ylab = "Mean Value")
    points(degrees, x@data[, N], pch = 16)
    title(sub = paste("kappa =", round(x@Kstat[N], 3)))
  }
  invisible(x)
})

setMethod("image", "LoopFeature", function(x, ...) {
  oo <- order(apply(x@fitted, 2, which.max))
  chock <- x@fitted[, oo]
  shock <- scale(chock)
  steps <- 360/nrow(x@data)
  image(seq(0, 345, steps), 1:ncol(shock), shock,
        xlab = "degrees", yaxt = "n", ylab = "", ...)
  mtext(colnames(shock), side = 4, at = 1:ncol(shock), line = 1, las = 2)
  invisible(x)
})


