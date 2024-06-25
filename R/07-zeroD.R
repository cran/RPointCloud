## zeroD.r
## Copyright (C) 2022-4 Kevin R. Coombes and Jake Reed
## LICENSE: Perl Artistic License 2.0

gammaSSE <- function(parm, X0, pdf) {
  alpha <- parm[1]
  beta <- parm[2]
  theo <- dgamma(X0, shape = alpha, rate = beta)
  SSE <- (theo - pdf)^2
  sum(SSE)
}

gammaFit <- function(edata, resn = 200) {
  ## Use the method-of-moments to get rough parameters
  mu <- mean(edata)
  sigma <- sd(edata)
  shape <- mu^2/sigma^2
  rate <- mu/sigma^2
  h0 <- hist(edata, breaks = 123, plot = FALSE)
  X0 <- h0$mids
  pdf <- h0$density
  ## Compute the crude recip parameter estimate
  oo <- optim(c(shape, rate), gammaSSE, X0 = X0, pdf = pdf)
  lambda <- oo$minimum
  val <- list(X0 = X0, pdf = pdf, mu = mu, h0 = h0,
              edata = edata, alpha = oo$par[1], beta = oo$par[2])
  val
}

