## significance.R
## Copyright (C) 2022-4 Kevin R. Coombes, RB McGee, and Jake Reed
## LICENSE: Perl Artistic License 2.0
## Need to attribute this code to: https://doi.org/10.1038/s41598-023-37842-2

## Weird 1D duration transformation in order to fit a left skewed gumbel distribution
process <- function(diag, A = 1, dimension = 1) {
  P <- diag[which(diag[, 1] == dimension), "Death"] / diag[which(diag[, 1] == dimension), "Birth"]
  P <- log(log(P))
  B <- digamma(1) - mean(P) * A
  P <- A * P + B
  return(P)
}

## Significance test for 1D durations
# the weird input is the 1D durations after transformation
test <- function(weird, alpha = 0.05) {
  x <- seq(0.000001, 0.999999, by = 0.000001)
  dist <- log(-log(x))
  X <- weird
  fn <- ecdf(dist)
  pvals <- 1 - fn(X)
  pv <- p.adjust(pvals, method = "bonferroni")
  res <- pv < alpha
  ret_ls <- list(sum(res == TRUE), pv, pvals)
  names(ret_ls) <- c("sig_loops", "bonferroni", "pvals")
  return(ret_ls)
}

## 95% Confidence Cutoff for 1D duration data
# This must be paired with log transformed births and deaths for the rips plot
con_band <- function(diag, alpha = 0.05, A = 1, dimension = 1) {
  P <- diag[which(diag[, 1] == dimension), "Death"] / diag[which(diag[, 1] == dimension), "Birth"]
  P <- log(log(P))
  M <- length(P)
  B <- digamma(1) - mean(P) * A
  ret <- exp((exp(-B) * log(M/alpha)) ^ (1/A))
  return(log(ret))
}
