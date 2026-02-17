## disentangle.R
## Copyright (C) 2022-4 Kevin R. Coombes, RB McGee, and Jake Reed
## LICENSE: Perl Artistic License 2.0

indexCycles <- function(cycle, dataset, digits = 6) {
  L <- length(dim(cycle))
##  if (L == 2) return(cycle)
  ## Now we have to convert the data-set coordinates to indices
  N1 <- ncol(dataset)
  N2 <- dim(cycle)[3]
  if (N1 != N2) {
    stop("Mismatched sizes. Dataset = ", N1, " and cycle = ", N2, "\n")
  }
  input <- apply(round(dataset, digits), 1, paste, collapse = "|")
  dex <- apply(cycle, c(1, 2), function(X) {
    check <- paste(round(X, digits), collapse = "|")
    w <- which(input == check)
    if (length(w) > 1) w <- w[1]
    w
  })
  return(dex)
}

lookup <- function(loc, dataset, digits = 6) {
  input <- apply(round(dataset, digits), 1, paste, collapse = "|")
  apply(loc, 1, function(L) {
    1.0*which(input == paste(round(L, digits), collapse = "|"))
  })
}


disentangle <- function(rips, dataset, digits = 6) {
  cl <- rips$cycleLocation
  fixed <- lapply(cl, indexCycles, dataset = dataset, digits = digits)
  rips$cycleLocation <- fixed
  rips$birthLocation <- lookup(rips$birthLocation, dataset, digits = digits)
  rips$deathLocation <- lookup(rips$deathLocation, dataset, digits = digits)
  return(rips)
}
