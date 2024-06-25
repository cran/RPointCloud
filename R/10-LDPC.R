takens <- function(r, dists) {
  k <- sum(dists < r)
  d <- -k / (sum(log(dists)[1:k]) - k*log(r))
  list(k=k, d=d)
}
LDPC <- function(CellID, dset, rg, quorum, samplesAreRows = TRUE) {
  if (!samplesAreRows) {
    return(LDPC(CellID, t(dset), rg, quorum))
  }
  X <- as.vector(as.matrix(dset[CellID, ]))
  SW <- sweep(dset, 2, X, "-")
  eucdist1cell <- sqrt(apply(SW^2, 1, sum))
  sortedDist <- eucdist1cell[order(eucdist1cell)][-1]
  cutter <- t(sapply(rg, takens, sortedDist))
  daft <- data.frame(R = rg, cutter)
  w <- which(daft$k > quorum)[1]
  if(is.na(w)) w <- length(rg)
  b <- daft[w,]
  b1 <- list(R = b$R, k = b$k[[1]], d = b$d[[1]])
  b1
}
