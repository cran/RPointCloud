library(RPointCloud)
data(cytof)
diag <- AML10.node287.rips[["diagram"]]
persistence <- diag[, "Death"] - diag[, "Birth"]
d1 <- persistence[diag[, "dimension"] == 1]
eb <- EBexpo(d1, 200)
hist(eb)
plot(eb, prior = 0.77)
abline(h=0, col= "gray")
cutoffSignificant(eb, 0.77, 0.8)
