library(RPointCloud)
library(TDA)
set.seed(54321)
theta <- seq(0, 359, by = 8) * pi/180
L <- length(theta)
x <- cos(theta) + rnorm(L, 0, 0.2)
y <- sin(theta) + rnorm(L, 0, 0.2)
z <- matrix(rnorm(2*L, 1, 0.5), ncol = 2)
joint <- rbind(as.matrix(data.frame(x, y)), z)
ripper <- ripsDiag(joint, maxdimension=1, maxscale = 3,
                   dist = "euclidean", library = "Dionysus",
                   location = TRUE)
RD <- ripper$diagram
gum <- RPointCloud:::process(RD)
tack <- RPointCloud:::test(gum)
cb <- RPointCloud:::con_band(RD)

set.seed(97531)
jag <- jitter(joint)
ripped <- ripsDiag(jag, maxdimension=1, maxscale = 3,
                   dist = "euclidean", library = "Dionysus",
                   location = TRUE)
RD2 <- ripped$diagram
gum <- RPointCloud:::process(RD)
tack <- RPointCloud:::test(gum)
cb2 <- RPointCloud:::con_band(RD2)

exp(cb)
exp(cb2)

opar <- par(mfrow=c(1,2))
plot(RD, band = cb)
plot(RD2, band = cb2)
par(opar)
