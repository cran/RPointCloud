library(RPointCloud)
library(TDA)
set.seed(54321)
theta <- seq(0, 359, by = 8) * pi/180
L <- length(theta)
x <- cos(theta) + rnorm(L, 0, 0.2)
y <- sin(theta) + rnorm(L, 0, 0.2)
z <- matrix(rnorm(2*L, 1, 0.5), ncol = 2)
joint <- rbind(as.matrix(data.frame(x, y)), z)
joint <- rbind(joint, joint[1:2,]) # deliberately duplicate data rows in longest loop
eucd <- dist(joint)

eucrip <- ripsDiag(eucd, maxdimension=1, maxscale = 3,
                   dist = "arbitrary", library = "Dionysus",
                   location = TRUE)
ripper <- ripsDiag(joint, maxdimension=1, maxscale = 3,
                   dist = "euclidean", library = "Dionysus",
                   location = TRUE)
RD <- ripper$diagram
oneD <- RD[,1] == 1
duration <- RD[oneD, 3] - RD[oneD,2]
w <- which.max(duration)

dexter <- disentangle(ripper, joint)
cyc1 <- Cycle(dexter, 1, 13, "blue") # fails with an error in RPointCloud 0.3.0
## possibly worse, it gives the wrong answer in 0.3.1 - 0.3.4
plot(cyc1, joint) # note the awful picture in 0.3.4!!

cyc2 <- Cycle(eucrip, 1, 13, "purple")
plot(cyc2, joint)

ifelse (all(cyc1@index == cyc2@index), # should be true
        "FIXED",
        "BROKEN")

