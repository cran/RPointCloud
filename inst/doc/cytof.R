## ----opts, echo=FALSE-------------------------------------------------------------------------
knitr::opts_chunk$set(fig.width=5, fig.height=5)
oopt <- options(width=96)
.format <- knitr::opts_knit$get("rmarkdown.pandoc.to")
.tag <- function(N, cap ) ifelse(.format == "html",
                                 paste("Figure", N, ":",  cap),
                                 cap)

## ----RPointCloud------------------------------------------------------------------------------
library(RPointCloud)

## ----libpack----------------------------------------------------------------------------------
suppressMessages( library(Mercator) )
library(ClassDiscovery)
library(Polychrome)
data(Dark24)
data(Light24)
suppressMessages( library(igraph) )
suppressMessages( library("ape") )
suppressPackageStartupMessages( library(circlize) )

## ----cytof------------------------------------------------------------------------------------
data(cytof)
ls()
dim(AML10.node287)
colnames(AML10.node287)
amldist <- dist(AML10.node287)

## ----fig01, fig.width = 9, fig.cap = .tag(1, "The Rips barcode diagram from TDA.")------------
diag <- AML10.node287.rips[["diagram"]]
opar <- par(mfrow = c(1,2))
plot(diag, barcode = TRUE, main = "Barcode")
plot(diag, main = "Rips Diagram")
par(opar)
rm(opar)

## ----merc-------------------------------------------------------------------------------------
mercury <- Mercator(amldist, metric = "euclidean", method = "hclust", K = 8)
mercury <- addVisualization(mercury, "mds")
mercury <- addVisualization(mercury, "tsne")
mercury <- addVisualization(mercury, "umap")
mercury <- addVisualization(mercury, "som")

## ----fig03, fig.width = 9, fig.height = 12, fig.cap = .tag(3, "Mercator Visualizations of the distance matrix.")----
opar <- par(mfrow = c(3,2), cex = 1.1)
plot(mercury, view = "hclust")
plot(mercury, view = "mds", main = "Mult-Dimensional Scaling")
plot(mercury, view = "tsne", main = "t-SNE")
plot(mercury, view = "umap", main = "UMAP")
barplot(mercury, main = "Silhouette Width")
plot(mercury, view = "som", main = "Self-Organizing Maps")
par(opar)
rm(opar)

## ----fig04, fig.cap = .tag(4, "Hierarchical connections between zero cycles.")----------------
nzero <- sum(diag[, "dimension"] == 0)
cycles <- AML10.node287.rips[["cycleLocation"]][1:nzero]
L <- sapply(cycles, length)
cycles <- cycles[L > 0]
W <- mercury@view[["umap"]]$layout
plot(W, main = "Connected Zero Cycles")
for (cyc in cycles) {
  points(W[cyc[1], , drop = FALSE], pch = 16,col = "red")
  X <- c(W[cyc[1], 1], W[cyc[2],1])
  Y <- c(W[cyc[1], 2], W[cyc[2],2])
  lines(X, Y)
}

## ----igraph-----------------------------------------------------------------------------------
edges <- t(do.call(cbind, cycles)) # this creates an "edgelist"
G <- graph_from_edgelist(edges)
G <- set_vertex_attr(G, "label", value = attr(amldist, "Labels"))

## ----layouts----------------------------------------------------------------------------------
set.seed(2734)
Lt <- layout_as_tree(G)
L <- layout_with_fr(G)

## ----fig05, fig.cap = .tag(5, "Two igraph depictions of the zero cycle structure."), fig.width = 10----
opar <- par(mfrow = c(1,2), mai = c(0.01, 0.01, 1.02, 0.01))
plot(G, layout = Lt, main = "As Tree")
plot(G, layout = L, main = "Fruchterman-Reingold")
par(opar)

## ----keg--------------------------------------------------------------------------------------
keg <- cluster_edge_betweenness(G)
table(membership(keg)) 
pal <- Dark24[membership(keg)]

## ----fig06, fig.width = 6, fig.height = 6, fig.cap = .tag(6, "Community structure, simplified.")----
is.hierarchical(keg)
H <- as.hclust(keg)
H$labels <- attr(amldist, "Labels")
K <-  10
colset <- Light24[cutree(H, k=K)]
G2 <- set_vertex_attr(G, "color", value = colset)
e <- 0.01
opar <- par(mai = c(e, e, e, e))
plot(G2, layout = L)
par(opar)

## ----fig08, fig.width=7, fig.height = 7, fig.cap = .tag(8, "Cladogram realization, from the ape package.")----
P <- as.phylo(H)
opar <- par(mai = c(0.01, 0.01, 1.0, 0.01))
plot(P, type = "u", tip.color = colset, cex = 1.2, main = "Ape/Cladogram")
par(opar)
rm(opar)

## ----views------------------------------------------------------------------------------------
U <- mercury@view[["mds"]]
V <- mercury@view[["tsne"]]$Y
W <- mercury@view[["umap"]]$layout

## ----fig10, fig.width = 9, fig.cap = .tag(10, "UMAP visualizations with AML10.node287 features.")----
featKi67 <- Feature(AML10.node287[,"Ki-67"], "Ki-67", c("cyan", "red"), c("Low", "High"))
featCD99 <- Feature(AML10.node287[,"CD99"], "CD99", c("green", "magenta"), c("Low", "High"))
opar <- par(mfrow = c(1,2))
plot(W, main = "UMAP; Ki-67", xlab = "U1", ylab = "U2")
points(featKi67, W, pch = 16, cex = 1.4)
plot(W, main = "UMAP; CD99", xlab = "U1", ylab = "U2")
points(featCD99, W, pch = 16, cex = 1.4)
par(opar)
rm(opar)

## ----persistence------------------------------------------------------------------------------
persistence <- diag[, "Death"] - diag[, "Birth"]

## ----d0---------------------------------------------------------------------------------------
d0 <- persistence[diag[, "dimension"] == 0]
d0 <- d0[d0 < 5]
summary(d0)
mu <- mean(d0)
nu <- median(d0)
sigma <- sd(d0)
shape <- mu^2/sigma^2
rate <- mu/sigma^2
xx <- seq(0, 4, length = 100)
yy <- dgamma(xx, shape = shape, rate = rate)
hist(d0, breaks = 123, freq = FALSE)
lines(xx, yy, col = "purple", lwd = 2)

## ----d1---------------------------------------------------------------------------------------
d1 <- persistence[diag[, "dimension"] == 1]
ef <- ExpoFit(d1) # should be close to log(2)/median? 
plot(ef)
eb <- EBexpo(d1, 200)
hist(eb)
plot(eb, prior = 0.56)
sum(d1 > cutoff(0.8, 0.56, eb)) # posterior 80%, prior 0.56
cutoff(0.8, 0.56, eb)
sum(d1 > 0.237) # post 80%
which(d1 > 0.237)
which.max(d1)

## ---------------------------------------------------------------------------------------------
cyc1 <- Cycle(AML10.node287.rips, 1, 103, "forestgreen")
cyc1@index

## ----fig.width = 12---------------------------------------------------------------------------
cyc2 <- Cycle(AML10.node287.rips, 1, 96, "red")
cyc3 <- Cycle(AML10.node287.rips, 1, 87, "purple")

opar <- par(mfrow = c(1, 3))
plot(cyc1, W, lwd = 2, pch = 16, col = "gray", xlab = "U1", ylab = "U2", main = "UMAP")
lines(cyc2, W, lwd=2)
lines(cyc3, W, lwd=2)

plot(U, pch = 16, col = "gray", main = "MDS")
lines(cyc1, U, lwd = 2)
lines(cyc2, U, lwd = 2)
lines(cyc3, U, lwd = 2)

plot(V, pch = 16, col = "gray", main = "t-SNE")
lines(cyc1, V, lwd = 2)
lines(cyc2, V, lwd = 2)
lines(cyc3, V, lwd = 2)
par(opar)
rm(opar)

## ----fig.width = 8, fig.height = 8------------------------------------------------------------
poof <- angleMeans(W, AML10.node287.rips, cyc3, AML10.node287)
poof[is.na(poof)] <- 0
angle.df <- poof[, c("Ki-67", "CD99", "pRb", "PCNA",
                     "CycA", "CycB")]
colorScheme <- list(c(M = "green", U = "magenta"),
                    c(Hi = "cyan", Lo ="red"),
                    c(Hi = "blue", Lo = "yellow"),
                    c(Hi = "#dddddd", Lo = "#111111"),
                    c(No = "#dddddd", Yes = "brown"),
                    c(No = "#dddddd", Yes = "purple"))
annote <- LoopCircos(cyc1, angle.df, colorScheme)
image(annote)

## ----echo = FALSE, eval = FALSE---------------------------------------------------------------
#  M <-  matrix(U <- unlist(colorScheme), ncol = 2, byrow = TRUE)
#  N <- matrix(1:12, nrow = 2)
#  opar <- par(mai = c(0, 2, 0, 2))
#  image(1:2, 1:6, N, col = U, xaxt = "n", yaxt = "n", xlab = "", ylab = "")
#  mtext(sapply(colorScheme, names)[1,], side = 2 , at = 1:6, las = 2, line = 1, adj = 1)
#  mtext(paste(sapply(colorScheme, names)[2,],
#              colnames(angle.df), sep = ", "),
#        side = 4 , at = 1:6, las = 2, line = 1, adj = 0)
#  par(opar)

## ----d2---------------------------------------------------------------------------------------
d2 <- persistence[diag[, "dimension"] == 2]
ef <- ExpoFit(d2) # should be close to log(2)/median? 
plot(ef)
eb <- EBexpo(d2, 200)
hist(eb)
plot(eb, prior = 0.75)
sum(d2 > cutoff(0.8, 0.75, eb)) # posterior 80%, prior 0.56
sum(d2 > cutoff(0.95, 0.75, eb)) # posterior 90%, prior 0.56
cutoff(0.95, 0.75, eb)
sum(d2 > 0.032) # post 90%
which(d2 > 0.032)

## ---------------------------------------------------------------------------------------------
vd <- getCycle(AML10.node287.rips, 2)
mds <- cmdscale(amldist, k = 3)
support <- cycleSupport(vd, mds)

## ----cleanup------------------------------------------------------------------
options(oopt)
#rm(list = ls())

