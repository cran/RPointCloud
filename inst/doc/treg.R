## ----opts, echo=FALSE-------------------------------------------------------------------------
knitr::opts_chunk$set(fig.width=8, fig.height=5)
oopt <- options(width=96)
.format <- knitr::opts_knit$get("rmarkdown.pandoc.to")
.tag <- function(N, cap ) ifelse(.format == "html",
                                 paste("Figure", N, ":",  cap),
                                 cap)

## ----RPointCloud------------------------------------------------------------------------------
library(RPointCloud)

## ----libpack----------------------------------------------------------------------------------
library(TDA)
library(Polychrome)
data(Dark24)
data(Light24)
library(Mercator)
library(igraph)
library(ClassDiscovery)
library(PCDimension)
library("ape")

## ----CLL--------------------------------------------------------------------------------------
data(treg)
ls()

## ----echo=FALSE, eval = FALSE-----------------------------------------------------------------
#  dmat <- distanceMatrix(dset, "pearson")
#  picked <- Mercator::downsample(target = 250,
#                                 distanceMat = as.matrix(dmat),
#                                 cutoff = 1E-6)
#  treg <- dset[, picked]
#  tmat <- distanceMatrix(treg, "pearson")
#  rm(picked)
#  set.seed(84263)
#  rip <- ripsDiag(tmat, 2, 0.7, "arbitrary", "Dionysus", TRUE)
#  save(rip, treg, tmat, file = "treg.rda")

## ----fig01, fig.width = 9, fig.cap = .tag(1, "The Rips barcode diagram from TDA.")------------
diag <- rip[["diagram"]]
opar <- par(mfrow = c(1,2))
plot(diag, barcode = TRUE, main = "Barcode")
plot(diag, main = "Rips Diagram")
par(opar)
rm(opar)

## ----fig02, fig.width = 12, fig.cap = .tag(2, "Landscape, silhouette, and lambda-cluster plots of the Rips diagram, from TDA.")----
L <- TDA::landscape(diag, KK = 1)
S <- TDA::silhouette(diag)
crt <- TDA::clusterTree(as.matrix(tmat), k = 5, dist = "arbitrary")

opar <- par(mfrow = c(1, 3))
plot(L, type = "l", main = "Landscape")
plot(S, type = "l", main = "Silhouette")
plot(crt, type = "lambda",main = "Lambda Cluster Tree")
par(opar)
rm(L, S, opar)

## ----Mercator---------------------------------------------------------------------------------
M <- Mercator(tmat, metric ="pearson", method = "mds", K = 8)
M <- addVisualization(M, "hclust")
M <- addVisualization(M, "tsne")
M <- addVisualization(M, "umap")
M <- addVisualization(M, "som")
M@palette <- Light24
set.seed(72345)
clue <- kmeans(t(treg), centers = 8, iter.max = 100, nstart = 20)
M <- setClusters(M, clue$cluster)

## ----fig03, fig.width = 9, fig.height = 12, fig.cap = .tag(3, "Mercator Visualizations of the distance matrix.")----
opar <- par(mfrow = c(3,2), cex = 1.1)
plot(M, view = "hclust")
plot(M, view = "mds", main = "Mult-Dimensional Scaling")
plot(M, view = "tsne", main = "t-SNE")
plot(M, view = "umap", main = "UMAP")
barplot(M, main = "Silhouette Width")
plot(M, view = "som", main = "Self-Organizing Maps")
par(opar)
rm(opar)

## ----fig04, fig.cap = .tag(4, "Hierarchical connections between zero cycles.")----------------
nzero <- sum(diag[, "dimension"] == 0)
cycles <- rip[["cycleLocation"]][2:nzero]
W <- M@view[["umap"]]$layout
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
G <- set_vertex_attr(G, "label", value = attr(tmat, "Labels"))

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
keg <- cluster_edge_betweenness(G) # 19
table(membership(keg)) 
pal <- Light24[membership(keg)]

## ----fig06, fig.width = 6, fig.height = 6, fig.cap = .tag(6, "Community structure, simplified.")----
is.hierarchical(keg)
H <- as.hclust(keg)
H$labels <- vertex_attr(G, "names")
K <-  12
colset <- Light24[cutree(H, k=K)]
G2 <- set_vertex_attr(G, "color", value = colset)
e <- 0.01
opar <- par(mai = c(e, e, e, e))
plot(G2, layout = L)
par(opar)

## ----fig08, fig.width=7, fig.height = 7, fig.cap = .tag(8, "Cladogram realization, from the ape package.")----
P <- as.phylo(H)
opar <- par(mai = c(0.01, 0.01, 1.0, 0.01))
plot(P, type = "u", tip.color = colset, cex = 0.8, main = "Ape/Cladogram")
par(opar)
rm(opar)

## ----views------------------------------------------------------------------------------------
U <- M@view[["mds"]]
V <- M@view[["tsne"]]$Y
W <- M@view[["umap"]]$layout

## ----fig10, fig.width = 9, fig.cap = .tag(10, "UMAP visualizations with clinical features.")----
FOXP3 <- Feature(treg["FOXP3",], "FOXP3", c("pink", "skyblue"), c("Low", "High"))
CTLA4 <- Feature(treg["CTLA4", ], "CTLA4", c("green", "magenta"), c("Low", "High"))
opar <- par(mfrow = c(1,2))
plot(W, main = "UMAP; FOXP3", xlab = "U1", ylab = "U2")
points(FOXP3, W, pch = 16, cex = 1.4)
plot(W, main = "UMAP; CTLA4", xlab = "U1", ylab = "U2")
points(CTLA4, W, pch = 16, cex = 1.4)
par(opar)
rm(opar)

## ----cleanup------------------------------------------------------------------
options(oopt)
#rm(list = ls())

