## 04-angles.R
## Copyright (C) 2022-4 Kevin R. Coombes, RB McGee, and Jake Reed
## LICENSE: Perl Artistic License 2.0

## helper function to get cycle-centroid
getCentroid <- function(cycle, view) {
  pts <- cycleSupport(cycle, view)
  apply(pts, 2, mean)
}

## helper function to compute angles of data points around centroid
getAngles <- function(dset, centroid) {
  recentered <- sweep(dset, 2, centroid, "-") # centroid now at (0,0)
  atan2(recentered[,1], recentered[,2])
}

## Function to calculate mean and SD of section of graph starting from
## centroid
angleMeans <- function(view, rips, cycle = NULL, dset, angleWidth = 20, incr = 15) {
  if (is.null(cycle)) {
    cycle <- new("Cycle")
    cycle@index <- getCycle(rips, 1) # longest cycle
  }
  ctr <- getCentroid(cycle@index, view)
  theta <- getAngles(view, ctr)
  degrees <- theta*180/pi
  deg <- c(degrees, degrees + 360)
  od <- order(deg)
  deg <- deg[od]
  magic <- rbind(dset, dset)[od,]
  partition <- seq(0, 360 - incr, incr)
  GM <- t(sapply(partition, function(center) {
    lb <- center - angleWidth/2
    ub <- center + angleWidth/2
    set <- subset(magic, deg > lb & deg < ub)
    m.gene <- apply(set, 2, mean, na.rm = TRUE)
    m.gene
  }))
  rownames(GM) <- partition
  GM
}

LoopCircos <- function(cycle, angles, colors) {
  new("LoopCircos",
      cycle,
      angles = angles,
      colors = colors)
}

setMethod("image", "LoopCircos", function(x, na.col = "grey", ...) {
  opar <- par(cex = 1.5, mai = c(0, 0, 0, 0))
  on.exit(par(opar))
  circos.clear()
  ## Should probably compute the parameters
  th <- 0.5 / dim(x@angles)[2]
  circos.par(track.height = th, start.degree = 90,
    circle.margin = 0.1)
  ## For each clinical feature/gene/what6ever
  for(i in 1:length(x@angles[1,])) {
#    cat(i, "\n", file = stderr())
    data <- as.data.frame(x@angles[,i])
    col <- colorRamp2(range(data, na.rm = TRUE), x@colors[[i]])
    args <- list(mat = data, cluster = FALSE, col = col, na.col = na.col)
    ## If statement to add angle designations on outside of 
    ## first track only
    if (i == 1) {
      args$rownames.side = "outside"
      args$rownames.cex = 1
    }
    do.call(circos.heatmap, args)
  }
  invisible(x)
})
