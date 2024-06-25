## 06-voids.R
## Copyright (C) 2022-4 Kevin R. Coombes, RB McGee, and Jake Reed
## LICENSE: Perl Artistic License 2.0

voidPlot <- function(cycle, view, feature = NULL, radius = 0.01, ...) {
  index <- cycle@index
  creak <- unique(rbind(index[, 1:2], index[, 2:3], index[, c(1,3)]))
  left <- t(sapply(creak[,1], getCoords, view = view))
  right <- t(sapply(creak[,2], getCoords, view = view))
  if (is.null(feature)) {
    colors <- "gray"
  } else {
    colors <- feature@colRamp(feature@values)
  }
  open3d(windowRect=c(100, 100, 900, 900))
  bg3d(color="white")
  axes3d(col='black')
  view3d(theta=-40, phi=35)
  spheres3d(view[,1], view[,2], view[,3], radius=radius,
            alpha=1, shininess=50, col = colors)
  for (J in 1:nrow(left)) {
    spheres3d(left[J,], color = cycle@color, radius=radius)
    spheres3d(right[J,], color = cycle@color, radius=radius)
    lines3d(rbind(left[J,], right[J,]), color = cycle@color, radius=0.5, lwd=3)
  }
  if (!is.null(feature)) {
    title3d(feature@name, line = 6, cex = 2)
  }
  invisible(cycle)
}

voidFeature <- function(feature, view, radius = 0.01, ...) {
  colors <- feature@colRamp(feature@values)
  spheres3d(view[,1], view[,2], view[,3], radius = radius,
            alpha=1, shininess=90, col = colors, ...)
  title3d(feature@name, line = 6, cex = 2)
  invisible(feature)
}

setFeature <- function(object, feature, span = 0.3) {
  E <- feature@values
  phi <- object@phi
  psi <- object@psi
  wicked <- loess(E ~ phi + psi, span = span)
  daft <- data.frame(phi = rep(object@displayphi,
                               times = length(object@displaypsi)),
                     psi = rep(object@displaypsi,
                               each = length(object@displayphi)))
  eek <- matrix(predict(wicked, daft),
                ncol = length(object@displaypsi))
  object@value <- eek
  object@feature <- feature
  object
}

Projection <- function(cycle, view, feature, span = 0.3, resn = 25) {
  if (ncol(view) < 3) {
    stop("Need a 3D view to compute a projection.\n")
  }
  if (cycle@dimension < 2) {
    stop("Need a void (dimension == 2) to compute a projection.\n")
  }
  centroid <- getCentroid(cycle@index, view)
  recent <- sweep(view, 2, centroid, "-")
  ## spherical coordinates
  r <- sqrt(apply(recent^2, 1, sum))
  obj <-  new("Projection",
              phi = atan2(recent[,2], recent[,1]),
              psi = acos(recent[,3]/r),
              displayphi = seq(-pi, pi, length = 1 + 2*resn),
              displaypsi = seq(0, pi, length = 1 + resn),
              feature = new("Feature"),
              value = matrix()
              )
  if (!missing(feature)) {
    obj <- setFeature(obj, feature, span)
  }
  obj
}

setMethod("plot", c("Projection", "missing"), function(x, y, pch = 16, ...) {
  if (length(x@feature@values) > 0) {
    ptcol <- x@feature@colRamp(x@feature@values)
    main <- x@feature@name
  } else {
    ptcol <- "black"
    main <- ""
  }
  plot(x@phi, x@psi, col = ptcol, pch = pch, main = main,
       xlab = "Longitude (phi)", ylab = "Latitude (psi)", ...)
  invisible(x)
})

setMethod("image", "Projection", function(x, ...) {
  image(x@displayphi, x@displaypsi, x@value,
        main = x@feature@name,
        xlab = "Longitude (phi)", ylab = "Latitude (psi)", ...)
  invisible(x)
})

