\name{Projection-class}
\docType{class}
\alias{Projection-class}
\alias{Projection}
\alias{plot,Projection,missing-method}
\alias{image,Projection-method}
\title{
\code{Projection} Objects For Visualizing 3D Topological Data Analysis
Results in 2D
}
\description{
  The \code{Projection} class is a tool for understanding and
  visualizing voids (topological spheres) and the features that can be
  used to define and interpret them. The core idea is to first project all
  data points in dimension (at least) 3 radially onto the surface of the
  sphere. We then use spherical coordinates on the sphere (the latitude,
  \code{psi} from 0 to 180 degrees and the longitude, \code{phi}, from -180 to 180
  degrees) to locate the points. The \code{plot} method shows these
  points in two dimensions, colored based on the expression of the
  chosen \code{Feature}. The \code{image} method uses a loess fit to
  smooth the observed data points across the entire surface of the
  sphere and displays the resulting image. Using latitude and longitude
  in this way is similar to a Mercator projection of the surface of the
  Earth onto a two-dimensional plot.

  It is worth noting that, for purposes of computing the smoothed
  version, we actually take two "trips" around the "equator" to fit the
  model, but only display the middle section of that fit. The point of
  this approach is to avoid edge effects in smoothing the data.
}
\usage{
Projection(cycle, view, feature, span = 0.3, resn = 25)
\S4method{plot}{Projection,missing}(x, y, pch = 16, ...)
\S4method{image}{Projection}(x, ...)
}
\arguments{
  \item{cycle}{A \code{Cycle} object representing a void, or 2-cycle.}
  \item{view}{ A (layout) matrix with coordinates showing where to plot
    each point in the feature. Must be at least three-dimensional.}
  \item{feature}{A \code{Feature} object, representing the expression
    levels of one of the features in the underlying data set.}
  \item{span}{Parameter used for a \code{loess} fit.}
  \item{resn}{The number of steps used to create the vectors \code{psi}
  and half-\code{phi}.}
  \item{x}{A \code{Projection} object.}
  \item{y}{ignored.}
  \item{pch}{A number; the graphical plotting character parameter \code{pch}.}
  \item{\dots}{The usual set of additional graphical parameters. }
}
\section{Slots}{
  \describe{
    \item{\code{phi}:}{A numeric vector}
    \item{\code{psi}:}{A numeric vector}
    \item{\code{displayphi}:}{A numeric vector}
    \item{\code{displaypsi}:}{A numeric vector}
    \item{\code{value}:}{A matrix.}
    \item{\code{feature}:}{A \code{Feature} object.}
  }
}
\section{Methods}{
  \describe{
    \item{plot(x, y, lwd = 2, ...):}{
      Produce a 2-dimensional scatter plot of a \code{Projection} object.
    }
    \item{image(x, ...)}{
      Produce a 2D heatmap image of a \code{Projection} object.
    }
  }
}
\value{
  The \code{Projection} function constructs and returns an object of the
  \code{Projection} class

  The \code{plot} and \code{image} methods return (invisibly) the
  Projection object that was their first argument.
}
\author{Kevin R. Coombes <krc@silicovore.com>}
\examples{
data(CLL)
cyc <- Cycle(ripdiag, 2, NULL, "black")
V <- cmdscale(daisydist, 3)
feat <- Feature(clinical[,"stat13"], "Deletion 13q",
               c("green", "magenta"), c("Abnormal", "Normal"))
ob <- Projection(cyc, V, feat, span = 0.2)
plot(ob)
image(ob, col = colorRampPalette(c("green", "gray", "magenta"))(64))

}
\keyword{ classes }
\keyword{ graphics }
