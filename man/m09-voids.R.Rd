\name{voids}
\alias{voids}
\alias{voidPlot}
\alias{voidFeature}
\title{
3D Scatter Plots of Cycles and Features
}
\description{
  Voids are two-dimensional cycles and should be thought of as the
  topological equivalents of (the surface) of a sphere. These are
  typically displayed as a point cloud in space, with wire-frame edges
  outlining the component (triangular) simplices.
}
\usage{
voidPlot(cycle, view, feature = NULL, radius = 0.01, ...)
voidFeature(feature, view, radius = 0.01, ...)
}
\arguments{
  \item{cycle}{An object of the \code{\link{Cycle}} class.}
  \item{view}{A matrix (with three columns, x, y, and, z) specifying the
    coordinates to be used to display the data.}
  \item{feature}{An object of the \code{\link{Feature}} class.}
  \item{radius}{The radius of spheres used to display the points in
    three dimensions.}
  \item{\dots}{Additional graphical parameters.}
}
\value{
  Both \code{voidPlot} and \code{voidFeature} invisibly return their
  first argument.
}
\details{
  The 3D displays are implemented using the \code{\link[rgl]{rgl}}
  package. The \code{voidPlot} function adds a wire-frame two-cycle
  to the 3d scatter plot of the underlying data. If a feature is
  present, then the point are colored by the expression lecvel of that
  feature. If not present, they are simply colored gray. However, you can
  always add the expression levels of a feature to an existing plot
  using the \code{voidFeature} function.

  In an interactive session, everything is displayed in an OpenGL window,
  where the graph can be manipulated with a mouse. Programmatically, you
  can use the \code{rgl::rglwidget} function to write out the code for
  an interactive HTML display, and you can save the widget to a file
  using the \code{htmlwidgets::saveWidget} function.
}
\author{Kevin R. Coombes <krc@silicovore.com>}
\examples{
data(CLL)
featRai <- Feature(clinical[,"CatRAI"], "Rai Stage", c("green", "magenta"), c("High", "Low"))
vd <- Cycle(ripdiag, 2, 95, "blue")
mds <- cmdscale(daisydist, k = 3)
\donttest{voidPlot(vd, mds)}
\donttest{voidFeature(featRai, mds, radius = 0.011)}
}
\keyword{ graphics }
