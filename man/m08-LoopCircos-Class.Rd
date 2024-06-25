\name{LoopCircos-class}
\docType{class}
\alias{LoopCircos-class}
\alias{LoopCircos}
\alias{image,LoopCircos-method}
\alias{angleMeans}
\title{
The \code{LoopCircos} Class
}
\description{
  A \code{LoopCircos} object represents a one-dimensional cycle (a
  loop, or topogical eequivalent of a circle), along with a set of
  features that vary around the loop and can be plotted in a "circos"
  plot to help explain why the loop exists. The \code{LoopCircos} class
  inherits from the \code{\link{Cycle}} class, so you can think of it as
  a cycle object with extra information (namely, the explanatory
  features).
}
\usage{
LoopCircos(cycle, angles, colors)
angleMeans(view, rips, cycle = NULL, dset, angleWidth = 20, incr = 15)
\S4method{image}{LoopCircos}(x, na.col = "grey", ...)
}
\arguments{
  \item{cycle}{In the \code{LoopCircos} function, an object of the
    \code{\link{Cycle}} class. In the \code{angleMeans} function, an 
    object of the \code{\link{Cycle}} class.} 
  \item{angles}{A matrix where columns are features, rows are angles,
    and the value is the mean expression of that feature in a sector
    around that central angle. These are almost always constructed using
    the \code{angleMeans} function.}
  \item{colors}{A list of character vectors, each of length two, containing
    the names or hexadecimal representations of colors used to create a
    color ramp to display the feature. The list should be the same
    length as the number of features to be displayed.}
  \item{view}{ A (layout) matrix with coordinates showing where to plot
    each point in the feature while computing means in sectors.}
  \item{rips}{A Rips \code{diagram} object.}
  \item{dset}{The numeric matrix from which Featuers are drawn.}
  \item{angleWidth}{A numeric vector of length one; the width (in
    degrees) of each sector.}
  \item{incr}{A numeric vector; the increment (in degeres)  between sector centers.}
  \item{x}{A \code{LoopCircos} object.}
  \item{na.col}{A character string; the color in which to plot missing
    (\code{NA}) data.}
  \item{\dots}{The usual set of additional graphical parameters. }
}
\section{Slots}{
  \describe{
    \item{\code{angles}:}{A matrix of angle means.}
    \item{\code{colors}:}{A list of color-pairs for displayign features.}
  }
}
\section{Methods}{
  \describe{
    \item{image(x, ...):}{
      Produce a circos plot of a \code{LoopCircos} object.
    }
  }
}
\value{
  The \code{angleMeans} function computes the mean expression of each
  numeric feature in a two-dimensional sector swept out around the
  centroid of the support of a cycle. It returns these values as a
  matrix, with each row corresponding to the angle (in degrees) around
  the centroid of the cycle.

  The \code{LoopCircos} function constructs and returns an object of the
  \code{LoopCircos} class

  The \code{image} method returns (invisibly) the LoopCircos
  object that was its first argument.

}
\seealso{
  \code{\link{Cycle}}, \code{\link{LoopFeature}}, \code{\link[TDA]{ripsDiag}}
}
\author{Kevin R. Coombes <krc@silicovore.com>}
\examples{
data(CLL)
cyc1 <- Cycle(ripdiag, 1, 236, "forestgreen")
V <- cmdscale(daisydist)
poof <- angleMeans(V, ripdiag, cyc1, clinical)
poof[is.na(poof)] <- 0
angle.df <- poof[, c("mutation.status", "CatCD38", "CatB2M", "CatRAI",
                     "Massive.Splenomegaly", "Hypogammaglobulinemia")]
colorScheme <- list(c(M = "green", U = "magenta"),
                    c(Hi = "cyan", Lo ="red"),
                    c(Hi = "blue", Lo = "yellow"),
                    c(Hi = "#dddddd", Lo = "#111111"),
                    c(No = "#dddddd", Yes = "brown"),
                    c(No = "#dddddd", Yes = "purple"))
annote <- LoopCircos(cyc1, angle.df, colorScheme)
image(annote)
}
\keyword{ classes }
\keyword{ graphics }
