\name{Cycle-class}
\docType{class}
\alias{Cycle-class}
\alias{Cycle}
\alias{lines,Cycle-method}
\alias{plot,Cycle,matrix-method}
\alias{getCycle}
\alias{cycleSupport}
\title{
\code{Cycle} Objects For Visualizing Topological Data Analysis Results
}
\description{
  The \code{Cycle} object represents a cycle found by performing
  topological data analysis (using the \code{TDA} package) on a data
  matrix or distance matrix.
}
\usage{
Cycle(rips, dimen, J, color)
getCycle(rips, dimension = 1, target = NULL)
cycleSupport(cycle, view)
\S4method{plot}{Cycle,matrix}(x, y, lwd = 2, ...)
\S4method{lines}{Cycle}(x, view, ...)
}
\arguments{
  \item{rips}{A Rips \code{diagram} from the \code{TDA} package.}
  \item{dimen}{An integer; the dimension of the cycle.}
  \item{J}{An integer; the index locating the cycle in the
    \code{cycleLocation} portion of the Rips diagram. If \code{NULL},
    then looks for the most persistent cycle.}
  \item{color}{A character vector of length one; the color in which to
    display the cycle.}
  \item{x}{A \code{Cycle} object.}
  \item{y}{ A (layout) matrix with coordinates showing where to plot
    each point supporting the cycle.}
  \item{view}{ A (layout) matrix with coordinates showing where to plot
    each point supporting the cycle.}
  \item{lwd}{A number; the graphical line width parameter \code{lwd}.}
  \item{\dots}{The usual set of additional graphical parameters. }
  \item{dimension}{An integer; the dimension of the cycle.}
  \item{target}{An integer indexing the desired cycle. Note that the index
    should be relative to other cuycles of the same dimension. If
    \code{NULL}, gets the longest persisting cycle.}
  \item{cycle}{A raw cycle, meaning a simple element of the
    \code{cycleLocation} part of a Rips diagram.}
}
\section{Slots}{
  \describe{
    \item{\code{index}:}{A matrix, containing the indices into the data
      matrix or distance matrix defining the simplices that realize the
      cycle.}
    \item{\code{dimension}:}{A matrix; the dimension of the cycle. }
    \item{\code{color}:}{A character vector; the color in which to
      display the cycle.}
  }
}
\section{Methods}{
  \describe{
    \item{plot(x, y, lwd = 2, ...):}{
      Produce a plot of a \code{Cycle} object.
    }
    \item{lines(x, view, ...)}{
      Add a depiction of a cycle to an existing plot. 
    }
  }
}
\value{
  The \code{Cycle} function constructs and returns an object of the
  \code{Cycle} class

  The \code{plot} and \code{lines} methods return (invisibly) the Cycle
  object that was their first argument.

  The \code{getCycle} function extracts a single raw cycle from a Rips
  diagram and returns it. The \code{cycleSupport} function combines such
  a cycle with a layout/view matrix to extract a list of the coordinates
  of the points contained in the cycle.
}
\seealso{
  \code{\link[TDA]{TDA}}, \code{\link[TDA]{ripsDiag}}
}
\author{Kevin R. Coombes <krc@silicovore.com>}
\examples{
data(CLL)
cyc1 <- Cycle(ripdiag, 1, 236, "forestgreen")
V <- cmdscale(daisydist)
plot(cyc1, V)
plot(V, pch = 16, col = "gray")
lines(cyc1, V)
}
\keyword{ classes }
\keyword{ graphics }
