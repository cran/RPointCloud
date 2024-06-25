\name{Feature-class}
\docType{class}
\alias{Feature-class}
\alias{Feature}
\alias{points,Feature-method}
\alias{plot,Feature,matrix-method}
\title{
\code{Feature} Objects For Visualizing Topological Data Anal;ysis Results
}
\description{
  The \code{Feature} object represents a feature from a data matrix used
  when performing topological data analysis with the \code{TDA}
  package. Features might be genes, proteins, clinical or demographic
  covariates, or any other item measured on a set of patient samples or
  cells.
}
\usage{
Feature(values, name, colors, meaning, ...)
\S4method{plot}{Feature,matrix}(x, y, pch = 16, ...)
\S4method{points}{Feature}(x, view, ...)
}
\arguments{
  \item{values}{A numeric vector.}
  \item{name}{A character vector of length one.}
  \item{colors}{A character vector of length at least two, containing
    the names or hexadecimal representations of colors used to create a
    color ramp to display the feature.}
  \item{meaning}{A character vector of length two containing the
    interpretations of the low and high extreme values of the
    feature. Note that this works perfectly well for binary factors
    represented as numeric 0-1 vectors.}
  \item{x}{A \code{Feature} object.}
  \item{y}{A (layout) matrix with coordinates showing where to plot
    each point in the feature.}
  \item{view}{ A (layout) matrix with coordinates showing where to plot
    each point in the feature.}
  \item{pch}{A number; the graphical plotting character parameter \code{pch}.}
  \item{\dots}{The usual set of additional graphical parameters. }
}
\section{Slots}{
  \describe{
    \item{\code{name}:}{A character vector of length one; the name of
      the feature..}
    \item{\code{values}:}{A numeric vector of the values of this feature. }
    \item{\code{meaning}:}{A character vector of length two containing the
    interpretations of the low and high extreme values of the
    feature. Note that this works perfectly well for binary factors
    represented as numeric 0-1 vectors. }
    \item{\code{colRamp}:}{A function created using the
      \code{\link[circlize]{colorRamp2}} function from the \code{\link[circlize]{circlize}}
      package.}
  }
}
\section{Methods}{
  \describe{
    \item{plot(x, y, lwd = 2, ...):}{
      Produce a plot of a \code{Feature} object.
    }
    \item{points(x, view, ...)}{
      Add a depiction of a feature to an existing plot.
    }
  }
}
\value{
  The \code{Feature} function constructs and returns an object of the
  \code{Feature} class

  The \code{plot} and \code{points} methods return (invisibly) the Feature
  object that was their first argument.
}
\author{Kevin R. Coombes <krc@silicovore.com>}
\examples{
data(CLL)
featSex <- Feature(clinical[,"Sex"], "Sex",
            c("pink", "skyblue"), c("Female", "Male"))
V <- cmdscale(daisydist)
plot(featSex, V)
plot(V)
points(featSex, V)
}
\keyword{ classes }
\keyword{ graphics }
