\name{LoopFeature-class}
\docType{class}
\alias{LoopFeature-class}
\alias{LoopFeature}
\alias{plot,LoopFeature,missing-method}
\alias{plot,LoopFeature,character-method}
\alias{image,LoopFeature-method}
\title{
\code{LoopFeature} Objects For Visualizing Features That Define Loops
}
\description{
  The \code{LoopFeature} class is a tool for understanding and
  visualizing loops (topological circles) and the features that can be
  used to define and interpret them. Having found a (statistically
  significant) loop, we investigate a feature by computing its mean
  expression in sectors of a fixed width (usually 20 degrees) at a grid
  of angles around the circle (usually multiples of 15 degrees from 0
  to 360). We model these data using the function
  \deqn{f(\theta) = A + B\sin(\theta) + C\cos(\theta).}{f(x) = A +
    B*sin(x) + C*cos(x).}
  We then compute the "fraction of unexplained variance" by dividing the
  residual sum of squares from this model by the total variance of the
  feature. Smaller values of this statistic are more likely to identify
  features that vary sytematically around the circle with a single peak
  and a single trough.
}
\usage{
LoopFeature(circMeans)
\S4method{plot}{LoopFeature,missing}(x, y, ...)
\S4method{plot}{LoopFeature,character}(x, y, ...)
\S4method{image}{LoopFeature}(x, ...)
}
\arguments{
  \item{circMeans}{A matrix, assumed to be the output from a call to the
    \code{angleMeans} function. Columns are features and rows are angles.}
  \item{x}{A \code{LoopFeature} object.}
  \item{y}{A character vector; the set of features to plot.}
  \item{\dots}{The usual set of additional graphical parameters. }
}
\section{Slots}{
  \describe{
    \item{\code{data}:}{The input \code{circMeans} data matrix.}
    \item{\code{fitted}:}{A matrix that is the same size as \code{data};
      the results of fitting a model for each feature as a linear
      combination of sine and cosine.}
    \item{\code{RSS}:}{A numeric vector; the residual sum of squares for
      each model.}
    \item{\code{V}:}{A numeric vector; the total variance for each feature.}
    \item{\code{Kstat}:}{A numeric vector, the unexplained variance
      statistic, RSS/V.}
  }
}
\section{Methods}{
  \describe{
    \item{plot(x, y, ...):}{
      For the selected features listed in \code{y} (which can be missing
      or "all" to plot all features), plots the fitted model as a curve
      along with the observed data.
    }
    \item{image(x, ...)}{
      Produce a 2D image of all the features, with each one scaled to
      the range [0,1] and with the rows ordered by where around the loop
      the maximum value occurs.
    }
  }
}
\value{
  The \code{LoopFeature} function constructs and returns an object of the
  \code{LoopFeature} class

  The \code{plot} and \code{image} methods return (invisibly) the
  LoopFeature object that was their first argument.
}
\author{Kevin R. Coombes <krc@silicovore.com>}
\examples{
data(CLL)
view <- cmdscale(daisydist)
circular <- angleMeans(view, ripdiag, NULL, clinical)
lf <- LoopFeature(circular)
sort(lf@Kstat)
plot(lf, "Serum.beta.2.microglobulin")
opar <- par(mai = c(0.82, 0.2, 0.82, 1.82))
image(lf, main = "Clinical Factors in CLL")
par(opar)

}
\keyword{ classes }
\keyword{ graphics }
