\name{ExpoFit-class}
\docType{class}
\alias{ExpoFit-class}
\alias{ExpoFit}
\alias{plot,ExpoFit,missing-method}
\title{
The \code{ExpoFit} Class
}
\description{
  An \code{ExpoFit} object represents a robust fit to an exponential
  distribution, in a form that can conveniently be used as part of an
  Empirical Bayes approach to decompose the distributions of cycle
  presistence or duration for a topological data analysis performed
  using the \code{TDA} package.
}
\usage{
ExpoFit(edata, resn = 200)
\S4method{plot}{ExpoFit,missing}(x, y, ...)
}
\arguments{
  \item{edata}{A numeric vector; the observed data that we think comes
    mainly from an exponential distribution.}
  \item{resn}{A numeric vector of length 1; the resolution (number of
    breaks) used to estimate a histogram.}
  \item{x}{An \code{ExpoFit} object.}
  \item{y}{Ignored.}
  \item{\dots}{The usual set of graphical parameters. }
}
\section{Slots}{
  \describe{
    \item{\code{edata}:}{A numeric vector; the observed data that we
      think comes from an exponential distribution.}
    \item{\code{h0}:}{A \code{histogram} object produced by the
      \code{\link{hist}} function applied to the supplied \code{edata}.}
    \item{\code{X0}:}{A numeric vector containing the midpoints of the
      breaks in the histogram object.}
    \item{\code{pdf}:}{The empirical density function extracted from the
      histogram object.}
    \item{\code{mu}:}{The observed mean of the putative exponential
      distribution.}
    \item{\code{lambda}:}{The robustly estimated parameter of the
      exponential distribution. Originally crudely repesented by the
      reciprocal of the mean.}
  }
}
\section{Methods}{
  \describe{
    \item{plot(x, y, lwd = 2, ...):}{
      Produce a plot of a \code{ExpoFit} object.
    }
  }
}
\value{
  The \code{ExpoFit} function constructs and returns an object of the
  \code{ExpoFit} class.

  The \code{plot} method returns (invisibly) the ExpoFit object that was
  its first argument.
}
\seealso{
  \code{\link{EBexpo}}
}
\author{Kevin R. Coombes <krc@silicovore.com>}
\examples{
data(cytof)
diag <- AML10.node287.rips[["diagram"]]
persistence <- diag[, "Death"] - diag[, "Birth"]
d1 <- persistence[diag[, "dimension"] == 1]
ef <- ExpoFit(d1) # should be close to log(2)/median? 
plot(ef)
}
\keyword{ classes }
\keyword{ graphics }
