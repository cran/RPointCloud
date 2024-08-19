\name{EBexpo-class}
\docType{class}
\alias{EBexpo-class}
\alias{EBexpo}
\alias{plot,EBexpo,missing-method}
\alias{hist,EBexpo-method}
\alias{cutoff}
\title{
The \code{EBexpo} Class
}
\description{
  The \code{EBexpo} object represents the results of an Empirical Bayes
  approach to estimate a distribution as a mixture of a (more or less)
  known exponential distribution along with a completely unknown
  "interesting" distribution. The basic method was described by Efron
  and Tibshirani with an application to differential expression in
  microarray data.
}
\usage{
EBexpo(edata, resn = 200)
cutoff(target, prior, object)
\S4method{plot}{EBexpo,missing}(x, prior = 1, significance = c(0.5, 0.8, 0.9),
                 ylim = c(-0.5, 1), xlab = "Duration",
                 ylab = "Probability(Interesting | Duration)", ...)
\S4method{hist}{EBexpo}(x, ...)
}
\arguments{
  \item{edata}{A numeric vector; the observed data that we think comes
    mainly from an exponential distribution.}
  \item{resn}{A numeric vector; the resolution used to estimate a histogream.}
  \item{x}{An \code{EBexpo} object.}
  \item{prior}{A numeric vector of length 1; the prior probability of an
    observed data point coming from the known exponential distribution.}
  \item{significance}{A numeric vector with values between 0 and 1; the
    target posterior probabiltiites.}
  \item{ylim}{A numeric vector of length two.}
  \item{xlab}{A character vector; the label for the x-axis.}
  \item{ylab}{A character vector; the label for the y-axis.}
  \item{\dots}{The usual set of graphical parameters. }
  \item{target}{The target posterior probability.}
  \item{object}{An \code{EBexpo} object.}
}
\section{Slots}{
  \describe{
    \item{\code{xvals}:}{Inherited from
      \code{\link[ClassComparison]{MultiWilcoxonTest}}} 
    \item{\code{statistics}:}{Inherited from
      \code{\link[ClassComparison]{MultiWilcoxonTest}}, Here, these are
      the same a the \code{edata} slot from an \code{link{ExpoFit}} object.}
    \item{\code{pdf}:}{Inherited from
      \code{\link[ClassComparison]{MultiWilcoxonTest}}} 
    \item{\code{theoretical.pdf}:}{Inherited from
      \code{\link[ClassComparison]{MultiWilcoxonTest}}} 
    \item{\code{unravel}:}{Inherited from
      \code{\link[ClassComparison]{MultiWilcoxonTest}}} 
    \item{\code{groups}:}{Inherited from
      \code{\link[ClassComparison]{MultiWilcoxonTest}}, but not used}
    \item{\code{call}:}{Inherited from
      \code{\link[ClassComparison]{MultiWilcoxonTest}}} 
    \item{\code{h0}:}{See \code{\link{ExpoFit}}}
    \item{\code{lambda}:}{See \code{\link{ExpoFit}}}
    \item{\code{mu}:}{See \code{\link{ExpoFit}}}
  }
}
\section{Methods}{
  \describe{
    \item{plot(x, prior, post = c(0.5, 0.8, 0.9), ...):}{
      Produce a plot of a \code{EBexpo} object.
    }
    \item{hist(x, ...):}{
      Produce a histogram of the observed distibution, with overlays.
    }
  }
}
\value{
  The \code{EBexpo} function constructs and returns an object of the
  \code{EBexpo} class

  The \code{plot} and \code{hist} methods return (invisibly) the EBexpo
  object that was their first argument.
}
\references{
  Efron B, Tibshirani R. Empirical bayes methods and false discovery
  rates for microarrays. Genet Epidemiol. 2002 Jun;23(1):70-86. doi:
  10.1002/gepi.1124.
}
\author{Kevin R. Coombes <krc@silicovore.com>}
\examples{
data(cytof)
diag <- AML10.node287.rips[["diagram"]]
persistence <- diag[, "Death"] - diag[, "Birth"]
d1 <- persistence[diag[, "dimension"] == 1]
eb <- EBexpo(d1, 200)
hist(eb)
plot(eb, prior = 0.56)
cutoff(0.8, 0.56, eb)
}
\keyword{ classes }
\keyword{ graphics }
