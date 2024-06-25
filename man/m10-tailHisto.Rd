\name{tailHisto}
\alias{tailHisto}
\title{
Function to extract the right tail of a histogram for display purposes.
}
\description{
  Sometimes, the most interesting part of a histogram lies in the tail,
  where details are obscured because of the scale of earlier peaks in
  the distribution. This function uses a user-defined target cutoff to
  extract the right tail of a histogram, preserving its structure as a
  histogram object.
}
\usage{
tailHisto(H, target)
}
\arguments{
  \item{H}{A \code{histogram} object.}
  \item{target}{A real number; the target cutoff defining the portion of
    the tail of the histogram to be extracted.}
}
\value{
  Returns another \code{histogram} object that only contains the portion
  of the histogram beyond the target cutoff.
}
\details{
  There is nothing special going on. The only sanity check is to ensure
  that the \code{target} is small enough that there is actually a part
  of the histogram that can be extracted. After that, we simply cut out
  the approprioate pieces, make sure they are structured properly, and
  return them.
}
\author{Kevin R. Coombes <krc@silicovore.com>}
\examples{
set.seed(12345)
fakeData <- rexp(2000, rate = 10)
H <- hist(fakeData, breaks = 123, plot = FALSE)
H2 <- tailHisto(H, 0.3)
opar <- par(mai=c(0.9, 0.9, 0.6, 0.2))
plot(H, freq = FALSE)
par(mai=c(3.4, 3.0, 0.6, 0.6), new=TRUE)
plot(H2, freq = FALSE, main = "", col = "skyblue")
par(opar)
}
\keyword{ manip }
