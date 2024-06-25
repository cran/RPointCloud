\name{disentangle}
\alias{disentangle}
\title{
Disentangling Rips Diagrams From Their Initial Data Coordinates
}
\description{
The \code{\link[TDA]{ripsDiag}} function in the \code{TDA} package produces
very different results depending on whether you invoke it on a data
matrix (expressed in terms of specific data coordinates) or a distance
matrix (expressed as abstract indices). This function converts the
specific coordinates into indices, allowing one to more easily plot
different views of the data structures.
}
\usage{
disentangle(rips, dataset)
}
\arguments{
  \item{rips}{A Rips \code{diagram} produced from a data matrix using
       the \code{TDA} package.}
  \item{dataset}{The original dataset used to create the \code{rips} object.}
}
\value{
  Returns a Rips diagram nearly identical to the one that would be
  produced if the \code{ripsDiag} function had been invoked instead on
  the Euclidean distance matrix.
}
\details{
  The core algorithm is, quite simply, to recombine the coordinates of a
  point in the Rips diagram in a manner consistent with their storage in
  the original data set, and find the index (i.e., row number) of that
  point.
}
\author{Kevin R. Coombes <krc@silicovore.com>}
\examples{
data(cytof)
fixed <- disentangle(Arip, AML10.node287)
}
\keyword{ manip }
