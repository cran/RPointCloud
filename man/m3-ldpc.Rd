\name{LDPC}
\alias{LDPC}
\alias{takens}
\title{
Local Dimension of Point Clouds
}
\description{
  Given a data set viewed as a point cloud in N-dimesnional space,
  compute the local estimate of the dimension of an underlying manifold,
  as defined by Ellis and McDermott, analogous to earlier related work
  by Takens.
}
\usage{
takens(r, dists)
LDPC(CellID, dset, rg, quorum, samplesAreRows = TRUE)
}
\arguments{
  \item{CellID}{An integer indexing one of the cells-samples in the data set.}
  \item{dset}{A data set. The usual orientation is that rows are cells
    and columns are features.}
  \item{rg}{A numerical vector of radial distances at which to compute
    Takens estimates of the local dimension.}
  \item{quorum}{The minimum number of neighboring cells required for
    the computation to be meaningful.}
  \item{samplesAreRows}{A logical value: do rows or columns represent
    samples at which to compute local dimensions.}
  \item{r}{A radial distances at which to compute the Takens estimate.}
  \item{dists}{A sorted vector of distances from one cell to all other cells.}
}
\value{
  The \code{takens} function returns a list with two items: the number
  of neighbors $k$ and the dimension estimate $d$ at each value
  of the radius from the input vector.

  The \code{LDPC} function returns a list containing vectors $R$, $k$,
  and $d$ values for each cell in the data set.
}
\details{
  "[T]he procedure is carried out as follows. A 'bin increment', $A$; a
  number, $m$, of bin increments; and a 'quorum', $q > 0$, are chosen
  and raw dimensions are calculated for $r=A$, $2A$, $\dots$, $m$. Next,
  for each observation, $x_i$, let $r_i$ be the smallest multiple of $A$
  not exceeding $mA$ such that the ball with radius $r_i$ centered at
  $x_i$ contains at least $q$ observations, providing that there are at
  least $q$ observations within $mA$ of $x_i$. Otherwise let $r_i = mA$."

  The \code{takens} function computes the Takens estimate of the local
  dimension of a point cloud at radius $r$ around a data point. For each
  cell-sample, we must compute and sort the distances from that cell to
  all other cells. (For the \code{takens} function, these distances are
  passed in as the second argument to the function.) Preliminary
  histograms of distance distributions may be used to inform a good set
  of radial distances. Note that the local dimension estimates are
  infinite if the radius is so small that there are no neighbors. The
  estmiates decrease as the radius increases os as the number of local
  neighbors increases. The reference paper by Ellis and McDermott says:

  The \code{LDPC} function iterates over all cells-samples int he data sets,
  computes and sorts their distance to all other cells, and invokes teh
  \code{takens} function to compute local estimates os dimension.
}
\author{Kevin R. Coombes <krc@silicovore.com>}
\references{
  Ellis and McDermott
}
\examples{
data(cytof)
localdim <- LDPC(1, AML10.node287, seq(1, 6, length=20), 30, TRUE)
}
\keyword{ graphics }
