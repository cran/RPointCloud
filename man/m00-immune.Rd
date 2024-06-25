\name{treg-data}
\alias{treg-data}
\alias{treg}
\alias{tmat}
\alias{rip}
\docType{data}
\title{Single Cell Data on T Regulatory (Treg) Cells}
\description{
  This data set contains mRNA and protein-antibody data on T-regulatory
  immune cells. It is a subset of a much larger data set collected from
  the peripheral blood of patients with a variety of health comnditions.
}
\usage{
data(treg)
}
\format{ Note that there are three distinct objects
  included in the data set: \code{treg}, \code{tmat}, and \code{rip}.
  \describe{
    \item{\code{treg}}{A numerical data matrix with 538 rows and 255
      columns. Each column represents a single cell from one of 61
      samples that were assayed by (mixed-omics) single cell
      sequencing. Each row represents one of the features that was
      measured in the assay. Of these, 51 are antibodies that were
      tagged with an RNA barcode to identify them; their names all end
      with the string \code{pAbO}. The remaining 487 features are mRNA
      measurements, named by their official gene symbol at the time the
      experiment was performed. Each column represents a different
      single cell. This matrix is a subset of a more complete data set
      of T regulatory cells (Tregs). It was produced using the
      \code{downsample} function from the \code{\link[Mercator]{Mercator}}
      package, which was in turn inspired by a similar routine used in
      the SPADE algorithm by Peng Qiu. A key point of the algorithm is
      to make sampling less likely from the densest part of the
      distribution in order to preserve rare cell types in the
      population.}
    \item{\code{tmat}}{A distance matrix, stored as a
      \code{\link{dist}} object, produced using Pearson correlation as a
      measure of distance between sigle-cell vectors in the \code{treg}
      data set.}
    \item{\code{rip}}{This object  is a "Rips diagram". It was produced
      by running the the \code{\link[TDA]{ripsDiag}} function from the
      \code{\link[TDA]{TDA}} R package on the \code{treg} subset of single cells.}
  }
}
\source{
  Data were kindly provided by Dr. Klaus Ley, director of the Georgia
  Immunolgy Center at Augusta University. Analysis to perform cell
  typing with Seurat and isolate the subset of 769 T regulatory cells
  found in \code{dset} was performed by Jake Reed, as was the
  downsampling to select the random subset of 255 cells in the
  \code{treg} subset and use them for topolgical data analysis to
  compute the \code{rip} object.
}
\author{Kevin R. Coombes \email{krc@silicovore.com},
  Jake Reed \email{hreed@augusta.edu}}

\references{
  Qiu P, Simonds EF, Bendall SC, Gibbs KD Jr, Bruggner RV, Linderman MD,
  Sachs K, Nolan GP, Plevritis SK. Extracting a cellular hierarchy from
  high-dimensional cytometry data with SPADE. Nat Biotechnol. 2011 Oct
  2;29(10):886-91. doi: 10.1038/nbt.1991.
}
\keyword{datasets}

