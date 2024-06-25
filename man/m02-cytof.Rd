\name{cytof-data}
\alias{cytof-data}
\alias{AML10.node287}
\alias{AML10.node287.rips}
\alias{Arip}
\docType{data}
\title{CyTOF Data}
\description{
  This data set contains 18 columns of protein expression values
  measured in 214 early monocyte cells using mass cytometry (CyTOF).
}
\usage{
data(cytof)
}
\format{ Note that there are three distinct objects included in the data
  set: \code{AML10.node287}, \code{AML10,node287.rips}, and \code{Arip}.
  \describe{
    \item{\code{AML10.node287}}{A numerical data matrix with 214 rows
      (single cells) amd 18 columns (antibody-protein markers). Values
      were obtained using mass cytometry applied to a sample collected
      from a patient with acute myeloid leukemia (AML).} 
    \item{\code{AML10.node287.rips}}{This object  is a
      "Rips diagram". It was produced by running the
      \code{\link[TDA]{ripsDiag}} function from the \code{\link[TDA]{TDA}} R
      package on the \code{AML10.node287} data matrix.}
    \item{\code{Arip}}{This object  is a "Rips diagram". It was produced
      by running the \code{\link[TDA]{ripsDiag}} function from the
      \code{\link[TDA]{TDA}} R package on the Euclidean distance matrix
      between the single cells.}
  }
}
\source{
  Data were collected by Dr. Greg Behbehani while working in the
  laboratory of Dr. Garry Nolan at Stanford University. This data set is
  a tiny subset of the full data set, which was described previously in
  the paper by Behbehani et al.
}
\author{Kevin R. Coombes \email{krc@silicovore.com},
  RB McGee \email{mcgee@sholycross.edu},
  Jake Reed \email{hreed@augusta.edu}}

\references{
  Behbehani GK, Samusik N, Bjornson ZB, Fantl WJ, Medeiros BC, Nolan
  GP. Mass Cytometric Functional Profiling of Acute Myeloid Leukemia
  Defines Cell-Cycle and Immunophenotypic Properties That Correlate with
  Known Responses to Therapy. Cancer Discov. 2015
  Sep;5(9):988-1003. doi: 10.1158/2159-8290.CD-15-0298. 
}
\keyword{datasets}

