\name{CLL-data}
\alias{CLL-data}
\alias{clinical}
\alias{daisydist}
\alias{ripdiag}
\docType{data}
\title{Chronic Lymphocytic Leukemia Clinical Data}
\description{
  Contains 29 columns of deidentified clinical data on 266 patients with
  chronic lymphocytic leukemia (CLL).
}
\usage{
data(CLL)
}
\format{ Note that there are three distinct objects included in the data
  set: \code{clinical}, \code{daisydist}, and \code{ripDiag}.
  \describe{
    \item{\code{clinical}}{A numerical data matrix with 266 rows
      (patients) amd 29 columns (clinical features). Patients were newly
      diagnosed with CLL and previously untreated at the time the
      clinical measurements were recorded.}
    \item{\code{daisydist}}{A distance matrix, stored as a \code{dist}
      object, recording pairwise distances between the 266 CLL
      patients. Distances were computed using the \code{\link[cluster]{daisy}}
      function of Kaufmann and Rooseeuw, as implemented in the
      \code{cluster} R package.}
    \item{\code{ripDiag}}{This object  is a "Rips diagram". It was produced
      by running the the \code{\link[TDA]{ripsDiag}} function from the
      \code{\link[TDA]{TDA}} R package on the \code{daisydist} distance matrix.}
  }
}
\source{
  Data were collected in the clinic of Dr. Michael Keating at the
  University of Texas M.D. Anderson Cancer Center to accompany patient
  samples sent to the laboratory of Dr. Lynne Abruzzo. Various subsets
  of the data have been reported previously; see the references for some
  publications. Computation of the daisy distance was performed by
  Caitlin E. Coombes, whose Master's Thesis also showed that this was 
  the most appropriate method when dealing with heterogeneous clinical
  data of mixed type. Computation of the topological data analysis (TDA) using the
  Rips diagram was performed by Jake Reed.
}
\author{Kevin R. Coombes \email{krc@silicovore.com},
  Caitlin E. Coombes \email{ccoombes@stanford.edu},
  Jake Reed \email{hreed@augusta.edu}}

\references{
  Schweighofer CD, Coombes KR, Barron LL, Diao L, Newman RJ, Ferrajoli
  A, O'Brien S, Wierda WG, Luthra R, Medeiros LJ, Keating MJ, Abruzzo
  LV.  A two-gene signature, SKI and SLAMF1, predicts time-to-treatment
  in previously untreated patients with chronic lymphocytic
  leukemia. PLoS One. 2011;6(12):e28277. doi:
  10.1371/journal.pone.0028277.

  Duzkale H, Schweighofer CD, Coombes KR, Barron LL, Ferrajoli A,
  O'Brien S, Wierda WG, Pfeifer J, Majewski T, Czerniak BA, Jorgensen
  JL, Medeiros LJ, Freireich EJ, Keating MJ, Abruzzo LV.  LDOC1 mRNA is
  differentially expressed in chronic lymphocytic leukemia and predicts
  overall survival in untreated patients. Blood. 2011 Apr
  14;117(15):4076-84. doi: 10.1182/blood-2010-09-304881.

  Schweighofer CD, Coombes KR, Majewski T, Barron LL, Lerner S, Sargent
  RL, O'Brien S, Ferrajoli A, Wierda WG, Czerniak BA, Medeiros LJ,
  Keating MJ, Abruzzo LV. Genomic variation by whole-genome SNP mapping
  arrays predicts time-to-event outcome in patients with chronic
  lymphocytic leukemia: a comparison of CLL and HapMap genotypes. J Mol
  Diagn. 2013 Mar;15(2):196-209. doi: 10.1016/j.jmoldx.2012.09.006.

  Herling CD, Coombes KR, Benner A, Bloehdorn J, Barron LL, Abrams ZB,
  Majewski T, Bondaruk JE, Bahlo J, Fischer K, Hallek M, Stilgenbauer S,
  Czerniak BA, Oakes CC, Ferrajoli A, Keating MJ, Abruzzo
  LV. Time-to-progression after front-line fludarabine,
  cyclophosphamide, and rituximab chemoimmunotherapy for chronic
  lymphocytic leukaemia: a retrospective, multicohort study. Lancet
  Oncol. 2019 Nov;20(11):1576-1586. doi: 10.1016/S1470-2045(19)30503-0.

  Coombes CE, Abrams ZB, Li S, Abruzzo LV, Coombes KR. Unsupervised
  machine learning and prognostic factors of survival in chronic
  lymphocytic leukemia. J Am Med Inform Assoc. 2020 Jul
  1;27(7):1019-1027. doi: 10.1093/jamia/ocaa060. 
}
\keyword{datasets}

