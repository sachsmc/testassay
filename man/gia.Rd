% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/testassay.R
\docType{data}
\name{gia}
\alias{gia}
\title{Growth Inhibition Assay}
\format{
A data frame with variables: \code{parasite}, \code{assay},
   \code{elisa}, \code{gia}, \code{sample}, \code{meanAAgia},
  and \code{varAAgia}
}
\usage{
gia
}
\description{
Data from a growth inhibition assay experiment. Samples were run repeatedly
on different assays, for two different strains of parasites (3d7 and FVO).
\code{elisa} is a measure of the amount of antibody and is measured once for each sample.
\code{sample} is a unique name for each sample and is defined as \code{paste(gia$parasite,gia$elisa,sep=".")}.
\code{gia} is the value
of interest, and the \code{meanAAgia} is the sample level mean, which is the best
estimate of the "true" gia level for that sample. \code{varAAgia} is the sample level
variance.
}
\keyword{datasets}
