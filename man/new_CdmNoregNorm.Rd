% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cdm-constructors.R
\name{new_CdmNoregNorm}
\alias{new_CdmNoregNorm}
\title{Create an object of class Rcpp_CdmNoregNorm}
\usage{
new_CdmNoregNorm(counts_data, ratio, sd)
}
\arguments{
\item{counts_data}{A numeric vector of non-negative whole numbers.
\code{NA}s are allowed.}

\item{ratio}{A numeric vector of non-negative numbers, the
same length as \code{counts_data}. \code{NA}s are allowed
at positions where \code{counts_data} is \code{NA}.}

\item{sd}{A numeric vector of positive numbers, the
same length as \code{counts_data}. \code{NA}s are allowed
at positions where \code{counts_data} is \code{NA}.}
}
\value{
An S4 object with class "Rcpp_CdmNoregNorm"
}
\description{
Create an object of (C++) class \code{Rcpp_CdmNoregNorm}.
Objects of class \code{Rcpp_CdmNoregNorm} are cohort-level
data models, plus cohort-level data, for a single dataset. The data
does not have a region dimension. The data model
is based on the normal distribution.
}
\details{
\code{new_CdmNoregNorm} is an internal function, and should
never be needed by end users.
}
