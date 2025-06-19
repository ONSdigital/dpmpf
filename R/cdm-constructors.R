
## Poisson binomial -----------------------------------------------------------

## HAS_TESTS
#' Create an object of class Rcpp_CdmNoregPoibin
#'
#' Create an object of (C++) class \code{Rcpp_CdmNoregPoibin}.
#' Objects of class \code{Rcpp_CdmNoregPoibin} are cohort-level
#' data models, plus cohort-level data, for a single dataset. The data
#' does not have a region dimension. The data model
#' is based on a Poisson-binomial mixture.
#'
#' \code{new_CdmNoregPoibin} is an internal function, and should
#' never be needed by end users.
#'
#' @param counts_data  A numeric vector of non-negative whole numbers.
#'   \code{NA}s are allowed.
#' @param prob A numeric scalar between 0 and 1 (exclusive).
#'
#' @return An S4 object with class "Rcpp_CdmNoregPoibin"
#'
#' @export
new_CdmNoregPoibin <- function(counts_data, prob) {
  check_counts_data(counts_data)
  check_prob(prob)
  methods::new(CdmNoregPoibin,
    counts_data = counts_data,
    prob = prob
  )
}

## HAS_TESTS
#' Create an object of class Rcpp_CdmWithregPoibin
#'
#' Create an object of (C++) class \code{Rcpp_CdmWithregPoibin}.
#' Objects of class \code{Rcpp_CdmWithregPoibin} are cohort-level
#' data models, plus cohort-level data, for a single dataset. The data
#' has a region dimension. The data model
#' is based on a Poisson-binomial mixture.
#'
#' \code{new_CdmWithregPoibin} is an internal function, and should
#' never be needed by end users.
#'
#' @inheritParams new_CdmNoregPoibin
#' @param counts_data  A numeric matrix of non-negative whole numbers.
#'   \code{NA}s are allowed.
#'
#' @return An S4 object with class "Rcpp_CdmWithregPoibin"
#'
#' @export
new_CdmWithregPoibin <- function(counts_data, prob) {
  check_counts_data(counts_data)
  if (!is.matrix(counts_data)) {
    stop("'counts_data' has class '", class(counts_data), "'")
  }
  check_prob(prob)
  methods::new(CdmWithregPoibin,
    counts_data = counts_data,
    prob = prob
  )
}



## Negative binomial ----------------------------------------------------------

## HAS_TESTS
#' Create an object of class Rcpp_CdmNoregNbinom
#'
#' Create an object of (C++) class \code{Rcpp_CdmNoregNbinom}.
#' Objects of class \code{Rcpp_CdmNoregNbinom} are cohort-level
#' data models, plus cohort-level data, for a single dataset. The data
#' does not have a region dimension. The data model
#' is based on the negative binomial distribution.
#'
#' \code{new_CdmNoregNbinom} is an internal function, and should
#' never be needed by end users.
#'
#' @param counts_data  A numeric vector of non-negative whole numbers.
#'   \code{NA}s are allowed.
#' @param ratio A numeric vector of non-negative numbers, the
#'   same length as \code{counts_data}. \code{NA}s are allowed
#'   at positions where \code{counts_data} is \code{NA}.
#' @param disp A numeric vector of non-negative numbers, the
#'   same length as \code{counts_data}. \code{NA}s are allowed
#'   at positions where \code{counts_data} is \code{NA}.
#'
#' @return An S4 object with class "Rcpp_CdmNoregNbinom"
#'
#' @export
new_CdmNoregNbinom <- function(counts_data, ratio, disp) {
  check_counts_data(counts_data)
  check_ratio_cdm(
    ratio = ratio,
    counts_data = counts_data
  )
  check_disp_cdm(
    disp = disp,
    counts_data = counts_data
  )
  methods::new(CdmNoregNbinom,
    counts_data = counts_data,
    ratio = ratio,
    disp = disp
  )
}

## HAS_TESTS
#' Create an object of class Rcpp_CdmWithregNbinom
#'
#' Create an object of (C++) class \code{Rcpp_CdmWithregNbinom}.
#' Objects of class \code{Rcpp_CdmWithregNbinom} are cohort-level
#' data models, plus cohort-level data, for a single dataset. The data
#' has a region dimension. The data model
#' is based on a negative binomial distribution.
#'
#' \code{new_CdmWithregNbinom} is an internal function, and should
#' never be needed by end users.
#'
#' @param counts_data  A numeric matrix of non-negative whole numbers.
#'   \code{NA}s are allowed.
#' @param ratio A numeric matrix of non-negative numbers, with the
#'   same dimensions as \code{counts_data}. \code{NA}s are allowed
#'   at positions where \code{counts_data} is \code{NA}.
#' @param disp A numeric matrix of non-negative numbers, with the
#'   same dimensions as \code{counts_data}. \code{NA}s are allowed
#'   at positions where \code{counts_data} is \code{NA}.
#'
#' @return An S4 object with class "Rcpp_CdmWithregNbinom"
#'
#' @export
new_CdmWithregNbinom <- function(counts_data, ratio, disp) {
  check_counts_data(counts_data)
  check_ratio_cdm(
    ratio = ratio,
    counts_data = counts_data
  )
  check_disp_cdm(
    disp = disp,
    counts_data = counts_data
  )
  for (name in c("counts_data", "ratio", "disp")) {
    val <- get(name)
    if (!is.matrix(val)) {
      stop("'", name, "' has class '", class(val), "'")
    }
  }
  methods::new(CdmWithregNbinom,
    counts_data = counts_data,
    ratio = ratio,
    disp = disp
  )
}




## Normal -----------------------------------------------------------

## HAS_TESTS
#' Create an object of class Rcpp_CdmNoregNorm
#'
#' Create an object of (C++) class \code{Rcpp_CdmNoregNorm}.
#' Objects of class \code{Rcpp_CdmNoregNorm} are cohort-level
#' data models, plus cohort-level data, for a single dataset. The data
#' does not have a region dimension. The data model
#' is based on the normal distribution.
#'
#' \code{new_CdmNoregNorm} is an internal function, and should
#' never be needed by end users.
#'
#' @param counts_data  A numeric vector of non-negative whole numbers.
#'   \code{NA}s are allowed.
#' @param ratio A numeric vector of non-negative numbers, the
#'   same length as \code{counts_data}. \code{NA}s are allowed
#'   at positions where \code{counts_data} is \code{NA}.
#' @param sd A numeric vector of positive numbers, the
#'   same length as \code{counts_data}. \code{NA}s are allowed
#'   at positions where \code{counts_data} is \code{NA}.
#'
#' @return An S4 object with class "Rcpp_CdmNoregNorm"
#'
#' @export
new_CdmNoregNorm <- function(counts_data, ratio, sd) {
  check_counts_data(counts_data)
  check_ratio_cdm(
    ratio = ratio,
    counts_data = counts_data
  )
  check_sd_cdm(
    sd,
    counts_data = counts_data
  )
  methods::new(CdmNoregNorm,
    counts_data = counts_data,
    ratio = ratio,
    sd = sd
  )
}
