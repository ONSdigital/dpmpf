
## HAS_TESTS
#' Create an object of class Rcpp_CdmsNoreg
#'
#' Objects of class \code{Rcpp_CdmsNoreg} hold
#' cohort no-region data models for a single demographic
#' series, eg population or deaths.
#'
#' @param children A list of objects of class "CdmNoreg"
#'
#' @return An object of class \code{Rcpp_CdmsNoreg}
#'
#' @seealso \code{\link{new_CdmNoregPoibin}},
#' \code{\link{new_CdmNoregNbinom}},
#' \code{\link{new_CdmsWithreg}}
#'
#'
#' @export
new_CdmsNoreg <- function(children = list()) {
  methods::new(CdmsNoreg,
    children = children
  )
}

## HAS_TESTS
#' Create an object of class Rcpp_CdmsNoreg
#'
#' Objects of class \code{Rcpp_CdmsNoreg} hold
#' cohort with-region data models for a single demographic
#' series, eg population or deaths.
#'
#' @param children A list of objects of class "CdmWithreg"
#'
#' @return An object of class \code{Rcpp_CdmsWithreg}
#'
#' @seealso \code{\link{new_CdmWithregPoibin}},
#' \code{\link{new_CdmWithregNbinom}},
#' \code{\link{new_CdmsNoreg}}
#'
#' @export
new_CdmsWithreg <- function(children = list()) {
  methods::new(CdmsWithreg,
    children = children
  )
}
