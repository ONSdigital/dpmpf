
## CdmsNoreg ------------------------------------------------------------------

test_that("'new_CdmsNoreg' creates expected object", {
  el1 <- new_CdmNoregPoibin(
    counts_data = as.numeric(1:4),
    prob = 0.98
  )
  el2 <- new_CdmNoregNbinom(
    counts_data = as.numeric(1:4),
    ratio = rep(1, 4),
    disp = c(0.5, 1, 1, 2.5)
  )
  x <- list(el1, el2)
  ans <- new_CdmsNoreg(x)
  expect_s4_class(ans, class = "Rcpp_CdmsNoreg")
})

## CdmsWithreg ----------------------------------------------------------------

test_that("'new_CdmsWithreg' creates expected object", {
  el1 <- new_CdmWithregPoibin(
    counts_data = matrix(as.numeric(1:8), nr = 2),
    prob = 0.98
  )
  el2 <- new_CdmWithregNbinom(
    counts_data = matrix(as.numeric(1:8), nr = 2),
    ratio = matrix(1.2, nr = 2, nc = 4),
    disp = matrix(2.2, nr = 2, nc = 4)
  )
  x <- list(el1, el2)
  ans <- new_CdmsWithreg(x)
  expect_s4_class(ans, class = "Rcpp_CdmsWithreg")
})
