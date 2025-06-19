
## Poisson binomial -----------------------------------------------------------

test_that("'new_CdmNoregPoibin creates expected object", {
  counts_data <- as.numeric(1:10)
  prob <- 0.99
  ans <- new_CdmNoregPoibin(
    counts_data = counts_data,
    prob = prob
  )
  expect_s4_class(ans, class = "Rcpp_CdmNoregPoibin")
})

test_that("'new_CdmWithregPoibin created expected object", {
  counts_data <- matrix(as.numeric(1:10),
    nr = 2
  )
  prob <- 0.99
  ans <- new_CdmWithregPoibin(
    counts_data = counts_data,
    prob = prob
  )
  expect_s4_class(ans, "Rcpp_CdmWithregPoibin")
})


## Negative binomial -----------------------------------------------------------

test_that("'new_CdmNoregNbinom creates expected object", {
  counts_data <- as.numeric(1:10)
  ratio <- rep(1.1, 10)
  disp <- rep(2, 10)
  ans <- new_CdmNoregNbinom(
    counts_data = counts_data,
    ratio = ratio,
    disp = disp
  )
  expect_s4_class(ans, class = "Rcpp_CdmNoregNbinom")
})

test_that("'new_CdmWithregPoibin created expected object", {
  counts_data <- matrix(as.numeric(1:10),
    nrow = 2
  )
  ratio <- matrix(1, nrow = 2, ncol = 5)
  disp <- matrix(1.3, nrow = 2, ncol = 5)
  ans <- new_CdmWithregNbinom(
    counts_data = counts_data,
    ratio = ratio,
    disp = disp
  )
  expect_s4_class(ans, "Rcpp_CdmWithregNbinom")
})




## Normal -----------------------------------------------------------

test_that("'new_CdmNoregNorm creates expected object", {
  counts_data <- as.numeric(1:10)
  ratio <- rep(1.1, 10)
  sd <- rep(2, 10)
  ans <- new_CdmNoregNorm(
    counts_data = counts_data,
    ratio = ratio,
    sd = sd
  )
  expect_s4_class(ans, class = "Rcpp_CdmNoregNorm")
})
