
## 'check_arg_cdm' ------------------------------------------------------------

test_that("'check_arg_cdm' returns TRUE when 'arg' is valid", {
  expect_true(check_arg_cdm(
    arg = as.numeric(c(1:10, NA)),
    nm_arg = "ratio",
    counts_data = as.numeric(c(10:1, NA)),
    neg_ok = FALSE,
    zero_ok = TRUE
  ))
})

test_that("'check_arg_cdm' throws correct error when 'arg' is non-numeric", {
  expect_error(
    check_arg_cdm(
      arg = as.character(c(1:10, NA)),
      nm_arg = "ratio",
      counts_data = as.numeric(c(10:1, NA)),
      neg_ok = FALSE,
      zero_ok = TRUE
    ),
    "'ratio' has class 'character'"
  )
})

test_that("'check_arg_cdm' throws correct error when 'arg' is integer", {
  expect_error(
    check_arg_cdm(
      arg = c(1:10, NA),
      nm_arg = "ratio",
      counts_data = as.numeric(c(10:1, NA)),
      neg_ok = FALSE,
      zero_ok = TRUE
    ),
    "'ratio' has class 'integer'"
  )
})

test_that("'check_arg_cdm' throws correct error when 'arg' is wrong length", {
  expect_error(
    check_arg_cdm(
      arg = as.numeric(1:10),
      nm_arg = "ratio",
      counts_data = as.numeric(c(10:1, NA)),
      neg_ok = FALSE,
      zero_ok = TRUE
    ),
    "'ratio' has length 10 but 'counts_data' has length 11"
  )
})

test_that("'check_arg_cdm' throws correct error when 'arg' has NAs in wrong places", {
  expect_error(
    check_arg_cdm(
      arg = as.numeric(c(1:10, NA)),
      nm_arg = "disp",
      counts_data = as.numeric(10:0),
      neg_ok = FALSE,
      zero_ok = TRUE
    ),
    "'disp' has NA at position 11 but 'counts_data' does not"
  )
})

test_that("'check_arg_cdm' throws correct error when 'arg' is negative", {
  expect_error(
    check_arg_cdm(
      arg = as.numeric(c(1:10, -1)),
      nm_arg = "disp",
      counts_data = as.numeric(10:0),
      neg_ok = FALSE,
      zero_ok = FALSE
    ),
    "'disp' has negative values"
  )
})

test_that("'check_arg_cdm' throws correct error when 'arg' is zero", {
  expect_error(
    check_arg_cdm(
      arg = as.numeric(c(1:10, 0)),
      nm_arg = "disp",
      counts_data = as.numeric(10:0),
      neg_ok = FALSE,
      zero_ok = FALSE
    ),
    "'disp' has zeros"
  )
})


## 'check_counts_data' --------------------------------------------------------

test_that("'check_counts_data' returns TRUE when 'counts_data' is valid", {
  expect_true(check_counts_data(c(as.numeric(1:10), NA)))
})

test_that("'check_counts_data' throws correct error when 'counts_data' is valid", {
  expect_error(
    check_counts_data("wrong"),
    "is\\.numeric\\(counts_data\\) is not TRUE"
  )
  expect_error(
    check_counts_data(1L),
    "is\\.integer\\(counts_data\\) is TRUE"
  )
  expect_error(
    check_counts_data(numeric()),
    "'counts_data' has length 0"
  )
  expect_error(
    check_counts_data(c(-1, NA)),
    "'counts_data' has negative values"
  )
  expect_error(
    check_counts_data(c(NA, 1.0001)),
    "'counts_data' has non-integer values"
  )
})


## 'check_disp_cdm' ------------------------------------------------------------

test_that("'check_disp_cdm' returns TRUE when 'disp' is valid", {
  expect_true(check_disp_cdm(
    disp = as.numeric(c(1:10, NA)),
    counts_data = as.numeric(c(10:1, NA))
  ))
})


## 'check_prob' ---------------------------------------------------------------

test_that("'check_prob' returns TRUE when 'prob' is valid", {
  expect_true(check_prob(0.8))
})

test_that("'check_prob' throws correct error when 'prob' is invalid", {
  expect_error(
    check_prob("a"),
    "'prob' is non-numeric"
  )
  expect_error(
    check_prob(c(0.8, 0.8)),
    "'prob' does not have length 1"
  )
  expect_error(
    check_prob(NA_real_),
    "'prob' is NA"
  )
  expect_error(
    check_prob(0),
    "'prob' is less than or equal to 0"
  )
  expect_error(
    check_prob(1),
    "'prob' is greater than or equal to 1"
  )
})


## 'check_ratio_cdm' ------------------------------------------------------------

test_that("'check_ratio_cdm' returns TRUE when 'ratio' is valid", {
  expect_true(check_ratio_cdm(
    ratio = as.numeric(c(0:9, NA)),
    counts_data = as.numeric(c(10:1, NA))
  ))
})



## 'check_sd_cdm' ------------------------------------------------------------

test_that("'check_sd_cdm' returns TRUE when 'sd' is valid", {
  expect_true(check_sd_cdm(
    sd = as.numeric(c(1:10, NA)),
    counts_data = as.numeric(c(10:1, NA))
  ))
})
