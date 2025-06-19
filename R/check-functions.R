
## HAS_TESTS
## Check a vector argument to a cohort data model constructor.
## 'arg' - the argument being checked, a numeric vector
## 'nm_arg' - the name of the argument being checked
## 'counts_data' - the cohort data being modelled
## 'neg_ok' - whether  'arg' can have negative values
## 'zero_ok' - whether 'arg' can have zeros
## return value - TRUE or an error
check_arg_cdm <- function(arg, nm_arg, counts_data, neg_ok, zero_ok) {
  ## 'arg' is numeric
  if (!is.numeric(arg)) {
    stop("'", nm_arg, "' has class '", class(arg), "'")
  }
  ## 'arg' is not integer
  if (is.integer(arg)) {
    stop("'", nm_arg, "' has class '", class(arg), "'")
  }
  ## 'arg' has the same length as 'counts_data'
  n_arg <- length(arg)
  n_counts <- length(counts_data)
  if (n_arg != n_counts) {
    stop(
      "'", nm_arg, "' has length ", n_arg,
      " but 'counts_data' has length ", n_counts
    )
  }
  ## 'arg' only has NA if 'counts_data' has NA
  is_na_clash <- is.na(arg) & !is.na(counts_data)
  i_na_clash <- match(TRUE, is_na_clash, nomatch = 0L)
  if (i_na_clash > 0L) {
    stop(
      "'", nm_arg, "' has NA at position ", i_na_clash,
      " but 'counts_data' does not"
    )
  }
  ## if 'neg_ok' is FALSE, no elements of 'arg' are negative
  if (!neg_ok && any(arg < 0, na.rm = TRUE)) {
    stop("'", nm_arg, "' has negative values")
  }
  ## if 'zero_ok' is FALSE, no elements of 'arg' are zero
  if (!zero_ok && any(arg == 0, na.rm = TRUE)) {
    stop("'", nm_arg, "' has zeros")
  }
  ## if all tests pass, return TRUE
  TRUE
}


## HAS_TESTS
## Check that 'counts_data' is a vector of non-negative
## whole numbers with non-zero length
## 'counts_data' A numeric vector (could be a matrix)
## return value - TRUE or an error
check_counts_data <- function(counts_data) {
  if (!is.numeric(counts_data)) {
    stop("is.numeric(counts_data) is not TRUE")
  }
  if (is.integer(counts_data)) {
    stop("is.integer(counts_data) is TRUE")
  }
  if (identical(length(counts_data), 0L)) {
    stop("'counts_data' has length 0")
  }
  counts_non_na <- counts_data[!is.na(counts_data)]
  if (any(counts_non_na < 0)) {
    stop("'counts_data' has negative values")
  }
  if (any(counts_non_na != round(counts_non_na))) {
    stop("'counts_data' has non-integer values")
  }
  TRUE
}


## HAS_TESTS
## Check that a 'disp' argument to a cdm constructor is valid
## 'disp' - the disp argument being checked, a numeric vector
## 'counts_data' - the cohort data being modelled
## return value - TRUE or an error
check_disp_cdm <- function(disp, counts_data) {
  check_arg_cdm(
    arg = disp,
    nm_arg = "disp",
    counts_data = counts_data,
    neg_ok = FALSE,
    zero_ok = TRUE
  )
}


## HAS_TESTS
## Check that 'prob' is a number between 0 and 1 (exclusive)
## 'prob' - A numeric scalar
## return value - TRUE or an error
check_prob <- function(prob) {
  if (!is.numeric(prob)) {
    stop("'prob' is non-numeric")
  }
  if (!identical(length(prob), 1L)) {
    stop("'prob' does not have length 1")
  }
  if (is.na(prob)) {
    stop("'prob' is NA")
  }
  if (prob <= 0) {
    stop("'prob' is less than or equal to 0")
  }
  if (prob >= 1) {
    stop("'prob' is greater than or equal to 1")
  }
  TRUE
}


## HAS_TESTS
## Check that a 'ratio' argument to a cdm constructor is valid
## 'ratio' - the ratio argument being checked, a numeric vector
## 'counts_data' - the cohort data being modelled
## return value - TRUE or an error
check_ratio_cdm <- function(ratio, counts_data) {
  check_arg_cdm(
    arg = ratio,
    nm_arg = "ratio",
    counts_data = counts_data,
    neg_ok = FALSE,
    zero_ok = TRUE
  )
}

## HAS_TESTS
## Check that a 'sd' argument to a cdm constructor is valid
## 'sd' - the sd argument being checked, a numeric vector
## 'counts_data' - the cohort data being modelled
## return value - TRUE or an error
check_sd_cdm <- function(sd, counts_data) {
  check_arg_cdm(
    arg = sd,
    nm_arg = "sd",
    counts_data = counts_data,
    neg_ok = FALSE,
    zero_ok = FALSE
  )
}
