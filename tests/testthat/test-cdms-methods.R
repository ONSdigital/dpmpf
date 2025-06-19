
## These tests are for internal C++ classes and methods

## calc_loglik ----------------------------------------------------------------

test_that("'calc_loglik' method for CdmsNoreg' works with non-empty models", {
    ## 10 intervals
    ## 5 particles
    counts_data_1 <- as.numeric(1:10) ## n_interval
    counts_data_2 <- as.numeric(3:12) ## n_interval
    counts_true <- as.numeric(3:7) ## n_particle
    prob <- 0.98
    ratio <- rep(1.1, 10)
    disp <- rep(0.1, 10)
    obs_zero <- 0.6
    i_interval <- 3L
    e1 <- new_CdmNoregPoibin(
        counts_data = counts_data_1,
        prob = prob
    )
    e2 <- new_CdmNoregNbinom(
        counts_data = counts_data_2,
        ratio = ratio,
        disp = disp
    )
    x <- new_CdmsNoreg(list(e1, e2))
    ans_obtained <- x$calc_loglik(
                          counts_true = counts_true,
                          i_interval = i_interval,
                          obs_zero = obs_zero
                      )
    ans_expected <- e1$calc_loglik(
                           counts_true = counts_true,
                           i_interval = i_interval,
                           obs_zero = obs_zero
                       ) +
        e2$calc_loglik(
               counts_true = counts_true,
               i_interval = i_interval,
               obs_zero = obs_zero
           )
    expect_identical(ans_obtained, ans_expected)
})

test_that("'calc_loglik' method for CdmsNoreg' works with empty models", {
  ## 10 intervals
  ## 5 particles
  counts_true <- as.numeric(3:7) ## n_particle
  i_interval <- 3L
  obs_zero <- 0.6
  x <- new_CdmsNoreg()
  ans_obtained <- x$calc_loglik(
    counts_true = counts_true,
    i_interval = i_interval,
    obs_zero = obs_zero
  )
  ans_expected <- rep(0, 5)
  expect_identical(ans_obtained, ans_expected)
})

test_that("'calc_loglik' method for 'CdmsWithreg' works with non-empty models", {
    ## 10 intervals
    ## 5 particles
    ## 2 regions
    counts_data_1 <- matrix(as.numeric(1:20), nrow = 2) ## n_interval
    counts_data_2 <- matrix(as.numeric(3:22), nrow = 2) ## n_interval
    counts_true <- matrix(as.numeric(3:12), nrow = 5) ## n_particle
    prob <- 0.98
    ratio <- matrix(1, nrow = 2, ncol = 10)
    disp <- matrix(1.1, nrow = 2, ncol = 10)
    obs_zero <- 0.6
    i_interval <- 3L
    e1 <- new_CdmWithregPoibin(
        counts_data = counts_data_1,
        prob = prob
    )
    e2 <- new_CdmWithregNbinom(
        counts_data = counts_data_2,
        ratio = ratio,
        disp = disp
    )
    x <- new_CdmsWithreg(list(e1, e2))
    ans_obtained <- x$calc_loglik(
                          counts_true = counts_true,
                          i_interval = i_interval,
                          obs_zero = obs_zero
                      )
    ans_expected <- e1$calc_loglik(
                           counts_true = counts_true,
                           i_interval = i_interval,
                           obs_zero = obs_zero
                       ) +
        e2$calc_loglik(
               counts_true = counts_true,
               i_interval = i_interval,
               obs_zero = obs_zero
           )
    expect_equal(ans_obtained, ans_expected)
})

test_that("'calc_loglik' method for 'CdmsWithreg' works with empty models", {
    ## 10 intervals
    ## 5 particles
    ## 2 regions
    counts_true <- matrix(as.numeric(3:12), nrow = 5) ## n_particle
    obs_zero <- 0.6
    i_interval <- 3L
    x <- new_CdmsWithreg()
    ans_obtained <- x$calc_loglik(
                          counts_true = counts_true,
                          i_interval = i_interval,
                          obs_zero = obs_zero
                      )
    ans_expected <- rep(0, 5)
    expect_equal(ans_obtained, ans_expected)
})


## draw_counts_true -----------------------------------------------------------

test_that("'draw_counts_true' method for CdmsNoreg' works - first cdm has data", {
  ## 10 intervals
  ## 5 particles
  counts_data_1 <- as.numeric(1:10) ## n_interval
  counts_data_2 <- as.numeric(3:12) ## n_interval
  prob <- 0.98
  ratio <- rep(1.1, 10)
  disp <- rep(1.2, 10)
  i_interval <- 3L
  e1 <- new_CdmNoregPoibin(
    counts_data = counts_data_1,
    prob = prob
  )
  e2 <- new_CdmNoregNbinom(
    counts_data = counts_data_2,
    ratio = ratio,
    disp = disp
  )
  x <- new_CdmsNoreg(list(e1, e2))
  for (seed in 1:20) {
    set.seed(seed)
    ans_obtained <- x$draw_counts_true(
      i_interval = 3L,
      n_particle = 5L
    )
    set.seed(seed)
    counts_true <- rep(NA_real_, times = 5)
    ## 'counts_true' modified in place
    e1$fill_counts_true(
      counts_true = counts_true,
      i_interval = i_interval
    )
    ans_expected <- counts_true
    expect_identical(ans_obtained, ans_expected)
  }
})

test_that("'draw_counts_true' method for 'CdmsNoreg' works - first cdm does not have data", {
  ## 10 intervals
  ## 5 particles
  counts_data_1 <- as.numeric(c(1:3, NA, 5:10)) ## n_interval
  counts_data_2 <- as.numeric(3:12) ## n_interval
  prob <- 0.98
  ratio <- rep(1.1, 10)
  disp <- rep(1.2, 10)
  i_interval <- 3L
  e1 <- new_CdmNoregPoibin(
    counts_data = counts_data_1,
    prob = prob
  )
  e2 <- new_CdmNoregNbinom(
    counts_data = counts_data_2,
    ratio = ratio,
    disp = disp
  )
  x <- new_CdmsNoreg(list(e1, e2))
  for (seed in 1:20) {
    set.seed(seed)
    ans_obtained <- x$draw_counts_true(
      i_interval = i_interval,
      n_particle = 5L
    )
    set.seed(seed)
    counts_true <- rep(NA_real_, times = 5)
    ## 'counts_true' modified in place
    e2$fill_counts_true(
      counts_true = counts_true,
      i_interval = i_interval
    )
    ans_expected <- counts_true
    expect_identical(ans_obtained, ans_expected)
  }
})

test_that("'draw_counts_true' method for 'CdmsNoreg' works - no data", {
  ## 10 intervals
  ## 5 particles
  i_interval <- 3L
  x <- new_CdmsNoreg()
  ans_obtained <- x$draw_counts_true(
    i_interval = i_interval,
    n_particle = 5L
  )
  ans_expected <- rep(NA_real_, times = 5)
  expect_identical(ans_obtained, ans_expected)
})

test_that("'draw_counts_true' method for 'CdmsWithreg' works - first cdm has complete data", {
  ## 10 intervals
  ## 5 particles
  ## 2 regions
  counts_data_1 <- matrix(as.numeric(1:20), nrow = 2) ## n_region x n_interval
  counts_data_2 <- matrix(as.numeric(3:22), nrow = 2) ## n_region x n_interval
  prob <- 0.98
  ratio <- matrix(1, nrow = 2, ncol = 10)
  disp <- matrix(1.1, nrow = 2, ncol = 10)
  i_interval <- 3L
  n_particle <- 5L
  n_region <- 2L
  e1 <- new_CdmWithregPoibin(
    counts_data = counts_data_1,
    prob = prob
  )
  e2 <- new_CdmWithregNbinom(
    counts_data = counts_data_2,
    ratio = ratio,
    disp = disp
  )
  x <- new_CdmsWithreg(list(e1, e2))
  for (seed in 1:20) {
    set.seed(seed)
    ans_obtained <- x$draw_counts_true(
      i_interval = i_interval,
      n_particle = n_particle,
      n_region = n_region
    )
    set.seed(seed)
    counts_true <- matrix(NA_real_, nrow = 5, ncol = 2L) ## n_particle x n_region
    ## 'counts_true' modified in place
    e1$fill_counts_true(
      counts_true = counts_true,
      i_interval = i_interval
    )
    ans_expected <- counts_true
    expect_identical(ans_obtained, ans_expected)
  }
})

test_that("'draw_counts_true' method for 'CdmsWithreg' works - first cdm has partial data", {
  ## 10 intervals
  ## 5 particles
  ## 2 regions
  i_interval <- 3L
  n_particle <- 5L
  n_region <- 2L
  counts_data_1 <- matrix(as.numeric(1:20), nrow = 2) ## n_region x n_interval
  counts_data_1[1, i_interval + 1] <- NA
  counts_data_2 <- matrix(as.numeric(3:22), nrow = 2) ## n_region x n_interval
  prob <- 0.98
  ratio <- matrix(1, nrow = 2, ncol = 10)
  disp <- matrix(1.1, nrow = 2, ncol = 10)
  e1 <- new_CdmWithregPoibin(
    counts_data = counts_data_1,
    prob = prob
  )
  e2 <- new_CdmWithregNbinom(
    counts_data = counts_data_2,
    ratio = ratio,
    disp = disp
  )
  x <- new_CdmsWithreg(list(e1, e2))
  for (seed in 1:20) {
    set.seed(seed)
    ans_obtained <- x$draw_counts_true(
      i_interval = i_interval,
      n_particle = n_particle,
      n_region = n_region
    )
    set.seed(seed)
    counts_true <- matrix(NA_real_, nrow = 5, ncol = 2L) ## n_particle x n_region
    ## 'counts_true' is modified in place
    e1$fill_counts_true(
      counts_true = counts_true,
      i_interval = i_interval
    )
    e2$fill_counts_true(
      counts_true = counts_true,
      i_interval = i_interval
    )
    ans_expected <- counts_true
    expect_identical(ans_obtained, ans_expected)
  }
})

test_that("'draw_counts_true' method for 'CdmsWithreg' works - no data", {
  ## 10 intervals
  ## 5 particles
  ## 2 regions
  i_interval <- 3L
  n_particle <- 5L
  n_region <- 2L
  x <- new_CdmsWithreg()
  ans_obtained <- x$draw_counts_true(
    i_interval = i_interval,
    n_particle = n_particle,
    n_region = n_region
  )
  ans_expected <- matrix(NA_real_,
    nrow = n_particle,
    ncol = n_region
  )
  expect_identical(ans_obtained, ans_expected)
})


## calc_logimp ---------------------------------------------------------------

test_that("'calc_logimp' method for CdmsNoreg' works - first cdm has data", {
  ## 10 intervals
  ## 5 particles
  i_interval <- 3L
  n_particle <- 5L
  prob <- 0.98
  counts_data_1 <- as.numeric(1:10) ## n_interval
  counts_data_2 <- as.numeric(3:12) ## n_interval
  counts_true <- as.numeric(0:4) ## n_particle
  ratio <- rep(1.1, 10)
  disp <- rep(0.1, 10)
  e1 <- new_CdmNoregPoibin(
    counts_data = counts_data_1,
    prob = prob
  )
  e2 <- new_CdmNoregNbinom(
    counts_data = counts_data_2,
    ratio = ratio,
    disp = disp
  )
  x <- new_CdmsNoreg(list(e1, e2))
  for (seed in 1:20) {
    set.seed(seed)
    ans_obtained <- x$calc_logimp(
      counts_true = counts_true,
      i_interval = i_interval
    )
    set.seed(seed)
    logimp <- rep(NA_real_, times = n_particle)
    ## 'logimp' is modified in place
    e1$fill_logimp(
      logimp = logimp,
      counts_true = counts_true,
      i_interval = i_interval
    )
    ans_expected <- logimp
    expect_identical(ans_obtained, ans_expected)
  }
})

test_that("'calc_logimp' method for CdmsNoreg' works - first cdm missing data", {
  ## 10 intervals
  ## 5 particles
  i_interval <- 3L
  n_particle <- 5L
  prob <- 0.98
  ratio <- rep(1.1, 10)
  disp <- rep(0.1, 10)
  counts_data_1 <- as.numeric(c(1:3, NA, 5:10)) ## n_interval
  counts_data_2 <- as.numeric(3:12) ## n_interval
  counts_true <- as.numeric(0:4) ## n_particle
  e1 <- new_CdmNoregPoibin(
    counts_data = counts_data_1,
    prob = prob
  )
  e2 <- new_CdmNoregNbinom(
    counts_data = counts_data_2,
    ratio = ratio,
    disp = disp
  )
  x <- new_CdmsNoreg(list(e1, e2))
  for (seed in 1:20) {
    set.seed(seed)
    ans_obtained <- x$calc_logimp(
      counts_true = counts_true,
      i_interval = i_interval
    )
    set.seed(seed)
    logimp <- rep(NA_real_, times = n_particle)
    ## 'logimp' is modified in place
    e2$fill_logimp(
      logimp = logimp,
      counts_true = counts_true,
      i_interval = i_interval
    )
    ans_expected <- logimp
    expect_identical(ans_obtained, ans_expected)
  }
})

test_that("'calc_logimp' method for CdmsNoreg' works - no data", {
  ## 10 intervals
  ## 5 particles
  i_interval <- 3L
  n_particle <- 5L
  counts_true <- as.numeric(0:4) ## n_particle
  x <- new_CdmsNoreg()
  ans_obtained <- x$calc_logimp(
    counts_true = counts_true,
    i_interval = i_interval
  )
  ans_expected <- rep(NA_real_, 5)
  expect_identical(ans_obtained, ans_expected)
})

test_that("'calc_logimp' method for 'CdmsWithreg' works - first cdm has complete data", {
  ## 10 intervals
  ## 5 particles
  ## 2 regions
  i_interval <- 3L
  n_particle <- 5L
  prob <- 0.98
  ratio <- matrix(1, nrow = 2, ncol = 10)
  disp <- matrix(1.1, nrow = 2, ncol = 10)
  counts_data_1 <- matrix(as.numeric(1:20), nrow = 2) ## n_region x n_interval
  counts_data_2 <- matrix(as.numeric(3:22), nrow = 2) ## n_region x n_interval
  counts_true <- matrix(as.numeric(0:9), nrow = 5) ## n_particle x n_region
  e1 <- new_CdmWithregPoibin(
    counts_data = counts_data_1,
    prob = prob
  )
  e2 <- new_CdmWithregNbinom(
    counts_data = counts_data_2,
    ratio = ratio,
    disp = disp
  )
  x <- new_CdmsWithreg(list(e1, e2))
  for (seed in 1:20) {
    set.seed(seed)
    ans_obtained <- x$calc_logimp(
      counts_true = counts_true,
      i_interval = i_interval
    )
    set.seed(seed)
    logimp <- matrix(NA_real_, nrow = 5, ncol = 2)
    ## 'logimp' is modified in place
    e1$fill_logimp(
      logimp = logimp,
      counts_true = counts_true,
      i_interval = i_interval
    )
    ans_expected <- logimp
    expect_identical(ans_obtained, ans_expected)
  }
})

test_that("'calc_logimp' method for 'CdmsWithreg' works - first cdm has partial data", {
  ## 10 intervals
  ## 5 particles
  ## 2 regions
  i_interval <- 3L
  n_particle <- 5L
  prob <- 0.98
  ratio <- matrix(1, nrow = 2, ncol = 10)
  disp <- matrix(1.1, nrow = 2, ncol = 10)
  counts_data_1 <- matrix(as.numeric(1:20), nrow = 2) ## n_region x n_interval
  counts_data_1[1, i_interval + 1] <- NA
  counts_data_2 <- matrix(as.numeric(3:22), nrow = 2) ## n_region x n_interval
  counts_true <- matrix(as.numeric(0:9), nrow = 5) ## n_particle x n_region
  e1 <- new_CdmWithregPoibin(
    counts_data = counts_data_1,
    prob = prob
  )
  e2 <- new_CdmWithregNbinom(
    counts_data = counts_data_2,
    ratio = ratio,
    disp = disp
  )
  x <- new_CdmsWithreg(list(e1, e2))
  for (seed in 1:20) {
    set.seed(seed)
    ans_obtained <- x$calc_logimp(
      counts_true = counts_true,
      i_interval = i_interval
    )
    set.seed(seed)
    logimp <- matrix(NA_real_, nrow = n_particle, ncol = 2)
    e1$fill_logimp(
      logimp = logimp,
      counts_true = counts_true,
      i_interval = i_interval
    )
    e2$fill_logimp(
      logimp = logimp,
      counts_true = counts_true,
      i_interval = i_interval
    )
    ans_expected <- logimp
    expect_identical(ans_obtained, ans_expected)
  }
})

test_that("'calc_logimp' method for 'CdmsWithreg' works - no data", {
  ## 10 intervals
  ## 5 particles
  ## 2 regions
  i_interval <- 3L
  n_particle <- 5L
  counts_true <- matrix(as.numeric(1:10), nrow = 5) ## n_particle x n_region
  x <- new_CdmsWithreg()
  ans_obtained <- x$calc_logimp(
    counts_true = counts_true,
    i_interval = i_interval
  )
  ans_expected <- matrix(NA_real_, nrow = n_particle, ncol = 2L)
  expect_identical(ans_obtained, ans_expected)
})
