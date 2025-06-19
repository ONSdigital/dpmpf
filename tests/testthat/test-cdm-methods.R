
## calc_loglik - Poisson-binomial no regions ----------------------------------

test_that("'calc_loglik method for CdmNoregPoibin' works - no NA", {
    ## 10 intervals
    ## 5 particles
    counts_true <- as.numeric(0:4) ## n_particle
    counts_data <- as.numeric(1:10) ## n_interval
    x <- new_CdmNoregPoibin(
        counts_data = counts_data,
        prob = 0.98
    )
    ans_obtained <- x$calc_loglik(
                          counts_true = counts_true,
                          i_interval = 3L,
                          obs_zero = 0.5
                      )
    ans_expected <- numeric(5)
    for (i in 1:5) {
        if (counts_true[[i]] >= 1)
            ans_expected[[i]] <-
                dpoibin1(counts_data[3L + 1L], ## i_interval is a C-style, not R-style, index
                         size = counts_true[[i]],
                         prob = 0.98,
                         use_log = TRUE
                         )
        else
            ans_expected[[i]] <-
                dpois(counts_data[3L + 1L], ## i_interval is a C-style, not R-style, index
                      lambda = 0.5,
                      log = TRUE
                      )
    }
    expect_identical(ans_obtained, ans_expected)
})

test_that("'calc_loglik method for CdmNoregPoibin' works - with NA", {
  ## 10 intervals
  ## 5 particles
  counts_true <- as.numeric(0:4) ## n_particle
  counts_data <- as.numeric(1:10) ## n_interval
  counts_data[4] <- NA
  x <- new_CdmNoregPoibin(
    counts_data = counts_data,
    prob = 0.98
  )
  ans_obtained <- x$calc_loglik(
    counts_true = counts_true,
    i_interval = 3L,
    obs_zero = 0.5
  )
  ans_expected <- rep(0, 5)
  expect_identical(ans_obtained, ans_expected)
})


## calc_loglik - negative-binomial no regions ---------------------------------

test_that("'calc_loglik method for CdmNoregNbinom' works - no NA", {
    ## 10 intervals
    ## 5 particles
    counts_true <- as.numeric(0:4) ## n_particle
    counts_data <- as.numeric(1:10) ## n_interval
    ratio <- seq(from = 0.9, to = 1.3, length.out = 10)
    disp <- seq(from = 1.5, to = 1.3, length.out = 10)
    x <- new_CdmNoregNbinom(
        counts_data = counts_data,
        ratio = ratio,
        disp = disp
    )
    ans_obtained <- x$calc_loglik(
                          counts_true = counts_true,
                          i_interval = 3L,
                          obs_zero = 0.7
                      )
    ans_expected <- numeric(5)
    for (i in 1:5) {
        if (counts_true[[i]] > 0)
            ans_expected[[i]] <- dnbinom(counts_data[3L + 1L], ## i_interval is a C-style, not R-style, index
                                         size = 1 / disp[[3L + 1L]],
                                         mu = ratio[3L + 1L] * counts_true[[i]],
                                         log = TRUE
                                         )
        else
            ans_expected[[i]] <- dpois(counts_data[3L + 1L], ## i_interval is a C-style, not R-style, index
                                       lambda = 0.7,
                                       log = TRUE
                                       )
    }
    expect_identical(ans_obtained, ans_expected)
})

test_that("'calc_loglik method for CdmNoregNbinom' works - with NA", {
  ## 10 intervals
  ## 5 particles
  counts_true <- as.numeric(0:4) ## n_particle
  counts_data <- as.numeric(1:10) ## n_interval
  counts_data[3 + 1] <- NA
  ratio <- seq(from = 0.9, to = 1.3, length.out = 10)
  disp <- seq(from = 1.5, to = 1.3, length.out = 10)
  x <- new_CdmNoregNbinom(
    counts_data = counts_data,
    ratio = ratio,
    disp = disp
  )
  ans_obtained <- x$calc_loglik(
    counts_true = counts_true,
    i_interval = 3L,
    obs_zero = 0.7
  )
  ans_expected <- rep(0, 5)
  expect_identical(ans_obtained, ans_expected)
})

test_that("'calc_loglik method for CdmNoregNbinom' works - dispersion equals 0", {
    ## 10 intervals
    ## 5 particles
    counts_true <- as.numeric(0:4) ## n_particle
    counts_data <- as.numeric(1:10) ## n_interval
    ratio <- seq(from = 0.9, to = 1.3, length.out = 10)
    disp <- rep(0, 10)
    x <- new_CdmNoregNbinom(
        counts_data = counts_data,
        ratio = ratio,
        disp = disp
    )
    ans_obtained <- x$calc_loglik(
                          counts_true = counts_true,
                          i_interval = 3L,
                          obs_zero = 0.7
                      )
    ans_expected <- numeric(5)
    for (i in 1:5) {
        if (counts_true[[i]] > 0)
            ans_expected[[i]] <- dpois(counts_data[3L + 1L], ## i_interval is a C-style, not R-style, index
                                       lambda = ratio[3L + 1L] * counts_true[[i]],
                                       log = TRUE
                                       )                
        else
            ans_expected[[i]] <- dpois(counts_data[3L + 1L], ## i_interval is a C-style, not R-style, index
                                       lambda = 0.7,
                                       log = TRUE
                                       )
    }
    expect_identical(ans_obtained, ans_expected)
})


## calc_loglik - normal no regions ---------------------------------

test_that("'calc_loglik method for CdmNoregNorm' works - no NA", {
    ## 10 intervals
    ## 5 particles
    counts_true <- as.numeric(0:4) ## n_particle
    counts_data <- as.numeric(1:10) ## n_interval
    ratio <- seq(from = 0.9, to = 1.3, length.out = 10)
    sd <- seq(from = 1.5, to = 1.3, length.out = 10)
    x <- new_CdmNoregNorm(
        counts_data = counts_data,
        ratio = ratio,
        sd = sd
    )
    ans_obtained <- x$calc_loglik(
                          counts_true = counts_true,
                          i_interval = 3L,
                          obs_zero = 0.7
                      )
    ans_expected <- numeric(5)
    for (i in 1:5) {
        if (counts_true[[i]] > 0)
            ans_expected[[i]] <- log(pnorm(counts_data[3L + 1L] + 0.5, ## i_interval is a C-style, not R-style, index
                                           mean = ratio[3L + 1L] * counts_true[[i]],
                                           sd = sd[3L + 1L])
                                     - pnorm(counts_data[3L + 1L] - 0.5, ## i_interval is a C-style, not R-style, index
                                             mean = ratio[3L + 1L] * counts_true[[i]],
                                             sd = sd[3L + 1L]))
        else
            ans_expected[[i]] <- dpois(counts_data[3L + 1L], ## i_interval is a C-style, not R-style, index
                                       lambda = 0.7,
                                       log = TRUE
                                       )
    }
    expect_equal(ans_obtained, ans_expected)
})

test_that("'calc_loglik method for CdmNoregNorm' works - with NA", {
  ## 10 intervals
  ## 5 particles
  counts_true <- as.numeric(0:4) ## n_particle
  counts_data <- as.numeric(1:10) ## n_interval
  counts_data[3 + 1] <- NA
  ratio <- seq(from = 0.9, to = 1.3, length.out = 10)
  sd <- seq(from = 1.5, to = 1.3, length.out = 10)
  x <- new_CdmNoregNorm(
    counts_data = counts_data,
    ratio = ratio,
    sd = sd
  )
  ans_obtained <- x$calc_loglik(
    counts_true = counts_true,
    i_interval = 3L,
    obs_zero = 0.7
  )
  ans_expected <- rep(0, 5)
  expect_identical(ans_obtained, ans_expected)
})


## calc_loglik - Poisson-binomial, with regions -------------------------------

test_that("'calc_loglik method for CdmWithregPoibin' works", {
    ## 10 intervals
    ## 2 regions
    ## 5 particles
    counts_true <- matrix(as.numeric(0:9), nr = 5) ## n_particle x n_region
    counts_data <- matrix(as.numeric(1:20), nr = 2) ## n_region x n_interval
    prob <- 0.98
    obs_zero <- 0.6
    i_interval <- 3L
    x <- new_CdmWithregPoibin(
        counts_data = counts_data,
        prob = 0.98
    )
    ans_obtained <- x$calc_loglik(
                          counts_true = counts_true,
                          i_interval = 3L,
                          obs_zero = obs_zero
                      )
    ans_expected <- numeric(5) ## n_particle
    for (i in 1:5) {
        if (counts_true[i, 1] > 0) ## only first row has zero
            ans_expected[[i]] <- dpoibin1(counts_data[1, 3 + 1], ## switch to R indexing
                                          size = counts_true[i, 1],
                                          prob = 0.98,
                                          use_log = TRUE
                                          ) +
                dpoibin1(counts_data[2, 3 + 1],
                         size = counts_true[i, 2],
                         prob = 0.98,
                         use_log = TRUE
                         )
        else
            ans_expected[[i]] <- dpois(counts_data[1, 3 + 1], ## switch to R indexing
                                       lambda = obs_zero,
                                       log = TRUE
                                       ) +
                dpoibin1(counts_data[2, 3 + 1],
                         size = counts_true[i, 2],
                         prob = 0.98,
                         use_log = TRUE
                         )
    }
    expect_identical(ans_obtained, ans_expected)
})


## calc_loglik - negative binomial, with regions ------------------------------

test_that("'calc_loglik method for CdmWithregNbinom' works - no NA", {
    ## 10 intervals
    ## 2 regions
    ## 5 particles
    counts_true <- matrix(as.numeric(0:9), nr = 5) ## n_particle x n_region
    counts_data <- matrix(as.numeric(1:20), nr = 2) ## n_region x n_interval
    ratio <- matrix(1.1, nrow = 2, ncol = 10) ## n_region x n_interval
    disp <- matrix(1.3, nrow = 2, ncol = 10) ## n_region x n_interval
    i_interval <- 3L
    obs_zero <- 0.6
    x <- new_CdmWithregNbinom(
        counts_data = counts_data,
        ratio = ratio,
        disp = disp
    )
    ans_obtained <- x$calc_loglik(
                          counts_true = counts_true,
                          i_interval = 3L,
                          obs_zero = obs_zero
                      )
    ans_expected <- numeric(5) ## n_particle
    for (i in 1:5) { ## zero only in row 1
        if (counts_true[i, 1] > 0) {
            ans_expected[[i]] <- dnbinom(counts_data[1, 3 + 1], ## switch to R indexing
                                         size = 1 / disp[1, 3 + 1],
                                         mu = ratio[1, 3 + 1] * counts_true[i, 1],
                                         log = TRUE
                                         ) +
                dnbinom(counts_data[2, 3 + 1],
                        size = 1 / disp[2, 3 + 1],
                        mu = ratio[2, 3 + 1] * counts_true[i, 2],
                        log = TRUE
                        )
        }
        else {
            ans_expected[[i]] <- dpois(counts_data[1, 3 + 1], ## switch to R indexing
                                       lambda = obs_zero,
                                       log = TRUE
                                       ) +
                dnbinom(counts_data[2, 3 + 1],
                        size = 1 / disp[2, 3 + 1],
                        mu = ratio[2, 3 + 1] * counts_true[i, 2],
                        log = TRUE
                        )
        }
    }
    expect_identical(ans_obtained, ans_expected)
})

test_that("'calc_loglik method for CdmWithregNbinom' works - with NA", {
  ## 10 intervals
  ## 2 regions
  ## 5 particles
  counts_true <- matrix(as.numeric(0:9), nr = 5) ## n_particle x n_region
  counts_data <- matrix(as.numeric(1:20), nr = 2) ## n_region x n_interval
  counts_data[2, 4] <- NA
  ratio <- matrix(1.1, nrow = 2, ncol = 10) ## n_region x n_interval
  disp <- matrix(1.3, nrow = 2, ncol = 10) ## n_region x n_interval
  obs_zero <- 0.8
  i_interval <- 3L
  x <- new_CdmWithregNbinom(
    counts_data = counts_data,
    ratio = ratio,
    disp = disp
  )
  ans_obtained <- x$calc_loglik(
    counts_true = counts_true,
    i_interval = 3L,
    obs_zero = obs_zero
  )
  ans_expected <- numeric(5) ## n_particle
  for (i in 1:5) {
      if (counts_true[i, 1] > 0)
          ans_expected[[i]] <- dnbinom(counts_data[1, 3 + 1], ## switch to R indexing
                                       size = 1 / disp[1, 3 + 1],
                                       mu = ratio[1, 3 + 1] * counts_true[i, 1],
                                       log = TRUE
                                       )
      else
          ans_expected[[i]] <- dpois(counts_data[1, 3 + 1], ## switch to R indexing
                                     lambda = obs_zero,
                                     log = TRUE
                                     )
  }
  expect_identical(ans_obtained, ans_expected)
})

test_that("'calc_loglik method for CdmWithregNbinom' works - disp is 0", {
    ## 10 intervals
    ## 2 regions
    ## 5 particles
    counts_true <- matrix(as.numeric(0:9), nr = 5) ## n_particle x n_region
    counts_data <- matrix(as.numeric(1:20), nr = 2) ## n_region x n_interval
    ratio <- matrix(1.1, nrow = 2, ncol = 10) ## n_region x n_interval
    disp <- matrix(0, nrow = 2, ncol = 10) ## n_region x n_interval
    i_interval <- 3L
    obs_zero <- 0.6
    x <- new_CdmWithregNbinom(
        counts_data = counts_data,
        ratio = ratio,
        disp = disp
    )
    ans_obtained <- x$calc_loglik(
                          counts_true = counts_true,
                          i_interval = 3L,
                          obs_zero = obs_zero
                      )
    ans_expected <- numeric(5) ## n_particle
    for (i in 1:5) { ## zero only in row 1
        if (counts_true[i, 1] > 0) {
            ans_expected[[i]] <- dpois(counts_data[1, 3 + 1], ## switch to R indexing
                                       lambda = ratio[1, 3 + 1] * counts_true[i, 1],
                                       log = TRUE
                                       ) +
                dpois(counts_data[2, 3 + 1],
                      lambda = ratio[2, 3 + 1] * counts_true[i, 2],
                      log = TRUE
                      )
        }
        else {
            ans_expected[[i]] <- dpois(counts_data[1, 3 + 1], ## switch to R indexing
                                       lambda = obs_zero,
                                       log = TRUE
                                       ) +
                dnbinom(counts_data[2, 3 + 1],
                        size = 1 / disp[2, 3 + 1],
                        mu = ratio[2, 3 + 1] * counts_true[i, 2],
                        log = TRUE
                        )
        }
    }
    expect_identical(ans_obtained, ans_expected)
})


## fill_counts_true - Poisson-binomial, no regions ----------------------------

test_that("'fill_counts_true' method for 'CdmNoregPoibin' works - counts_data not NA", {
  ## 10 intervals
  ## 5 particles
  counts_true <- rep(NA_real_, times = 5) ## n_particle
  counts_data <- as.numeric(1:10) ## n_interval
  prob <- 0.98
  x <- new_CdmNoregPoibin(
    counts_data = counts_data,
    prob = prob
  )
  for (seed in 1:5) {
    set.seed(seed)
    ## modify 'counts_true' in place
    x$fill_counts_true(
      counts_true = counts_true,
      i_interval = 3L
    )
    ans_obtained <- counts_true
    ans_expected <- numeric(5) ## n_particle
    set.seed(seed)
    for (i in 1:5) {
      ans_expected[[i]] <- rpoibin1(
        size = counts_data[3 + 1],
        prob = prob
      )
    }
    expect_identical(ans_obtained, ans_expected)
  }
})

test_that("'fill_counts_true' method for 'CdmNoregPoibin' works - counts_data is NA", {
  ## 10 intervals
  ## 5 particles
  counts_data <- as.numeric(c(1:3, NA_real_, 5:10)) ## n_interval
  counts_true <- rep(NA_real_, times = 5) ## n_particle
  prob <- 0.98
  i_interval <- 3L
  x <- new_CdmNoregPoibin(
    counts_data = counts_data,
    prob = prob
  )
  ## modify 'counts_true' in place - but nothing should happen because data is NA
  x$fill_counts_true(
    counts_true = counts_true,
    i_interval = i_interval
  )
  ans_obtained <- counts_true
  ans_expected <- counts_true <- rep(NA_real_, times = 5)
  expect_identical(ans_obtained, ans_expected)
})


## fill_counts_true - negative binomial, no regions ---------------------------

test_that("'fill_counts_true' method for 'CdmNoregNbinom' works - counts_data not NA", {
  ## 10 intervals
  ## 5 particles
  set.seed(0)
  counts_true <- rep(NA_real_, times = 5) ## n_particle
  counts_data <- as.numeric(1:10) ## n_interval
  ratio <- seq(from = 0.9, to = 1.3, length.out = 10)
  disp <- seq(from = 1.5, to = 1.3, length.out = 10)
  x <- new_CdmNoregNbinom(
    counts_data = counts_data,
    ratio = ratio,
    disp = disp
  )
  for (seed in 1:5) {
    set.seed(seed)
    ## modify 'counts_true' in place
    x$fill_counts_true(
      counts_true = counts_true,
      i_interval = 3L
    )
    ans_obtained <- counts_true
    ans_expected <- numeric(5) ## n_particle
    set.seed(seed)
    for (i in 1:5) {
      ans_expected[[i]] <- rnbinom(
        n = 1L,
        size = 1 / disp[3 + 1],
        mu = counts_data[3 + 1] / ratio[3 + 1]
      )
    }
    expect_identical(ans_obtained, ans_expected)
  }
})

test_that("'fill_counts_true' method for 'CdmNoregNbinom' works - counts_data is NA", {
  ## 10 intervals
  ## 5 particles
  counts_true <- rep(NA_real_, times = 5) ## n_particle
  counts_data <- as.numeric(1:10) ## n_interval
  counts_data[4] <- NA
  ratio <- seq(from = 0.9, to = 1.3, length.out = 10)
  disp <- seq(from = 1.5, to = 1.3, length.out = 10)
  x <- new_CdmNoregNbinom(
    counts_data = counts_data,
    ratio = ratio,
    disp = disp
  )
  for (seed in 1:5) {
    set.seed(seed)
    ## modify 'counts_true' in place
    x$fill_counts_true(
      counts_true = counts_true,
      i_interval = 3L
    )
    ans_obtained <- counts_true
    ans_expected <- rep(NA_real_, 5)
    expect_identical(ans_obtained, ans_expected)
  }
})

test_that("'fill_counts_true' method for 'CdmNoregNbinom' works - disp is 0", {
  ## 10 intervals
  ## 5 particles
  set.seed(0)
  counts_true <- rep(NA_real_, times = 5) ## n_particle
  counts_data <- as.numeric(1:10) ## n_interval
  ratio <- seq(from = 0.9, to = 1.3, length.out = 10)
  disp <- rep(0, 10)
  x <- new_CdmNoregNbinom(
    counts_data = counts_data,
    ratio = ratio,
    disp = disp
  )
  for (seed in 1:5) {
    set.seed(seed)
    ## modify 'counts_true' in place
    x$fill_counts_true(
      counts_true = counts_true,
      i_interval = 3L
    )
    ans_obtained <- counts_true
    ans_expected <- numeric(5) ## n_particle
    set.seed(seed)
    for (i in 1:5) {
      ans_expected[[i]] <- rpois(
        n = 1L,
        lambda = counts_data[3 + 1] / ratio[3 + 1]
      )
    }
    expect_identical(ans_obtained, ans_expected)
  }
})


## fill_counts_true - normal, no regions ---------------------------

test_that("'fill_counts_true' method for 'CdmNoregNorm' works - counts_data not NA", {
  ## 10 intervals
  ## 5 particles
  set.seed(0)
  counts_true <- rep(NA_real_, times = 5) ## n_particle
  counts_data <- as.numeric(1:10) ## n_interval
  ratio <- seq(from = 0.9, to = 1.3, length.out = 10)
  sd <- seq(from = 1.5, to = 1.3, length.out = 10)
  x <- new_CdmNoregNorm(
    counts_data = counts_data,
    ratio = ratio,
    sd = sd
  )
  for (seed in 1:5) {
    set.seed(seed)
    ## modify 'counts_true' in place
    x$fill_counts_true(
      counts_true = counts_true,
      i_interval = 3L
    )
    ans_obtained <- counts_true
    ans_expected <- numeric(5) ## n_particle
    set.seed(seed)
    for (i in 1:5) {
      ans_expected[[i]] <- as.integer(rnorm(
        n = 1L,
        mean = counts_data[3 + 1] / ratio[3 + 1],
        sd = sd[3 + 1] 
      ) + 0.5)
    }
    expect_identical(ans_obtained, ans_expected)
  }
})

test_that("'fill_counts_true' method for 'CdmNoregNorm' works - counts_data is NA", {
  ## 10 intervals
  ## 5 particles
  counts_true <- rep(NA_real_, times = 5) ## n_particle
  counts_data <- as.numeric(1:10) ## n_interval
  counts_data[4] <- NA
  ratio <- seq(from = 0.9, to = 1.3, length.out = 10)
  sd <- seq(from = 1.5, to = 1.3, length.out = 10)
  x <- new_CdmNoregNorm(
    counts_data = counts_data,
    ratio = ratio,
    sd = sd
  )
  for (seed in 1:5) {
    set.seed(seed)
    ## modify 'counts_true' in place
    x$fill_counts_true(
      counts_true = counts_true,
      i_interval = 3L
    )
    ans_obtained <- counts_true
    ans_expected <- rep(NA_real_, 5)
    expect_identical(ans_obtained, ans_expected)
  }
})


## fill_counts_true - Poisson-binomial, with regions --------------------------

test_that("'fill_counts_true' method for 'CdmWithregPoibin' works - counts_data not NA", {
  ## 10 intervals
  ## 5 particles
  ## 2 regions
  counts_data <- matrix(as.numeric(1:20), nr = 2) ## n_region x n_interval
  prob <- 0.98
  i_interval <- 3L
  x <- new_CdmWithregPoibin(
    counts_data = counts_data,
    prob = prob
  )
  for (seed in 1:5) {
    counts_true <- matrix(NA_real_, nrow = 5, ncol = 2) ## n_particle x n_region
    set.seed(seed)
    ## modify 'counts_true' in place
    x$fill_counts_true(
      counts_true = counts_true,
      i_interval = i_interval
    )
    ans_obtained <- counts_true
    ans_expected <- matrix(nrow = 5, ncol = 2) ## n_particle x n_region
    set.seed(seed)
    for (i_region in 1:2) {
      for (i_particle in 1:5) {
        ans_expected[i_particle, i_region] <- rpoibin1(
          size = counts_data[i_region, i_interval + 1],
          prob = prob
        )
      }
    }
    expect_identical(ans_obtained, ans_expected)
  }
})

test_that("'fill_counts_true' method for CdmWithregPoibin' works - counts_true partly full", {
  ## 10 intervals
  ## 5 particles
  ## 2 regions
  counts_data <- matrix(as.numeric(1:20), nr = 2) ## n_region x n_interval
  counts_data[1, 3] <- NA
  prob <- 0.98
  i_interval <- 3L
  x <- new_CdmWithregPoibin(
    counts_data = counts_data,
    prob = prob
  )
  for (seed in 1:5) {
    counts_true <- matrix(NA_real_, nrow = 5, ncol = 2) ## n_particle x n_region
    counts_true[, 1] <- 1:5
    set.seed(seed)
    ## modify 'counts_true' in place
    x$fill_counts_true(
      counts_true = counts_true,
      i_interval = i_interval
    )
    ans_obtained <- counts_true
    ans_expected <- matrix(nrow = 5, ncol = 2) ## n_particle x n_region
    ans_expected[, 1] <- as.numeric(1:5)
    set.seed(seed)
    for (i_particle in 1:5) {
      ans_expected[i_particle, 2] <- rpoibin1(
        size = counts_data[2, i_interval + 1],
        prob = prob
      )
    }
    expect_identical(ans_obtained, ans_expected)
  }
})


## fill_counts_true - negative binomial, with regions -------------------------

test_that("'fill_counts_true' method for 'CdmWithregNbinom' works - counts_data not NA", {
  ## 10 intervals
  ## 5 particles
  ## 2 regions
  set.seed(0)
  counts_data <- matrix(as.numeric(1:20), nr = 2) ## n_region x n_interval
  ratio <- matrix(runif(n = 20, min = 0.5, max = 1.5), nr = 2)
  disp <- matrix(runif(n = 20, min = 0.5, max = 1.5), nr = 2)
  i_interval <- 3L
  x <- new_CdmWithregNbinom(
    counts_data = counts_data,
    ratio = ratio,
    disp = disp
  )
  for (seed in 1:5) {
    counts_true <- matrix(NA_real_, nrow = 5, ncol = 2) ## n_particle x n_region
    set.seed(seed)
    ## modify 'counts_true' in place
    x$fill_counts_true(
      counts_true = counts_true,
      i_interval = i_interval
    )
    ans_obtained <- counts_true
    ans_expected <- matrix(nrow = 5, ncol = 2) ## n_particle x n_region
    set.seed(seed)
    for (i_region in 1:2) {
      for (i_particle in 1:5) {
        size <- 1 / disp[i_region, i_interval + 1]
        mu <- counts_data[i_region, i_interval + 1] / ratio[i_region, i_interval + 1]
        ans_expected[i_particle, i_region] <- rnbinom(
          n = 1L,
          size = size,
          mu = mu
        )
      }
    }
    expect_identical(ans_obtained, ans_expected)
  }
})

test_that("'fill_counts_true' method for CdmWithregNbinom' works - counts_true partly full", {
  ## 10 intervals
  ## 5 particles
  ## 2 regions
  set.seed(0)
  counts_data <- matrix(as.numeric(1:20), nr = 2) ## n_region x n_interval
  counts_data[1, 3] <- NA
  ratio <- matrix(runif(n = 20, min = 0.5, max = 1.5), nr = 2)
  disp <- matrix(runif(n = 20, min = 0.5, max = 1.5), nr = 2)
  i_interval <- 3L
  x <- new_CdmWithregNbinom(
    counts_data = counts_data,
    ratio = ratio,
    disp = disp
  )
  for (seed in 1:5) {
    counts_true <- matrix(NA_real_, nrow = 5, ncol = 2) ## n_particle x n_region
    counts_true[, 1] <- 1:5
    set.seed(seed)
    ## modify 'counts_true' in place
    x$fill_counts_true(
      counts_true = counts_true,
      i_interval = i_interval
    )
    ans_obtained <- counts_true
    ans_expected <- matrix(nrow = 5, ncol = 2) ## n_particle x n_region
    ans_expected[, 1] <- as.numeric(1:5)
    set.seed(seed)
    for (i_particle in 1:5) {
      ans_expected[i_particle, 2] <- rnbinom(
        n = 1L,
        size = 1 / disp[2, i_interval + 1],
        mu = counts_data[2, i_interval + 1] / ratio[2, i_interval + 1]
      )
    }
    expect_identical(ans_obtained, ans_expected)
  }
})

test_that("'fill_counts_true' method for 'CdmWithregNbinom' works - disp is 0", {
  ## 10 intervals
  ## 5 particles
  ## 2 regions
  set.seed(0)
  counts_data <- matrix(as.numeric(1:20), nr = 2) ## n_region x n_interval
  ratio <- matrix(runif(n = 20, min = 0.5, max = 1.5), nr = 2)
  disp <- matrix(0, nr = 2, nc = 10)
  i_interval <- 3L
  x <- new_CdmWithregNbinom(
    counts_data = counts_data,
    ratio = ratio,
    disp = disp
  )
  for (seed in 1:5) {
    counts_true <- matrix(NA_real_, nrow = 5, ncol = 2) ## n_particle x n_region
    set.seed(seed)
    ## modify 'counts_true' in place
    x$fill_counts_true(
      counts_true = counts_true,
      i_interval = i_interval
    )
    ans_obtained <- counts_true
    ans_expected <- matrix(nrow = 5, ncol = 2) ## n_particle x n_region
    set.seed(seed)
    for (i_region in 1:2) {
      for (i_particle in 1:5) {
        lambda <- counts_data[i_region, i_interval + 1] / ratio[i_region, i_interval + 1]
        ans_expected[i_particle, i_region] <- as.double(rpois(
          n = 1L,
          lambda = lambda
        ))
      }
    }
    expect_identical(ans_obtained, ans_expected)
  }
})


## fill_logimp - Poisson-binomial, no regions ---------------------------------

test_that("'fill_logimp' method for 'CdmNoregPoibin' works - counts_data not NA", {
    ## 10 intervals
    ## 5 particles
    counts_data <- as.numeric(1:10) ## n_interval
    counts_true <- as.numeric(0:4) ## n_particle
    prob <- 0.98
    i_interval <- 3L
    x <- new_CdmNoregPoibin(
        counts_data = counts_data,
        prob = prob
    )
    for (seed in 1:5) {
        logimp <- rep(NA_real_, times = 5) ## n_particle
        set.seed(seed)
        ## modify 'logimp' in place
        x$fill_logimp(
              logimp = logimp,
              counts_true = counts_true,
              i_interval = i_interval
          )
        ans_obtained <- logimp
        ans_expected <- numeric(5) ## n_particle
        set.seed(seed)
        for (i in 1:5) {
            ans_expected[[i]] <- dpoibin1(
                x = counts_true[[i]],
                size = counts_data[i_interval + 1],
                prob = prob,
                use_log = TRUE
            )
        }
        expect_identical(ans_obtained, ans_expected)
    }
})

test_that("'fill_logimp' method for 'CdmNoregPoibin' works - counts_data is NA", {
    ## 10 intervals
    ## 5 particles
    counts_data <- as.numeric(c(1:3, NA, 5:10)) ## n_interval
    logimp <- rep(NA_real_, times = 5) ## n_particle
    counts_true <- as.numeric(6:10) ## n_particle
    prob <- 0.98
    i_interval <- 3L
    x <- new_CdmNoregPoibin(
        counts_data = counts_data,
        prob = 0.98
    )
    ## modify 'logimp' in place
    x$fill_logimp(
          logimp = logimp,
          counts_true = counts_true,
          i_interval = i_interval
      )
    ans_obtained <- logimp
    ans_expected <- rep(NA_real_, times = 5)
    expect_identical(ans_obtained, ans_expected)
})


## fill_logimp - negative binomial, no regions --------------------------------

test_that("'fill_logimp' method for 'CdmNoregNbinom' works - counts_data not NA", {
    ## 10 intervals
    ## 5 particles
    set.seed(0)
    counts_data <- as.numeric(1:10) ## n_interval
    ratio <- runif(n = 10, min = 0.5, max = 1.5)
    disp <- runif(n = 10, min = 0.5, max = 1.5)
    counts_true <- as.numeric(0:4) ## n_particle
    i_interval <- 3L
    x <- new_CdmNoregNbinom(
        counts_data = counts_data,
        ratio = ratio,
        disp = disp
    )
    for (seed in 1:5) {
        logimp <- rep(NA_real_, times = 5) ## n_particle
        set.seed(seed)
        ## modify 'logimp' in place
        x$fill_logimp(
              logimp = logimp,
              counts_true = counts_true,
              i_interval = i_interval
          )
        ans_obtained <- logimp
        ans_expected <- numeric(5) ## n_particle
        set.seed(seed)
        for (i in 1:5) {
            ans_expected[[i]] <- dnbinom(
                x = counts_true[[i]],
                size = 1 / disp[i_interval + 1],
                mu = counts_data[i_interval + 1] / ratio[i_interval + 1],
                log = TRUE
            )
        }
        expect_identical(ans_obtained, ans_expected)
    }
})

test_that("'fill_logimp' method for 'CdmNoregNbinom' works - counts_data has NA", {
  ## 10 intervals
  ## 5 particles
  set.seed(0)
  counts_data <- as.numeric(1:10) ## n_interval
  counts_data[4] <- NA
  ratio <- runif(n = 10, min = 0.5, max = 1.5)
  disp <- runif(n = 10, min = 0.5, max = 1.5)
  counts_true <- as.numeric(0:4) ## n_particle
  i_interval <- 3L
  x <- new_CdmNoregNbinom(
    counts_data = counts_data,
    ratio = ratio,
    disp = disp
  )
  logimp <- rep(NA_real_, times = 5) ## n_particle
  ## modify 'logimp' in place
  x$fill_logimp(
    logimp = logimp,
    counts_true = counts_true,
    i_interval = i_interval
  )
  ans_obtained <- logimp
  ans_expected <- rep(NA_real_, 5)
  expect_identical(ans_obtained, ans_expected)
})

test_that("'fill_logimp' method for 'CdmNoregNbinom' works - disp is 0", {
    ## 10 intervals
    ## 5 particles
    set.seed(0)
    counts_data <- as.numeric(1:10) ## n_interval
    ratio <- runif(n = 10, min = 0.5, max = 1.5)
    disp <- rep(0, 10)
    counts_true <- as.numeric(0:4) ## n_particle
    i_interval <- 3L
    x <- new_CdmNoregNbinom(
        counts_data = counts_data,
        ratio = ratio,
        disp = disp
    )
    for (seed in 1:5) {
        logimp <- rep(NA_real_, times = 5) ## n_particle
        set.seed(seed)
        ## modify 'logimp' in place
        x$fill_logimp(
              logimp = logimp,
              counts_true = counts_true,
              i_interval = i_interval
          )
        ans_obtained <- logimp
        ans_expected <- numeric(5) ## n_particle
        set.seed(seed)
        for (i in 1:5) {
            ans_expected[[i]] <- dpois(
                x = counts_true[[i]],
                lambda = counts_data[i_interval + 1] / ratio[i_interval + 1],
                log = TRUE
            )
        }
        expect_identical(ans_obtained, ans_expected)
    }
})


## fill_logimp - normal, no regions --------------------------------

test_that("'fill_logimp' method for 'CdmNoregNorm' works - counts_data not NA", {
    ## 10 intervals
    ## 5 particles
    set.seed(0)
    counts_data <- as.numeric(1:10) ## n_interval
    ratio <- runif(n = 10, min = 0.5, max = 1.5)
    sd <- runif(n = 10, min = 0.5, max = 1.5)
    counts_true <- as.numeric(0:4) ## n_particle
    i_interval <- 3L
    x <- new_CdmNoregNorm(
        counts_data = counts_data,
        ratio = ratio,
        sd = sd
    )
    for (seed in 1:5) {
        logimp <- rep(NA_real_, times = 5) ## n_particle
        set.seed(seed)
        ## modify 'logimp' in place
        x$fill_logimp(
              logimp = logimp,
              counts_true = counts_true,
              i_interval = i_interval
          )
        ans_obtained <- logimp
        ans_expected <- numeric(5) ## n_particle
        set.seed(seed)
        for (i in 1:5) {
            ans_expected[[i]] <- log(pnorm(
                q = counts_true[[i]] + 0.5,
                mean = counts_data[i_interval + 1] / ratio[i_interval + 1],
                sd = sd[i_interval + 1]) - pnorm(
                                               q = counts_true[[i]] - 0.5,
                                               mean = counts_data[i_interval + 1] / ratio[i_interval + 1],
                                               sd = sd[i_interval + 1]))
        }
        expect_equal(ans_obtained, ans_expected)
    }
})

test_that("'fill_logimp' method for 'CdmNoregNorm' works - counts_data has NA", {
  ## 10 intervals
  ## 5 particles
  set.seed(0)
  counts_data <- as.numeric(1:10) ## n_interval
  counts_data[4] <- NA
  ratio <- runif(n = 10, min = 0.5, max = 1.5)
  sd <- runif(n = 10, min = 0.5, max = 1.5)
  counts_true <- as.numeric(0:4) ## n_particle
  i_interval <- 3L
  x <- new_CdmNoregNorm(
    counts_data = counts_data,
    ratio = ratio,
    sd = sd
  )
  logimp <- rep(NA_real_, times = 5) ## n_particle
  ## modify 'logimp' in place
  x$fill_logimp(
    logimp = logimp,
    counts_true = counts_true,
    i_interval = i_interval
  )
  ans_obtained <- logimp
  ans_expected <- rep(NA_real_, 5)
  expect_identical(ans_obtained, ans_expected)
})



## fill_logimp - Poisson-binomial, with regions -------------------------------

test_that("'fill_logimp' method for CdmWithregPoibin' works - counts_data complete", {
    ## 10 intervals
    ## 5 particles
    ## 2 regions
    counts_data <- matrix(as.numeric(1:20), nr = 2) ## n_region x n_interval
    prob <- 0.98
    i_interval <- 3L
    counts_true <- matrix(as.numeric(0:9), nrow = 5, ncol = 2) ## n_particle x n_region
    x <- new_CdmWithregPoibin(
        counts_data = counts_data,
        prob = prob
    )
    for (seed in 1:5) {
        logimp <- matrix(NA_real_, nrow = 5, ncol = 2) ## n_particle x n_region
        set.seed(seed)
        ## modify 'logimp' in place
        x$fill_logimp(
              logimp = logimp,
              counts_true = counts_true,
              i_interval = i_interval
          )
        ans_obtained <- logimp
        ans_expected <- matrix(NA_real_,
                               nrow = 5,
                               ncol = 2
                               ) ## n_particle x n_region
        set.seed(seed)
        for (i_region in 1:2) {
            for (i_particle in 1:5) {
                ans_expected[i_particle, i_region] <- dpoibin1(
                    x = counts_true[i_particle, i_region],
                    size = counts_data[i_region, i_interval + 1],
                    prob = prob,
                    use_log = TRUE
                )
            }
        }
        expect_identical(ans_obtained, ans_expected)
    }
})

test_that("'fill_logimp' method for CdmWithregPoibin' works - counts_true partly full", {
    ## 10 intervals
    ## 5 particles
    ## 2 regions
    counts_data <- matrix(as.numeric(1:20), nr = 2) ## n_region x n_interval
    counts_data[1, 4] <- NA
    prob <- 0.98
    i_interval <- 3L
    counts_true <- matrix(as.numeric(9:0), nrow = 5, ncol = 2) ## n_particle x n_region
    x <- new_CdmWithregPoibin(
        counts_data = counts_data,
        prob = prob
    )
    for (seed in 1:5) {
        logimp <- matrix(NA_real_, nrow = 5, ncol = 2) ## n_particle x n_region
        logimp[, 1] <- 0.2
        set.seed(seed)
        ## modify 'logimp' in place
        x$fill_logimp(
              logimp = logimp,
              counts_true = counts_true,
              i_interval = i_interval
          )
        ans_obtained <- logimp
        ans_expected <- matrix(nrow = 5, ncol = 2) ## n_particle x n_region
        ans_expected[, 1] <- 0.2
        set.seed(seed)
        for (i_particle in 1:5) {
            ans_expected[i_particle, 2] <- dpoibin1(
                x = counts_true[i_particle, 2],
                size = counts_data[2, i_interval + 1],
                prob = prob,
                use_log = TRUE
            )
        }
        expect_identical(ans_obtained, ans_expected)
    }
})


## fill_logimp - negative binomial, with regions ------------------------------

test_that("'fill_logimp' method for CdmWithregPoibin' works - counts_data is NA", {
    ## 10 intervals
    ## 5 particles
    ## 2 regions
    set.seed(0)
    counts_data <- matrix(as.numeric(1:20), nr = 2) ## n_region x n_interval
    ratio <- matrix(runif(n = 20, min = 0.5, max = 1.5), nr = 2)
    disp <- matrix(runif(n = 20, min = 0.5, max = 1.5), nr = 2)
    i_interval <- 3L
    counts_true <- matrix(as.numeric(0:9), nrow = 5, ncol = 2) ## n_particle x n_region
    x <- new_CdmWithregNbinom(
        counts_data = counts_data,
        ratio = ratio,
        disp = disp
    )
    for (seed in 1:5) {
        logimp <- matrix(NA_real_, nrow = 5, ncol = 2) ## n_particle x n_region
        set.seed(seed)
        ## modify 'logimp' in place
        x$fill_logimp(
              logimp = logimp,
              counts_true = counts_true,
              i_interval = i_interval
          )
        ans_obtained <- logimp
        ans_expected <- matrix(NA_real_,
                               nrow = 5,
                               ncol = 2
                               ) ## n_particle x n_region
        set.seed(seed)
        for (i_region in 1:2) {
            for (i_particle in 1:5) {
                size <- 1 / disp[i_region, i_interval + 1]
                mu <- counts_data[i_region, i_interval + 1] / ratio[i_region, i_interval + 1]
                ans_expected[i_particle, i_region] <- dnbinom(
                    x = counts_true[i_particle, i_region],
                    size = size,
                    mu = mu,
                    log = TRUE
                )
            }
        }
        expect_identical(ans_obtained, ans_expected)
    }
})

test_that("'fill_logimp' method for CdmWithregNbinom' works - counts_true partly full", {
    ## 10 intervals
    ## 5 particles
    ## 2 regions
    set.seed(0)
    counts_data <- matrix(as.numeric(1:20), nr = 2) ## n_region x n_interval
    counts_data[1, 4] <- NA
    ratio <- matrix(runif(n = 20, min = 0.5, max = 1.5), nr = 2)
    disp <- matrix(runif(n = 20, min = 0.5, max = 1.5), nr = 2)
    i_interval <- 3L
    counts_true <- matrix(as.numeric(9:0), nrow = 5, ncol = 2) ## n_particle x n_region
    x <- new_CdmWithregNbinom(
        counts_data = counts_data,
        ratio = ratio,
        disp = disp
    )
    for (seed in 1:5) {
        logimp <- matrix(NA_real_, nrow = 5, ncol = 2) ## n_particle x n_region
        logimp[, 1] <- 0.2
        set.seed(seed)
        ## modify 'logimp' in place
        x$fill_logimp(
              logimp = logimp,
              counts_true = counts_true,
              i_interval = i_interval
          )
        ans_obtained <- logimp
        ans_expected <- matrix(nrow = 5, ncol = 2) ## n_particle x n_region
        ans_expected[, 1] <- 0.2
        set.seed(seed)
        for (i_particle in 1:5) {
            ans_expected[i_particle, 2] <- dnbinom(
                x = counts_true[i_particle, 2],
                size = 1 / disp[2, i_interval + 1],
                mu = counts_data[2, i_interval + 1] / ratio[2, i_interval + 1],
                log = TRUE
            )
        }
        expect_identical(ans_obtained, ans_expected)
    }
})

test_that("'fill_logimp' method for CdmWithregPoibin' works - disp is 0", {
    ## 10 intervals
    ## 5 particles
    ## 2 regions
    set.seed(0)
    counts_data <- matrix(as.numeric(1:20), nr = 2) ## n_region x n_interval
    ratio <- matrix(runif(n = 20, min = 0.5, max = 1.5), nr = 2)
    disp <- matrix(0, nr = 2, nc = 10)
    i_interval <- 3L
    counts_true <- matrix(as.numeric(0:9), nrow = 5, ncol = 2) ## n_particle x n_region
    x <- new_CdmWithregNbinom(
        counts_data = counts_data,
        ratio = ratio,
        disp = disp
    )
    for (seed in 1:5) {
        logimp <- matrix(NA_real_, nrow = 5, ncol = 2) ## n_particle x n_region
        set.seed(seed)
        ## modify 'logimp' in place
        x$fill_logimp(
              logimp = logimp,
              counts_true = counts_true,
              i_interval = i_interval
          )
        ans_obtained <- logimp
        ans_expected <- matrix(NA_real_,
                               nrow = 5,
                               ncol = 2
                               ) ## n_particle x n_region
        set.seed(seed)
        for (i_region in 1:2) {
            for (i_particle in 1:5) {
                lambda <- counts_data[i_region, i_interval + 1] / ratio[i_region, i_interval + 1]
                ans_expected[i_particle, i_region] <- dpois(
                    x = counts_true[i_particle, i_region],
                    lambda = lambda,
                    log = TRUE
                )
            }
        }
        expect_identical(ans_obtained, ans_expected)
    }
})

