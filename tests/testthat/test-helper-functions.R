

## draw_gross_mig_inner -------------------------------------------------------

test_that("'draw_gross_mig_inner' give correct answer with valid inputs", {
    rates_im <- runif(n = 5, max = 3)
    rates_em <- runif(n = 5, max = 0.5)
    exposure <- runif(n = 5, min = 3, max = 10)
    net_mig <- round(rnorm(n = 5, sd = 5))
    set.seed(10)
    ans_obtained <- draw_gross_mig_inner(rates_im = rates_im,
                                         rates_em = rates_em,
                                         exposure = exposure,
                                         net_mig = net_mig)
    set.seed(10)
    ans_expected_gross_mig <- rpoistr(
        n = 5,
        lambda = rates_im + rates_em * exposure,
        lower = abs(net_mig)
    )
    different_parity <- (ans_expected_gross_mig %% 2L) != (abs(net_mig) %% 2L)
    ans_expected_gross_mig[different_parity] <- ans_expected_gross_mig[different_parity] + 1
    logprob <- rbind( dpoistr(ans_expected_gross_mig,
                              lambda = rates_im + rates_em * exposure,
                              lower = abs(net_mig),
                              use_log = TRUE
                              ),
                     ifelse(ans_expected_gross_mig > abs(net_mig),
                            dpoistr(ans_expected_gross_mig - 1,
                                    lambda = rates_im + rates_em * exposure,
                                    lower = abs(net_mig),
                                    use_log = TRUE
                                    ),
                            -Inf
                            ))
    ans_expected_logimp <- sapply(seq_len(5),
                                  function(i) log_sum_exp_2(logprob[1,i], logprob[2,i]))
    ans_expected <- list(ans_expected_gross_mig,
                         ans_expected_logimp)
    expect_identical(ans_obtained, ans_expected)
})    


## draw_index_parent ----------------------------------------------------------

test_that("'draw_index_parent' selects elements in proportion to weights", {
  set.seed(0)
  wt <- (1:4) / sum(1:4)
  fast <- TRUE
  if (fast) {
    n <- 1000
    tolerance <- 0.01
  } else {
    n <- 1000000
    tolerance <- 0.001
  }
  x <- replicate(n = n, draw_index_parent(wt))
  obs_propn <- as.numeric(prop.table(table(as.integer(x))))
  expect_equal(obs_propn, wt, tolerance = tolerance)
})

## draw_rates_inner ----------------------------------------------------------

test_that("'draw_rates_inner' gives results with expeced mean and variance", {
    for (seed in 1:2) {
        set.seed(seed)
        n <- 100000L
        mean <- runif(1, max = 10)
        disp <- runif(1, max = 2)
        ans <- draw_rates_inner(n = n, mean = mean, disp = disp)
        expect_equal(mean(ans), mean, tolerance = 0.01)
        expect_equal(var(ans), mean * disp, tolerance = 0.01)
        expect_identical(length(ans), n)
    }
    ans <- draw_rates_inner(n = n, mean = 30, disp = 0)
    expect_true(all(ans == 30))
    expect_identical(length(ans), n)
    ans <- draw_rates_inner(n = n, mean = 0, disp = 1)
    expect_true(all(ans == 0))
    expect_identical(length(ans), n)
    ans <- draw_rates_inner(n = n, mean = 0, disp = 0)
    expect_true(all(ans == 0))
    expect_identical(length(ans), n)
})


## has_region -----------------------------------------------------------------

test_that("'has_region' works with numeric vector", {
  df <- data.frame(x = 1, n_region = 2)
  expect_true(has_region(df))
  df <- data.frame(x = 1, n_interval = 2)
  expect_false(has_region(df))
})


## log_sum_exp_2 -----------------------------------------------------------------

test_that("'log_sum_exp_2' gives same answer as calculating directly (with moderate sized numbers)", {
    for (seed in 1:10) {
        x <- rnorm(10)
        y <- rnorm(10)
        expect_equal(log_sum_exp_2(x = x, y = y),
                     log(exp(x) + exp(y)))
    }
})

test_that("'log_sum_exp_2' returns -Inf when both arguments are -Inf", {
        expect_equal(log_sum_exp_2(x = -Inf, y = -Inf),
                     -Inf)
})


## 'logprob_trans' ------------------------------------------------------------

test_that("'logprob_trans' gives expected answer", {
    for (seed in seq_len(5)) {
        set.seed(seed)
        stock_start <- as.numeric(rpois(n = 1, lambda = 100))
        counts_bth <- as.numeric(rpois(n = 1, lambda = 10))
        counts_dth <- as.numeric(rpois(n = 1, lambda = 10))
        counts_im <- as.numeric(rpois(n = 1, lambda = 10))
        counts_em <- as.numeric(rpois(n = 1, lambda = 10))
        rates_bth <- runif(n = 1, max = 2)
        rates_dth <- runif(n = 1, max = 2)
        rates_im <- runif(n = 1, max = 2)
        rates_em <- runif(n = 1, max = 2)
        exposure <- runif(n = 1, min = 0.5, max = 20)
        ans_obtained <- logprob_trans(stock_start = stock_start,
                                      counts_bth = counts_bth,
                                      counts_dth = counts_dth,
                                      counts_im = counts_im,
                                      counts_em = counts_em,
                                      rates_bth = rates_bth,
                                      rates_dth = rates_dth,
                                      rates_im = rates_im,
                                      rates_em = rates_em,
                                      exposure = exposure,
                                      is_dominant = TRUE)
        logtrans_im <- dpois(counts_im, lambda = rates_im, log = TRUE)
        prob_exit <- 1 - exp(-0.5 * (rates_dth + rates_em))
        prob_dth <- rates_dth / (rates_dth + rates_em)
        count_im <- counts_im
        if (count_im %% 2 == 0) {
            size <- stock_start + count_im / 2
            logtrans_dth_em <- dmultinom(x = c(counts_dth,
                                               counts_em,
                                               size - counts_dth - counts_em),
                                         size = size,
                                         prob = c(prob_dth * prob_exit, (1-prob_dth) * prob_exit, (1-prob_exit)),
                                         log = TRUE)
        } else {
            size1 <- stock_start + count_im / 2 + 0.5
            logtrans_dth_em <- log(0.5 * dmultinom(x = c(counts_dth,
                                                         counts_em,
                                                         size1 - counts_dth - counts_em),
                                                   size = size1,
                                                   prob = c(prob_dth * prob_exit, (1-prob_dth) * prob_exit, (1-prob_exit)))
                                   + 0.5 * dmultinom(x = c(counts_dth,
                                                           counts_em,
                                                           size1 - 1 - counts_dth - counts_em),
                                                     size = size1 - 1,
                                                     prob = c(prob_dth * prob_exit, (1-prob_dth) * prob_exit, (1-prob_exit))))
        }
        logtrans_bth <- dpois(counts_bth, lambda = rates_bth * exposure, log = TRUE)
        ans_expected <- logtrans_im + logtrans_dth_em + logtrans_bth
        expect_equal(ans_obtained, ans_expected)
    }
})

    

## make_pf --------------------------------------------------------------------

test_that("'make_pf' works with for models without region", {
  ## constructing a valid 'df' is a lot of work, so test
  ## via function 'simulated_pf_noreg', which calls 'make_pf'
  ans <- simulated_pf_noreg()
  expect_s3_class(ans, "PFilterNoReg")
})

test_that("'make_pf' works with for models with region", {
  ## constructing a valid 'df' is a lot of work, so test
  ## via function 'simulated_pf_withreg', which calls 'make_pf'
  ans <- simulated_pf_withreg()
  expect_s3_class(ans, "PFilterWithReg")
})


## show_start_message ---------------------------------------------------------

test_that("show_start_message works", {
    df <- list(cohort = 2000L, sexgender = "Female")
    expect_message(show_start_message(df),
                   "starting cohort :")
})


## draw_index_parent ----------------------------------------------------------

test_that("'draw_index_parent' selects elements in proportion to weights", {
  set.seed(0)
  wt <- (1:4) / sum(1:4)
  fast <- TRUE
  if (fast) {
    n <- 1000
    tolerance <- 0.01
  } else {
    n <- 1000000
    tolerance <- 0.001
  }
  x <- replicate(n = n, draw_index_parent(wt))
  obs_propn <- as.numeric(prop.table(table(as.integer(x))))
  expect_equal(obs_propn, wt, tolerance = tolerance)
})


## has_region -----------------------------------------------------------------

test_that("'has_region' works with numeric vector", {
  df <- data.frame(x = 1, n_region = 2)
  expect_true(has_region(df))
  df <- data.frame(x = 1, n_interval = 2)
  expect_false(has_region(df))
})


## make_pf --------------------------------------------------------------------

test_that("'make_pf' works with for models without region", {
  ## constructing a valid 'df' is a lot of work, so test
  ## via function 'simulated_pf_noreg', which calls 'make_pf'
  ans <- simulated_pf_noreg()
  expect_s3_class(ans, "PFilterNoReg")
})

test_that("'make_pf' works with for models with region", {
  ## constructing a valid 'df' is a lot of work, so test
  ## via function 'simulated_pf_withreg', which calls 'make_pf'
  ans <- simulated_pf_withreg()
  expect_s3_class(ans, "PFilterWithReg")
})


## softmax --------------------------------------------------------------------

test_that("'softmax' works with valid 'x'", {
  x <- (-5):5
  expect_equal(
    softmax(x),
    exp(x) / sum(exp(x))
  )
  expect_equal(
    softmax(0.5),
    1
  )
})


## sum_by ---------------------------------------------------------------------

test_that("'sum_by' works with valid 'x' and 'g'", {
  set.seed(0)
  x <- rnorm(n = 100)
  g <- sample(10, size = 100, replace = TRUE)
  ans_obtained <- sum_by(x = x, g = g)
  ans_expected <- as.numeric(tapply(x, factor(g, levels = unique(g)), sum))
  expect_equal(ans_obtained, ans_expected)
})


## softmax --------------------------------------------------------------------

test_that("'softmax' works with valid 'x'", {
  x <- (-5):5
  expect_equal(
    softmax(x),
    exp(x) / sum(exp(x))
  )
  expect_equal(
    softmax(0.5),
    1
  )
})


## sum_by ---------------------------------------------------------------------

test_that("'sum_by' works with valid 'x' and 'g'", {
  set.seed(0)
  x <- rnorm(n = 100)
  g <- sample(10, size = 100, replace = TRUE)
  ans_obtained <- sum_by(x = x, g = g)
  ans_expected <- as.numeric(tapply(x, factor(g, levels = unique(g)), sum))
  expect_equal(ans_obtained, ans_expected)
})


## write_diagnostics ----------------------------------------------------------

test_that("'write_diagnostics' works", {
    diagnostics <- data.frame(x = 1:2)
    work_dir <- tempdir()
    write_diagnostics(
        diagnostics = diagnostics,
        cohort = 2000,
        sexgender = "Female",
        work_dir = work_dir
    )
    fname <- file.path(work_dir, "tmp-diagnostics-2000-Female.rds")
    ans_obtained <- readRDS(fname)
    ans_expected <- diagnostics
    expect_identical(ans_obtained, ans_expected)
})


## write_outputs --------------------------------------------------------------

test_that("'write_outputs' works", {
    meta <- serialize(data.frame(x = 1:2), connection = NULL)
    n_meta <- length(meta)
    outputs <- list(
        a = list(n_meta, 2L, meta, c(1, 2)),
        b = list(n_meta, 2L, meta, c(5, 6))
        )
    work_dir <- tempdir()
    write_outputs(
        outputs = outputs,
        work_dir = work_dir
    )
    fname1 <- file.path(work_dir, "tmp-a.bin")
    con <- file(fname1, "rb")
    expect_identical(list(readBin(con, what = "integer", n = 1L),
                          readBin(con, what = "integer", n = 1L),
                          readBin(con, what = "raw", n = n_meta),
                          readBin(con, what = "double", n = 2L)),
                     outputs[[1L]])
    close(con)
    fname2 <- file.path(work_dir, "tmp-b.bin")
    con <- file(fname2, "rb")
    expect_identical(list(readBin(con, what = "integer", n = 1L),
                          readBin(con, what = "integer", n = 1L),
                          readBin(con, what = "raw", n = n_meta),
                          readBin(con, what = "double", n = 2L)),
                     outputs[[2L]])
    close(con)
})

