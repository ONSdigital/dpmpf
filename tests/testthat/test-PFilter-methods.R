
## counts_bth ------------------------------------------------------------

test_that("'calc_counts_bth' gives expected answer - no regions", {
  set.seed(0)
  i_interval <- 2L
  pf <- simulated_pf_noreg()
  pf$calc_counts_bth(i_interval = i_interval)
  expect_identical(pf$counts_bth, pf$counts_births[, i_interval])
})

test_that("'calc_counts_bth' gives expected answer - parallelogram", {
  set.seed(0)
  i_interval <- 1L
  pf <- simulated_pf_noreg(parallelogram = TRUE)
  pf$calc_counts_bth(i_interval = i_interval)
  expect_identical(pf$counts_bth_second, pf$counts_births[, i_interval + 1])
})

test_that("'calc_counts_bth' gives expected answer - with regions", {
  set.seed(0)
  i_interval <- 2L
  pf <- simulated_pf_withreg()
  pf$calc_counts_bth(i_interval = i_interval)
  expect_identical(pf$counts_bth, pf$counts_births[, , i_interval])
})


## calc_counts_dth ------------------------------------------------------------

test_that("'calc_counts_dth' gives expected answer - no regions", {
  set.seed(0)
  i_interval <- 2L
  pf <- simulated_pf_noreg()
  pf$calc_counts_dth(i_interval = i_interval)
  expect_identical(pf$counts_dth, pf$counts_deaths[, i_interval])
})

test_that("'calc_counts_dth' gives expected answer - parallelogram", {
  set.seed(0)
  i_interval <- 1L
  pf <- simulated_pf_noreg(parallelogram = TRUE)
  pf$calc_counts_dth(i_interval = i_interval)
  expect_identical(pf$counts_dth_second, pf$counts_deaths[, i_interval + 1])
})

test_that("'calc_counts_dth' gives expected answer - with regions", {
  set.seed(0)
  i_interval <- 2L
  pf <- simulated_pf_withreg()
  pf$calc_counts_dth(i_interval = i_interval)
  expect_identical(pf$counts_dth, pf$counts_deaths[, , i_interval])
})


## calc_exposure --------------------------------------------------------------

test_that("'calc_exposure' gives expected answer - no regions", {
    for (seed in 1:5) {
        set.seed(seed)
        pf <- simulated_pf_noreg()
        pf$stock_start <- rpois(n = pf$n_particle, lambda = 1)
        pf$stock_end <- rpois(n = pf$n_particle, lambda = 1)
        pf$counts_bth <- rep(0, pf$n_particle)
        pf$counts_im1 <- rpois(n = pf$n_particle, lambda = 0.5)
        pf$counts_im2 <- rpois(n = pf$n_particle, lambda = 0.5)
        pf$calc_exposure(1L)
        ans_obtained <- pf$exposure
        ans_expected <- 0.25 * pf$stock_start + 0.25 * pf$stock_end
        count <- pf$counts_im1 + pf$counts_im2
        ans_expected[ans_expected == 0 & count > 0] <- 0.25
        expect_identical(ans_obtained, ans_expected)
    }
})

test_that("'calc_exposure' gives expected answer - parallelogram", {
    for (seed in 1:5) {
        set.seed(seed)
        pf <- simulated_pf_noreg(parallelogram = TRUE)
        pf$stock_end <- rpois(n = pf$n_particle, lambda = 1)
        pf$stock_end_second <- rpois(n = pf$n_particle, lambda = 1)
        pf$counts_bth_second <- rpois(n = pf$n_particle, lambda = 0.5)
        pf$counts_im1_second <- rpois(n = pf$n_particle, lambda = 0.5)
        pf$calc_exposure(1L)
        ans_obtained <- pf$exposure_second
        ans_expected <- 0.25 * pf$stock_end + 0.25 * pf$stock_end_second
        count <- pf$counts_im1_second + pf$counts_bth_second
        ans_expected[ans_expected == 0 & count > 0] <- 0.25
        expect_identical(ans_obtained, ans_expected)
    }
})


test_that("'calc_exposure' gives expected answer - with regions", {
    for (seed in seq_len(5)) {
        set.seed(seed)
        pf <- simulated_pf_withreg()
        pf$stock_start <- matrix(rpois(n = pf$n_particle * pf$n_region, lambda = 10),
                                 nrow = pf$n_particle
                                 )
        pf$stock_end <- matrix(rpois(n = pf$n_particle * pf$n_region, lambda = 10),
                               nrow = pf$n_particle
                               )
        pf$counts_bth <- matrix(0, nrow = pf$n_particle, ncol = pf$n_region)
        pf$counts_dth <- matrix(0, nrow = pf$n_particle, ncol = pf$n_region)
        pf$counts_in <- matrix(rpois(n = pf$n_particle * pf$n_region, lambda = 0.5),
                               nrow = pf$n_particle)
        pf$counts_out <- matrix(rpois(n = pf$n_particle * pf$n_region, lambda = 0.5),
                                nrow = pf$n_particle)
        pf$counts_im1 <- matrix(rpois(n = pf$n_particle * pf$n_region, lambda = 0.5),
                                nrow = pf$n_particle)
        pf$counts_em1 <- matrix(rpois(n = pf$n_particle * pf$n_region, lambda = 0.5),
                                nrow = pf$n_particle)
        pf$counts_im2 <- matrix(rpois(n = pf$n_particle * pf$n_region, lambda = 0.5),
                                nrow = pf$n_particle)
        pf$counts_em2 <- matrix(rpois(n = pf$n_particle * pf$n_region, lambda = 0.5),
                                nrow = pf$n_particle)
        pf$calc_exposure()
        ans_obtained <- pf$exposure
        ans_expected <- 0.25 * pf$stock_start + 0.25 * pf$stock_end
        gross_mig <- pf$counts_in + pf$counts_out + pf$counts_im1 +
            pf$counts_em1 + pf$counts_im2 + pf$counts_em2
        ans_expected[ans_expected == 0 & gross_mig > 0] <- 0.25
        expect_identical(ans_obtained, ans_expected)
    }
})


## calc_exposure_approx1 ------------------------------------------------------

test_that("'calc_exposure_approx1' gives expected answer - no regions", {
  set.seed(0)
  i_interval <- 1L
  pf <- simulated_pf_noreg(has_stock_init = TRUE)
  pf$index_parent[, i_interval] <- sample(pf$n_particle)
  pf$calc_stock_start(i_interval = i_interval)
  pf$calc_rates(i_interval = i_interval)
  pf$calc_exposure_approx1()
  ans_obtained <- pf$exposure_approx1
  ans_expected <- 0.5 * pf$stock_start
  ans_cim <- 0.25 * pf$rates_cim
  ans_expected[ans_cim > ans_expected] <- ans_cim[ans_cim > ans_expected]
  expect_identical(ans_obtained, ans_expected)
})

test_that("'calc_exposure_approx1' gives expected answer - with regions", {
  set.seed(0)
  i_interval <- 1L
  pf <- simulated_pf_withreg(has_stock_init = TRUE)
  pf$index_parent[, i_interval] <- sample(pf$n_particle)
  pf$calc_stock_start(i_interval = i_interval)
  pf$calc_rates(i_interval = i_interval)
  pf$calc_exposure_approx1()
  ans_obtained <- pf$exposure_approx1
  ans_expected <- 0.5 * pf$stock_start
  ans_cin <- 0.25 * pf$rates_cin
  ans_expected[ans_cin > ans_expected] <- ans_cin[ans_cin > ans_expected]
  expect_identical(ans_obtained, ans_expected)
})


## calc_exposure_approx2 ------------------------------------------------------

test_that("'calc_exposure_approx2' gives expected answer - no regions", {
  set.seed(0)
  i_interval <- 1L
  pf <- simulated_pf_noreg(has_stock_init = TRUE)
  pf$index_parent[, i_interval] <- sample(pf$n_particle)
  pf$calc_stock_start(i_interval = i_interval)
  pf$draw_stock_end_net_mig(i_interval = i_interval, is_forecast = FALSE)
  pf$calc_rates(i_interval = i_interval)
  pf$calc_exposure_approx2(i_interval)
  ans_obtained <- pf$exposure_approx2
  ans_expected <- 0.25 * (pf$stock_start + pf$stock_end)
  ans_cim <- 0.25 * pf$rates_cim
  ans_expected[ans_cim > ans_expected] <- ans_cim[ans_cim > ans_expected]
  expect_identical(ans_obtained, ans_expected)
})

test_that("'calc_exposure_approx2' gives expected answer - parallelogram", {
  set.seed(0)
  i_interval <- 1L
  pf <- simulated_pf_noreg(parallelogram = TRUE)
  pf$stock_start <- rpois(pf$n_particle, lambda = 5)
  pf$stock_end <- rpois(pf$n_particle, lambda = 5)
  pf$stock_end_second <- rpois(pf$n_particle, lambda = 5)
  pf$rates_cim <- runif(pf$n_particle, max = 2)
  pf$rates_im1_second <- runif(pf$n_particle, max = 2)
  pf$calc_exposure_approx2(i_interval)
  ans_obtained_first <- pf$exposure_approx2
  ans_expected_first <- 0.25 * (pf$stock_start + pf$stock_end)
  ans_cim <- 0.25 * pf$rates_cim
  ans_expected_first[ans_cim > ans_expected_first] <- ans_cim[ans_cim > ans_expected_first]
  expect_identical(ans_obtained_first, ans_expected_first)
  ans_obtained_second <- pf$exposure_approx2_second
  ans_expected_second <- 0.25 * (pf$stock_end + pf$stock_end_second)
  ans_im <- 0.25 * pf$rates_im_second
  ans_expected_second[ans_im > ans_expected_second] <- ans_im[ans_im > ans_expected_second]
  expect_identical(ans_obtained_second, ans_expected_second)
})


test_that("'calc_exposure_approx2' gives expected answer - with regions", {
  set.seed(0)
  i_interval <- 1L
  pf <- simulated_pf_withreg(has_stock_init = TRUE)
  pf$index_parent[, i_interval] <- sample(pf$n_particle)
  pf$calc_stock_start(i_interval = i_interval)
  pf$calc_rates(i_interval = i_interval)
  pf$calc_counts_dth(i_interval = i_interval)
  pf$calc_exposure_approx1()
  pf$draw_stock_end_net_mig(i_interval = i_interval, is_forecast = FALSE)
  pf$calc_rates(i_interval = i_interval)
  pf$calc_exposure_approx2()
  ans_obtained <- pf$exposure_approx2
  ans_expected <- 0.25 * (pf$stock_start + pf$stock_end)
  ans_cin <- 0.25 * pf$rates_cin
  ans_expected[ans_cin > ans_expected] <- ans_cin[ans_cin > ans_expected]
  expect_identical(ans_obtained, ans_expected)
})


## calc_exposure_approx_parallelogram -----------------------------------------

test_that("'calc_exposure_parallelogram' gives expected answer", {
    for (seed in 1:5) {
        set.seed(seed)
        pf <- simulated_pf_noreg(parallelogram = TRUE)
        pf$stock_start <- rpois(n = pf$n_particle, lambda = 1)
        pf$stock_end_second <- rpois(n = pf$n_particle, lambda = 1)
        pf$counts_dth <- rpois(n = pf$n_particle, lambda = 0.5)
        pf$counts_dth_second <- rpois(n = pf$n_particle, lambda = 0.5)
        pf$rates_im1 <- rpois(n = pf$n_particle, lambda = 0.5)
        pf$rates_im1_second <- rpois(n = pf$n_particle, lambda = 0.5)
        pf$rates_em1 <- rpois(n = pf$n_particle, lambda = 0.5)
        pf$rates_em1_second <- rpois(n = pf$n_particle, lambda = 0.5)
        pf$calc_exposure_approx_parallelogram()
        ans_obtained_first <- pf$exposure_approx_first
        ans_obtained_second <- pf$exposure_approx_second
        stock_mid <- 0.5 * (pf$stock_start + pf$rates_im1 - pf$counts_dth - 0.5 * pf$rates_em1 * pf$stock_start) +
            0.5 * (pf$stock_end_second - pf$rates_im1_second + pf$counts_dth_second
                + 0.5 * pf$rates_em1_second * pf$stock_end_second)
        ans_expected_first <- pmax(0, 0.25 * pf$stock_start + 0.25 * stock_mid)
        ans_expected_second <- pmax(0, 0.25 * pf$stock_end_second + 0.25 * stock_mid)
        expect_identical(ans_obtained_first, ans_expected_first)
        expect_identical(ans_obtained_second, ans_expected_second)
    }
})


## calc_index_ancestor --------------------------------------------------------

test_that("'calc_index_ancestor' gives expected answer", {
  pf <- simulated_pf_noreg()
  pf$index_parent <- cbind(
    c(1L, 1L, 2L, 2L, 3L),
    c(1L, 2L, 3L, 3L, 4L),
    c(1L, 2L, 3L, 4L, 5L),
    c(1L, 2L, 2L, 3L, 3L),
    c(1L, 1L, 2L, 3L, 4L)
  )
  pf$calc_index_ancestor()
  ans_obtained <- pf$index_ancestor
  ans_expected <- cbind(
    c(1L, 1L, 1L, 1L, 2L),
    c(1L, 1L, 2L, 2L, 3L),
    c(1L, 1L, 2L, 2L, 3L),
    c(1L, 1L, 2L, 2L, 3L),
    c(1L, 1L, 2L, 3L, 4L)
  )
  expect_identical(ans_obtained, ans_expected)
})


## calc_logimp -----------------------------------------------------------

test_that("'calc_logimp' gives expected answer - no regions", {
  pf <- simulated_pf_noreg()
  pf$logimp_counts_dth <- rep(1 / 10, times = pf$n_particle)
  pf$logimp_counts_im_em <- rep(1 / 11, times = pf$n_particle)
  pf$logimp_gross_mig <- rep(1 / 12, times = pf$n_particle)
  pf$logimp_stock_end_init <- rep(1 / 13, times = pf$n_particle)
  pf$logimp_stock_end_net_mig <- rep(1 / 14, times = pf$n_particle)
  expect_equal(
    pf$calc_logimp(is_forecast = FALSE),
    pf$logimp_counts_im_em +
      pf$logimp_gross_mig +
      pf$logimp_stock_end_net_mig
  )
  expect_equal(
    pf$calc_logimp(is_forecast = TRUE),
    pf$logimp_counts_dth +
      pf$logimp_counts_im_em +
      pf$logimp_gross_mig +
      pf$logimp_stock_end_net_mig
  )
})

test_that("'calc_logimp' gives expected answer - with regions", {
  pf <- simulated_pf_withreg()
  pf$logimp_counts_dth <- rep(1 / 10, times = pf$n_particle)
  pf$logimp_counts_im_em <- rep(1 / 11, times = pf$n_particle)
  pf$logimp_gross_mig <- rep(1 / 12, times = pf$n_particle)
  pf$logimp_stock_end_init <- rep(1 / 13, times = pf$n_particle)
  pf$logimp_stock_end_net_mig <- rep(1 / 14, times = pf$n_particle)
  pf$logimp_counts_in_out <- rep(1 / 15, times = pf$n_particle)
  expect_equal(
    pf$calc_logimp(is_forecast = FALSE),
    pf$logimp_counts_im_em +
      pf$logimp_gross_mig +
      pf$logimp_stock_end_net_mig +
      pf$logimp_counts_in_out
  )
  expect_equal(
    pf$calc_logimp(is_forecast = TRUE),
    pf$logimp_counts_dth +
      pf$logimp_counts_im_em +
      pf$logimp_gross_mig +
      pf$logimp_stock_end_net_mig +
      pf$logimp_counts_in_out
  )
})


## calc_logimp_init -----------------------------------------------------------

test_that("'calc_logimp_init' gives expected answer - no regions", {
  pf <- simulated_pf_noreg()
  pf$logimp_stock_end_init <- rep(0.3, pf$n_particle)
  expect_identical(pf$calc_logimp_init(), rep(0.3, pf$n_particle))
})

test_that("'calc_logimp_init' gives expected answer - with regions", {
  pf <- simulated_pf_withreg()
  pf$logimp_stock_end_init <- rep(0.3, pf$n_particle)
  expect_identical(pf$calc_logimp_init(), rep(0.3, pf$n_particle))
})


## calc_loglik ----------------------------------------------------------------

test_that("'calc_loglik' gives expected answer - no regions", {
  set.seed(0)
  pf <- simulated_pf_noreg()
  n <- pf$n_particle
  pf$stock_end <- as.numeric(rpois(n = n, lambda = 20))
  pf$counts_im1 <- as.numeric(rpois(n = n, lambda = 10))
  pf$counts_em1 <- as.numeric(rpois(n = n, lambda = 10))
  pf$counts_im2 <- as.numeric(rpois(n = n, lambda = 10))
  pf$counts_em2 <- as.numeric(rpois(n = n, lambda = 10))
  cdm_im1 <- new_CdmNoregPoibin(as.numeric(rpois(n = n, lambda = 10)), prob = 0.9)
  cdm_em1 <- new_CdmNoregPoibin(as.numeric(rpois(n = n, lambda = 10)), prob = 0.9)
  cdm_im2 <- new_CdmNoregPoibin(as.numeric(rpois(n = n, lambda = 10)), prob = 0.9)
  cdm_em2 <- new_CdmNoregPoibin(as.numeric(rpois(n = n, lambda = 10)), prob = 0.9)
  pf$cdms_immigration1 <- new_CdmsNoreg(list(cdm_im1))
  pf$cdms_emigration1 <- new_CdmsNoreg(list(cdm_em1))
  pf$cdms_immigration2 <- new_CdmsNoreg(list(cdm_im2))
  pf$cdms_emigration2 <- new_CdmsNoreg(list(cdm_em2))
  ans_obtained <- pf$calc_loglik(2L)
  ans_expected <- pf$cdms_stock$calc_loglik(pf$stock_end, 3L, pf$obs_zero) +
    pf$cdms_immigration1$calc_loglik(pf$counts_im1, 2L, pf$obs_zero) +
    pf$cdms_emigration1$calc_loglik(pf$counts_em1, 2L, pf$obs_zero) +
    pf$cdms_immigration2$calc_loglik(pf$counts_im2, 2L, pf$obs_zero) +
    pf$cdms_emigration2$calc_loglik(pf$counts_em2, 2L, pf$obs_zero)
  expect_equal(ans_obtained, ans_expected)
})

test_that("'calc_loglik' gives expected answer - parallel", {
  set.seed(0)
  pf <- simulated_pf_noreg(parallel = TRUE, has_one_imem = TRUE)
  n <- pf$n_particle
  pf$stock_end <- as.numeric(rpois(n = n, lambda = 20))
  pf$stock_end_second <- as.numeric(rpois(n = n, lambda = 20))
  pf$counts_im1 <- as.numeric(rpois(n = n, lambda = 10))
  pf$counts_em1 <- as.numeric(rpois(n = n, lambda = 10))
  pf$counts_im2 <- as.numeric(rpois(n = n, lambda = 10))
  pf$counts_em2 <- as.numeric(rpois(n = n, lambda = 10))
  ans_obtained <- pf$calc_loglik(0L)
  ans_expected <- pf$cdms_stock$calc_loglik(pf$stock_end, 1L, pf$obs_zero) +
      pf$cdms_stock$calc_loglik(pf$stock_end_second, 2L, pf$obs_zero)
  expect_equal(ans_obtained, ans_expected)
})

test_that("'calc_loglik' gives expected answer - with regions", {
  set.seed(0)
  pf <- simulated_pf_withreg()
  n <- pf$n_particle * pf$n_region
  nr <- pf$n_particle
  pf$stock_end <- matrix(as.numeric(rpois(n = n, lambda = 20)), nr = nr)
  pf$counts_in <- matrix(as.numeric(rpois(n = n, lambda = 10)), nr = nr)
  pf$counts_out <- matrix(as.numeric(rpois(n = n, lambda = 10)), nr = nr)
  pf$counts_im1 <- matrix(as.numeric(rpois(n = n, lambda = 10)), nr = nr)
  pf$counts_em1 <- matrix(as.numeric(rpois(n = n, lambda = 10)), nr = nr)
  pf$counts_im2 <- matrix(as.numeric(rpois(n = n, lambda = 10)), nr = nr)
  pf$counts_em2 <- matrix(as.numeric(rpois(n = n, lambda = 10)), nr = nr)
  cdm_in <- new_CdmNoregPoibin(matrix(as.numeric(rpois(n = n, lambda = 10)), nr = nr),
    prob = 0.9
  )
  cdm_out <- new_CdmNoregPoibin(matrix(as.numeric(rpois(n = n, lambda = 10)), nr = nr),
    prob = 0.9
  )
  cdm_im1 <- new_CdmNoregPoibin(matrix(as.numeric(rpois(n = n, lambda = 10)), nr = nr),
    prob = 0.9
  )
  cdm_em1 <- new_CdmNoregPoibin(matrix(as.numeric(rpois(n = n, lambda = 10)), nr = nr),
    prob = 0.9
  )
  cdm_im2 <- new_CdmNoregPoibin(matrix(as.numeric(rpois(n = n, lambda = 10)), nr = nr),
    prob = 0.9
  )
  cdm_em2 <- new_CdmNoregPoibin(matrix(as.numeric(rpois(n = n, lambda = 10)), nr = nr),
    prob = 0.9
  )
  pf$cdms_internal_in <- new_CdmsNoreg(list(cdm_in))
  pf$cdms_internal_out <- new_CdmsNoreg(list(cdm_out))
  pf$cdms_immigration1 <- new_CdmsNoreg(list(cdm_im1))
  pf$cdms_emigration1 <- new_CdmsNoreg(list(cdm_em1))
  pf$cdms_immigration2 <- new_CdmsNoreg(list(cdm_im2))
  pf$cdms_emigration2 <- new_CdmsNoreg(list(cdm_em2))
  ans_obtained <- pf$calc_loglik(2L)
  ans_expected <- pf$cdms_stock$calc_loglik(pf$stock_end, 3L, pf$obs_zero) +
    pf$cdms_internal_in$calc_loglik(pf$counts_in, 2L, pf$obs_zero) +
    pf$cdms_internal_out$calc_loglik(pf$counts_out, 2L, pf$obs_zero) +
    pf$cdms_immigration1$calc_loglik(pf$counts_im1, 2L, pf$obs_zero) +
    pf$cdms_emigration1$calc_loglik(pf$counts_em1, 2L, pf$obs_zero) +
    pf$cdms_immigration2$calc_loglik(pf$counts_im2, 2L, pf$obs_zero) +
    pf$cdms_emigration2$calc_loglik(pf$counts_em2, 2L, pf$obs_zero)
  expect_equal(ans_obtained, ans_expected)
})


## calc_loglik_init -----------------------------------------------------------

test_that("'calc_loglik' gives expected answer - no regions", {
  set.seed(0)
  pf <- simulated_pf_noreg()
  pf$stock_end <- rpois(n = pf$n_particle, lambda = 10)
  ans_obtained <- pf$calc_loglik_init()
  ans_expected <- pf$cdms_stock$calc_loglik(pf$stock_end, 0L, pf$obs_zero)
  expect_equal(ans_obtained, ans_expected)
})

test_that("'calc_loglik' gives expected answer - with regions", {
  set.seed(0)
  pf <- simulated_pf_withreg()
  pf$stock_end <- matrix(rpois(n = pf$n_particle * pf$n_region, lambda = 10),
    nrow = pf$n_particle
 ) 
  ans_obtained <- pf$calc_loglik_init()
  ans_expected <- pf$cdms_stock$calc_loglik(pf$stock_end, 0L, pf$obs_zero)
  expect_equal(ans_obtained, ans_expected)
})


## calc_logtrans --------------------------------------------------------------

test_that("'calc_logtrans' gives expected answer - no regions", {
    for (seed in seq_len(5)) {
        set.seed(seed)
        pf <- simulated_pf_noreg()
        n <- pf$n_particle
        nr <- pf$n_particle
        pf$exposure <- as.numeric(runif(n = n, min = 0.5, max = 20))
        pf$stock_start <- as.numeric(rpois(n = n, lambda = 100))
        pf$counts_bth <- as.numeric(rpois(n = n, lambda = 10))
        pf$counts_dth <- as.numeric(rpois(n = n, lambda = 10))
        pf$counts_im1 <- as.numeric(rpois(n = n, lambda = 10))
        pf$counts_em1 <- as.numeric(rpois(n = n, lambda = 10))
        pf$rates_bth <- as.numeric(runif(n = n, max = 2))
        pf$rates_dth <- as.numeric(runif(n = n, max = 2))
        pf$rates_im1 <- as.numeric(runif(n = n, max = 2))
        pf$rates_em1 <- as.numeric(runif(n = n, max = 2))
        ans_obtained <- pf$calc_logtrans(1L)[2] ## first particle only
        logtrans_im <- dpois(pf$counts_im1[2], lambda = pf$rates_im1[2], log = TRUE)
        prob_exit <- 1 - exp(-0.5 * (pf$rates_dth[2] + pf$rates_em1[2]))
        prob_dth <- pf$rates_dth[2] / (pf$rates_dth[2] + pf$rates_em1[2])
        count_im <- pf$counts_im1[2]
        if (count_im %% 2 == 0) {
            size <- pf$stock_start[2] + count_im / 2
            logtrans_dth_em <- dmultinom(x = c(pf$counts_dth[2],
                                               pf$counts_em1[2],
                                               size - pf$counts_dth[2] - pf$counts_em1[2]),
                                         size = size,
                                         prob = c(prob_dth * prob_exit, (1-prob_dth) * prob_exit, (1-prob_exit)),
                                         log = TRUE)
        } else {
            size1 <- pf$stock_start[2] + count_im / 2 + 0.5
            logtrans_dth_em <- log(0.5 * dmultinom(x = c(pf$counts_dth[2],
                                                         pf$counts_em1[2],
                                                         size1 - pf$counts_dth[2] - pf$counts_em1[2]),
                                                   size = size1,
                                                   prob = c(prob_dth * prob_exit, (1-prob_dth) * prob_exit, (1-prob_exit)))
                                   + 0.5 * dmultinom(x = c(pf$counts_dth[2],
                                                           pf$counts_em1[2],
                                                           size1 - 1 - pf$counts_dth[2] - pf$counts_em1[2]),
                                                     size = size1 - 1,
                                                     prob = c(prob_dth * prob_exit, (1-prob_dth) * prob_exit, (1-prob_exit))))
        }
        logtrans_bth <- dpois(pf$counts_bth[2], lambda = pf$rates_bth[2] * pf$exposure[2], log = TRUE)
        ans_expected <- logtrans_im + logtrans_dth_em + logtrans_bth
        expect_equal(ans_obtained, ans_expected)
    }
})

test_that("'calc_logtrans' gives expected answer - parallel", {
    for (seed in seq_len(5)) {
        set.seed(seed)
        pf <- simulated_pf_noreg(parallelogram = TRUE)
        n <- pf$n_particle
        nr <- pf$n_particle
        pf$exposure <- as.numeric(runif(n = n, min = 0.5, max = 20))
        pf$stock_start <- as.numeric(rpois(n = n, lambda = 100))
        pf$counts_bth <- as.numeric(rpois(n = n, lambda = 10))
        pf$counts_dth <- as.numeric(rpois(n = n, lambda = 10))
        pf$counts_im1 <- as.numeric(rpois(n = n, lambda = 10))
        pf$counts_em1 <- as.numeric(rpois(n = n, lambda = 10))
        pf$rates_bth <- as.numeric(runif(n = n, max = 2))
        pf$rates_dth <- as.numeric(runif(n = n, max = 2))
        pf$rates_im1 <- as.numeric(runif(n = n, max = 2))
        pf$rates_em1 <- as.numeric(runif(n = n, max = 2))
        pf$stock_end <- as.numeric(rpois(n = n, lambda = 100))
        pf$counts_bth_second <- as.numeric(rpois(n = n, lambda = 10))
        pf$counts_dth_second <- as.numeric(rpois(n = n, lambda = 10))
        pf$counts_im1_second <- as.numeric(rpois(n = n, lambda = 10))
        pf$counts_em1_second <- as.numeric(rpois(n = n, lambda = 10))
        pf$rates_bth_second <- as.numeric(runif(n = n, max = 2))
        pf$rates_dth_second <- as.numeric(runif(n = n, max = 2))
        pf$rates_im1_second <- as.numeric(runif(n = n, max = 2))
        pf$rates_em1_second <- as.numeric(runif(n = n, max = 2))
        ans_obtained <- pf$calc_logtrans(1L)[2] ## first particle only
        logtrans_im <- dpois(pf$counts_im1[2], lambda = pf$rates_im1[2], log = TRUE)
        prob_exit <- 1 - exp(-0.5 * (pf$rates_dth[2] + pf$rates_em1[2]))
        prob_dth <- pf$rates_dth[2] / (pf$rates_dth[2] + pf$rates_em1[2])
        count_im <- pf$counts_im1[2]
        if (count_im %% 2 == 0) {
            size <- pf$stock_start[2] + count_im / 2
            logtrans_dth_em <- dmultinom(x = c(pf$counts_dth[2],
                                               pf$counts_em1[2],
                                               size - pf$counts_dth[2] - pf$counts_em1[2]),
                                         size = size,
                                         prob = c(prob_dth * prob_exit, (1-prob_dth) * prob_exit, (1-prob_exit)),
                                         log = TRUE)
        } else {
            size1 <- pf$stock_start[2] + count_im / 2 + 0.5
            logtrans_dth_em <- log(0.5 * dmultinom(x = c(pf$counts_dth[2],
                                                         pf$counts_em1[2],
                                                         size1 - pf$counts_dth[2] - pf$counts_em1[2]),
                                                   size = size1,
                                                   prob = c(prob_dth * prob_exit, (1-prob_dth) * prob_exit, (1-prob_exit)))
                                   + 0.5 * dmultinom(x = c(pf$counts_dth[2],
                                                           pf$counts_em1[2],
                                                           size1 - 1 - pf$counts_dth[2] - pf$counts_em1[2]),
                                                     size = size1 - 1,
                                                     prob = c(prob_dth * prob_exit, (1-prob_dth) * prob_exit, (1-prob_exit))))
        }
        logtrans_bth <- dpois(pf$counts_bth[2], lambda = pf$rates_bth[2] * pf$exposure[2], log = TRUE)
        logtrans_im_second <- dpois(pf$counts_im1_second[2], lambda = pf$rates_im1_second[2], log = TRUE)
        prob_exit_second <- 1 - exp(-0.5 * (pf$rates_dth_second[2] + pf$rates_em1_second[2]))
        prob_dth_second <- pf$rates_dth_second[2] / (pf$rates_dth_second[2] + pf$rates_em1_second[2])
        count_im_second <- pf$counts_im1_second[2]
        if (count_im_second %% 2 == 0) {
            size <- pf$stock_end[2] + count_im_second / 2
            logtrans_dth_em_second <- dmultinom(x = c(pf$counts_dth_second[2],
                                               pf$counts_em1_second[2],
                                               size - pf$counts_dth_second[2] - pf$counts_em1_second[2]),
                                         size = size,
                                         prob = c(prob_dth_second * prob_exit_second,
                                         (1-prob_dth_second) * prob_exit_second, (1-prob_exit_second)),
                                         log = TRUE)
        } else {
            size1 <- pf$stock_end[2] + count_im_second / 2 + 0.5
            logtrans_dth_em_second <- log(0.5 * dmultinom(x = c(pf$counts_dth_second[2],
                                                         pf$counts_em1_second[2],
                                                         size1 - pf$counts_dth_second[2] - pf$counts_em1_second[2]),
                                                   size = size1,
                                                   prob = c(prob_dth_second * prob_exit_second,
                                                   (1-prob_dth_second) * prob_exit_second, (1-prob_exit_second)))
                                   + 0.5 * dmultinom(x = c(pf$counts_dth_second[2],
                                                           pf$counts_em1_second[2],
                                                           size1 - 1 - pf$counts_dth_second[2] - pf$counts_em1_second[2]),
                                                     size = size1 - 1,
                                                     prob = c(prob_dth_second * prob_exit_second,
                                                     (1-prob_dth_second) * prob_exit_second, (1-prob_exit_second))))
        }
        logtrans_bth_second <- dpois(pf$counts_bth_second[2], lambda = pf$rates_bth_second[2] * pf$exposure_second[2], log = TRUE)
        ans_expected <- logtrans_im + logtrans_dth_em + logtrans_bth +
            logtrans_im_second + logtrans_dth_em_second + logtrans_bth_second
        expect_equal(ans_obtained, ans_expected)
    }
})

test_that("'calc_logtrans' gives expected answer - with regions", {
    set.seed(0)
    pf <- simulated_pf_withreg()
    n <- pf$n_particle * pf$n_region
    nr <- pf$n_particle
    pf$exposure <- matrix(as.numeric(runif(n = n, min = 0.5, max = 20)), nr = nr)
    pf$stock_start <- matrix(as.numeric(rpois(n = n, lambda = 100)), nr = nr)
    pf$counts_bth <- matrix(as.numeric(rpois(n = n, lambda = 10)), nr = nr)
    pf$counts_dth <- matrix(as.numeric(rpois(n = n, lambda = 10)), nr = nr)
    pf$counts_in <- matrix(as.numeric(rpois(n = n, lambda = 10)), nr = nr)
    pf$counts_out <- matrix(as.numeric(rpois(n = n, lambda = 10)), nr = nr)
    pf$counts_im1 <- matrix(as.numeric(rpois(n = n, lambda = 10)), nr = nr)
    pf$counts_em1 <- matrix(as.numeric(rpois(n = n, lambda = 10)), nr = nr)
    pf$rates_bth <- matrix(as.numeric(runif(n = n, max = 2)), nr = nr)
    pf$rates_dth <- matrix(as.numeric(runif(n = n, max = 2)), nr = nr)
    pf$rates_in <- matrix(as.numeric(runif(n = n, max = 2)), nr = nr)
    pf$rates_out <- matrix(as.numeric(runif(n = n, max = 2)), nr = nr)
    pf$rates_im1 <- matrix(as.numeric(runif(n = n, max = 2)), nr = nr)
    pf$rates_em1 <- matrix(as.numeric(runif(n = n, max = 2)), nr = nr)
    ans_obtained <- pf$calc_logtrans()[2] ## single particle in single region
    logtrans_im1 <- dpois(pf$counts_im1[2,], lambda = pf$rates_im1[2,], log = TRUE)
    prob_exit <- 1 - exp(-0.5 * (pf$rates_dth[2,] + pf$rates_em1[2,]))
    prob_dth <- pf$rates_dth[2,] / (pf$rates_dth[2,] + pf$rates_em1[2,])
    logtrans_dth_em <- numeric(length = pf$n_region)
    for (i_region in seq_len(pf$n_region)) {
        if (pf$counts_im1[2,i_region] %% 2 == 0) {
            size <- pf$stock_start[2,i_region] + pf$counts_im1[2,i_region] / 2
            logtrans_dth_em[i_region] <- dmultinom(x = c(pf$counts_dth[2,i_region],
                                                         pf$counts_em1[2,i_region],
                                                         size - pf$counts_dth[2,i_region] - pf$counts_em1[2,i_region]),
                                                   size = size,
                                                   prob = c(prob_dth[i_region] * prob_exit[i_region],
                                                   (1-prob_dth[i_region]) * prob_exit[i_region], (1-prob_exit[i_region])),
                                                   log = TRUE)
        } else {
            size1 <- pf$stock_start[2,i_region] + pf$counts_im1[2,i_region] / 2 + 0.5
            logtrans_dth_em[i_region] <- log(0.5 * dmultinom(x = c(pf$counts_dth[2,i_region],
                                                                   pf$counts_em1[2,i_region],
                                                                   size1 - pf$counts_dth[2,i_region] - pf$counts_em1[2,i_region]),
                                                             size = size1,
                                                             prob = c(prob_dth[i_region] * prob_exit[i_region],
                                                             (1-prob_dth[i_region]) * prob_exit[i_region], (1-prob_exit[i_region])))
                                             + 0.5 * dmultinom(x = c(pf$counts_dth[2,i_region],
                                                                     pf$counts_em1[2,i_region],
                                                                     size1 - 1 - pf$counts_dth[2,i_region] - pf$counts_em1[2,i_region]),
                                                               size = size1 - 1,
                                                               prob = c(prob_dth[i_region] * prob_exit[i_region],
                                                               (1-prob_dth[i_region]) * prob_exit[i_region],
                                                               (1-prob_exit[i_region]))))
        }
    }
    logtrans_bth <- dpois(pf$counts_bth[2,], lambda = pf$rates_bth[2,] * pf$exposure[2,], log = TRUE)
    ans_expected <- sum(logtrans_im1 + logtrans_dth_em + logtrans_bth)
    expect_equal(ans_obtained, ans_expected)
}) 


## calc_logwt_unnorm ----------------------------------------------------------

test_that("'calc_logwt_unnorm' gives expected answer - no regions", {
  set.seed(0)
  pf <- simulated_pf_noreg(has_stock_init = TRUE)
  pf$skip_sampling_init()
  pf$draw_values(i_interval = 1L, is_forecast = FALSE) ## R-style index
  pf$calc_logwt_unnorm(i_interval = 1L, is_forecast = FALSE) ## R-style index
  ans_obtained <- pf$logwt_unnorm[, 2L]
  ans_expected <- pf$calc_loglik(i_interval = 0L) + ## C-style index
    pf$calc_logtrans(1L) -
      pf$calc_logimp(is_forecast = FALSE) +
      pf$calc_logwt_unnorm_prev(i_interval = 1L)
  expect_identical(ans_obtained, ans_expected)
  expect_identical(pf$sum_loglik[[2]], sum(pf$calc_loglik(i_interval = 0L)))
  expect_identical(pf$sum_logtrans[[1]], sum(pf$calc_logtrans(i_interval = 1L)))
  expect_identical(pf$sum_logimp[[2]], sum(pf$calc_logimp(is_forecast = TRUE)))
  expect_identical(pf$sum_logwt_unnorm[[2]], sum(pf$logwt_unnorm[,2]))
})

## calc_logwt_unnorm_init ----------------------------------------------------

test_that("'calc_logwt_unnorm' gives expected answer - no regions", {
  set.seed(0)
  pf <- simulated_pf_noreg()
  pf$draw_values_init()
  pf$calc_logwt_unnorm_init()
  ans_obtained <- pf$logwt_unnorm[, 1L]
  ans_expected <- pf$calc_loglik_init() - pf$calc_logimp_init()
  expect_identical(ans_obtained, ans_expected)
})


## calc_logwt_unnorm_prev -----------------------------------------------------

test_that("'calc_logwt_unnorm' gives expected answer - did not resample", {
  set.seed(0)
  pf <- simulated_pf_noreg()
  pf$draw_values_init()
  pf$calc_logwt_unnorm_init()
  pf$resampled[[1L]] <- FALSE
  ans_obtained <- pf$calc_logwt_unnorm_prev(i_interval = 1L)
  ans_expected <- pf$logwt_unnorm[, 1]
  expect_identical(ans_obtained, ans_expected)
})

test_that("'calc_logwt_unnorm' gives expected answer - resampled previously", {
  set.seed(0)
  pf <- simulated_pf_noreg()
  pf$threshold <- 0.9
  pf$logwt_unnorm[, 3L] <- log(abs(rt(n = pf$n_particle, df = 1)))
  pf$resample(i_interval = 2L)
  expect_true(pf$resampled[3L])
  ans_obtained <- pf$calc_logwt_unnorm_prev(i_interval = 3L)
  ans_expected <- rep(log(1/pf$n_particle), times = pf$n_particle)
  expect_identical(ans_obtained, ans_expected)
})


## calc_n_unique --------------------------------------------------------------

test_that("'calc_n_unique' gives expected answer", {
  set.seed(0)
  pf <- simulated_pf_noreg(n_particle = 6L, n_thin = 2L)
  pf$logwt_unnorm[] <- runif(n = pf$n_particle * (pf$n_interval + 1L), max = 10)
  pf$index_ancestor <- cbind(
    c(1L, 1L, 1L, 1L, 2L, 1L),
    c(1L, 1L, 2L, 2L, 3L, 1L),
    c(1L, 1L, 2L, 2L, 3L, 1L),
    c(1L, 1L, 2L, 2L, 3L, 1L),
    c(1L, 1L, 2L, 3L, 4L, 1L)
  )
  pf$index_output <- c(1L, 3L, 5L)
  pf$calc_n_unique()
  ans_obtained <- pf$n_unique
  ans_expected <- c(2L, 3L, 3L, 3L, 3L)
  expect_equal(ans_obtained, ans_expected)
})
    

## calc_rates -----------------------------------------------------------------

test_that("'calc_rates' gives expected answer - no regions", {
  set.seed(0)
  i_interval <- 2L
  n_particle <- 5L
  pf <- simulated_pf_noreg(n_particle = n_particle)
  pf$calc_rates(i_interval = i_interval)
  expect_identical(
    pf$rates_bth,
    rep(pf$rates_births[[i_interval]], times = n_particle)
  )
  expect_identical(
    pf$rates_dth,
    rep(pf$rates_deaths[[i_interval]], times = n_particle)
  )
  expect_identical(
    pf$rates_im1,
    rep(pf$rates_immigration1[[i_interval]], times = n_particle)
  )
  expect_identical(
    pf$rates_em1,
    rep(pf$rates_emigration1[[i_interval]], times = n_particle)
  )
  expect_identical(
    pf$rates_im2,
    rep(pf$rates_immigration2[[i_interval]], times = n_particle)
  )
  expect_identical(
    pf$rates_em2,
    rep(pf$rates_emigration2[[i_interval]], times = n_particle)
  )
  expect_identical(pf$rates_cim, pf$rates_im1 + pf$rates_im2)
})

test_that("'calc_rates' gives expected answer - parallelogram", {
  set.seed(0)
  i_interval <- 1L
  n_particle <- 5L
  pf <- simulated_pf_noreg(n_particle = n_particle, parallelogram = TRUE)
  pf$calc_rates(i_interval = i_interval)
  expect_identical(
    pf$rates_bth_second,
    rep(pf$rates_births[[i_interval + 1]], times = n_particle)
  )
  expect_identical(
    pf$rates_dth_second,
    rep(pf$rates_deaths[[i_interval + 1]], times = n_particle)
  )
  expect_identical(
    pf$rates_im1_second,
    rep(pf$rates_immigration1[[i_interval + 1]], times = n_particle)
  )
  expect_identical(
    pf$rates_em1_second,
    rep(pf$rates_emigration1[[i_interval + 1]], times = n_particle)
  )
})


test_that("'calc_rates' gives expected answer - with regions", {
  set.seed(0)
  i_interval <- 2L
  n_particle <- 5L
  pf <- simulated_pf_withreg(n_particle = n_particle)
  pf$calc_rates(i_interval = i_interval)
  expect_identical(
    pf$rates_bth,
    matrix(rep(pf$rates_births[, i_interval], each = n_particle),
      nrow = n_particle
    )
  )
  expect_identical(
    pf$rates_dth,
    matrix(rep(pf$rates_deaths[, i_interval], each = n_particle),
      nrow = n_particle
    )
  )
  expect_identical(
    pf$rates_in,
    matrix(rep(pf$rates_internal_in[, i_interval], each = n_particle),
      nrow = n_particle
    )
  )
  expect_identical(
    pf$rates_out,
    matrix(rep(pf$rates_internal_out[, i_interval], each = n_particle),
      nrow = n_particle
    )
  )
  expect_identical(
    pf$rates_im1,
    matrix(rep(pf$rates_immigration1[, i_interval], each = n_particle),
      nrow = n_particle
    )
  )
  expect_identical(
    pf$rates_em1,
    matrix(rep(pf$rates_emigration1[, i_interval], each = n_particle),
      nrow = n_particle
    )
  )
  expect_identical(
    pf$rates_im2,
    matrix(rep(pf$rates_immigration2[, i_interval], each = n_particle),
      nrow = n_particle
    )
  )
  expect_identical(
    pf$rates_em2,
    matrix(rep(pf$rates_emigration2[, i_interval], each = n_particle),
      nrow = n_particle
    )
  )
  expect_identical(pf$rates_cim, pf$rates_im1 + pf$rates_im2)
  expect_identical(pf$rates_cem, pf$rates_em1 + pf$rates_em2)
  expect_identical(pf$rates_cin, pf$rates_in + pf$rates_cim)
  expect_identical(pf$rates_cout, pf$rates_out + pf$rates_cem)
})


## calc_stock_start -----------------------------------------------------------

test_that("'calc_stock_start' gives expected answer - no regions", {
  set.seed(0)
  i_interval <- 2L
  n_particle <- 10L
  pf <- simulated_pf_noreg(n_particle = n_particle, has_stock_init = TRUE)
  pf$index_parent[, i_interval - 1L] <- sample(n_particle)
  pf$calc_stock_start(i_interval = i_interval)
  ans_obtained <- pf$stock_start
  ans_expected <- pf$counts_stock[, i_interval][pf$index_parent[, i_interval - 1L]]
  expect_identical(ans_obtained, ans_expected)
})

test_that("'calc_stock_start' gives expected answer - with regions", {
  set.seed(0)
  i_interval <- 2L
  n_particle <- 10L
  n_region <- 4L
  pf <- simulated_pf_withreg(
    n_particle = n_particle,
    has_stock_init = TRUE,
    n_region = n_region
  )
  pf$index_parent[, i_interval - 1L] <- sample(n_particle)
  pf$counts_stock[, , i_interval - 1L] <- rpois(n = n_particle * n_region, lambda = 10)
  pf$calc_stock_start(i_interval = i_interval)
  ans_obtained <- pf$stock_start
  i <- rep(pf$index_parent[, i_interval - 1L], times = n_region)
  val <- pf$counts_stock[, , i_interval]
  ans_expected <- matrix(val[i],
    nrow = n_particle,
    ncol = n_region
  )
  expect_identical(ans_obtained, ans_expected)
})


## calc_trajectories ----------------------------------------------------------

test_that("'calc_trajectories' gives expected answer - no regions", {
  set.seed(0)
  pf <- simulated_pf_noreg()
  pf$counts_stock[] <- as.numeric(rpois(n = pf$n_particle * (pf$n_interval + 1L), lambda = 100))
  pf$counts_births[] <- as.numeric(rpois(n = pf$n_particle * pf$n_interval, lambda = 100))
  pf$counts_deaths[] <- as.numeric(rpois(n = pf$n_particle * pf$n_interval, lambda = 100))
  pf$counts_immigration1[] <- as.numeric(rpois(n = pf$n_particle * pf$n_interval, lambda = 100))
  pf$counts_emigration1[] <- as.numeric(rpois(n = pf$n_particle * pf$n_interval, lambda = 100))
  pf$counts_immigration2[] <- as.numeric(rpois(n = pf$n_particle * pf$n_interval, lambda = 100))
  pf$counts_emigration2[] <- as.numeric(rpois(n = pf$n_particle * pf$n_interval, lambda = 100))
  pf$index_ancestor <- cbind(
    c(1L, 1L, 1L, 1L, 2L),
    c(1L, 1L, 2L, 2L, 3L),
    c(1L, 1L, 2L, 2L, 3L),
    c(1L, 1L, 2L, 2L, 3L),
    c(1L, 1L, 2L, 3L, 4L)
  )
  pf0 <- pf$clone()
  pf$calc_trajectories()
  i_stock <- pf$index_ancestor + rep(c(0, 5, 10, 15, 20), each = 5)
  i_events <- pf$index_ancestor[, -1] + rep(c(0, 5, 10, 15), each = 5)
  expect_identical(as.numeric(pf$counts_stock), pf0$counts_stock[i_stock])
  expect_identical(as.numeric(pf$counts_births), pf0$counts_births[i_events])
  expect_identical(as.numeric(pf$counts_deaths), pf0$counts_deaths[i_events])
  expect_identical(as.numeric(pf$counts_immigration1), pf0$counts_immigration1[i_events])
  expect_identical(as.numeric(pf$counts_emigration1), pf0$counts_emigration1[i_events])
  expect_identical(as.numeric(pf$counts_immigration2), pf0$counts_immigration2[i_events])
  expect_identical(as.numeric(pf$counts_emigration2), pf0$counts_emigration2[i_events])
})

test_that("'calc_trajectories' gives expected answer - with regions", {
  set.seed(0)
  pf <- simulated_pf_withreg()
  pf$counts_stock[] <- as.numeric(rpois(n = pf$n_particle * pf$n_region * (pf$n_interval + 1L), lambda = 100))
  pf$counts_births[] <- as.numeric(rpois(n = pf$n_particle * pf$n_region * pf$n_interval, lambda = 100))
  pf$counts_deaths[] <- as.numeric(rpois(n = pf$n_particle * pf$n_region * pf$n_interval, lambda = 100))
  pf$counts_internal_in[] <- as.numeric(rpois(n = pf$n_particle * pf$n_region * pf$n_interval, lambda = 100))
  pf$counts_internal_out[] <- as.numeric(rpois(n = pf$n_particle * pf$n_region * pf$n_interval, lambda = 100))
  pf$counts_immigration1[] <- as.numeric(rpois(n = pf$n_particle * pf$n_region * pf$n_interval, lambda = 100))
  pf$counts_emigration1[] <- as.numeric(rpois(n = pf$n_particle * pf$n_region * pf$n_interval, lambda = 100))
  pf$counts_immigration2[] <- as.numeric(rpois(n = pf$n_particle * pf$n_region * pf$n_interval, lambda = 100))
  pf$counts_emigration2[] <- as.numeric(rpois(n = pf$n_particle * pf$n_region * pf$n_interval, lambda = 100))
  pf$index_ancestor <- cbind(
    c(1L, 1L, 1L, 1L, 2L),
    c(1L, 1L, 2L, 2L, 3L),
    c(1L, 1L, 2L, 2L, 3L),
    c(1L, 1L, 2L, 2L, 3L),
    c(1L, 1L, 2L, 3L, 4L)
  )
  pf0 <- pf$clone()
  pf$calc_trajectories()
  i_stock <- pf$index_ancestor + rep(c(0, 5, 10, 15, 20), each = 5)
  i_events <- pf$index_ancestor[, -1] + rep(c(0, 5, 10, 15), each = 5)
  for (i_region in seq_len(pf$n_region)) {
    expect_identical(as.numeric(pf$counts_stock[, i_region, ]), pf0$counts_stock[, i_region, ][i_stock])
    expect_identical(as.numeric(pf$counts_births[, i_region, ]), pf0$counts_births[, i_region, ][i_events])
    expect_identical(as.numeric(pf$counts_deaths[, i_region, ]), pf0$counts_deaths[, i_region, ][i_events])
    expect_identical(as.numeric(pf$counts_internal_in[, i_region, ]), pf0$counts_internal_in[, i_region, ][i_events])
    expect_identical(as.numeric(pf$counts_internal_out[, i_region, ]), pf0$counts_internal_out[, i_region, ][i_events])
    expect_identical(as.numeric(pf$counts_immigration1[, i_region, ]), pf0$counts_immigration1[, i_region, ][i_events])
    expect_identical(as.numeric(pf$counts_emigration1[, i_region, ]), pf0$counts_emigration1[, i_region, ][i_events])
    expect_identical(as.numeric(pf$counts_immigration2[, i_region, ]), pf0$counts_immigration2[, i_region, ][i_events])
    expect_identical(as.numeric(pf$counts_emigration2[, i_region, ]), pf0$counts_emigration2[, i_region, ][i_events])
  }
})


## draw_counts_bth ------------------------------------------------------------

test_that("'draw_counts_bth' gives expected answer - no regions, is dominant", {
  set.seed(0)
  pf <- simulated_pf_noreg()
  pf$calc_rates(1L)
  pf$exposure_approx2 <- runif(n = pf$n_particle, max = 10)
  set.seed(0)
  pf$draw_counts_bth()
  ans_obtained <- pf$counts_bth
  set.seed(0)
  ans_expected <- rpois(
    n = pf$n_particle,
    lambda = pf$exposure_approx2 * pf$rates_bth
  )
  expect_identical(ans_obtained, ans_expected)
})

test_that("'draw_counts_bth' gives expected answer - no regions, not dominant", {
  set.seed(0)
  pf <- simulated_pf_noreg(is_dominant = FALSE)
  pf$calc_rates(1L)
  pf$exposure_approx2 <- runif(n = pf$n_particle, max = 10)
  set.seed(0)
  pf$draw_counts_bth()
  ans_obtained <- pf$counts_bth
  set.seed(0)
  ans_expected <- rep(0, pf$n_particle)
  expect_equal(ans_obtained, ans_expected)
})

test_that("'draw_counts_bth' gives expected answer - with regions, is dominant", {
  set.seed(0)
  pf <- simulated_pf_withreg()
  pf$calc_rates(1L)
  pf$exposure_approx2 <- matrix(runif(n = pf$n_particle * pf$n_region, max = 10),
    nrow = pf$n_particle
  )
  set.seed(0)
  pf$draw_counts_bth()
  ans_obtained <- pf$counts_bth
  set.seed(0)
  ans_expected <- matrix(rpois(
    n = pf$n_particle * pf$n_region,
    lambda = pf$exposure_approx2 * pf$rates_bth
  ),
  nrow = pf$n_particle
  )
  expect_identical(ans_obtained, ans_expected)
})

test_that("'draw_counts_bth' gives expected answer - with regions, not dominant", {
  set.seed(0)
  pf <- simulated_pf_withreg(is_dominant = FALSE)
  pf$calc_rates(1L)
  pf$exposure_approx2 <- matrix(runif(n = pf$n_particle * pf$n_region, max = 10),
    nrow = pf$n_particle
  )
  set.seed(0)
  pf$draw_counts_bth()
  ans_obtained <- pf$counts_bth
  set.seed(0)
  ans_expected <- matrix(0, nrow = pf$n_particle, ncol = pf$n_region)
  expect_equal(ans_obtained, ans_expected)
})


## draw_counts_dth ------------------------------------------------------------

test_that("'draw_counts_dth' gives expected answer - no regions", {
  set.seed(0)
  i_interval <- 1L
  pf <- simulated_pf_noreg(has_stock_init = TRUE)
  pf$index_parent[, i_interval] <- sample(pf$n_particle)
  pf$calc_stock_start(i_interval = i_interval)
  pf$calc_rates(i_interval = i_interval)
  pf$calc_exposure_approx1()
  set.seed(0)
  pf$draw_counts_dth()
  ans_obtained <- pf$counts_dth
  set.seed(0)
  ans_expected <- rpois(
    n = pf$n_particle,
    lambda = pf$rates_dth * pf$exposure_approx1
  )
  expect_identical(ans_obtained, ans_expected)
  ans_obtained <- pf$logimp_counts_dth
  ans_expected <- dpois(
    x = pf$counts_dth,
    lambda = pf$rates_dth * pf$exposure_approx1,
    log = TRUE
  )
  expect_identical(ans_obtained, ans_expected)
})

test_that("'draw_counts_dth' gives expected answer - with regions", {
  set.seed(0)
  i_interval <- 1L
  pf <- simulated_pf_withreg(has_stock_init = TRUE)
  pf$index_parent[, i_interval] <- sample(pf$n_particle)
  pf$calc_stock_start(i_interval = i_interval)
  pf$calc_rates(i_interval = i_interval)
  pf$calc_exposure_approx1()
  set.seed(0)
  pf$draw_counts_dth()
  ans_obtained <- pf$counts_dth
  set.seed(0)
  ans_expected <- matrix(rpois(
    n = pf$n_particle * pf$n_region,
    lambda = pf$rates_dth * pf$exposure_approx1
  ),
  nrow = pf$n_particle
  )
  expect_identical(ans_obtained, ans_expected)
  ans_obtained <- pf$logimp_counts_dth
  ans_expected <- rowSums(matrix(dpois(
    x = pf$counts_dth,
    lambda = pf$rates_dth * pf$exposure_approx1,
    log = TRUE
  ),
  nrow = pf$n_particle
  ))
  expect_identical(ans_obtained, ans_expected)
})


## draw_counts_forecast -------------------------------------------------------

test_that("'draw_counts_forecast' gives expected answer - no regions", {
    for (seed in seq_len(5)) {
        set.seed(seed)
        pf <- simulated_pf_noreg()
        n <- pf$n_particle
        nr <- pf$n_particle
        pf$stock_start <- as.numeric(rpois(n = n, lambda = 100))
        pf$rates_bth <- as.numeric(runif(n = n, max = 2))
        pf$rates_dth <- as.numeric(runif(n = n, max = 2))
        pf$rates_cim <- as.numeric(runif(n = n, max = 2))
        pf$rates_cem <- as.numeric(runif(n = n, max = 2))
        set.seed(seed)
        pf$draw_counts_forecast()
        set.seed(seed)
        expect_equal(pf$counts_im1, rpois(n = pf$n_particle, pf$rates_cim))
        expect_equal(pf$stock_end, pf$stock_start - pf$counts_dth + pf$counts_im1 - pf$counts_em1)
        expect_true(all(pf$counts_bth >= 0))
    }
})


## draw_counts_im_em ----------------------------------------------------------

test_that("'draw_counts_im_em' gives expected answer - no regions, has_one_imem is FALSE", {
  set.seed(0)
  pf <- simulated_pf_noreg()
  pf$calc_rates(1L)
  pf$net_mig <- round(rnorm(n = pf$n_particle, sd = 10))
  pf$gross_mig <- abs(pf$net_mig) + 2 * rpois(n = pf$n_particle, lambda = 5)
  set.seed(0)
  pf$draw_counts_im_em(1L)
  ans_obtained_counts_im_em <- cbind(pf$counts_im1, pf$counts_em1, pf$counts_im2, pf$counts_em2)
  ans_obtained_logimp <- pf$logimp_counts_im_em
  set.seed(0)
  counts_cim <- (pf$gross_mig + pf$net_mig) / 2
  counts_cem <- (pf$gross_mig - pf$net_mig) / 2
  counts_im1 <- rbinom(n = pf$n_particle, size = counts_cim, prob = pf$rates_im1 / pf$rates_cim)
  counts_em1 <- rbinom(n = pf$n_particle, size = counts_cem, prob = pf$rates_em1 / pf$rates_cem)
  ans_expected_counts_im_em <- cbind(
    counts_im1,
    counts_em1,
    counts_cim - counts_im1,
    counts_cem - counts_em1
  )
  dimnames(ans_expected_counts_im_em) <- NULL
  ans_expected_logimp <- dbinom(counts_im1, size = counts_cim, prob = pf$rates_im1 / pf$rates_cim, log = TRUE) +
    dbinom(counts_em1, size = counts_cem, prob = pf$rates_em1 / pf$rates_cem, log = TRUE)
  expect_identical(ans_obtained_counts_im_em, ans_expected_counts_im_em)
  expect_identical(ans_obtained_logimp, ans_expected_logimp)
})

test_that("'draw_counts_im_em' gives expected answer - no regions, has_one_imem is TRUE", {
  set.seed(0)
  pf <- simulated_pf_noreg(has_one_imem = TRUE)
  pf$calc_rates(1L)
  pf$net_mig <- round(rnorm(n = pf$n_particle, sd = 10))
  pf$gross_mig <- abs(pf$net_mig) + 2 * rpois(n = pf$n_particle, lambda = 5)
  set.seed(0)
  pf$draw_counts_im_em(1L)
  ans_obtained_counts_im_em <- cbind(pf$counts_im1, pf$counts_em1, pf$counts_im2, pf$counts_em2)
  ans_obtained_logimp <- pf$logimp_counts_im_em
  set.seed(0)
  counts_cim <- (pf$gross_mig + pf$net_mig) / 2
  counts_cem <- (pf$gross_mig - pf$net_mig) / 2
  counts_im1 <- counts_cim
  counts_em1 <- counts_cem
  ans_expected_counts_im_em <- cbind(
    counts_im1,
    counts_em1,
    counts_cim - counts_im1,
    counts_cem - counts_em1
  )
  dimnames(ans_expected_counts_im_em) <- NULL
  ans_expected_logimp <- rep(0, times = pf$n_particle)
  expect_identical(ans_obtained_counts_im_em, ans_expected_counts_im_em)
  expect_identical(ans_obtained_logimp, ans_expected_logimp)
})

test_that("'draw_counts_im_em' gives expected answer - parallelogram", {
  set.seed(0)
  pf <- simulated_pf_noreg(has_one_imem = TRUE, parallelogram = TRUE)
  pf$calc_rates(1L)
  pf$net_mig <- round(rnorm(n = pf$n_particle, sd = 10))
  pf$gross_mig <- abs(pf$net_mig) + 2 * rpois(n = pf$n_particle, lambda = 5)
  pf$net_mig_second <- round(rnorm(n = pf$n_particle, sd = 10))
  pf$gross_mig_second <- abs(pf$net_mig_second) + 2 * rpois(n = pf$n_particle, lambda = 5)
  set.seed(0)
  pf$draw_counts_im_em(1L)
  ans_obtained_counts_im_em <- cbind(pf$counts_im1, pf$counts_em1, pf$counts_im1_second, pf$counts_em1_second)
  ans_obtained_logimp <- pf$logimp_counts_im_em
  set.seed(0)
  counts_im1 <- (pf$gross_mig + pf$net_mig) / 2
  counts_em1 <- (pf$gross_mig - pf$net_mig) / 2
  counts_im1_second <- (pf$gross_mig_second + pf$net_mig_second) / 2
  counts_em1_second <- (pf$gross_mig_second - pf$net_mig_second) / 2
  ans_expected_counts_im_em <- cbind(
    counts_im1,
    counts_em1,
    counts_im1_second,
    counts_em1_second
  )
  dimnames(ans_expected_counts_im_em) <- NULL
  ans_expected_logimp <- rep(0, times = pf$n_particle)
  expect_identical(ans_obtained_counts_im_em, ans_expected_counts_im_em)
  expect_identical(ans_obtained_logimp, ans_expected_logimp)
})


test_that("'draw_counts_im_em' gives expected answer - with regions, has_one_imem is FALSE", {
  set.seed(0)
  pf <- simulated_pf_withreg(has_one_imem = FALSE)
  pf$calc_rates(1L)
  pf$net_mig <- matrix(round(rnorm(n = pf$n_particle * pf$n_region, sd = 10)),
    nrow = pf$n_particle
  )
  pf$gross_mig <- matrix(abs(pf$net_mig) + 2 * rpois(n = pf$n_particle * pf$n_region, lambda = 5),
    nrow = pf$n_particle
  )
  pf$counts_in <- matrix(0, nrow = pf$n_particle, ncol = pf$n_region)
  pf$counts_out <- matrix(0, nrow = pf$n_particle, ncol = pf$n_region)
  set.seed(0)
  pf$draw_counts_im_em(1L)
  ans_obtained_counts_im_em <- cbind(pf$counts_im1, pf$counts_em1, pf$counts_im2, pf$counts_em2)
  ans_obtained_logimp <- pf$logimp_counts_im_em
  set.seed(0)
  counts_cim <- (pf$gross_mig + pf$net_mig) / 2
  counts_cem <- (pf$gross_mig - pf$net_mig) / 2
  counts_im1 <- matrix(rbinom(n = pf$n_particle * pf$n_region,
                              size = counts_cim,
                              prob = pf$rates_im1 / pf$rates_cim),
                       nrow = pf$n_particle)
  counts_em1 <- matrix(rbinom(n = pf$n_particle * pf$n_region,
                              size = counts_cem,
                              prob = pf$rates_em1 / pf$rates_cem),
                       nrow = pf$n_particle)
  ans_expected_counts_im_em <- cbind(
    counts_im1,
    counts_em1,
    counts_cim - counts_im1,
    counts_cem - counts_em1
  )
  dimnames(ans_expected_counts_im_em) <- NULL
  ans_expected_logimp <- dbinom(counts_im1, size = counts_cim, prob = pf$rates_im1 / pf$rates_cim, log = TRUE) +
    dbinom(counts_em1, size = counts_cem, prob = pf$rates_em1 / pf$rates_cem, log = TRUE)
  ans_expected_logimp <- rowSums(ans_expected_logimp)
  expect_identical(ans_obtained_counts_im_em, ans_expected_counts_im_em)
  expect_identical(ans_obtained_logimp, ans_expected_logimp)
})

test_that("'draw_counts_im_em' gives expected answer - with regions, has_one_imem is TRUE", {
  set.seed(0)
  pf <- simulated_pf_withreg(has_one_imem = TRUE)
  pf$calc_rates(1L)
  pf$net_mig <- matrix(round(rnorm(n = pf$n_particle * pf$n_region, sd = 10)),
                       nrow = pf$n_particle
  )
  pf$gross_mig <- matrix(abs(pf$net_mig) + 2 * rpois(n = pf$n_particle * pf$n_region, lambda = 5),
                         nrow = pf$n_particle
                         )
  pf$counts_in <- matrix(0, nrow = pf$n_particle, ncol = pf$n_region)
  pf$counts_out <- matrix(0, nrow = pf$n_particle, ncol = pf$n_region)
  set.seed(0)
  pf$draw_counts_im_em(1L)
  ans_obtained_counts_im_em <- cbind(pf$counts_im1, pf$counts_em1, pf$counts_im2, pf$counts_em2)
  ans_obtained_logimp <- pf$logimp_counts_im_em
  set.seed(0)
  counts_cim <- (pf$gross_mig + pf$net_mig) / 2
  counts_cem <- (pf$gross_mig - pf$net_mig) / 2
  counts_im1 <- matrix(counts_cim, nrow = pf$n_particle)
  counts_em1 <- matrix(counts_cem, nrow = pf$n_particle)
  ans_expected_counts_im_em <- cbind(
    counts_im1,
    counts_em1,
    counts_cim - counts_im1,
    counts_cem - counts_em1
  )
  dimnames(ans_expected_counts_im_em) <- NULL
  ans_expected_logimp <- matrix(0, nrow = pf$n_particle, ncol = pf$n_region)
  ans_expected_logimp <- rowSums(ans_expected_logimp)
  expect_identical(ans_obtained_counts_im_em, ans_expected_counts_im_em)
  expect_identical(ans_obtained_logimp, ans_expected_logimp)
})


test_that("'draw_counts_im_em' gives expected answer - with regions, is_split_mig is FALSE", {
  set.seed(0)
  pf <- simulated_pf_withreg(is_split_mig = FALSE)
  pf$calc_rates(1L)
  pf$net_mig <- matrix(round(rnorm(n = pf$n_particle * pf$n_region, sd = 10)),
                       nrow = pf$n_particle
  )
  pf$gross_mig <- matrix(abs(pf$net_mig) + 2 * rpois(n = pf$n_particle * pf$n_region, lambda = 5),
                         nrow = pf$n_particle
                         )
  pf$counts_in <- matrix(0, nrow = pf$n_particle, ncol = pf$n_region)
  pf$counts_out <- matrix(0, nrow = pf$n_particle, ncol = pf$n_region)
  set.seed(0)
  pf$draw_counts_im_em(1L)
  ans_obtained_counts_im_em <- cbind(pf$counts_im1, pf$counts_em1, pf$counts_im2, pf$counts_em2)
  ans_obtained_logimp <- pf$logimp_counts_im_em
  set.seed(0)
  counts_cim <- (pf$gross_mig + pf$net_mig) / 2
  counts_cem <- (pf$gross_mig - pf$net_mig) / 2
  counts_im1 <- matrix(counts_cim, nrow = pf$n_particle)
  counts_em1 <- matrix(counts_cem, nrow = pf$n_particle)
  ans_expected_counts_im_em <- cbind(
    counts_im1,
    counts_em1,
    counts_cim - counts_im1,
    counts_cem - counts_em1
  )
  dimnames(ans_expected_counts_im_em) <- NULL
  ans_expected_logimp <- matrix(0, nrow = pf$n_particle, ncol = pf$n_region)
  ans_expected_logimp <- rowSums(ans_expected_logimp)
  expect_identical(ans_obtained_counts_im_em, ans_expected_counts_im_em)
  expect_identical(ans_obtained_logimp, ans_expected_logimp)
})


## draw_counts_in_out ---------------------------------------------------------

test_that("'draw_counts_in_out' gives expected answer - is_split_mig is TRUE", {
  set.seed(0)
  pf <- simulated_pf_withreg()
  pf$calc_rates(1L)
  pf$exposure_approx2 <- runif(n = pf$n_particle, max = 10)
  set.seed(0)
  pf$draw_counts_in_out()
  ans_obtained_counts_in <- pf$counts_in
  ans_obtained_counts_out <- pf$counts_out
  ans_obtained_logimp <- pf$logimp_counts_in_out
  set.seed(0)
  ans_expected_counts_out <- matrix(rpois(
    n = pf$n_particle * pf$n_region,
    lambda = pf$exposure_approx2 * pf$rates_out
  ),
  nrow = pf$n_particle
  )
  totals <- rowSums(pf$counts_out)
  ans_expected_counts_in <- t(sapply(
    seq_len(pf$n_particle),
    function(i) rmultinom(1, size = totals[[i]], prob = pf$rates_in[i, ])
  ))
  ans_expected_logimp <- rowSums(matrix(dpois(pf$counts_out,
    lambda = pf$exposure_approx2 * pf$rates_out,
    log = TRUE
  ),
  nrow = pf$n_particle
  )) +
    sapply(
      seq_len(pf$n_particle),
      function(i) {
        dmultinom(pf$counts_in[i, ],
          size = totals[[i]],
          prob = pf$rates_in[i, ],
          log = TRUE
        )
      }
    )
  expect_identical(ans_obtained_counts_in, ans_expected_counts_in)
  expect_identical(ans_obtained_counts_out, ans_expected_counts_out)
  expect_identical(ans_obtained_logimp, ans_expected_logimp)
})

test_that("'draw_counts_in_out' gives expected answer - is_split_mig is FALSE", {
    set.seed(0)
    pf <- simulated_pf_withreg(is_split_mig = FALSE)
    pf$calc_rates(1L)
    pf$exposure_approx2 <- runif(n = pf$n_particle, max = 10)
    set.seed(0)
    pf$draw_counts_in_out()
    ans_obtained_counts_in <- pf$counts_in
    ans_obtained_counts_out <- pf$counts_out
    ans_obtained_logimp <- pf$logimp_counts_in_out
    set.seed(0)
    ans_expected_counts_out <- matrix(0,
                                      nrow = pf$n_particle,
                                      ncol = pf$n_region
                                      )
    ans_expected_counts_in <- matrix(0,
                                     nrow = pf$n_particle,
                                     ncol = pf$n_region
                                     )
    ans_expected_logimp <- rep(0, times = pf$n_particle)
    expect_identical(ans_obtained_counts_in, ans_expected_counts_in)
    expect_identical(ans_obtained_counts_out, ans_expected_counts_out)
    expect_identical(ans_obtained_logimp, ans_expected_logimp)
})


## draw_index_output ----------------------------------------------------------

test_that("'draw_index_output' gives expected answer", {
  set.seed(0)
  pf <- simulated_pf_noreg(n_particle = 10L, n_thin = 2L)
  pf$draw_index_output()
  expect_identical(length(pf$index_output), 5L)
  expect_true(all(pf$index_output %in% 1:10))
  expect_identical(length(unique(pf$index_output)), 5L)
})  


## draw_gross_mig -------------------------------------------------------------

test_that("'draw_gross_mig' gives expected answer - no regions", {
  set.seed(10)
  pf <- simulated_pf_noreg()
  pf$calc_rates(1L)
  pf$exposure_approx2 <- runif(n = pf$n_particle, max = 10)
  pf$net_mig <- as.numeric(round(rnorm(n = pf$n_particle, sd = 10)))
  set.seed(0)
  pf$draw_gross_mig(1L)
  ans_obtained_gross_mig <- pf$gross_mig
  ans_obtained_logimp <- pf$logimp_gross_mig
  set.seed(0)
  ans_expected_gross_mig <- rpoistr(
    n = pf$n_particle,
    lambda = pf$rates_cim + pf$rates_cem * pf$exposure_approx2,
    lower = abs(pf$net_mig)
  )
  different_parity <- (ans_expected_gross_mig %% 2L) != (abs(pf$net_mig) %% 2L)
  ans_expected_gross_mig[different_parity] <- ans_expected_gross_mig[different_parity] + 1
  logprob <- rbind( dpoistr(ans_expected_gross_mig,
    lambda = pf$rates_cim + pf$rates_cem * pf$exposure_approx2,
    lower = abs(pf$net_mig),
    use_log = TRUE
  ),
    ifelse(ans_expected_gross_mig > abs(pf$net_mig),
      dpoistr(ans_expected_gross_mig - 1,
        lambda = pf$rates_cim + pf$rates_cem * pf$exposure_approx2,
        lower = abs(pf$net_mig),
        use_log = TRUE
      ),
      -Inf
      ))
  ans_expected_logimp <- sapply(seq_len(pf$n_particle),
                                function(i) log_sum_exp_2(logprob[1,i], logprob[2,i]))
  expect_identical(ans_obtained_gross_mig, ans_expected_gross_mig)
  expect_identical(ans_obtained_logimp, ans_expected_logimp)
})

test_that("'draw_gross_mig' gives expected answer - parallelogram", {
    set.seed(10)
    pf <- simulated_pf_noreg(parallelogram = TRUE)
    pf$calc_rates(1L)
    pf$exposure_approx2 <- runif(n = pf$n_particle, max = 10)
    pf$net_mig <- as.numeric(round(rnorm(n = pf$n_particle, sd = 10)))
    pf$net_mig_second <- as.numeric(round(rnorm(n = pf$n_particle, sd = 10)))
    set.seed(0)
    pf$draw_gross_mig(1L)
    ans_obtained_gross_mig <- c(pf$gross_mig, pf$gross_mig_second)
    ans_obtained_logimp <- pf$logimp_gross_mig
    set.seed(0)
    ans_expected_gross_mig <- rpoistr(
        n = pf$n_particle,
        lambda = pf$rates_cim + pf$rates_cem * pf$exposure_approx2,
        lower = abs(pf$net_mig)
    )
    different_parity <- (ans_expected_gross_mig %% 2L) != (abs(pf$net_mig) %% 2L)
    ans_expected_gross_mig[different_parity] <- ans_expected_gross_mig[different_parity] + 1
    logprob <- rbind( dpoistr(ans_expected_gross_mig,
                              lambda = pf$rates_cim + pf$rates_cem * pf$exposure_approx2,
                              lower = abs(pf$net_mig),
                              use_log = TRUE
                              ),
                     ifelse(ans_expected_gross_mig > abs(pf$net_mig),
                            dpoistr(ans_expected_gross_mig - 1,
                                    lambda = pf$rates_cim + pf$rates_cem * pf$exposure_approx2,
                                    lower = abs(pf$net_mig),
                                    use_log = TRUE
                                    ),
                            -Inf
                            ))
    ans_expected_logimp <- sapply(seq_len(pf$n_particle),
                                  function(i) log_sum_exp_2(logprob[1,i], logprob[2,i]))
    ans_expected_gross_mig_second <- rpoistr(
        n = pf$n_particle,
        lambda = pf$rates_im1_second + pf$rates_em1_second * pf$exposure_approx2_second,
        lower = abs(pf$net_mig_second)
    )
    different_parity <- (ans_expected_gross_mig_second %% 2L) != (abs(pf$net_mig_second) %% 2L)
    ans_expected_gross_mig_second[different_parity] <- ans_expected_gross_mig_second[different_parity] + 1
    logprob <- rbind( dpoistr(ans_expected_gross_mig_second,
                              lambda = pf$rates_im1_second + pf$rates_em1_second * pf$exposure_approx2_second,
                              lower = abs(pf$net_mig_second),
                              use_log = TRUE
                              ),
                     ifelse(ans_expected_gross_mig_second > abs(pf$net_mig_second),
                            dpoistr(ans_expected_gross_mig_second - 1,
                                    lambda = pf$rates_im1_second + pf$rates_em1_second * pf$exposure_approx2_second,
                                    lower = abs(pf$net_mig_second),
                                    use_log = TRUE
                                    ),
                            -Inf
                            ))
    ans_expected_logimp_second <- sapply(seq_len(pf$n_particle),
                                         function(i) log_sum_exp_2(logprob[1,i], logprob[2,i]))
    expect_identical(ans_obtained_gross_mig, c(ans_expected_gross_mig, ans_expected_gross_mig_second))
    expect_identical(ans_obtained_logimp, ans_expected_logimp + ans_expected_logimp_second)
})


test_that("'draw_gross_mig' gives expected answer - with regions, is_split_mig is TRUE", {
    set.seed(0)
    pf <- simulated_pf_withreg()
    pf$calc_rates(1L)
    pf$exposure_approx2 <- runif(n = pf$n_particle, max = 10)
    pf$net_mig <- matrix(as.numeric(round(rnorm(n = pf$n_particle * pf$n_region, sd = 10))),
                         nrow = pf$n_particle
                         )
    pf$counts_in <- matrix(as.numeric(rpois(n = pf$n_particle * pf$n_region, lambda = 10)),
                           nrow = pf$n_particle
                           )
    pf$counts_out <- matrix(as.numeric(rpois(n = pf$n_particle * pf$n_region, lambda = 10)),
                            nrow = pf$n_particle
                            )
    set.seed(0)
    pf$draw_gross_mig(1L)
    ans_obtained_gross_mig <- pf$gross_mig
    ans_obtained_logimp <- pf$logimp_gross_mig
    set.seed(0)
    net_mig_int <- pf$net_mig - pf$counts_in + pf$counts_out
    ans_expected_gross_mig <- rpoistr(
        n = pf$n_particle * pf$n_region,
        lambda = pf$rates_cim + pf$rates_cem * pf$exposure_approx2,
        lower = abs(net_mig_int)
    )
    different_parity <- (ans_expected_gross_mig %% 2L) != (abs(net_mig_int) %% 2L)
    ans_expected_gross_mig[different_parity] <- ans_expected_gross_mig[different_parity] + 1
    ans_expected_gross_mig <- matrix(ans_expected_gross_mig, nrow = pf$n_particle)
    logprob <- rbind(as.numeric(dpoistr(ans_expected_gross_mig,
                                        lambda = pf$rates_cim + pf$rates_cem * pf$exposure_approx2,
                                        lower = abs(net_mig_int),
                                        use_log = TRUE
                                        )),
                     as.numeric(ifelse(ans_expected_gross_mig > abs(net_mig_int),
                                       dpoistr(ans_expected_gross_mig - 1,
                                               lambda = pf$rates_cim + pf$rates_cem * pf$exposure_approx2,
                                               lower = abs(net_mig_int),
                                               use_log = TRUE
                                               ),
                                       -Inf
                                       ))
                     )
    ans_expected_logimp <- sapply(seq_len(pf$n_particle * pf$n_region),
                                  function(i) log_sum_exp_2(logprob[1,i], logprob[2,i]))
    ans_expected_logimp <- rowSums(matrix(ans_expected_logimp, nrow = pf$n_particle))
    expect_identical(ans_obtained_gross_mig, ans_expected_gross_mig)
    expect_identical(ans_obtained_logimp, ans_expected_logimp)
})

test_that("'draw_gross_mig' gives expected answer - with regions, is_split_mig is FALSE", {
    set.seed(0)
    pf <- simulated_pf_withreg(is_split_mig = FALSE)
    pf$calc_rates(1L)
    pf$exposure_approx2 <- runif(n = pf$n_particle, max = 10)
    pf$net_mig <- matrix(as.numeric(round(rnorm(n = pf$n_particle * pf$n_region, sd = 10))),
                         nrow = pf$n_particle
                         )
    set.seed(0)
    pf$draw_gross_mig(1L)
    ans_obtained_gross_mig <- pf$gross_mig
    ans_obtained_logimp <- pf$logimp_gross_mig
    set.seed(0)
    ans_expected_gross_mig <- rpoistr(
        n = pf$n_particle * pf$n_region,
        lambda = pf$rates_cin + pf$rates_cout * pf$exposure_approx2,
        lower = abs(pf$net_mig)
    )
    different_parity <- (ans_expected_gross_mig %% 2L) != (abs(pf$net_mig) %% 2L)
    ans_expected_gross_mig[different_parity] <- ans_expected_gross_mig[different_parity] + 1
    ans_expected_gross_mig <- matrix(ans_expected_gross_mig, nrow = pf$n_particle)
    logprob <- rbind(as.numeric(dpoistr(ans_expected_gross_mig,
                                        lambda = pf$rates_cin + pf$rates_cout * pf$exposure_approx2,
                                        lower = abs(pf$net_mig),
                                        use_log = TRUE
                                        )),
                     as.numeric(ifelse(ans_expected_gross_mig > abs(pf$net_mig),
                                       dpoistr(ans_expected_gross_mig - 1,
                                               lambda = pf$rates_cin + pf$rates_cout * pf$exposure_approx2,
                                               lower = abs(pf$net_mig),
                                               use_log = TRUE
                                               ),
                                       -Inf
                                       ))
                     )
    ans_expected_logimp <- sapply(seq_len(pf$n_particle * pf$n_region),
                                  function(i) log_sum_exp_2(logprob[1,i], logprob[2,i]))
    ans_expected_logimp <- rowSums(matrix(ans_expected_logimp, nrow = pf$n_particle))
    expect_identical(ans_obtained_gross_mig, ans_expected_gross_mig)
    expect_identical(ans_obtained_logimp, ans_expected_logimp)
})


## draw_net_mig_parallelogram -------------------------------------------------

test_that("'draw_net_mig_parallelogram' gives expected answer", {
    for (seed in 1:5) {
        set.seed(seed)
        pf <- simulated_pf_noreg(parallelogram = TRUE)
        pf$net_mig_combined <- round(rnorm(n = pf$n_particle, sd = 10))
        pf$stock_start <- rpois(pf$n_particle, lambda = 5)
        pf$counts_dth <- rpois(pf$n_particle, lambda = 3)
        pf$counts_dth_second <- rpois(pf$n_particle, lambda = 3)
        pf$stock_end_second <- pf$stock_start - pf$counts_dth - pf$counts_dth_second + pf$net_mig_combined
        pf$rates_im1 <- runif(pf$n_particle, max = 3)
        pf$rates_em1 <- runif(pf$n_particle, max = 0.5)
        pf$rates_im1_second <- runif(pf$n_particle, max = 3)
        pf$rates_em1_second <- runif(pf$n_particle, max = 0.5)
        logimp <- runif(pf$n_particle, min = -5, max = -1)
        pf$logimp_stock_end_net_mig <- logimp
        pf$calc_exposure_approx_parallelogram()
        set.seed(seed)
        pf$draw_net_mig_parallelogram()
        ans_obtained <- list(pf$net_mig, pf$net_mig_second, pf$stock_end, pf$logimp_stock_end_net_mig)
        set.seed(seed)
        diff <- rskeltr(n = pf$n_particle,
                        mu1 = pf$rates_im1 + pf$rates_em1_second * pf$exposure_approx_second,
                        mu2 = pf$rates_im1_second + pf$rates_em1 * pf$exposure_approx_first,
                        lower = -pf$net_mig_combined + 2 * (pf$counts_dth - pf$stock_start))
        different_parity <- (diff %% 2L) != (pf$net_mig_combined %% 2L)
        diff[different_parity] <- diff[different_parity] + 1
        net_mig <- (pf$net_mig_combined + diff) / 2
        net_mig_second <- pf$net_mig_combined - net_mig
        stock_end <- pf$stock_start - pf$counts_dth + net_mig
        logprob <- rbind(dskeltr(diff,
                                 mu1 = pf$rates_im1 + pf$rates_em1_second * pf$exposure_approx_second,
                                 mu2 = pf$rates_im1_second + pf$rates_em1 * pf$exposure_approx_first,
                                 lower = -pf$net_mig_combined + 2 * (pf$counts_dth - pf$stock_start),
                                 use_log = TRUE
                                 ),
                         ifelse(diff > -pf$net_mig_combined + 2 * (pf$counts_dth - pf$stock_start),
                                dskeltr(diff - 1,
                                        mu1 = pf$rates_im1 + pf$rates_em1_second * pf$exposure_approx_second,
                                        mu2 = pf$rates_im1_second + pf$rates_em1 * pf$exposure_approx_first,
                                        lower = -pf$net_mig_combined + 2 * (pf$counts_dth - pf$stock_start),
                                        use_log = TRUE
                                        ),
                                -Inf
                                ))
        diff_logimp <- sapply(seq_len(pf$n_particle),
                              function(i) log_sum_exp_2(logprob[1,i], logprob[2,i]))
        ans_expected <- list(net_mig, net_mig_second, stock_end, logimp + diff_logimp)
        expect_equal(ans_obtained, ans_expected)        
    }
})
 

## draw_rates -----------------------------------------------------------------

test_that("'draw_rates' gives expected answer - no regions, disp defaults to 0", {
  set.seed(0)
  i_interval <- 2L
  n_particle <- 5L
  pf <- simulated_pf_noreg(n_particle = n_particle)
  pf$draw_rates(i_interval = i_interval)
  expect_identical(
    pf$rates_bth,
    rep(pf$rates_births[[i_interval]], times = n_particle)
  )
  expect_identical(
    pf$rates_dth,
    rep(pf$rates_deaths[[i_interval]], times = n_particle)
  )
  expect_identical(
    pf$rates_im1,
    rep(pf$rates_immigration1[[i_interval]], times = n_particle)
  )
  expect_identical(
    pf$rates_em1,
    rep(pf$rates_emigration1[[i_interval]], times = n_particle)
  )
  expect_identical(
    pf$rates_im2,
    rep(pf$rates_immigration2[[i_interval]], times = n_particle)
  )
  expect_identical(
    pf$rates_em2,
    rep(pf$rates_emigration2[[i_interval]], times = n_particle)
  )
  expect_identical(pf$rates_cim, pf$rates_im1 + pf$rates_im2)
})


test_that("'draw_rates' gives expected answer - no regions, disp non-default", {
  set.seed(0)
  i_interval <- 2L
  n_particle <- 5L
  pf <- simulated_pf_noreg(n_particle = n_particle)
  pf$disp_births <- 0.1
  pf$disp_deaths <- 0.1
  pf$disp_immigration1 <- 0.1
  pf$disp_immigration2 <- 0
  pf$disp_emigration1 <- 0.1
  pf$disp_emigration2 <- 0
  set.seed(100)
  pf$draw_rates(i_interval = i_interval)
  set.seed(100)
  expect_identical(
      pf$rates_bth,
      draw_rates_inner(n = n_particle,
                       mean = pf$rates_births[[i_interval]],
                       disp = pf$disp_births)
  )
  expect_identical(
    pf$rates_dth,
      draw_rates_inner(n = n_particle,
                       mean = pf$rates_deaths[[i_interval]],
                       disp = pf$disp_deaths)
  )
  expect_identical(
    pf$rates_im1,
      draw_rates_inner(n = n_particle,
                       mean = pf$rates_immigration1[[i_interval]],
                       disp = pf$disp_immigration1)
  )
  expect_identical(
    pf$rates_em1,
      draw_rates_inner(n = n_particle,
                       mean = pf$rates_emigration1[[i_interval]],
                       disp = pf$disp_emigration1)
  )
  expect_identical(
    pf$rates_im2,
      draw_rates_inner(n = n_particle,
                       mean = pf$rates_immigration2[[i_interval]],
                       disp = pf$disp_immigration2)
  )
  expect_identical(
    pf$rates_em2,
      draw_rates_inner(n = n_particle,
                       mean = pf$rates_emigration2[[i_interval]],
                       disp = pf$disp_emigration2)
  )
  expect_identical(pf$rates_cim, pf$rates_im1 + pf$rates_im2)
  expect_identical(pf$rates_cem, pf$rates_em1 + pf$rates_em2)
})


test_that("'draw_rates' gives expected answer - parallelogram", {
  set.seed(0)
  i_interval <- 1L
  n_particle <- 5L
  pf <- simulated_pf_noreg(n_particle = n_particle, parallelogram = TRUE)
  pf$draw_rates(i_interval = i_interval)
  expect_identical(
    pf$rates_bth,
    rep(pf$rates_births[[i_interval]], times = n_particle)
  )
  expect_identical(
    pf$rates_dth,
    rep(pf$rates_deaths[[i_interval]], times = n_particle)
  )
  expect_identical(
    pf$rates_im1,
    rep(pf$rates_immigration1[[i_interval]], times = n_particle)
  )
  expect_identical(
    pf$rates_em1,
    rep(pf$rates_emigration1[[i_interval]], times = n_particle)
  )
  expect_identical(
    pf$rates_im2,
    rep(pf$rates_immigration2[[i_interval]], times = n_particle)
  )
  expect_identical(
    pf$rates_em2,
    rep(pf$rates_emigration2[[i_interval]], times = n_particle)
  )
  expect_identical(pf$rates_cim, pf$rates_im1 + pf$rates_im2)
  expect_identical(
    pf$rates_bth_second,
    rep(pf$rates_births[[i_interval + 1L]], times = n_particle)
  )
  expect_identical(
    pf$rates_dth_second,
    rep(pf$rates_deaths[[i_interval + 1L]], times = n_particle)
  )
  expect_identical(
    pf$rates_im1_second,
    rep(pf$rates_immigration1[[i_interval + 1]], times = n_particle)
  )
  expect_identical(
    pf$rates_em1_second,
    rep(pf$rates_emigration1[[i_interval + 1L]], times = n_particle)
  )
})



## draw_stock_end_init --------------------------------------------------------

test_that("'draw_stock_end_init' gives expected answer - no regions", {
  set.seed(0)
  pf <- simulated_pf_noreg()
  set.seed(0)
  pf$draw_stock_end_init()
  ans_obtained <- pf$stock_end
  set.seed(0)
  ans_expected <- pf$cdms_stock$draw_counts_true(0L, pf$n_particle)
  expect_identical(ans_obtained, ans_expected)
})

test_that("'draw_stock_end_init' gives expected answer - with regions", {
  set.seed(0)
  pf <- simulated_pf_withreg()
  set.seed(0)
  pf$draw_stock_end_init()
  ans_obtained <- pf$stock_end
  set.seed(0)
  ans_expected <- pf$cdms_stock$draw_counts_true(0L, pf$n_particle, pf$n_region)
  expect_identical(ans_obtained, ans_expected)
})

test_that("'draw_stock_end_init' throws expected error when missing population data - no regions", {
  pf <- simulated_pf_noreg()
  pf$cdms_stock <- new_CdmsNoreg()
  expect_error(
    pf$draw_stock_end_init(),
    "problem with cohort '2000' and sex/gender 'Female' : no data to estimate initial stock"
  )
})

test_that("'draw_stock_end_init' throws expected error when missing population data - with regions", {
  pf <- simulated_pf_withreg()
  pf$cdms_stock <- new_CdmsWithreg()
  expect_error(
    pf$draw_stock_end_init(),
    "problem with cohort '2000' and sex/gender 'Female' : no data to estimate initial stock"
  )
})

test_that("'draw_stock_end_init' throws expected error with partial population data - with regions", {
  pf <- simulated_pf_withreg()
  counts_data_stock <- as.numeric(stats::rpois(
    n = pf$n_region * (pf$n_interval + 1L),
    lambda = 20
  ))
  counts_data_stock <- matrix(counts_data_stock,
    nrow = pf$n_region
  )
  counts_data_stock[1, 1] <- NA
  cdm_stock <- new_CdmWithregPoibin(
    counts_data = counts_data_stock,
    prob = 0.95
  )
  pf$cdms_stock <- new_CdmsWithreg(list(cdm_stock))
  expect_error(
    pf$draw_stock_end_init(),
    "problem with cohort '2000' and sex/gender 'Female' : no data to estimate initial stock"
  )
})


## draw_stock_end_net_mig_noreg -----------------------------------------------

test_that("'draw_stock_end_net_mig' gives expected answer - no regions, have stock data", {
  set.seed(0)
  pf <- simulated_pf_noreg(has_stock_init = TRUE)
  pf$index_parent[, 1L] <- sample(pf$n_particle)
  pf$calc_stock_start(1L)
  pf$calc_counts_dth(1L)
  pf$calc_rates(1L)
  pf$calc_exposure_approx1()
  set.seed(0)
  pf$draw_stock_end_net_mig(i_interval = 1L, is_forecast = FALSE)
  ans_obtained_stock_end <- pf$stock_end
  ans_obtained_net_mig <- pf$net_mig
  ans_obtained_logimp <- pf$logimp_stock_end_net_mig
  set.seed(0)
  ans_expected_stock_end <- pf$cdms_stock$draw_counts_true(
    i_interval = 1L,
    n_particle = pf$n_particle
  )
  ans_expected_net_mig <- ans_expected_stock_end - pf$stock_start + pf$counts_dth
  ans_expected_logimp <- pf$cdms_stock$calc_logimp(
    counts_true = ans_expected_stock_end,
    i_interval = 1L
  )
  expect_identical(ans_obtained_stock_end, ans_expected_stock_end)
  expect_identical(ans_obtained_net_mig, ans_expected_net_mig)
  expect_identical(ans_obtained_logimp, ans_expected_logimp)
})

test_that("'draw_stock_end_net_mig' gives expected answer - no regions, no stock data", {
  set.seed(0)
  pf <- simulated_pf_noreg(has_stock_init = TRUE)
  pf$cdms_stock <- new_CdmsNoreg()
  pf$index_parent[, 1L] <- sample(pf$n_particle)
  pf$calc_stock_start(1L)
  pf$calc_counts_dth(1L)
  pf$calc_rates(1L)
  pf$calc_exposure_approx1()
  set.seed(0)
  pf$draw_stock_end_net_mig(i_interval = 1L, is_forecast = FALSE)
  ans_obtained_stock_end <- pf$stock_end
  ans_obtained_net_mig <- pf$net_mig
  ans_obtained_logimp <- pf$logimp_stock_end_net_mig
  set.seed(0)
  lower <- pf$counts_dth - pf$stock_start
  ans_expected_net_mig <- rskeltr(
    n = pf$n_particle,
    mu1 = pf$rates_cim,
    mu2 = pf$rates_cem * pf$exposure_approx1,
    lower = lower
  )
  ans_expected_stock_end <- ans_expected_net_mig + pf$stock_start - pf$counts_dth
  ans_expected_logimp <- dskeltr(
    x = ans_expected_net_mig,
    mu1 = pf$rates_cim,
    mu2 = pf$rates_cem * pf$exposure_approx1,
    lower = lower,
    use_log = TRUE
  )
  expect_identical(ans_obtained_stock_end, ans_expected_stock_end)
  expect_identical(ans_obtained_net_mig, ans_expected_net_mig)
  expect_identical(ans_obtained_logimp, ans_expected_logimp)
})

test_that("'draw_stock_end_net_mig' gives expected answer - no regions, forecasting", {
  set.seed(0)
  pf <- simulated_pf_noreg(has_stock_init = TRUE)
  pf$index_parent[, 1L] <- sample(pf$n_particle)
  pf$calc_stock_start(1L)
  pf$calc_counts_dth(1L)
  pf$calc_rates(1L)
  pf$calc_exposure_approx1()
  set.seed(0)
  pf$draw_stock_end_net_mig(i_interval = 1L, is_forecast = TRUE)
  ans_obtained_stock_end <- pf$stock_end
  ans_obtained_net_mig <- pf$net_mig
  ans_obtained_logimp <- pf$logimp_stock_end_net_mig
  set.seed(0)
  lower <- pf$counts_dth - pf$stock_start
  ans_expected_net_mig <- rskeltr(
    n = pf$n_particle,
    mu1 = pf$rates_cim,
    mu2 = pf$rates_cem * pf$exposure_approx1,
    lower = lower
  )
  ans_expected_stock_end <- ans_expected_net_mig + pf$stock_start - pf$counts_dth
  ans_expected_logimp <- dskeltr(
    x = ans_expected_net_mig,
    mu1 = pf$rates_cim,
    mu2 = pf$rates_cem * pf$exposure_approx1,
    lower = lower,
    use_log = TRUE
  )
  expect_identical(ans_obtained_stock_end, ans_expected_stock_end)
  expect_identical(ans_obtained_net_mig, ans_expected_net_mig)
  expect_identical(ans_obtained_logimp, ans_expected_logimp)
})

test_that("'draw_stock_end_net_mig' gives expected answer - with regions, have complete stock data", {
  set.seed(0)
  pf <- simulated_pf_withreg(has_stock_init = TRUE)
  pf$index_parent[, 1L] <- sample(pf$n_particle)
  pf$calc_stock_start(1L)
  pf$calc_exposure_approx1()
  pf$calc_counts_dth(1L)
  pf$calc_rates(1L)
  set.seed(0)
  pf$draw_stock_end_net_mig(i_interval = 1L, is_forecast = FALSE)
  ans_obtained_stock_end <- pf$stock_end
  ans_obtained_net_mig <- pf$net_mig
  ans_obtained_logimp <- pf$logimp_stock_end_net_mig
  set.seed(0)
  ans_expected_stock_end <- pf$cdms_stock$draw_counts_true(
    i_interval = 1L,
    n_particle = pf$n_particle,
    n_region = pf$n_region
  )
  ans_expected_net_mig <- ans_expected_stock_end - pf$stock_start + pf$counts_dth
  ans_expected_logimp <- rowSums(pf$cdms_stock$calc_logimp(
    counts_true = ans_expected_stock_end,
    i_interval = 1L
  ))
  expect_identical(ans_obtained_stock_end, ans_expected_stock_end)
  expect_identical(ans_obtained_net_mig, ans_expected_net_mig)
  expect_identical(ans_obtained_logimp, ans_expected_logimp)
})

test_that("'draw_stock_end_net_mig' gives expected answer - with regions, have no stock data", {
  set.seed(0)
  pf <- simulated_pf_withreg(has_stock_init = TRUE)
  pf$cdms_stock <- new_CdmsWithreg()
  pf$index_parent[, 1L] <- sample(pf$n_particle)
  pf$calc_stock_start(1L)
  pf$calc_counts_dth(1L)
  pf$calc_rates(1L)
  pf$calc_exposure_approx1()
  set.seed(0)
  pf$draw_stock_end_net_mig(i_interval = 1L, is_forecast = FALSE)
  ans_obtained_stock_end <- pf$stock_end
  ans_obtained_net_mig <- pf$net_mig
  ans_obtained_logimp <- pf$logimp_stock_end_net_mig
  set.seed(0)
  lower <- pf$counts_dth - pf$stock_start
  ans_expected_net_mig <- matrix(rskeltr(
    n = pf$n_particle * pf$n_region,
    mu1 = pf$rates_cin,
    mu2 = pf$rates_cout * pf$exposure_approx1,
    lower = lower
  ),
  nrow = pf$n_particle
  )
  ans_expected_stock_end <- ans_expected_net_mig + pf$stock_start - pf$counts_dth
  ans_expected_logimp <- rowSums(matrix(dskeltr(
    x = ans_expected_net_mig,
    mu1 = pf$rates_cin,
    mu2 = pf$rates_cout * pf$exposure_approx1,
    lower = lower,
    use_log = TRUE
  ),
  nrow = pf$n_particle
  ))
  expect_identical(ans_obtained_stock_end, ans_expected_stock_end)
  expect_identical(ans_obtained_net_mig, ans_expected_net_mig)
  expect_identical(ans_obtained_logimp, ans_expected_logimp)
})

test_that("'draw_stock_end_net_mig' gives expected answer - with regions, forecasting", {
  set.seed(0)
  pf <- simulated_pf_withreg(has_stock_init = TRUE)
  pf$cdms_stock <- new_CdmsWithreg()
  pf$index_parent[, 1L] <- sample(pf$n_particle)
  pf$calc_stock_start(1L)
  pf$calc_counts_dth(1L)
  pf$calc_rates(1L)
  pf$calc_exposure_approx1()
  set.seed(0)
  pf$draw_stock_end_net_mig(i_interval = 1L, is_forecast = TRUE)
  ans_obtained_stock_end <- pf$stock_end
  ans_obtained_net_mig <- pf$net_mig
  ans_obtained_logimp <- pf$logimp_stock_end_net_mig
  set.seed(0)
  lower <- pf$counts_dth - pf$stock_start
  ans_expected_net_mig <- matrix(rskeltr(
    n = pf$n_particle * pf$n_region,
    mu1 = pf$rates_cin,
    mu2 = pf$rates_cout * pf$exposure_approx1,
    lower = lower
  ),
  nrow = pf$n_particle
  )
  ans_expected_stock_end <- ans_expected_net_mig + pf$stock_start - pf$counts_dth
  ans_expected_logimp <- rowSums(matrix(dskeltr(
    x = ans_expected_net_mig,
    mu1 = pf$rates_cin,
    mu2 = pf$rates_cout * pf$exposure_approx1,
    lower = lower,
    use_log = TRUE
  ),
  nrow = pf$n_particle
  ))
  expect_identical(ans_obtained_stock_end, ans_expected_stock_end)
  expect_identical(ans_obtained_net_mig, ans_expected_net_mig)
  expect_identical(ans_obtained_logimp, ans_expected_logimp)
})

test_that("'draw_stock_end_net_mig' gives expected answer - with regions, have partial stock data", {
  set.seed(0)
  pf <- simulated_pf_withreg(has_stock_init = TRUE)
  counts_data_stock <- as.numeric(stats::rpois(
    n = pf$n_region * (pf$n_interval + 1L),
    lambda = 20
  ))
  counts_data_stock <- matrix(counts_data_stock,
    nrow = pf$n_region
  )
  counts_data_stock[1, 2] <- NA
  cdm_stock <- new_CdmWithregPoibin(
    counts_data = counts_data_stock,
    prob = 0.95
  )
  pf$cdms_stock <- new_CdmsWithreg(list(cdm_stock))
  pf$index_parent[, 1L] <- sample(pf$n_particle)
  pf$calc_stock_start(1L)
  pf$calc_counts_dth(1L)
  pf$calc_rates(1L)
  pf$calc_exposure_approx1()
  set.seed(0)
  pf$draw_stock_end_net_mig(i_interval = 1L, is_forecast = FALSE)
  ans_obtained_stock_end <- pf$stock_end
  ans_obtained_net_mig <- pf$net_mig
  ans_obtained_logimp <- pf$logimp_stock_end_net_mig
  set.seed(0)
  ans_expected_stock_end <- pf$cdms_stock$draw_counts_true(
    i_interval = 1L,
    n_particle = pf$n_particle,
    n_region = pf$n_region
  )
  ans_expected_net_mig <- ans_expected_stock_end - pf$stock_start + pf$counts_dth
  ans_expected_net_mig[, 1] <- rskeltr(
    n = pf$n_particle,
    mu1 = pf$rates_cin[, 1],
    mu2 = pf$rates_cout[, 1] * pf$exposure_approx1[, 1],
    lower = pf$counts_dth[, 1] - pf$stock_start[, 1]
  )
  ans_expected_stock_end[, 1] <- pf$stock_start[, 1] - pf$counts_dth[, 1] + ans_expected_net_mig[, 1]
  ans_expected_logimp <- pf$cdms_stock$calc_logimp(ans_expected_stock_end, i_interval = 1L)
  ans_expected_logimp[, 1] <- dskeltr(ans_expected_net_mig[, 1],
    mu1 = pf$rates_cin[, 1],
    mu2 = pf$rates_cout[, 1] * pf$exposure_approx1[, 1],
    lower = pf$counts_dth[, 1] - pf$stock_start[, 1],
    use_log = TRUE
  )
  ans_expected_logimp <- rowSums(ans_expected_logimp)
  expect_identical(ans_obtained_stock_end, ans_expected_stock_end)
  expect_identical(ans_obtained_net_mig, ans_expected_net_mig)
  expect_identical(ans_obtained_logimp, ans_expected_logimp)
})


## draw_stock_end_net_mig_parallelogram -----------------------------------------------

test_that("'draw_stock_end_net_mig_parallogram' gives expected answer", {
  set.seed(0)
  pf <- simulated_pf_noreg(has_stock_init = TRUE, parallelogram = TRUE)
  pf$index_parent[, 1L] <- sample(pf$n_particle)
  pf$calc_stock_start(1L)
  pf$calc_counts_dth(1L)
  set.seed(0)
  pf$draw_stock_end_net_mig_parallelogram(i_interval = 1L)
  ans_obtained_stock_end_second <- pf$stock_end_second
  ans_obtained_net_mig_combined <- pf$net_mig_combined
  ans_obtained_logimp <- pf$logimp_stock_end_net_mig
  set.seed(0)
  ans_expected_stock_end_second <- pf$cdms_stock$draw_counts_true(
    i_interval = 2L,
    n_particle = pf$n_particle
  )
  ans_expected_net_mig_combined <- ans_expected_stock_end_second - pf$stock_start +
      pf$counts_dth + pf$counts_dth_second
  ans_expected_logimp <- pf$cdms_stock$calc_logimp(
    counts_true = ans_expected_stock_end_second,
    i_interval = 2L
  )
  expect_identical(ans_obtained_stock_end_second, ans_expected_stock_end_second)
  expect_identical(ans_obtained_net_mig_combined, ans_expected_net_mig_combined)
  expect_identical(ans_obtained_logimp, ans_expected_logimp)
})


## draw_values ----------------------------------------------------------------

test_that("'draw_values' gives expected answer - no regions, not forecast", {
  set.seed(0)
  pf <- simulated_pf_noreg(has_stock_init = TRUE)
  pf$skip_sampling_init()
  pf$draw_values(i_interval = 1L, is_forecast = FALSE)
  expect_true(!anyNA(pf$counts_stock[, 2L]))
  expect_true(!anyNA(pf$counts_births[, 1L]))
  expect_true(!anyNA(pf$counts_deaths[, 1L]))
  expect_true(!anyNA(pf$counts_immigration1[, 1L]))
  expect_true(!anyNA(pf$counts_emigration1[, 1L]))
  expect_true(!anyNA(pf$counts_immigration2[, 1L]))
  expect_true(!anyNA(pf$counts_emigration2[, 1L]))
  expect_equal(
    pf$counts_stock[, 2L] - pf$counts_stock[, 1L],
    pf$counts_immigration1[, 1L] + pf$counts_immigration2[, 1L] -
      (pf$counts_deaths[, 1L] + pf$counts_emigration1[, 1L] + pf$counts_emigration2[, 1L])
  )
})

test_that("'draw_values' gives expected answer - parallelogram", {
  set.seed(0)
  pf <- simulated_pf_noreg(has_stock_init = TRUE, parallelogram = TRUE, has_one_imem = TRUE)
  pf$skip_sampling_init()
  pf$draw_values(i_interval = 1L, is_forecast = FALSE)
  expect_true(!anyNA(pf$counts_stock[, 2L]))
  expect_true(!anyNA(pf$counts_births[, 1:2]))
  expect_true(!anyNA(pf$counts_deaths[, 1:2]))
  expect_true(!anyNA(pf$counts_immigration1[, 1:2]))
  expect_true(!anyNA(pf$counts_emigration1[, 1:2]))
  expect_equal(
    pf$counts_stock[, 2L] - pf$counts_stock[, 1L],
    pf$counts_immigration1[, 1L] + pf$counts_immigration2[, 1L] -
      (pf$counts_deaths[, 1L] + pf$counts_emigration1[, 1L] + pf$counts_emigration2[, 1L])
  )
  expect_equal(
    pf$counts_stock[, 3L] - pf$counts_stock[, 2L],
    pf$counts_immigration1[, 2L] - pf$counts_deaths[, 2L] - pf$counts_emigration1[, 2L]
  )
})

test_that("'draw_values' gives expected answer - no regions, is forecast", {
  set.seed(0)
  pf <- simulated_pf_noreg(has_stock_init = TRUE)
  pf$skip_sampling_init()
  pf$draw_values(i_interval = 1L, is_forecast = TRUE)
  expect_true(!anyNA(pf$counts_stock[, 2L]))
  expect_true(!anyNA(pf$counts_births[, 1L]))
  expect_true(!anyNA(pf$counts_deaths[, 1L]))
  expect_true(!anyNA(pf$counts_immigration1[, 1L]))
  expect_true(!anyNA(pf$counts_emigration1[, 1L]))
  expect_true(!anyNA(pf$counts_immigration2[, 1L]))
  expect_true(!anyNA(pf$counts_emigration2[, 1L]))
  expect_equal(
    pf$counts_stock[, 2L] - pf$counts_stock[, 1L],
    pf$counts_immigration1[, 1L] + pf$counts_immigration2[, 1L] -
      (pf$counts_deaths[, 1L] + pf$counts_emigration1[, 1L] + pf$counts_emigration2[, 1L])
  )
})

test_that("'draw_values' gives expected answer - with regions, not forecast", {
  set.seed(0)
  pf <- simulated_pf_withreg(has_stock_init = TRUE)
  pf$skip_sampling_init()
  pf$draw_values(i_interval = 1L, is_forecast = FALSE)
  expect_true(!anyNA(pf$counts_stock[, , 2L]))
  expect_true(!anyNA(pf$counts_births[, , 1L]))
  expect_true(!anyNA(pf$counts_deaths[, , 1L]))
  expect_true(!anyNA(pf$counts_internal_in[, , 1L]))
  expect_true(!anyNA(pf$counts_internal_out[, , 1L]))
  expect_true(!anyNA(pf$counts_immigration1[, , 1L]))
  expect_true(!anyNA(pf$counts_emigration1[, , 1L]))
  expect_true(!anyNA(pf$counts_immigration2[, , 1L]))
  expect_true(!anyNA(pf$counts_emigration2[, , 1L]))
  expect_equal(
    pf$counts_stock[, , 2L] - pf$counts_stock[, , 1L],
    pf$counts_internal_in[, , 1L] +
      pf$counts_immigration1[, , 1L] +
      pf$counts_immigration2[, , 1L] -
      (pf$counts_internal_out[, , 1L]
      + pf$counts_deaths[, , 1L]
        + pf$counts_emigration1[, , 1L]
        + pf$counts_emigration2[, , 1L])
  )
})

test_that("'draw_values' gives expected answer - with regions, is forecast", {
  set.seed(0)
  pf <- simulated_pf_withreg(has_stock_init = TRUE)
  pf$skip_sampling_init()
  pf$draw_values(i_interval = 1L, is_forecast = TRUE)
  expect_true(!anyNA(pf$counts_stock[, , 2L]))
  expect_true(!anyNA(pf$counts_births[, , 1L]))
  expect_true(!anyNA(pf$counts_deaths[, , 1L]))
  expect_true(!anyNA(pf$counts_internal_in[, , 1L]))
  expect_true(!anyNA(pf$counts_internal_out[, , 1L]))
  expect_true(!anyNA(pf$counts_immigration1[, , 1L]))
  expect_true(!anyNA(pf$counts_emigration1[, , 1L]))
  expect_true(!anyNA(pf$counts_immigration2[, , 1L]))
  expect_true(!anyNA(pf$counts_emigration2[, , 1L]))
  expect_equal(
    pf$counts_stock[, , 2L] - pf$counts_stock[, , 1L],
    pf$counts_internal_in[, , 1L] +
      pf$counts_immigration1[, , 1L] +
      pf$counts_immigration2[, , 1L] -
      (pf$counts_internal_out[, , 1L]
      + pf$counts_deaths[, , 1L]
        + pf$counts_emigration1[, , 1L]
        + pf$counts_emigration2[, , 1L])
  )
})


## draw_values_init -----------------------------------------------------------

test_that("'draw_values_init' gives expected answer - no regions", {
  set.seed(0)
  pf <- simulated_pf_noreg()
  pf$draw_values_init()
  expect_true(!any(is.na(pf$stock_end)))
  expect_identical(pf$stock_end, pf$counts_stock[, 1L])
})

test_that("'draw_values_init' gives expected answer - with regions", {
  set.seed(0)
  pf <- simulated_pf_withreg()
  pf$draw_values_init()
  expect_true(!any(is.na(pf$stock_end)))
  expect_identical(pf$stock_end, pf$counts_stock[, , 1L])
})


## forecast_cohort -----------------------------------------------------------

test_that("'forecast_cohort' gives expected answer", {
    for (seed in seq_len(5)) {
        set.seed(seed)
        pf <- simulated_pf_noreg(has_stock_init = TRUE)
        pf$forecast_cohort()
        expect_equal(pf$counts_stock[,-1], pf$counts_stock[,-(pf$n_interval+1)] - pf$counts_deaths -
                                           pf$counts_emigration1 + pf$counts_immigration1)
    }
})



## make_diagnostics -----------------------------------------------------------

test_that("'make_diagnostics' gives expected answer - existing cohort", {
  pf <- simulated_pf_noreg()
  pf$ess <- seq_len(pf$n_interval + 1L)
  pf$resampled <- rep(c(TRUE, FALSE), length.out = pf$n_interval + 1L)
  pf$n_unique <- seq_len(pf$n_interval + 1L)
  ans_obtained <- pf$make_diagnostics()
  ans_expected <- data.frame(
    cohort = rep(pf$cohort, times = pf$n_interval + 1L),
    sexgender = rep(pf$sexgender, times = pf$n_interval + 1L),
    time = pf$time_levels_stock,
    age = c(pf$age_levels_events[[1]], pf$age_levels_events),
    triangle = c("<none>", rep(c("Upper", "Lower"), length.out = pf$n_interval)),
    sum_loglik = pf$sum_loglik,
    sum_logtrans = c(0, pf$sum_logtrans),
    sum_logimp = pf$sum_logimp,
    sum_logwt_unnorm = pf$sum_logwt_unnorm,
    ess = pf$ess,
    resampled = pf$resampled,
    n_particle = pf$n_particle,
    n_thin = pf$n_thin,
    n_unique = pf$n_unique
  )
  expect_identical(ans_obtained, ans_expected)
})

test_that("'make_diagnostics' gives expected answer - new cohort", {
  pf <- simulated_pf_noreg(cohort_type = "new")
  pf$ess <- seq_len(pf$n_interval + 1L)
  pf$resampled <- rep(c(TRUE, FALSE), length.out = pf$n_interval + 1L)
  pf$n_unique <- seq_len(pf$n_interval + 1L)
  ans_obtained <- pf$make_diagnostics()
  ans_expected <- data.frame(
    cohort = rep(pf$cohort, times = pf$n_interval + 1L),
    sexgender = rep(pf$sexgender, times = pf$n_interval + 1L),
    time = pf$time_levels_stock,
    age = c(-1L, pf$age_levels_events),
    triangle = c("<none>", rep(c("Lower", "Upper"), length.out = pf$n_interval)),
    sum_loglik = pf$sum_loglik,
    sum_logtrans = c(0, pf$sum_logtrans),
    sum_logimp = pf$sum_logimp,
    sum_logwt_unnorm = pf$sum_logwt_unnorm,
    ess = pf$ess,
    resampled = pf$resampled,
    n_particle = pf$n_particle,
    n_thin = pf$n_thin,
    n_unique = pf$n_unique
  )
  expect_identical(ans_obtained, ans_expected)
})



## make_output_events ---------------------------------------------------------

test_that("'make_output_events' gives expected answer - no regions, is_dominant and is_forecast both TRUE", {
  set.seed(0)
  n_particle <- 10L
  n_thin <- 2L
  pf <- simulated_pf_noreg(n_particle = n_particle, n_thin = n_thin)
  pf$counts_births[] <- as.numeric(seq_len(pf$n_particle * pf$n_interval))
  pf$counts_deaths[] <- as.numeric(seq_len(pf$n_particle * pf$n_interval))
  pf$counts_immigration1[] <- as.numeric(seq_len(pf$n_particle * pf$n_interval))
  pf$counts_emigration1[] <- as.numeric(seq_len(pf$n_particle * pf$n_interval))
  pf$counts_immigration2[] <- as.numeric(seq_len(pf$n_particle * pf$n_interval))
  pf$counts_emigration2[] <- as.numeric(seq_len(pf$n_particle * pf$n_interval))
  pf$index_ancestor[] <- 1:5
  pf$draw_index_output()
  ans_obtained <- pf$make_output_events(is_forecast = TRUE)
  meta <- data.frame(
    cohort = rep(pf$cohort, times = pf$n_interval),
    sexgender = rep(pf$sexgender, times = pf$n_interval),
    time = pf$time_levels_events,
    age = pf$age_levels_events
  )
  meta <- serialize(meta, connection = NULL)
  n_meta <- length(meta)
  n_count <- length(pf$counts_deaths[pf$index_output, ])
  ans_expected <- list(
    births = list(n_meta, n_count, meta,
                  as.numeric(pf$counts_births[pf$index_output, ])
    ),
    deaths = list(n_meta, n_count, meta,
      as.numeric(pf$counts_deaths[pf$index_output, ])
    ),
    immigration1 = list(n_meta, n_count, meta,
      as.numeric(pf$counts_immigration1[pf$index_output, ])
    ),
    emigration1 = list(n_meta, n_count, meta,
      as.numeric(pf$counts_emigration1[pf$index_output, ])
    ),
    immigration2 = list(n_meta, n_count, meta,
      as.numeric(pf$counts_immigration2[pf$index_output, ])
    ),
    emigration2 = list(n_meta, n_count, meta,
      as.numeric(pf$counts_emigration2[pf$index_output, ])
    )
  )
  expect_identical(ans_obtained, ans_expected)
})

test_that("'make_output_events' gives expected answer - no regions, is_dominant is FALSE", {
  set.seed(0)
  n_particle <- 5L
  pf <- simulated_pf_noreg(n_particle = n_particle, is_dominant = FALSE)
  pf$counts_deaths[] <- as.numeric(seq_len(pf$n_particle * pf$n_interval))
  pf$counts_immigration1[] <- as.numeric(seq_len(pf$n_particle * pf$n_interval))
  pf$counts_emigration1[] <- as.numeric(seq_len(pf$n_particle * pf$n_interval))
  pf$counts_immigration2[] <- as.numeric(seq_len(pf$n_particle * pf$n_interval))
  pf$counts_emigration2[] <- as.numeric(seq_len(pf$n_particle * pf$n_interval))
  pf$index_ancestor[] <- 1:5
  pf$draw_index_output()
  ans_obtained <- pf$make_output_events(is_forecast = TRUE)
  meta <- data.frame(
    cohort = rep(pf$cohort, times = pf$n_interval),
    sexgender = rep(pf$sexgender, times = pf$n_interval),
    time = pf$time_levels_events, 
    age = pf$age_levels_events
  )
  meta <- serialize(meta, connection = NULL)
  n_meta <- length(meta)
  n_count <- length(pf$counts_deaths[pf$index_output, ])
  ans_expected <- list(
    deaths = list(n_meta, n_count, meta,
      as.numeric(pf$counts_deaths[pf$index_output, ])
    ),
    immigration1 = list(n_meta, n_count, meta,
      as.numeric(pf$counts_immigration1[pf$index_output, ])
    ),
    emigration1 = list(n_meta, n_count, meta,
      as.numeric(pf$counts_emigration1[pf$index_output, ])
    ),
    immigration2 = list(n_meta, n_count, meta,
      as.numeric(pf$counts_immigration2[pf$index_output, ])
    ),
    emigration2 = list(n_meta, n_count, meta,
      as.numeric(pf$counts_emigration2[pf$index_output, ])
    )
  )
  expect_identical(ans_obtained, ans_expected)
})

test_that("'make_output_events' gives expected answer - with regions, is_dominant, is_forecast both TRUE", {
  set.seed(0)
  n_particle <- 10L
  n_thin <- 2L
  pf <- simulated_pf_withreg(n_particle = n_particle, n_thin = n_thin)
  pf$counts_births[] <- as.numeric(seq_len((pf$n_particle / pf$n_thin) * pf$n_region * pf$n_interval))
  pf$counts_deaths[] <- as.numeric(seq_len((pf$n_particle / pf$n_thin) * pf$n_region * pf$n_interval))
  pf$counts_internal_in[] <- as.numeric(seq_len((pf$n_particle / pf$n_thin) * pf$n_region * pf$n_interval))
  pf$counts_internal_out[] <- as.numeric(seq_len((pf$n_particle / pf$n_thin) * pf$n_region * pf$n_interval))
  pf$counts_immigration1[] <- as.numeric(seq_len((pf$n_particle / pf$n_thin) * pf$n_region * pf$n_interval))
  pf$counts_emigration1[] <- as.numeric(seq_len((pf$n_particle / pf$n_thin) * pf$n_region * pf$n_interval))
  pf$counts_immigration2[] <- as.numeric(seq_len((pf$n_particle / pf$n_thin) * pf$n_region * pf$n_interval))
  pf$counts_emigration2[] <- as.numeric(seq_len((pf$n_particle / pf$n_thin) * pf$n_region * pf$n_interval))
  pf$index_ancestor[] <- 1:5
  pf$draw_index_output()
  ans_obtained <- pf$make_output_events(is_forecast = TRUE)
  meta <- data.frame(
    cohort = rep(pf$cohort, times = pf$n_region * pf$n_interval),
    sexgender = rep(pf$sexgender, times = pf$n_region * pf$n_interval),
    time = rep(pf$time_levels_events, each = pf$n_region),
    age = rep(pf$age_levels_events, each = pf$n_region),
    region = rep(pf$region_levels, times = pf$n_interval)
  )
  meta <- serialize(meta, connection = NULL)
  n_meta <- length(meta)
  n_count <- length(pf$counts_births[pf$index_output , , ])
  ans_expected <- list(
    births = list(n_meta, n_count, meta, 
      as.numeric(pf$counts_births[pf$index_output , , ])
    ),
    deaths = list(n_meta, n_count, meta, 
      as.numeric(pf$counts_deaths[pf$index_output , , ])
    ),
    internal_in = list(n_meta, n_count, meta, 
      as.numeric(pf$counts_internal_in[pf$index_output , , ])
    ),
    internal_out = list(n_meta, n_count, meta, 
      as.numeric(pf$counts_internal_out[pf$index_output , , ])
    ),
    immigration1 = list(n_meta, n_count, meta, 
      as.numeric(pf$counts_immigration1[pf$index_output , , ])
    ),
    emigration1 = list(n_meta, n_count, meta, 
      as.numeric(pf$counts_emigration1[pf$index_output , , ])
    ),
    immigration2 = list(n_meta, n_count, meta, 
      as.numeric(pf$counts_immigration2[pf$index_output , , ])
    ),
    emigration2 = list(n_meta, n_count, meta, 
      as.numeric(pf$counts_emigration2[pf$index_output , , ])
    )
  )
  expect_identical(ans_obtained, ans_expected)
})

test_that("'make_output_events' gives expected answer - with regions, is_dominant is FALSE", {
  set.seed(0)
  n_particle <- 10L
  n_thin <- 2L
  pf <- simulated_pf_withreg(n_particle = n_particle, n_thin = n_thin, is_dominant = FALSE)
  pf$counts_deaths[] <- as.numeric(seq_len((pf$n_particle / pf$n_thin) * pf$n_region * pf$n_interval))
  pf$counts_internal_in[] <- as.numeric(seq_len((pf$n_particle / pf$n_thin) * pf$n_region * pf$n_interval))
  pf$counts_internal_out[] <- as.numeric(seq_len((pf$n_particle / pf$n_thin) * pf$n_region * pf$n_interval))
  pf$counts_immigration1[] <- as.numeric(seq_len((pf$n_particle / pf$n_thin) * pf$n_region * pf$n_interval))
  pf$counts_emigration1[] <- as.numeric(seq_len((pf$n_particle / pf$n_thin) * pf$n_region * pf$n_interval))
  pf$counts_immigration2[] <- as.numeric(seq_len((pf$n_particle / pf$n_thin) * pf$n_region * pf$n_interval))
  pf$counts_emigration2[] <- as.numeric(seq_len((pf$n_particle / pf$n_thin) * pf$n_region * pf$n_interval))
  pf$index_ancestor[] <- 1:5
  pf$draw_index_output()
  ans_obtained <- pf$make_output_events(is_forecast = TRUE)
  meta <- data.frame(
    cohort = rep(pf$cohort, times = pf$n_region * pf$n_interval),
    sexgender = rep(pf$sexgender, times = pf$n_region * pf$n_interval),
    time = rep(pf$time_levels_events, each = pf$n_region),
    age = rep(pf$age_levels_events, each = pf$n_region),
    region = rep(pf$region_levels, times = pf$n_interval)
  )
  meta <- serialize(meta, connection = NULL)
  n_meta <- length(meta)
  n_count <- length(pf$counts_deaths[pf$index_output , , ])
  ans_expected <- list(
    deaths = list(n_meta, n_count, meta,
      as.numeric(pf$counts_deaths[pf$index_output, ,])
    ),
    internal_in = list(n_meta, n_count, meta,
      as.numeric(pf$counts_internal_in[pf$index_output, ,])
    ),
    internal_out = list(n_meta, n_count, meta,
      as.numeric(pf$counts_internal_out[pf$index_output, ,])
    ),
    immigration1 = list(n_meta, n_count, meta,
      as.numeric(pf$counts_immigration1[pf$index_output, ,])
    ),
    emigration1 = list(n_meta, n_count, meta,
      as.numeric(pf$counts_emigration1[pf$index_output, ,])
    ),
    immigration2 = list(n_meta, n_count, meta,
      as.numeric(pf$counts_immigration2[pf$index_output, ,])
    ),
    emigration2 = list(n_meta, n_count, meta,
      as.numeric(pf$counts_emigration2[pf$index_output, ,])
    )
  )
  expect_identical(ans_obtained, ans_expected)
})


## make_output_population --------------------------------------------------

test_that("'make_output_population' gives expected answer - no regions", {
  set.seed(0)
  n_particle <- 10L
  n_thin <- 2L
  pf <- simulated_pf_noreg(n_particle = n_particle, n_thin = n_thin)
  pf$index_ancestor[] <- 1:5
  pf$draw_index_output()
  pf$counts_stock[] <- seq_len(pf$n_particle * (pf$n_interval + 1))
  ans_obtained <- pf$make_output_population()
  meta <- data.frame(
    cohort = rep(pf$cohort, times = pf$n_popn),
    sexgender = rep(pf$sexgender, times = pf$n_popn),
    time = pf$time_levels_stock[pf$is_popn],
    age = pf$age_levels_stock[pf$is_popn]
  )
  meta <- serialize(meta, connection = NULL)
  count = as.numeric(pf$counts_stock[pf$index_output, pf$is_popn])
  ans_expected <- list(length(meta), length(count), meta, count)
  expect_identical(ans_obtained, ans_expected)
})

test_that("'make_output_population' gives expected answer - with regions", {
  set.seed(0)
  n_particle <- 10L
  n_thin <- 2L
  pf <- simulated_pf_withreg(n_particle = n_particle, n_thin = n_thin)
  pf$counts_stock[] <- seq_len(pf$n_particle * pf$n_region * (pf$n_interval + 1))
  pf$index_ancestor[] <- 1:5
  pf$draw_index_output()
  ans_obtained <- pf$make_output_population()
  meta <- data.frame(
    cohort = rep(pf$cohort, times = pf$n_popn * pf$n_region),
    sexgender = rep(pf$sexgender, times = pf$n_popn * pf$n_region),
    time = rep(pf$time_levels_stock[pf$is_popn], each = pf$n_region),
    age = rep(pf$age_levels_stock[pf$is_popn], each = pf$n_region),
    region = rep(pf$region_levels, times = pf$n_popn)
  )
  count <- as.double(pf$counts_stock[pf$index_output, , pf$is_popn])
  meta <- serialize(meta, connection = NULL)
  ans_expected <- list(length(meta), length(count), meta, count)
  expect_identical(ans_obtained, ans_expected)
})


## resample -------------------------------------------------------------------

test_that("'resample' gives expected answer - do resampling", {
  set.seed(0)
  pf <- simulated_pf_noreg()
  pf$threshold <- 0.9
  i_interval <- 2L
  pf$logwt_unnorm[, i_interval + 1L] <- log(abs(rt(n = pf$n_particle, df = 1)))
  pf$resample(i_interval)
  expect_true(pf$resampled[i_interval + 1L])
  expect_equal(pf$ess[i_interval + 1L], 1 / sum(softmax(pf$logwt_unnorm[, i_interval + 1])^2))
})

test_that("'resample' gives expected answer - do not do resampling because ess low", {
  set.seed(0)
  pf <- simulated_pf_noreg()
  i_interval <- 2L
  pf$logwt_unnorm[, i_interval + 1L] <- rep(1, times = pf$n_particle)
  pf$resample(i_interval)
  expect_false(pf$resampled[i_interval + 1L])
  expect_equal(pf$ess[i_interval + 1L], 1 / sum(softmax(pf$logwt_unnorm[, i_interval + 1L])^2))
  expect_equal(pf$index_parent[, i_interval + 1L], seq_len(pf$n_particle))
})


## resample_init --------------------------------------------------------------

test_that("'resample_init' gives expected answer - do resampling", {
  set.seed(0)
  pf <- simulated_pf_noreg()
  pf$threshold <- 0.9
  pf$logwt_unnorm[, 1L] <- log(abs(rt(n = pf$n_particle, df = 1)))
  pf$resample_init()
  expect_true(pf$resampled[1L])
  expect_equal(pf$ess[1L], 1 / sum(softmax(pf$logwt_unnorm[, 1])^2))
})

test_that("'resample_init' gives expected answer - do not do resampling", {
  set.seed(0)
  pf <- simulated_pf_noreg()
  pf$logwt_unnorm[, 1L] <- rep(1, times = pf$n_particle)
  pf$resample_init()
  expect_false(pf$resampled[1L])
  expect_equal(pf$ess[1L], 1 / sum(softmax(pf$logwt_unnorm[, 1])^2))
  expect_equal(pf$index_parent[, 1L], seq_len(pf$n_particle))
})


## run ------------------------------------------------------------------------

test_that("'run' gives valid results with no regions, not forecast", {
  set.seed(0)
  is_forecast <- FALSE
  pf <- simulated_pf_noreg(is_forecast = is_forecast)
  pf$run(is_forecast = is_forecast)
  expect_true(!anyNA(pf$counts_stock))
  expect_true(!anyNA(pf$counts_births))
  expect_true(!anyNA(pf$counts_deaths))
  expect_true(!anyNA(pf$counts_immigration1))
  expect_true(!anyNA(pf$counts_emigration1))
  expect_true(!anyNA(pf$counts_immigration2))
  expect_true(!anyNA(pf$counts_emigration2))
  expect_equal(
    pf$counts_stock[, -1] - pf$counts_stock[, -(pf$n_interval + 1)],
    pf$counts_immigration1 + pf$counts_immigration2 -
      pf$counts_deaths - pf$counts_emigration1 - pf$counts_emigration2
  )
  expect_true(!anyNA(pf$ess))
  expect_true(!anyNA(pf$resampled))
  expect_true(!anyNA(pf$n_unique))
})

test_that("'run' gives valid results with no regions, is forecast", {
  set.seed(0)
  is_forecast <- TRUE
  pf <- simulated_pf_noreg(is_forecast = is_forecast)
  pf$run(is_forecast = is_forecast)
  expect_true(!anyNA(pf$counts_stock))
  expect_true(!anyNA(pf$counts_births))
  expect_true(!anyNA(pf$counts_deaths))
  expect_true(!anyNA(pf$counts_immigration1))
  expect_true(!anyNA(pf$counts_emigration1))
  expect_true(!anyNA(pf$counts_immigration2))
  expect_true(!anyNA(pf$counts_emigration2))
  expect_equal(
    pf$counts_stock[, -1] - pf$counts_stock[, -(pf$n_interval + 1)],
    pf$counts_immigration1 + pf$counts_immigration2 -
      pf$counts_deaths - pf$counts_emigration1 - pf$counts_emigration2
  )
  expect_true(!anyNA(pf$ess))
  expect_true(!anyNA(pf$resampled))
  expect_true(!anyNA(pf$n_unique))
})


## skip_sampling_init ---------------------------------------------------------

test_that("'skip_sampling_init' gives expected answer", {
  set.seed(0)
  pf <- simulated_pf_noreg(has_stock_init = TRUE)
  pf$skip_sampling_init()
  expect_identical(pf$logwt_unnorm[, 1L], rep(log(1/pf$n_particle), times = pf$n_particle))
  expect_identical(pf$index_parent[, 1L], seq_len(pf$n_particle))
  expect_false(pf$resampled[1L])
  expect_identical(pf$ess[1L], as.numeric(pf$n_particle))
})


## update_counts --------------------------------------------------------------

test_that("'update_counts' gives expected answer with no regions", {
  set.seed(0)
  pf <- simulated_pf_noreg()
  pf$stock_end <- as.double(rpois(n = pf$n_particle, lambda = 10))
  pf$counts_bth <- as.double(rpois(n = pf$n_particle, lambda = 10))
  pf$counts_dth <- as.double(rpois(n = pf$n_particle, lambda = 10))
  pf$counts_im1 <- as.double(rpois(n = pf$n_particle, lambda = 10))
  pf$counts_em1 <- as.double(rpois(n = pf$n_particle, lambda = 10))
  pf$counts_im2 <- as.double(rpois(n = pf$n_particle, lambda = 10))
  pf$counts_em2 <- as.double(rpois(n = pf$n_particle, lambda = 10))
  i_interval <- 2L
  pf$update_counts(i_interval = i_interval)
  expect_identical(pf$counts_stock[, i_interval + 1L], pf$stock_end)
  expect_identical(pf$counts_births[, i_interval], pf$counts_bth)
  expect_identical(pf$counts_deaths[, i_interval], pf$counts_dth)
  expect_identical(pf$counts_immigration1[, i_interval], pf$counts_im1)
  expect_identical(pf$counts_emigration1[, i_interval], pf$counts_em1)
  expect_identical(pf$counts_immigration2[, i_interval], pf$counts_im2)
  expect_identical(pf$counts_emigration2[, i_interval], pf$counts_em2)
})

test_that("'update_counts' gives expected answer with parallelogram", {
  set.seed(0)
  pf <- simulated_pf_noreg(parallelogram = TRUE)
  pf$stock_end <- as.double(rpois(n = pf$n_particle, lambda = 10))
  pf$counts_bth <- as.double(rpois(n = pf$n_particle, lambda = 10))
  pf$counts_dth <- as.double(rpois(n = pf$n_particle, lambda = 10))
  pf$counts_im1 <- as.double(rpois(n = pf$n_particle, lambda = 10))
  pf$counts_em1 <- as.double(rpois(n = pf$n_particle, lambda = 10))
  pf$counts_im2 <- rep(0, times = pf$n_particle)
  pf$counts_em2 <- rep(0, times = pf$n_particle)
  pf$stock_end_second <- as.double(rpois(n = pf$n_particle, lambda = 10))
  pf$counts_bth_second <- as.double(rpois(n = pf$n_particle, lambda = 10))
  pf$counts_dth_second <- as.double(rpois(n = pf$n_particle, lambda = 10))
  pf$counts_im1_second <- as.double(rpois(n = pf$n_particle, lambda = 10))
  pf$counts_em1_second <- as.double(rpois(n = pf$n_particle, lambda = 10))
  i_interval <- 1L
  pf$update_counts(i_interval = i_interval)
  expect_identical(pf$counts_stock[, i_interval + 1L], pf$stock_end)
  expect_identical(pf$counts_births[, i_interval], pf$counts_bth)
  expect_identical(pf$counts_deaths[, i_interval], pf$counts_dth)
  expect_identical(pf$counts_immigration1[, i_interval], pf$counts_im1)
  expect_identical(pf$counts_emigration1[, i_interval], pf$counts_em1)
  expect_identical(pf$counts_immigration2[, i_interval], pf$counts_im2)
  expect_identical(pf$counts_emigration2[, i_interval], pf$counts_em2)
  expect_identical(pf$counts_stock[, i_interval + 2L], pf$stock_end_second)
  expect_identical(pf$counts_births[, i_interval + 1L], pf$counts_bth_second)
  expect_identical(pf$counts_deaths[, i_interval + 1L], pf$counts_dth_second)
  expect_identical(pf$counts_immigration1[, i_interval + 1L], pf$counts_im1_second)
  expect_identical(pf$counts_emigration1[, i_interval + 1L], pf$counts_em1_second)
})

test_that("'update_counts' gives expected answer with regions", {
  set.seed(0)
  pf <- simulated_pf_withreg()
  pf$stock_end <- matrix(as.double(rpois(n = pf$n_particle * pf$n_region, lambda = 10)), nr = pf$n_particle)
  pf$counts_bth <- matrix(as.double(rpois(n = pf$n_particle * pf$n_region, lambda = 10)), nr = pf$n_particle)
  pf$counts_dth <- matrix(as.double(rpois(n = pf$n_particle * pf$n_region, lambda = 10)), nr = pf$n_particle)
  pf$counts_in <- matrix(as.double(rpois(n = pf$n_particle * pf$n_region, lambda = 10)), nr = pf$n_particle)
  pf$counts_out <- matrix(as.double(rpois(n = pf$n_particle * pf$n_region, lambda = 10)), nr = pf$n_particle)
  pf$counts_im1 <- matrix(as.double(rpois(n = pf$n_particle * pf$n_region, lambda = 10)), nr = pf$n_particle)
  pf$counts_em1 <- matrix(as.double(rpois(n = pf$n_particle * pf$n_region, lambda = 10)), nr = pf$n_particle)
  pf$counts_im2 <- matrix(as.double(rpois(n = pf$n_particle * pf$n_region, lambda = 10)), nr = pf$n_particle)
  pf$counts_em2 <- matrix(as.double(rpois(n = pf$n_particle * pf$n_region, lambda = 10)), nr = pf$n_particle)
  i_interval <- 2L
  pf$update_counts(i_interval = i_interval)
  expect_identical(pf$counts_stock[, , i_interval + 1L], pf$stock_end)
  expect_identical(pf$counts_births[, , i_interval], pf$counts_bth)
  expect_identical(pf$counts_deaths[, , i_interval], pf$counts_dth)
  expect_identical(pf$counts_internal_in[, , i_interval], pf$counts_in)
  expect_identical(pf$counts_internal_out[, , i_interval], pf$counts_out)
  expect_identical(pf$counts_immigration1[, , i_interval], pf$counts_im1)
  expect_identical(pf$counts_emigration1[, , i_interval], pf$counts_em1)
  expect_identical(pf$counts_immigration2[, , i_interval], pf$counts_im2)
  expect_identical(pf$counts_emigration2[, , i_interval], pf$counts_em2)
})


## update_counts_init --------------------------------------------------------------

test_that("'update_counts_init' gives expected answer with no regions", {
  set.seed(0)
  pf <- simulated_pf_noreg()
  pf$stock_end <- as.double(rpois(n = pf$n_particle, lambda = 10))
  pf$update_counts_init()
  expect_identical(pf$counts_stock[, 1L], pf$stock_end)
})

test_that("'update_counts_init' gives expected answer with regions", {
  set.seed(0)
  pf <- simulated_pf_withreg()
  pf$stock_end <- matrix(as.double(rpois(n = pf$n_particle * pf$n_region, lambda = 10)),
                         nr = pf$n_particle)
  pf$update_counts_init()
  expect_identical(pf$counts_stock[, , 1L], pf$stock_end)
})


## write_results --------------------------------------------------------------

test_that("'write_results' gives expected answer - is_forecast and is_dominant are TRUE", {
    set.seed(0)
    pf <- simulated_pf_noreg()
    pf$counts_stock[] <- as.numeric(seq_len(pf$n_particle * (pf$n_interval + 1L)))
    pf$counts_births[] <- as.numeric(seq_len(pf$n_particle * pf$n_interval))
    pf$counts_deaths[] <- as.numeric(seq_len(pf$n_particle * pf$n_interval))
    pf$counts_immigration1[] <- as.numeric(seq_len(pf$n_particle * pf$n_interval))
    pf$counts_emigration1[] <- as.numeric(seq_len(pf$n_particle * pf$n_interval))
    pf$counts_immigration2[] <- as.numeric(seq_len(pf$n_particle * pf$n_interval))
    pf$counts_emigration2[] <- as.numeric(seq_len(pf$n_particle * pf$n_interval))
    pf$draw_index_output()
    work_dir <- tempfile()
    dir.create(work_dir)
    pf$write_results(is_forecast = TRUE, work_dir = work_dir)
    expect_setequal(list.files(work_dir),
                    c(paste0("tmp-",
                             c("population", "births", "deaths",
                               "immigration1", "emigration1",
                               "immigration2", "emigration2"),
                             "-",
                             pf$cohort,
                             "-",
                             pf$sexgender,
                             ".bin"),
                      paste0("tmp-diagnostics-", pf$cohort, "-", pf$sexgender, ".rds")))
    unlink(work_dir, recursive = TRUE)
})

test_that("'write_results' gives expected answer - is_forecast is FALSE", {
    set.seed(0)
    pf <- simulated_pf_noreg()
    pf$counts_stock[] <- as.numeric(seq_len(pf$n_particle * (pf$n_interval + 1L)))
    pf$counts_births[] <- as.numeric(seq_len(pf$n_particle * pf$n_interval))
    pf$counts_deaths[] <- as.numeric(seq_len(pf$n_particle * pf$n_interval))
    pf$counts_immigration1[] <- as.numeric(seq_len(pf$n_particle * pf$n_interval))
    pf$counts_emigration1[] <- as.numeric(seq_len(pf$n_particle * pf$n_interval))
    pf$counts_immigration2[] <- as.numeric(seq_len(pf$n_particle * pf$n_interval))
    pf$counts_emigration2[] <- as.numeric(seq_len(pf$n_particle * pf$n_interval))
    pf$draw_index_output()
    work_dir <- tempfile()
    dir.create(work_dir)
    pf$write_results(is_forecast = FALSE, work_dir = work_dir)
    expect_setequal(list.files(work_dir),
                    c(paste0("tmp-",
                             c("population", "deaths",
                               "immigration1", "emigration1",
                               "immigration2", "emigration2"),
                             "-",
                             pf$cohort,
                             "-",
                             pf$sexgender,
                             ".bin"),
                      paste0("tmp-diagnostics-", pf$cohort, "-", pf$sexgender, ".rds")))
    unlink(work_dir, recursive = TRUE)
})


