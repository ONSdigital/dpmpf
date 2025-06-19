
## WARNING - BE CAREFUL WITH i_interval INDICES. METHODS FOR cdms USE
## C-STYLE 0-BASE INDICES, WHILE EVERYTHING ELSE USES R-STYLE 1-BASED
## INDICES. WE'LL NEED TO BE CAREFUL WHEN TRANSLATING FROM R TO C.


## counts_bth ------------------------------------------------------------

## Extract counts of births for interval 'i_interval'.
## (since births treated as known in estimates.)

## HAS_TESTS
calc_counts_bth_noreg <- function(i_interval) {
    self$counts_bth <- self$counts_births[, i_interval]
    if (self$is_parallelogram_first[i_interval])
        self$counts_bth_second <- self$counts_births[, i_interval + 1L]
}

PFilterNoReg$set(
  which = "public",
  name = "calc_counts_bth",
  value = calc_counts_bth_noreg
)

## HAS_TESTS
calc_counts_bth_withreg <- function(i_interval) {
  self$counts_bth <- self$counts_births[, , i_interval]
}

PFilterWithReg$set(
  which = "public",
  name = "calc_counts_bth",
  value = calc_counts_bth_withreg
)


## calc_counts_dth ------------------------------------------------------------

## Extract counts of deaths for interval 'i_interval'.
## Used for doing estimates rather than forecasts,
## (since deaths treated as known in estimates.)

## HAS_TESTS
calc_counts_dth_noreg <- function(i_interval) {
    self$counts_dth <- self$counts_deaths[, i_interval]
    if (self$is_parallelogram_first[i_interval])
        self$counts_dth_second <- self$counts_deaths[, i_interval + 1L]
}

PFilterNoReg$set(
  which = "public",
  name = "calc_counts_dth",
  value = calc_counts_dth_noreg
)

## HAS_TESTS
calc_counts_dth_withreg <- function(i_interval) {
  self$counts_dth <- self$counts_deaths[, , i_interval]
}

PFilterWithReg$set(
  which = "public",
  name = "calc_counts_dth",
  value = calc_counts_dth_withreg
)


## calc_exposure --------------------------------------------------------------

## Calculate exposure for current interval. Assumes that
## stock at end of interval is known.

## HAS_TESTS
calc_exposure_noreg <- function(i_interval) {
    stock_start <- self$stock_start ## dbl vector n_particle
    stock_end <- self$stock_end ## dbl vector n_particle
    counts_events <- (self$counts_bth +
                      self$counts_im1 +
                      self$counts_im2)
    exposure <- pmax(
        0.25 * (stock_start + stock_end),
        0.25 * (counts_events > 0)
    )
    self$exposure <- exposure
    if (self$is_parallelogram_first[[i_interval]]) {
        stock_end_second <- self$stock_end_second
        counts_events_second <- (self$counts_bth_second +
                                 self$counts_im1_second)
        exposure_second <- pmax(
            0.25 * (stock_end + stock_end_second),
            0.25 * (counts_events_second > 0)
        )
        self$exposure_second <- exposure_second
    }
}

PFilterNoReg$set(
  which = "public",
  name = "calc_exposure",
  value = calc_exposure_noreg
)


## HAS_TESTS
calc_exposure_withreg <- function(i_interval) {
    n_particle <- self$n_particle
    n_region <- self$n_region
    stock_start <- self$stock_start
    stock_end <- self$stock_end
    counts_events <- (self$counts_bth +
                      self$counts_dth +
                      self$counts_in +
                      self$counts_out +
                      self$counts_im1 +
                      self$counts_em1 +
                      self$counts_im2 +
                      self$counts_em2)
    exposure <- pmax(
        0.25 * (stock_start + stock_end),
        0.25 * (counts_events > 0)
    )
    exposure <- matrix(exposure,
                       nrow = n_particle,
                       ncol = n_region
                       )
    self$exposure <- exposure
}

PFilterWithReg$set(
  which = "public",
  name = "calc_exposure",
  value = calc_exposure_withreg
)


## calc_exposure_approx1 ------------------------------------------------------

## Calculate first version of approximate exposure
## for current interval.
## Does not assume that stock at end of interval is known.
## Tries to allow for the situation where the starting
## population is small or zero, and there is lots
## of in-migration.

## HAVE_TESTS
calc_exposure_approx1_noreg <- function() {
    stock_start <- self$stock_start
    rates_cim <- self$rates_cim
    exposure <- pmax(
        0.5 * stock_start,
        0.25 * rates_cim
    )
    self$exposure_approx1 <- exposure
}

PFilterNoReg$set(
  which = "public",
  name = "calc_exposure_approx1",
  value = calc_exposure_approx1_noreg
)


## HAVE_TESTS
calc_exposure_approx1_withreg <- function() {
  n_particle <- self$n_particle
  n_region <- self$n_region
  stock_start <- self$stock_start
  rates_cin <- self$rates_cin ## 'rates_cin' not 'rates_cim'
  exposure <- pmax(
    0.5 * stock_start,
    0.25 * rates_cin
  )
  exposure <- matrix(exposure,
    nrow = n_particle,
    ncol = n_region
  )
  self$exposure_approx1 <- exposure
}

PFilterWithReg$set(
  which = "public",
  name = "calc_exposure_approx1",
  value = calc_exposure_approx1_withreg
)


## calc_exposure_approx2 ------------------------------------------------------

## Calculate exposure for current interval once we
## know stock at end of interval but do not know
## events during interval.

## HAS_TESTS
calc_exposure_approx2_noreg <- function(i_interval) {
    stock_start <- self$stock_start ## dbl vector n_particle
    stock_end <- self$stock_end ## dbl vector n_particle
    rates_cim <- self$rates_cim
    exposure <- pmax(
        0.25 * (stock_start + stock_end),
        0.25 * rates_cim
    )
    self$exposure_approx2 <- exposure
    if (self$is_parallelogram_first[[i_interval]]) {
        stock_end_second <- self$stock_end_second
        rates_cim_second <- self$rates_im1_second
        exposure_second <- pmax(
            0.25 * (stock_end + stock_end_second),
            0.25 * rates_cim_second
        )
        self$exposure_approx2_second <- exposure_second
    }
}

PFilterNoReg$set(
  which = "public",
  name = "calc_exposure_approx2",
  value = calc_exposure_approx2_noreg
)


## HAS_TESTS
calc_exposure_approx2_withreg <- function(i_interval) {
  n_particle <- self$n_particle
  n_region <- self$n_region
  rates_cin <- self$rates_cin
  stock_start <- self$stock_start
  stock_end <- self$stock_end
  exposure <- pmax(
      0.25 * (stock_start + stock_end),
      0.25 * rates_cin
  )
  exposure <- matrix(exposure,
    nrow = n_particle,
    ncol = n_region
  )
  self$exposure_approx2 <- exposure
}

PFilterWithReg$set(
  which = "public",
  name = "calc_exposure_approx2",
  value = calc_exposure_approx2_withreg
)



## calc_exposure_approx_parallelogram -----------------------------------------

## Calculate approximate exposure used in parallelogram
## importance function. 

## HAVE_TESTS
calc_exposure_approx_parallelogram <- function() {
    stock_start <- self$stock_start
    stock_end_second <- self$stock_end_second
    counts_dth_first <- self$counts_dth
    counts_dth_second <- self$counts_dth_second
    rates_im_first <- self$rates_im1
    rates_im_second <- self$rates_im1_second
    rates_em_first <- self$rates_em1
    rates_em_second <- self$rates_em1_second
    stock_end_forward <- (stock_start
        - counts_dth_first
        + rates_im_first
        - 0.5 * rates_em_first * stock_start)
    stock_end_backward <- (stock_end_second
        + counts_dth_second
        - rates_im_second
        + 0.5 * rates_em_second * stock_end_second)
    stock_end_av <- 0.5 * stock_end_forward + 0.5 * stock_end_backward
    self$exposure_approx_first <- pmax(0.25 * stock_start + 0.25 * stock_end_av,
                                       0)
    self$exposure_approx_second <- pmax(0.25 * stock_end_av + 0.25 * stock_end_second,
                                        0)
}

PFilter$set(
  which = "public",
  name = "calc_exposure_approx_parallelogram",
  value = calc_exposure_approx_parallelogram
)

## calc_index_ancestor --------------------------------------------------------

## Obtain the index of the ancestors of
## the particles that are present at the completion
## of filtering. We do this recursively, by
## identifying the parents of each generation.
## The calculations are the same for the no-region
## and with-region models.

## HAS_TESTS
calc_index_ancestor <- function() {
  n_interval <- self$n_interval
  index_parent <- self$index_parent
  index_ancestor <- index_parent
  ## the final column of 'index_parent' and 'index_ancestor' is identical
  for (i_interval in seq.int(from = n_interval, to = 1L)) {
    index_ancestor_prev <- index_ancestor[, i_interval + 1L]
    index_ancestor[, i_interval] <- index_parent[index_ancestor_prev, i_interval]
  }
  self$index_ancestor <- index_ancestor
}

PFilter$set(
  which = "public",
  name = "calc_index_ancestor",
  value = calc_index_ancestor
)


## calc_logimp ----------------------------------------------------------------

## Calculate total importance (aka proposal)
## probabilities for proposals.

## HAS_TESTS
calc_logimp_noreg <- function(is_forecast) {
  ans <- (self$logimp_stock_end_net_mig
    + self$logimp_counts_im_em
    + self$logimp_gross_mig)
  if (is_forecast) {
    ans <- ans + self$logimp_counts_dth
  }
  ans
}

PFilterNoReg$set(
  which = "public",
  name = "calc_logimp",
  value = calc_logimp_noreg
)


## HAS_TESTS
calc_logimp_withreg <- function(is_forecast) {
  ans <- (self$logimp_stock_end_net_mig
    + self$logimp_counts_im_em
    + self$logimp_gross_mig
    + self$logimp_counts_in_out)
  if (is_forecast) {
    ans <- ans + self$logimp_counts_dth
  }
  ans
}

PFilterWithReg$set(
  which = "public",
  name = "calc_logimp",
  value = calc_logimp_withreg
)


## calc_logimp_init -----------------------------------------------------------

## Calculate total importance (aka proposal)
## probabilities for proposals. This consists
## of a single probability: the probability
## for the initial stock (ie the stock at
## the end of interval 0).

## HAS_TESTS
calc_logimp_init_noreg <- function() {
  self$logimp_stock_end_init
}

PFilterNoReg$set(
  which = "public",
  name = "calc_logimp_init",
  value = calc_logimp_init_noreg
)


## HAS_TESTS
calc_logimp_init_withreg <- function() {
  self$logimp_stock_end_init
}

PFilterWithReg$set(
  which = "public",
  name = "calc_logimp_init",
  value = calc_logimp_init_withreg
)


## calc_loglik ----------------------------------------------------------------

## Calculate log-likelihood for counts for stock and migration.
## Births and deaths are not included since they are being
## treated as known.

## HAS_TESTS
calc_loglik_noreg <- function(i_interval) {
    obs_zero <- self$obs_zero
    calc_stock <- self$cdms_stock$calc_loglik
    calc_im1 <- self$cdms_immigration1$calc_loglik
    calc_em1 <- self$cdms_emigration1$calc_loglik
    calc_im2 <- self$cdms_immigration2$calc_loglik
    calc_em2 <- self$cdms_emigration2$calc_loglik
    stock_end <- self$stock_end
    counts_im1 <- self$counts_im1
    counts_em1 <- self$counts_em1
    counts_im2 <- self$counts_im2
    counts_em2 <- self$counts_em2
    ans <- (calc_stock(
        counts_true = stock_end,
        i_interval = i_interval + 1L,
        obs_zero = obs_zero
    )
        + calc_im1(
              counts_true = counts_im1,
              i_interval = i_interval,
              obs_zero = obs_zero
          )
        + calc_em1(
              counts_true = counts_em1,
              i_interval = i_interval,
              obs_zero = obs_zero
          )
        + calc_im2(
              counts_true = counts_im2,
              i_interval = i_interval,
              obs_zero = obs_zero
          )
        + calc_em2(
              counts_true = counts_em2,
              i_interval = i_interval,
              obs_zero = obs_zero
          ))
    if (self$is_parallelogram_first[[i_interval + 1L]]) {
        stock_end_second <- self$stock_end_second
        ans <- ans +
            calc_stock(counts_true = stock_end_second,
                       i_interval = i_interval + 2L,
                       obs_zero = obs_zero)
    }
    ans
}


PFilterNoReg$set(
  which = "public",
  name = "calc_loglik",
  value = calc_loglik_noreg
)


## HAS_TESTS
calc_loglik_withreg <- function(i_interval) {
    obs_zero <- self$obs_zero
    calc_stock <- self$cdms_stock$calc_loglik
    calc_in <- self$cdms_internal_in$calc_loglik
    calc_out <- self$cdms_internal_out$calc_loglik
    calc_im1 <- self$cdms_immigration1$calc_loglik
    calc_em1 <- self$cdms_emigration1$calc_loglik
    calc_im2 <- self$cdms_immigration2$calc_loglik
    calc_em2 <- self$cdms_emigration2$calc_loglik
    stock_end <- self$stock_end
    counts_in <- self$counts_in
    counts_out <- self$counts_out
    counts_im1 <- self$counts_im1
    counts_em1 <- self$counts_em1
    counts_im2 <- self$counts_im2
    counts_em2 <- self$counts_em2
    ans <- (calc_stock(
        counts_true = stock_end,
        i_interval = i_interval + 1L,
        obs_zero = obs_zero
    )
        + calc_in(
              counts_true = counts_in,
              i_interval = i_interval,
              obs_zero = obs_zero
          )
        + calc_out(
              counts_true = counts_out,
              i_interval = i_interval,
              obs_zero = obs_zero
          )
        + calc_im1(
              counts_true = counts_im1,
              i_interval = i_interval,
              obs_zero = obs_zero
          )
        + calc_em1(
              counts_true = counts_em1,
              i_interval = i_interval,
              obs_zero = obs_zero
          )
        + calc_im2(
              counts_true = counts_im2,
              i_interval = i_interval,
              obs_zero = obs_zero
          )
        + calc_em2(
              counts_true = counts_em2,
              i_interval = i_interval,
              obs_zero = obs_zero
          ))
    ans
}

PFilterWithReg$set(
  which = "public",
  name = "calc_loglik",
  value = calc_loglik_withreg
)


## calc_loglik_init -----------------------------------------------------------

## Calculate log-likelihood of values for initial stock

## HAS_TESTS
calc_loglik_init_noreg <- function() {
    obs_zero <- self$obs_zero
    cdms <- self$cdms_stock
    stock_end <- self$stock_end
    cdms$calc_loglik(
             counts_true = stock_end,
             i_interval = 0L,
             obs_zero = obs_zero
         )
}

PFilterNoReg$set(
  which = "public",
  name = "calc_loglik_init",
  value = calc_loglik_init_noreg
)


## HAS_TESTS
calc_loglik_init_withreg <- function() {
    obs_zero <- self$obs_zero
    cdms <- self$cdms_stock
    stock_end <- self$stock_end
    cdms$calc_loglik(
             counts_true = stock_end,
             i_interval = 0L,
             obs_zero = obs_zero
         )
}

PFilterWithReg$set(
  which = "public",
  name = "calc_loglik_init",
  value = calc_loglik_init_withreg
)


## calc_logtrans --------------------------------------------------------------

## Calculate the log of the transition probability, ie p(x_k | x_{k-1}).

## HAS_TESTS
calc_logtrans_noreg <- function(i_interval) {
    ans <- logprob_trans(stock_start = self$stock_start,
                         counts_bth = self$counts_bth,
                         counts_dth = self$counts_dth,
                         counts_im = self$counts_im1,
                         counts_em = self$counts_em1,
                         rates_bth = self$rates_bth,
                         rates_dth = self$rates_dth,
                         rates_im = self$rates_im1,
                         rates_em = self$rates_em1,
                         exposure = self$exposure,
                         is_dominant = self$is_dominant)
    if (self$is_parallelogram_first[[i_interval]])
        ans <- ans + logprob_trans(stock_start = self$stock_end,
                                   counts_bth = self$counts_bth_second,
                                   counts_dth = self$counts_dth_second,
                                   counts_im = self$counts_im1_second,
                                   counts_em = self$counts_em1_second,
                                   rates_bth = self$rates_bth_second,
                                   rates_dth = self$rates_dth_second,
                                   rates_im = self$rates_im1_second,
                                   rates_em = self$rates_em1_second,
                                   exposure = self$exposure_second,
                                   is_dominant = self$is_dominant)
    ans
}

PFilterNoReg$set(
  which = "public",
  name = "calc_logtrans",
  value = calc_logtrans_noreg
  )


## HAS_TESTS
calc_logtrans_withreg <- function(i_interval) {
    n_particle <- self$n_particle
    n_region <- self$n_region
    stock_start <- self$stock_start
    counts_bth <- self$counts_bth
    counts_dth <- self$counts_dth
    counts_im1 <- self$counts_im1
    counts_em1 <- self$counts_em1
    rates_bth <- self$rates_bth
    rates_dth <- self$rates_dth
    rates_im1 <- self$rates_im1
    rates_em1 <- self$rates_em1
    exposure <- self$exposure
    is_dominant <- self$is_dominant
    ## immigration
    ans_im1 <- stats::dpois(
                          x = counts_im1,
                          lambda = rates_im1,
                          log = TRUE
                      )
    ans_im1 <- matrix(ans_im1, nrow = n_particle, ncol = n_region)
    ## deaths and emigration
    ## dmultinom not vectorised, so do one particle at a time
    ans_dth_em1 <- matrix(nrow = n_particle, ncol = n_region)
    for (i_particle in seq_len(n_particle)) {
        for (i_region in seq_len(n_region)) {
            stk_start <- stock_start[i_particle, i_region]
            cnt_dth <- counts_dth[i_particle, i_region]
            cnt_im1 <- counts_im1[i_particle, i_region]
            cnt_em1 <- counts_em1[i_particle, i_region]
            rt_dth <- rates_dth[i_particle, i_region]
            rt_em1 <- rates_em1[i_particle, i_region]
            prob_exit <- 1 - exp(-0.5 * (rt_dth + rt_em1))
            prob_dth <- rt_dth / (rt_dth + rt_em1)
            prob <- c(prob_dth * prob_exit, (1 - prob_dth) * prob_exit, 1 - prob_exit)
            is_im1_even <- (cnt_im1 %% 2L) == 0L
            if (is_im1_even) {
                size <- stk_start + 0.5 * cnt_im1
                if (size >= cnt_dth + cnt_em1) {
                    x <- c(cnt_dth, cnt_em1, size - cnt_dth - cnt_em1)
                    ans_particle <- stats::dmultinom(x = x, size = size, prob = prob, log = TRUE)
                } else {
                    ans_particle <- -1e308
                }
            } else {
                size1 <- stk_start + 0.5 * cnt_im1 + 0.5
                size2 <- size1 - 1
                x1 <- c(cnt_dth, cnt_em1, size1 - cnt_dth - cnt_em1)
                x2 <- c(cnt_dth, cnt_em1, size2 - cnt_dth - cnt_em1)
                if (size1 >= cnt_dth + cnt_em1) {
                    ans_part1 <- stats::dmultinom(x = x1, size = size1, prob = prob, log = TRUE)
                } else {
                    ans_part1 <- -1e308
                }
                if (size2 >= cnt_dth + cnt_em1) {
                    ans_part2 <- stats::dmultinom(x = x2, size = size2, prob = prob, log = TRUE)
                } else {
                    ans_part2 <- -1e308
                }
                ## numerically safe version of log(0.5(exp(ans_part1) + exp(ans_part2)))
                if (ans_part1 > ans_part2)
                    ans_particle <- log(0.5) + (ans_part1 + log1p(exp(ans_part2 - ans_part1)))
                else
                    ans_particle <- log(0.5) + (ans_part2 + log1p(exp(ans_part1 - ans_part2)))
            }
            ans_dth_em1[i_particle, i_region] <- ans_particle
        }
    }
    ## births
    if (is_dominant) {
        ans_bth <- stats::dpois(
                              x = counts_bth,
                              lambda = rates_bth * exposure,
                              log = TRUE
                          )
    } else {
        ans_bth <- 0
    }
    ans_bth <- matrix(ans_bth, nrow = n_particle, ncol = n_region)
    ## combine
    ans <- ans_im1 + ans_dth_em1 + ans_bth
    rowSums(ans)
}


PFilterWithReg$set(
  which = "public",
  name = "calc_logtrans",
  value = calc_logtrans_withreg
  )


## calc_logwt_unnorm ----------------------------------------------------------

## Calculate unnormalised weights for proposals.
## The 'NoReg' and 'WithReg' methods for 'calc_loglik!',
## 'calc_logtrans', and 'calc_logimp' all return
## vectors of length 'n_particle', so we can
## use one method for 'calc_wt_unnorm'
## for both cases. Also record the sums of the
## components of the weights, as diagnostics.

## HAS_TESTS
calc_logwt_unnorm <- function(i_interval, is_forecast) {
    if (self$is_parallelogram_second[[i_interval]]) {
        self$logwt_unnorm[, i_interval + 1L] <- self$logwt_unnorm[, i_interval]
        self$sum_loglik[[i_interval + 1L]] <- self$sum_loglik[[i_interval]]
        self$sum_logtrans[[i_interval]] <- self$sum_logtrans[[i_interval - 1L]]
        self$sum_logimp[[i_interval + 1L]] <- self$sum_logimp[[i_interval]]
        self$sum_logwt_unnorm[[i_interval + 1L]] <- self$sum_logwt_unnorm[[i_interval]]
    }
    else {
        loglik <- self$calc_loglik(i_interval - 1L) ## delete -1L when switching to C
        logtrans <- self$calc_logtrans(i_interval = i_interval)
        logimp <- self$calc_logimp(is_forecast)
        logwt_unnorm_prev <- self$calc_logwt_unnorm_prev(i_interval)
        logwt_unnorm <- loglik + logtrans - logimp + logwt_unnorm_prev
        self$logwt_unnorm[, i_interval + 1L] <- logwt_unnorm
        self$sum_loglik[[i_interval + 1L]] <- sum(loglik)
        self$sum_logtrans[[i_interval]] <- sum(logtrans)
        self$sum_logimp[[i_interval + 1L]] <- sum(logimp)
        self$sum_logwt_unnorm[[i_interval + 1L]] <- sum(logwt_unnorm)
    }
}

PFilter$set(
  which = "public",
  name = "calc_logwt_unnorm",
  value = calc_logwt_unnorm
)


## calc_logwt_unnorm_init -----------------------------------------------------

## Calculate unnormalised weights for the initial stock.
## We use a flat prior for the initial stock, so
## there is no equivalent of the 'logtrans' term
## found in 'calc_logwt_unnorm'.
## The 'NoReg' and 'WithReg' methods for 'calc_loglik_init'
## and 'calc_logimp_init' both return vectors of length
## 'n_particle', so we can use one method for 'calc_wt_unnorm_init'
## for both cases. Also record the sums of the
## components of the weights, as diagnostics.

## HAS_TESTS
calc_logwt_unnorm_init <- function() {
  loglik <- self$calc_loglik_init()
  logimp <- self$calc_logimp_init()
  logwt_unnorm <- loglik - logimp
  self$logwt_unnorm[, 1L] <- logwt_unnorm
  self$sum_loglik[[1L]] <- sum(loglik)
  self$sum_logimp[[1L]] <- sum(logimp)
  self$sum_logwt_unnorm[[1L]] <- sum(logwt_unnorm)
}

PFilter$set(
  which = "public",
  name = "calc_logwt_unnorm_init",
  value = calc_logwt_unnorm_init
  )


## calc_logwt_unnorm_prev -----------------------------------------------------

## Calculate unnormalised weights from previous interval.
## These differ depending on whether resampling was carried out.
## Note that 'resampled' has length 'n_interval + 1' and
## 'logwt_unnorm' has 'n_interval + 1' columns, so using
## 'i_interval' pulls out previous interval.

## NO_TESTS
calc_logwt_unnorm_prev <- function(i_interval) {
    n_particle <- self$n_particle
    logwt_unnorm <- self$logwt_unnorm
    resampled <- self$resampled
    resampled_prev <- resampled[[i_interval]]
    if (resampled_prev)
        ans <- rep(-1 * log(n_particle), times = n_particle)
    else
        ans <- logwt_unnorm[, i_interval]
    ans
}

PFilter$set(
  which = "public",
  name = "calc_logwt_unnorm_prev",
  value = calc_logwt_unnorm_prev
  )


## calc_n_unique --------------------------------------------------------------

## Calculate the number of unique particles
## at each index, after resampling, at the
## end of the the estimation period.
## Algorithm same for no-region and with-region
## cases.

## HAS_TESTS
calc_n_unique <- function() {
  n_interval <- self$n_interval
  index_ancestor <- self$index_ancestor
  index_output <- self$index_output
  n_unique <- numeric(length = n_interval + 1L)
  for (i_interval in seq_len(n_interval + 1L)) {
      i_ancestor <- index_ancestor[index_output, i_interval]
      n_unique[[i_interval]] <- length(unique(i_ancestor))
  }
  self$n_unique <- n_unique
}

PFilter$set(
  which = "public",
  name = "calc_n_unique",
  value = calc_n_unique
)


## calc_rates -----------------------------------------------------------------

## Generate reformated rates for interval 'i_interval'.

## HAS_TESTS

calc_rates_noreg <- function(i_interval) {
    n_particle <- self$n_particle
    self$rates_bth <- rep(self$rates_births[[i_interval]],
                          times = n_particle
                          )
    self$rates_dth <- rep(self$rates_deaths[[i_interval]],
                          times = n_particle
                          )
    self$rates_im1 <- rep(self$rates_immigration1[[i_interval]],
                          times = n_particle
                          )
    self$rates_em1 <- rep(self$rates_emigration1[[i_interval]],
                          times = n_particle
                          )
    self$rates_im2 <- rep(self$rates_immigration2[[i_interval]],
                          times = n_particle
                          )
    self$rates_em2 <- rep(self$rates_emigration2[[i_interval]],
                          times = n_particle
                          )
    self$rates_cim <- self$rates_im1 + self$rates_im2
    self$rates_cem <- self$rates_em1 + self$rates_em2
    if (self$is_parallelogram_first[i_interval]) {
        self$rates_bth_second <- rep(self$rates_births[[i_interval + 1L]],
                                     times = n_particle
                                     )
        self$rates_dth_second <- rep(self$rates_deaths[[i_interval + 1L]],
                                     times = n_particle
                                     )
        self$rates_im1_second <- rep(self$rates_immigration1[[i_interval + 1L]],
                                     times = n_particle
                                     )
        self$rates_em1_second <- rep(self$rates_emigration1[[i_interval + 1L]],
                                     times = n_particle
                                     )
    }
}

PFilterNoReg$set(
  which = "public",
  name = "calc_rates",
  value = calc_rates_noreg
)


## HAS_TESTS
calc_rates_withreg <- function(i_interval) {
  n_particle <- self$n_particle
  n_region <- self$n_region
  rates_bth <- rep(self$rates_births[, i_interval],
    each = n_particle
  )
  rates_dth <- rep(self$rates_deaths[, i_interval],
    each = n_particle
  )
  rates_in <- rep(self$rates_internal_in[, i_interval],
    each = n_particle
  )
  rates_out <- rep(self$rates_internal_out[, i_interval],
    each = n_particle
  )
  rates_im1 <- rep(self$rates_immigration1[, i_interval],
    each = n_particle
  )
  rates_em1 <- rep(self$rates_emigration1[, i_interval],
    each = n_particle
  )
  rates_im2 <- rep(self$rates_immigration2[, i_interval],
    each = n_particle
  )
  rates_em2 <- rep(self$rates_emigration2[, i_interval],
    each = n_particle
  )
  rates_cim <- rates_im1 + rates_im2
  rates_cem <- rates_em1 + rates_em2
  rates_cin <- rates_in + rates_cim
  rates_cout <- rates_out + rates_cem
  self$rates_bth <- matrix(rates_bth,
    nrow = n_particle,
    ncol = n_region
  )
  self$rates_dth <- matrix(rates_dth,
    nrow = n_particle,
    ncol = n_region
  )
  self$rates_in <- matrix(rates_in,
    nrow = n_particle,
    ncol = n_region
  )
  self$rates_out <- matrix(rates_out,
    nrow = n_particle,
    ncol = n_region
  )
  self$rates_im1 <- matrix(rates_im1,
    nrow = n_particle,
    ncol = n_region
  )
  self$rates_em1 <- matrix(rates_em1,
    nrow = n_particle,
    ncol = n_region
  )
  self$rates_im2 <- matrix(rates_im2,
    nrow = n_particle,
    ncol = n_region
  )
  self$rates_em2 <- matrix(rates_em2,
    nrow = n_particle,
    ncol = n_region
  )
  self$rates_cim <- matrix(rates_cim,
    nrow = n_particle,
    ncol = n_region
  )
  self$rates_cem <- matrix(rates_cem,
    nrow = n_particle,
    ncol = n_region
  )
  self$rates_cin <- matrix(rates_cin,
    nrow = n_particle,
    ncol = n_region
  )
  self$rates_cout <- matrix(rates_cout,
    nrow = n_particle,
    ncol = n_region
  )
}

PFilterWithReg$set(
  which = "public",
  name = "calc_rates",
  value = calc_rates_withreg
)


## calc_stock_start -----------------------------------------------------------

## Stock at start of this interval equals stock at end of last
## interval. We need to use 'index_parent' to account for
## the effect of resampling.

## HAS_TESTS
calc_stock_start_noreg <- function(i_interval) {
  stock <- self$counts_stock[, i_interval] ## dbl vector n_particle
  i_parent <- self$index_parent[, i_interval] ## int vector n_particle
  self$stock_start <- stock[i_parent] ## dble vector n_particle
}

PFilterNoReg$set(
  which = "public",
  name = "calc_stock_start",
  value = calc_stock_start_noreg
)

## HAS_TESTS
calc_stock_start_withreg <- function(i_interval) {
  stock <- self$counts_stock[, , i_interval] ## dbl matrix n_particle x n_region
  i_parent <- self$index_parent[, i_interval] ## int vector n_particle
  self$stock_start <- stock[i_parent, ] ## dbl matrix n_particle x n_region
}

PFilterWithReg$set(
  which = "public",
  name = "calc_stock_start",
  value = calc_stock_start_withreg
)


## calc_trajectories ----------------------------------------------------------

## Rearrange counts estimates so that each row corresponds
## one trajectory, as viewed from the perspective of
## particles at the end of filtering.

## HAS_TESTS
calc_trajectories_noreg <- function() {
  n_interval <- self$n_interval
  index_ancestor <- self$index_ancestor
  counts_stock <- self$counts_stock
  counts_births <- self$counts_births
  counts_deaths <- self$counts_deaths
  counts_immigration1 <- self$counts_immigration1
  counts_emigration1 <- self$counts_emigration1
  counts_immigration2 <- self$counts_immigration2
  counts_emigration2 <- self$counts_emigration2
  for (i_interval in seq.int(from = n_interval, to = 1L)) {
    i_ancestor <- index_ancestor[, i_interval + 1L]
    counts_stock[, i_interval + 1L] <- counts_stock[i_ancestor, i_interval + 1L]
    counts_births[, i_interval] <- counts_births[i_ancestor, i_interval]
    counts_deaths[, i_interval] <- counts_deaths[i_ancestor, i_interval]
    counts_immigration1[, i_interval] <- counts_immigration1[i_ancestor, i_interval]
    counts_emigration1[, i_interval] <- counts_emigration1[i_ancestor, i_interval]
    counts_immigration2[, i_interval] <- counts_immigration2[i_ancestor, i_interval]
    counts_emigration2[, i_interval] <- counts_emigration2[i_ancestor, i_interval]
  }
  i_ancestor <- index_ancestor[, 1L]
  counts_stock[, 1L] <- counts_stock[i_ancestor, 1L]
  self$counts_stock <- counts_stock
  self$counts_births <- counts_births
  self$counts_deaths <- counts_deaths
  self$counts_immigration1 <- counts_immigration1
  self$counts_emigration1 <- counts_emigration1
  self$counts_immigration2 <- counts_immigration2
  self$counts_emigration2 <- counts_emigration2
}

PFilterNoReg$set(
  which = "public",
  name = "calc_trajectories",
  value = calc_trajectories_noreg
)


## HAS_TESTS
calc_trajectories_withreg <- function() {
  n_region <- self$n_region
  n_interval <- self$n_interval
  index_ancestor <- self$index_ancestor
  counts_stock <- self$counts_stock
  counts_births <- self$counts_births
  counts_deaths <- self$counts_deaths
  counts_internal_in <- self$counts_internal_in
  counts_internal_out <- self$counts_internal_out
  counts_immigration1 <- self$counts_immigration1
  counts_emigration1 <- self$counts_emigration1
  counts_immigration2 <- self$counts_immigration2
  counts_emigration2 <- self$counts_emigration2
  for (i_interval in seq.int(from = n_interval, to = 1L)) {
    i_ancestor <- index_ancestor[, i_interval + 1L]
    counts_stock[, , i_interval + 1L] <- counts_stock[i_ancestor, , i_interval + 1L]
    counts_births[, , i_interval] <- counts_births[i_ancestor, , i_interval]
    counts_deaths[, , i_interval] <- counts_deaths[i_ancestor, , i_interval]
    counts_internal_in[, , i_interval] <- counts_internal_in[i_ancestor, , i_interval]
    counts_internal_out[, , i_interval] <- counts_internal_out[i_ancestor, , i_interval]
    counts_immigration1[, , i_interval] <- counts_immigration1[i_ancestor, , i_interval]
    counts_emigration1[, , i_interval] <- counts_emigration1[i_ancestor, , i_interval]
    counts_immigration2[, , i_interval] <- counts_immigration2[i_ancestor, , i_interval]
    counts_emigration2[, , i_interval] <- counts_emigration2[i_ancestor, , i_interval]
  }
  i_ancestor <- index_ancestor[, 1L]
  counts_stock[, , 1L] <- counts_stock[i_ancestor, , 1L]
  self$counts_stock <- counts_stock
  self$counts_births <- counts_births
  self$counts_deaths <- counts_deaths
  self$counts_internal_in <- counts_internal_in
  self$counts_internal_out <- counts_internal_out
  self$counts_immigration1 <- counts_immigration1
  self$counts_emigration1 <- counts_emigration1
  self$counts_immigration2 <- counts_immigration2
  self$counts_emigration2 <- counts_emigration2
}

PFilterWithReg$set(
  which = "public",
  name = "calc_trajectories",
  value = calc_trajectories_withreg
)



## draw_counts_bth ------------------------------------------------------------

## Generate births to members of cohort,
## given stock and births rates.
## Used in 'forecast_account'.
## Births are drawn straight from the
## posterior distribution, so no need
## to calculate 'logimp'.

## HAS_TESTS
draw_counts_bth_noreg <- function() {
  n_particle <- self$n_particle
  is_dominant <- self$is_dominant
  if (is_dominant) {
    rates_bth <- self$rates_bth
    exposure <- self$exposure_approx2
    lambda <- rates_bth * exposure
    self$counts_bth <- stats::rpois(
      n = n_particle,
      lambda = lambda
    )
  } else {
    self$counts_bth <- rep(0, times = n_particle)
  }
}

PFilterNoReg$set(
  which = "public",
  name = "draw_counts_bth",
  value = draw_counts_bth_noreg
)


## HAS_TESTS
draw_counts_bth_withreg <- function() {
  n_particle <- self$n_particle
  n_region <- self$n_region
  is_dominant <- self$is_dominant
  if (is_dominant) {
    rates_bth <- self$rates_bth
    exposure <- self$exposure_approx2
    lambda <- rates_bth * exposure
    counts_bth <- stats::rpois(
      n = n_particle * n_region,
      lambda = lambda
    )
    self$counts_bth <- matrix(counts_bth,
      nrow = n_particle,
      ncol = n_region
    )
  } else {
    self$counts_bth <- matrix(0,
      nrow = n_particle,
      ncol = n_region
    )
  }
}

PFilterWithReg$set(
  which = "public",
  name = "draw_counts_bth",
  value = draw_counts_bth_withreg
)


## draw_counts_dth_noreg ------------------------------------------------------

## Randomly generate deaths, using approximate exposure.
## Only used in forecasts. We need to adjust for the
## fact we are using approximate exposure, so we
## calculate 'logimp'.

## HAS_TESTS
draw_counts_dth_noreg <- function() {
  n_particle <- self$n_particle
  rates_dth <- self$rates_dth
  exposure <- self$exposure_approx1
  lambda <- rates_dth * exposure
  counts_dth <- stats::rpois(
    n = n_particle,
    lambda = lambda
  )
  logimp_counts_dth <- stats::dpois(
    x = counts_dth,
    lambda = lambda,
    log = TRUE
  )
  self$counts_dth <- counts_dth
  self$logimp_counts_dth <- logimp_counts_dth
}

PFilterNoReg$set(
  which = "public",
  name = "draw_counts_dth",
  value = draw_counts_dth_noreg
)


## HAS_TESTS
draw_counts_dth_withreg <- function() {
  n_particle <- self$n_particle
  n_region <- self$n_region
  rates_dth <- self$rates_dth
  exposure <- self$exposure_approx1
  lambda <- rates_dth * exposure
  counts_dth <- stats::rpois(
    n = n_particle * n_region,
    lambda = lambda
  )
  logimp_counts_dth <- stats::dpois(
    x = counts_dth,
    lambda = lambda,
    log = TRUE
  )
  counts_dth <- matrix(counts_dth,
    nrow = n_particle,
    ncol = n_region
  )
  logimp_counts_dth <- matrix(logimp_counts_dth,
    nrow = n_particle,
    ncol = n_region
  )
  self$counts_dth <- counts_dth
  self$logimp_counts_dth <- rowSums(logimp_counts_dth)
}

PFilterWithReg$set(
  which = "public",
  name = "draw_counts_dth",
  value = draw_counts_dth_withreg
)


## draw_counts_forecast_noreg -------------------------------------------------

## Draw counts for stock, deaths, immigration, emigration, and,
## in the case of the dominant sex/gender, births

## HAS_TESTS
draw_counts_forecast_noreg <- function() {
    n_particle <- self$n_particle
    is_dominant <- self$is_dominant
    stock_start <- self$stock_start
    rates_bth <- self$rates_bth
    rates_dth <- self$rates_dth
    rates_cim <- self$rates_cim
    rates_cem <- self$rates_cem
    counts_cim <- stats::rpois(n = n_particle, lambda = rates_cim)
    prob_exit <- 1 - exp(-0.5 * (rates_dth + rates_cem))
    prob_dth <- (rates_dth / (rates_dth + rates_cem))
    im_is_odd  <- (counts_cim %% 2L) != 0L
    adj_odd <- sample(c(-0.5, 0.5), size = n_particle, replace = TRUE)
    cim_start  <- 0.5 * counts_cim + im_is_odd * adj_odd
    cim_end <- counts_cim - cim_start
    stock_start_im <- stock_start + cim_start
    counts_exit <- stats::rbinom(n = n_particle, prob = prob_exit, size = stock_start_im)
    counts_dth <- stats::rbinom(n = n_particle, prob = prob_dth, size = counts_exit)
    counts_cem <- counts_exit - counts_dth
    stock_end <- stock_start_im - counts_exit + cim_end
    self$stock_end <- stock_end
    self$counts_dth <- counts_dth
    self$counts_im1 <- counts_cim
    self$counts_im2 <- rep(0, times = n_particle)
    self$counts_em1 <- counts_cem
    self$counts_em2 <- rep(0, times = n_particle)
    if (is_dominant) {
        exposure <- 0.25 * (stock_start + stock_end)
        counts_bth <- stats::rpois(n = n_particle, lambda = rates_bth * exposure)
        self$counts_bth <- counts_bth
    }
}

PFilterNoReg$set(
  which = "public",
  name = "draw_counts_forecast",
  value = draw_counts_forecast_noreg
  )


## draw_counts_im_em ----------------------------------------------------------

## Given values for international net and gross migration,
## and for migration rates, generate counts for immigration
## and emigration, for current interval


## HAS_TESTS
draw_counts_im_em_noreg <- function(i_interval) {
    has_one_imem <- self$has_one_imem
    n_particle <- self$n_particle
    rates_im1 <- self$rates_im1
    rates_em1 <- self$rates_em1
    rates_cim <- self$rates_cim
    rates_cem <- self$rates_cem
    net_mig <- self$net_mig
    gross_mig <- self$gross_mig
    counts_cim <- (gross_mig + net_mig) / 2
    counts_cem <- (gross_mig - net_mig) / 2
    if (has_one_imem) {
        counts_im1 <- counts_cim
        counts_em1 <- counts_cem
        logimp_counts_im_em <- rep(0, times = n_particle)
    } else {
        prob_im1 <- ifelse(rates_cim > 0, rates_im1 / rates_cim, 0)
        prob_em1 <- ifelse(rates_cem > 0, rates_em1 / rates_cem, 0)
        counts_im1 <- stats::rbinom(
                                 n = n_particle,
                                 size = counts_cim,
                                 prob = prob_im1
                             )
        counts_em1 <- stats::rbinom(
                                 n = n_particle,
                                 size = counts_cem,
                                 prob = prob_em1
                             )
        logimp_counts_im1 <- stats::dbinom(
                                        x = counts_im1,
                                        size = counts_cim,
                                        prob = prob_im1,
                                        log = TRUE
                                    )
        logimp_counts_em1 <- stats::dbinom(
                                        x = counts_em1,
                                        size = counts_cem,
                                        prob = prob_em1,
                                        log = TRUE
                                    )
        logimp_counts_im_em <- logimp_counts_im1 + logimp_counts_em1
    }
    counts_im2 <- counts_cim - counts_im1
    counts_em2 <- counts_cem - counts_em1
    self$counts_im1 <- counts_im1
    self$counts_em1 <- counts_em1
    self$counts_im2 <- counts_im2
    self$counts_em2 <- counts_em2
    self$logimp_counts_im_em <- logimp_counts_im_em
    if (self$is_parallelogram_first[[i_interval]]) {
        net_mig_second <- self$net_mig_second
        gross_mig_second <- self$gross_mig_second
        counts_im1_second <- (gross_mig_second + net_mig_second) / 2
        counts_em1_second <- (gross_mig_second - net_mig_second) / 2
        self$counts_im1_second <- counts_im1_second
        self$counts_em1_second <- counts_em1_second
    }
}

PFilterNoReg$set(
  which = "public",
  name = "draw_counts_im_em",
  value = draw_counts_im_em_noreg
)


## HAS_TESTS
draw_counts_im_em_withreg <- function(i_interval) {
    has_one_imem <- self$has_one_imem
    is_split_mig <- self$is_split_mig
    n_particle <- self$n_particle
    n_region <- self$n_region
    net_mig <- self$net_mig ## includes internal migration
    if (is_split_mig) {
        counts_in <- self$counts_in
        counts_out <- self$counts_out
        net_mig <- net_mig - counts_in + counts_out ## excludes internal migration
    }
    gross_mig <- self$gross_mig
    counts_cim <- (gross_mig + net_mig) / 2
    counts_cem <- (gross_mig - net_mig) / 2
    if (has_one_imem || !is_split_mig) {
        counts_im1 <- counts_cim
        counts_em1 <- counts_cem
        logimp_counts_im_em <- rep(0, times = n_particle * n_region)
    } else {
        rates_im1 <- self$rates_im1
        rates_em1 <- self$rates_em1
        rates_cim <- self$rates_cim
        rates_cem <- self$rates_cem
        prob_im1 <- ifelse(rates_cim > 0, rates_im1 / rates_cim, 0)
        prob_em1 <- ifelse(rates_cem > 0, rates_em1 / rates_cem, 0)
        counts_im1 <- stats::rbinom(
                                 n = n_particle * n_region,
                                 size = counts_cim,
                                 prob = prob_im1
                             )
        counts_em1 <- stats::rbinom(
                                 n = n_particle * n_region,
                                 size = counts_cem,
                                 prob = prob_em1
                             )
        logimp_im1 <- stats::dbinom(
                                 x = counts_im1,
                                 size = counts_cim,
                                 prob = prob_im1,
                                 log = TRUE
                             )
        logimp_em1 <- stats::dbinom(
                                 x = counts_em1,
                                 size = counts_cem,
                                 prob = prob_em1,
                                 log = TRUE
                             )
        logimp_counts_im_em <- logimp_im1 + logimp_em1
    }
    counts_cim <- matrix(counts_cim,
                         nrow = n_particle,
                         ncol = n_region
                         )
    counts_cem <- matrix(counts_cem,
                         nrow = n_particle,
                         ncol = n_region
                         )
    counts_im1 <- matrix(counts_im1,
                         nrow = n_particle,
                         ncol = n_region
                         )
    counts_em1 <- matrix(counts_em1,
                         nrow = n_particle,
                         ncol = n_region
                         )
    logimp_counts_im_em <- matrix(logimp_counts_im_em,
                                  nrow = n_particle,
                                  ncol = n_region
                                  )
    counts_im2 <- counts_cim - counts_im1
    counts_em2 <- counts_cem - counts_em1
    self$counts_im1 <- counts_im1
    self$counts_em1 <- counts_em1
    self$counts_im2 <- counts_im2
    self$counts_em2 <- counts_em2
    self$logimp_counts_im_em <- rowSums(logimp_counts_im_em) ## sum over regions
}

PFilterWithReg$set(
  which = "public",
  name = "draw_counts_im_em",
  value = draw_counts_im_em_withreg
)


## draw_counts_in_out ---------------------------------------------------------

## Draw values for 'counts_internal_in' and 'counts_internal_out',
## by region. The sum of values for 'counts_internal_in' must equal
## the sum of values for 'counts_internal_out'.
## There is only a method for 'WithReg', since the "NoReg" models
## do not include internal migration.

## Once we have data for migration we could use
## that to help determine proposals, alongside rates.

## HAS_TESTS
draw_counts_in_out <- function() {
    n_particle <- self$n_particle
    n_region <- self$n_region
    is_split_mig <- self$is_split_mig
    if (is_split_mig) {
        rates_out <- self$rates_out
        rates_in <- self$rates_in
        exposure <- self$exposure_approx2
        lambda <- rates_out * exposure
        counts_out <- stats::rpois(
                                 n = n_particle * n_region,
                                 lambda = lambda
                             )
        logimp_out <- stats::dpois(
                                 x = counts_out,
                                 lambda = lambda,
                                 log = TRUE
                             )
        counts_out <- matrix(counts_out,
                             nrow = n_particle,
                             ncol = n_region
                             )
        logimp_out <- matrix(logimp_out,
                             nrow = n_particle,
                             ncol = n_region
                             )
        counts_out_totals <- rowSums(counts_out)
        logimp_out_totals <- rowSums(logimp_out)
        counts_in <- matrix(
            nrow = n_particle,
            ncol = n_region
        )
        logimp_in_totals <- numeric(length = n_particle)
        for (i_particle in seq_len(n_particle)) {
            size <- counts_out_totals[[i_particle]]
            prob <- rates_in[i_particle, ]
            counts <- stats::rmultinom(
                                 n = 1L,
                                 size = size,
                                 prob = prob
                             )
            logimp_in_total <- stats::dmultinom(
                                          x = counts,
                                          size = size,
                                          prob = prob,
                                          log = TRUE
                                      )
            counts_in[i_particle, ] <- counts
            logimp_in_totals[[i_particle]] <- logimp_in_total
        }
        logimp_counts_in_out <- logimp_out_totals + logimp_in_totals
    } else {
        counts_in <- matrix(0,
                            nrow = n_particle,
                            ncol = n_region
                            )
        counts_out <- matrix(0,
                             nrow = n_particle,
                             ncol = n_region
                             )
        logimp_counts_in_out <- rep(0, times = n_particle)
    }
    self$counts_in <- counts_in
    self$counts_out <- counts_out
    self$logimp_counts_in_out <- logimp_counts_in_out
}

PFilterWithReg$set(
  which = "public",
  name = "draw_counts_in_out",
  value = draw_counts_in_out
  )



## draw_indices_output --------------------------------------------------------

## Randomly select index for trajectories that will be
## reported in final outputs

## HAS_TESTS
draw_index_output <- function() {
    n_particle <- self$n_particle
    n_thin <- self$n_thin
    size <- n_particle %/% n_thin # integer division
    index_output <- sample(x = n_particle,
                           size = size,
                           replace = FALSE)
    self$index_output <- index_output
}

PFilter$set(
  which = "public",
  name = "draw_index_output",
  value = draw_index_output
  )


## draw_gross_mig -------------------------------------------------------------

## Randomly generate 'gross_mig', the sum of inward and outward
## international migration. Note that in the regional model,
## 'gross_mig' does not include internal migration.

## HAS_TESTS
draw_gross_mig_noreg <- function(i_interval) {
    rates_im <- self$rates_cim
    rates_em <- self$rates_cem
    exposure <- self$exposure_approx2
    net_mig <- self$net_mig
    ans <- draw_gross_mig_inner(rates_im = rates_im,
                                rates_em = rates_em,
                                exposure = exposure,
                                net_mig = net_mig)
    self$gross_mig <- ans[[1L]]
    self$logimp_gross_mig <- ans[[2]]
    if (self$is_parallelogram_first[[i_interval]]) {
        rates_im <- self$rates_im1_second
        rates_em <- self$rates_em1_second
        exposure <- self$exposure_approx2_second
        net_mig <- self$net_mig_second
        ans <- draw_gross_mig_inner(rates_im = rates_im,
                                    rates_em = rates_em,
                                    exposure = exposure,
                                    net_mig = net_mig)
        self$gross_mig_second <- ans[[1L]]
        self$logimp_gross_mig <- self$logimp_gross_mig + ans[[2]]
    }
}

PFilterNoReg$set(
  which = "public",
  name = "draw_gross_mig",
  value = draw_gross_mig_noreg
)


## HAS_TESTS
draw_gross_mig_withreg <- function(i_interval) {
  n_particle <- self$n_particle
  n_region <- self$n_region
  is_split_mig <- self$is_split_mig
  if (is_split_mig) {
      rates_cim <- self$rates_cim
      rates_cem <- self$rates_cem
      counts_in <- self$counts_in
      counts_out <- self$counts_out
  } else {
      rates_cin <- self$rates_cin
      rates_cout <- self$rates_cout
  }
  exposure <- self$exposure_approx2
  net_mig <- self$net_mig ## includes internal migration
  if (is_split_mig) {
      net_mig <- net_mig - counts_in + counts_out ## excludes internal migration
      lambda <- rates_cim + rates_cem * exposure
  } else {
      lambda <- rates_cin + rates_cout * exposure
  }
  lower <- abs(net_mig)
  gross_mig <- rpoistr(
    n = n_particle * n_region,
    lambda = lambda,
    lower = lower
  )
  is_diff_parity <- (gross_mig %% 2L) != (lower %% 2L)
  gross_mig[is_diff_parity] <- gross_mig[is_diff_parity] + 1L
  ## when calculate the log importance probability,
  ## we need to account for the fact that there are two ways
  ## getting to the final value, except when the final
  ## value equals 'lower'
  logprob_high <- dpoistr(gross_mig,
    lambda = lambda,
    lower = lower,
    use_log = TRUE
  )
  logprob_low <- dpoistr(gross_mig - 1L,
    lambda = lambda,
    lower = lower,
    use_log = TRUE
  )
  logprob_low[gross_mig == lower] <- -Inf
  logimp_gross_mig <- log_sum_exp_2(x = logprob_high,
                                    y = logprob_low)
  gross_mig <- matrix(gross_mig,
    nrow = n_particle,
    ncol = n_region
  )
  logimp_gross_mig <- matrix(logimp_gross_mig,
    nrow = n_particle,
    ncol = n_region
  )
  self$gross_mig <- gross_mig
  self$logimp_gross_mig <- rowSums(logimp_gross_mig)
}

PFilterWithReg$set(
  which = "public",
  name = "draw_gross_mig",
  value = draw_gross_mig_withreg
  )


## draw_net_mig_parallelogram --------------------------------------------------

## HAS_TESTS
draw_net_mig_parallelogram <- function() {
    n_particle <- self$n_particle
    net_mig_combined <- self$net_mig_combined
    stock_start <- self$stock_start
    counts_dth <- self$counts_dth
    rates_im_first <- self$rates_im1
    rates_em_first <- self$rates_em1
    rates_im_second <- self$rates_im1_second
    rates_em_second <- self$rates_em1_second
    exposure_approx_first <- self$exposure_approx_first
    exposure_approx_second <- self$exposure_approx_second
    lower <- -1 * net_mig_combined + 2 * (counts_dth - stock_start)
    mu1 <- rates_im_first + rates_em_second * exposure_approx_second
    mu2 <- rates_im_second + rates_em_first * exposure_approx_first
    diff_net_mig <- rskeltr(n = n_particle,
                            mu1 = mu1,
                            mu2 = mu2,
                            lower = lower)
    is_diff_parity <- (net_mig_combined %% 2L) != (diff_net_mig %% 2L)
    diff_net_mig[is_diff_parity] <- diff_net_mig[is_diff_parity] + 1
    ## when we calculate the log importance probability,
    ## we need to account for the fact that there are two ways
    ## of getting to the final value, except when the final
    ## value equals 'lower'
    logprob_high <- dskeltr(diff_net_mig,
                            mu1 = mu1,
                            mu2 = mu2,
                            lower = lower,
                            use_log = TRUE
                            )
    logprob_low <- dskeltr(diff_net_mig - 1,
                           mu1 = mu1,
                           mu2 = mu2,
                           lower = lower,
                           use_log = TRUE
                           )
    logprob_low[diff_net_mig == lower] <- -Inf
    logimp_stock_end_net_mig = log_sum_exp_2(x = logprob_high,
                                             y = logprob_low)
    net_mig <- (net_mig_combined + diff_net_mig) / 2
    net_mig_second <- (net_mig_combined - diff_net_mig) / 2
    self$net_mig <- net_mig
    self$net_mig_second <- net_mig_second
    self$stock_end <- stock_start - counts_dth + net_mig
    self$logimp_stock_end_net_mig <- self$logimp_stock_end_net_mig + logimp_stock_end_net_mig
}


PFilterNoReg$set(
  which = "public",
  name = "draw_net_mig_parallelogram",
  value = draw_net_mig_parallelogram
  )


## draw_rates -----------------------------------------------------------------

## Generate rates for interval 'i_interval' by drawing
## from gamma distributions with mean governed by
## 'rates_*' and dispersion governed by 'disp_*'.

## HAS_TESTS
draw_rates_noreg <- function(i_interval) {
    n_particle <- self$n_particle
    self$rates_bth <- draw_rates_inner(n = n_particle,
                                       mean = self$rates_births[[i_interval]],
                                       disp = self$disp_births)
    self$rates_dth <- draw_rates_inner(n = n_particle,
                                       mean = self$rates_deaths[[i_interval]],
                                       disp = self$disp_deaths)
    self$rates_im1 <- draw_rates_inner(n = n_particle,
                                       mean = self$rates_immigration1[[i_interval]],
                                       disp = self$disp_immigration1)
    self$rates_em1 <- draw_rates_inner(n = n_particle,
                                       mean = self$rates_emigration1[[i_interval]],
                                       disp = self$disp_emigration1)
    self$rates_im2 <- draw_rates_inner(n = n_particle,
                                       mean = self$rates_immigration2[[i_interval]],
                                       disp = self$disp_immigration2)
    self$rates_em2 <- draw_rates_inner(n = n_particle,
                                       mean = self$rates_emigration2[[i_interval]],
                                       disp = self$disp_emigration2)
    self$rates_cim <- self$rates_im1 + self$rates_im2
    self$rates_cem <- self$rates_em1 + self$rates_em2
    if (self$is_parallelogram_first[i_interval]) {
        self$rates_bth_second <- draw_rates_inner(n = n_particle,
                                                  mean = self$rates_births[[i_interval + 1L]],
                                                  disp = self$disp_births)
        self$rates_dth_second <- draw_rates_inner(n = n_particle,
                                                  mean = self$rates_deaths[[i_interval + 1L]],
                                                  disp = self$disp_deaths)
        self$rates_im1_second <- draw_rates_inner(n = n_particle,
                                                  mean = self$rates_immigration1[[i_interval + 1L]],
                                                  disp = self$disp_immigration1)
        self$rates_em1_second <- draw_rates_inner(n = n_particle,
                                                  mean = self$rates_emigration1[[i_interval + 1L]],
                                                  disp = self$disp_emigration1)
    }
}

PFilterNoReg$set(
  which = "public",
  name = "draw_rates",
  value = draw_rates_noreg
)




## draw_stock_end_init --------------------------------------------------------

## Draw value for stock at end of initial period.
## 'x_0' consists entirely of stock at the end of
## the period, so drawing stock at the end of
## period 0 is equivalent to drawing 'x_0'.
## We have a non-informative prior for stock at
## then end of 'x_0', and are relying entirely on
## the data, so it is an error to not have data.

## HAS_TESTS
draw_stock_end_init_noreg <- function() {
    n_particle <- self$n_particle
  cdms <- self$cdms_stock
  stock_end <- cdms$draw_counts_true(
    i_interval = 0L, ## C function, so use C indexing
    n_particle = n_particle
  )
  logimp_stock_end_init <- cdms$calc_logimp(
    counts_true = stock_end,
    i_interval = 0L
  )
  val_first_particle <- stock_end[[1L]] ## in C need to change 1 to 0
  have_stock_data <- !is.na(val_first_particle)
  if (!have_stock_data) {
    stop(
      "problem with cohort '", self$cohort, "' and sex/gender '",
      self$sexgender, "' : no data to estimate initial stock"
    )
  }
  self$stock_end <- stock_end
  self$logimp_stock_end_init <- logimp_stock_end_init
}

PFilterNoReg$set(
  which = "public",
  name = "draw_stock_end_init",
  value = draw_stock_end_init_noreg
)


## HAS_TESTS
draw_stock_end_init_withreg <- function() {
    n_particle <- self$n_particle
    n_region <- self$n_region
    cdms <- self$cdms_stock
    stock_end <- cdms$draw_counts_true(
                          i_interval = 0L,
                          n_particle = n_particle,
                          n_region = n_region
                      )
    logimp_stock_end_init <- cdms$calc_logimp(
                                      counts_true = stock_end,
                                      i_interval = 0L
                                  )
    vals_first_particles <- stock_end[1L, ] ## in C need to change 1 to 0
    all_have_stock_data <- !anyNA(vals_first_particles)
    if (!all_have_stock_data) {
        stop(
            "problem with cohort '", self$cohort, "' and sex/gender '",
            self$sexgender, "' : no data to estimate initial stock in one or more regions"
        )
    }
    self$stock_end <- stock_end
    self$logimp_stock_end_init <- rowSums(logimp_stock_end_init)
}

PFilterWithReg$set(
  which = "public",
  name = "draw_stock_end_init",
  value = draw_stock_end_init_withreg
)


## draw_stock_end_net_mig -----------------------------------------------------

## Jointly update stock at end of interval 'stock_end',
## and net migration during interval 'net_mig' as well
## as flag 'have_stock_data'. While doing this,
## we also calculate the log importance probability.

## HAS_TESTS
draw_stock_end_net_mig_noreg <- function(i_interval, is_forecast) {
  n_particle <- self$n_particle
  stock_start <- self$stock_start
  counts_dth <- self$counts_dth
  rates_cim <- self$rates_cim
  rates_cem <- self$rates_cem
  exposure <- self$exposure_approx1
  cdms <- self$cdms_stock
  ## when forecasting, we never have stock data
  ## for the end of the interval
  if (is_forecast) {
    have_stock_data <- FALSE
  } ## when estimating, where there might be data,
  ## we try to derive data-based estimates
  else {
    stock_end <- cdms$draw_counts_true(
      i_interval = i_interval, ## change to 'i_interval + 1' when switching to C
      n_particle = n_particle
    )
    val_first_particle <- stock_end[[1L]] ## in C need to change 1 to 0
    have_stock_data <- !is.na(val_first_particle)
  }
  ## if we succeed, we obtain the implied net migration
  ## and calculate the log importance probability
  if (have_stock_data) {
    net_mig <- stock_end - stock_start + counts_dth
    logimp_stock_end_net_mig <- cdms$calc_logimp(
      counts_true = stock_end,
      i_interval = i_interval
    ) ## change to 'i_interval + 1' when switching to C
  }
  ## otherwise we randomly generate migration, subject
  ## to the constraint that the implied stock at the
  ## end of the interval must be non-negative
  else {
    lower <- counts_dth - stock_start
    expected_cem <- rates_cem * exposure
    net_mig <- rskeltr(
      n = n_particle,
      mu1 = rates_cim,
      mu2 = expected_cem,
      lower = lower
    )
    stock_end <- stock_start - counts_dth + net_mig
    ## calculate the implied log importance probability
    logimp_stock_end_net_mig <- dskeltr(
      x = net_mig,
      mu1 = rates_cim,
      mu2 = expected_cem,
      lower = lower,
      use_log = TRUE
    )
    logimp_stock_end_net_mig[is.infinite(logimp_stock_end_net_mig)] <- -1e300
  }
  ## record stock, net migration, and log importance probabilities
  self$stock_end <- stock_end
  self$net_mig <- net_mig
  self$logimp_stock_end_net_mig <- logimp_stock_end_net_mig
}

PFilterNoReg$set(
  which = "public",
  name = "draw_stock_end_net_mig",
  value = draw_stock_end_net_mig_noreg
)


## HAS_TESTS
draw_stock_end_net_mig_withreg <- function(i_interval, is_forecast) {
  n_particle <- self$n_particle
  n_region <- self$n_region
  stock_start <- self$stock_start
  counts_dth <- self$counts_dth
  rates_cin <- self$rates_cin
  rates_cout <- self$rates_cout
  exposure <- self$exposure_approx1
  cdms <- self$cdms_stock
  ## when forecasting, we never have stock data
  ## for the end of the interval
  if (is_forecast) {
    have_stock_data <- rep(FALSE, times = n_region)
    stock_end <- matrix(NA_real_,
      nrow = n_particle,
      ncol = n_region
    )
    logimp_stock_end_net_mig <- matrix(NA_real_,
      nrow = n_particle,
      ncol = n_region
    )
  }
  ## When estimating, where there might be data,
  ## so we try to derive data-based estimates,
  ## and associated log importance probabilities.
  ## If, and only if, a region does not have data,
  ## the corresponding columns of 'stock_end' and
  ## 'logimp_stock_end_net_mig' will consist entirely of NAs.
  else {
    stock_end <- cdms$draw_counts_true(
      i_interval = i_interval, ## add '+1' when switching to C
      n_particle = n_particle,
      n_region = n_region
    )
    logimp_stock_end_net_mig <- cdms$calc_logimp(
      counts_true = stock_end,
      i_interval = i_interval
    ) ## add '+1' switching to C
    val_first_particle <- stock_end[1L, ] ## change 1L to 0 in C
    have_stock_data <- !is.na(val_first_particle)
  }
  ## Go through region by region (ie column by column)
  ## and (i) if a data was available for the region, so that
  ## the regional stock has beeen estimated, calculate the implied
  ## regional net migration, or (ii) if no data was available,
  ## use regional migration rates to generate regional net
  ## migration, and the associated values for stock and
  ## for log importance probabilities.
  net_mig <- matrix(nrow = n_particle, ncol = n_region)
  for (i_region in seq_len(n_region)) {
    stock_start_reg <- stock_start[, i_region]
    counts_dth_reg <- counts_dth[, i_region]
    have_stock_data_reg <- have_stock_data[[i_region]]
    if (have_stock_data_reg) {
      stock_end_reg <- stock_end[, i_region]
      net_mig_reg <- stock_end_reg - stock_start_reg + counts_dth_reg
    } else {
      lower_reg <- counts_dth_reg - stock_start_reg
      rates_cin_reg <- rates_cin[, i_region]
      rates_cout_reg <- rates_cout[, i_region]
      exposure_reg <- exposure[, i_region]
      expected_cout_reg <- rates_cout_reg * exposure_reg
      net_mig_reg <- rskeltr(
        n = n_particle,
        mu1 = rates_cin_reg,
        mu2 = expected_cout_reg,
        lower = lower_reg
      )
      stock_end_reg <- stock_start_reg - counts_dth_reg + net_mig_reg
      logimp_reg <- dskeltr(
        x = net_mig_reg,
        mu1 = rates_cin_reg,
        mu2 = expected_cout_reg,
        lower = lower_reg,
        use_log = TRUE
      )
    }
    net_mig[, i_region] <- net_mig_reg
    stock_end[, i_region] <- stock_end_reg
    if (!have_stock_data_reg) {
      logimp_stock_end_net_mig[, i_region] <- logimp_reg
    }
  }
  ## collect results
  self$stock_end <- stock_end
  self$net_mig <- net_mig
  self$logimp_stock_end_net_mig <- rowSums(logimp_stock_end_net_mig) ## sum over regions
}

PFilterWithReg$set(
  which = "public",
  name = "draw_stock_end_net_mig",
  value = draw_stock_end_net_mig_withreg
  )



## draw_stock_end_net_mig_parallelogram ---------------------------------------

## Jointly update stock at end of interval 'stock_end_second',
## and net migration during upper and lower triangles.
## While doing this, we also calculate the log importance probability.

## HAS_TESTS
draw_stock_end_net_mig_parallelogram <- function(i_interval) {
    n_particle <- self$n_particle
    stock_start <- self$stock_start
    counts_dth <- self$counts_dth
    counts_dth_second <- self$counts_dth_second
    cdms <- self$cdms_stock
    stock_end_second <- cdms$draw_counts_true(
                                 i_interval = i_interval + 1L,
                                 n_particle = n_particle
                             )
    net_mig_combined <- stock_end_second - stock_start + counts_dth + counts_dth_second
    logimp_stock_end_net_mig <- cdms$calc_logimp(
                                         counts_true = stock_end_second,
                                         i_interval = i_interval + 1L
                                     ) 
    self$stock_end_second <- stock_end_second
    self$net_mig_combined <- net_mig_combined
    self$logimp_stock_end_net_mig <- logimp_stock_end_net_mig
}

PFilter$set(
  which = "public",
  name = "draw_stock_end_net_mig_parallelogram",
  value = draw_stock_end_net_mig_parallelogram
)


## draw_values ----------------------------------------------------------------

## Generate values for migration, population and
## (if forecasting) for births and deaths,
## at interval 'i_interval'.
## We calculate log importance probabilities ('logimp')
## whenever we draw a random value. Doing random
## number generation and calculation of probabilities
## together allows us to re-use intermediate
## quantities, and to ensure that the
## calculation of probabilities
## is consistent and complete. The main disadvantage
## of taking this approach is that we forgo a
## few opportunties for cancelling importance
## probabilities and prior probabilities.

## HAS_TESTS
draw_values_noreg <- function(i_interval, is_forecast) {
  self$calc_stock_start(i_interval = i_interval)
  self$draw_rates(i_interval = i_interval)
  ## in next version, just draw from transition function
  if (is_forecast) {
      self$calc_exposure_approx1()
      self$draw_counts_dth()
      self$draw_stock_end_net_mig(i_interval = i_interval,
                                  is_forecast = is_forecast)
      self$calc_exposure_approx2(i_interval = i_interval)
      self$draw_gross_mig(i_interval = i_interval)
      self$draw_counts_im_em(i_interval = i_interval)
      self$draw_counts_bth()
      self$calc_exposure(i_interval = i_interval)
      self$update_counts(i_interval = i_interval)
  }
  else if (self$is_parallelogram_first[[i_interval]]) {
      self$calc_counts_dth(i_interval = i_interval)
      self$draw_stock_end_net_mig_parallelogram(i_interval = i_interval)
      self$calc_exposure_approx_parallelogram()
      self$draw_net_mig_parallelogram()
      self$calc_exposure_approx2(i_interval = i_interval)
      self$draw_gross_mig(i_interval = i_interval)
      self$draw_counts_im_em(i_interval = i_interval)
      self$calc_counts_bth(i_interval = i_interval)
      self$calc_exposure(i_interval = i_interval)
      self$update_counts(i_interval = i_interval)
  }
  else if (self$is_parallelogram_second[[i_interval]]) {
      NULL
  }
  else {
      self$calc_exposure_approx1()
      self$calc_counts_dth(i_interval = i_interval)
      self$draw_stock_end_net_mig(i_interval = i_interval,
                                  is_forecast = is_forecast)
      self$calc_exposure_approx2(i_interval = i_interval)
      self$draw_gross_mig(i_interval = i_interval)
      self$draw_counts_im_em(i_interval = i_interval)
      self$calc_counts_bth(i_interval = i_interval)
      self$calc_exposure(i_interval = i_interval)
      self$update_counts(i_interval = i_interval)
  }
}

PFilterNoReg$set(
  which = "public",
  name = "draw_values",
  value = draw_values_noreg
)


## HAS_TESTS
draw_values_withreg <- function(i_interval, is_forecast) {
  self$calc_stock_start(i_interval = i_interval)
  self$calc_rates(i_interval = i_interval)
  self$calc_exposure_approx1()
  if (is_forecast) {
    self$draw_counts_dth()
  } else {
    self$calc_counts_dth(i_interval = i_interval)
  }
  self$draw_stock_end_net_mig(
    i_interval = i_interval,
    is_forecast = is_forecast
  )
  self$calc_exposure_approx2(i_interval = i_interval)
  self$draw_counts_in_out() ## behaviour depends in 'is_split_mig'
  self$draw_gross_mig(i_interval = i_interval)     ## behaviour depends in 'is_split_mig'
  self$draw_counts_im_em(i_interval = i_interval)  ## behaviour depends in 'is_split_mig'
  if (is_forecast) {
    self$draw_counts_bth()
  } else {
    self$calc_counts_bth(i_interval = i_interval)
  }
  self$calc_exposure(i_interval = i_interval)
  self$update_counts(i_interval = i_interval)
}

PFilterWithReg$set(
  which = "public",
  name = "draw_values",
  value = draw_values_withreg
)


## draw_values_init ---------------------------------------------------------

## Draw value for initial stock, for existing cohorts.
## Value obtained by perturbing stock data. Note
## that 'x0' consists only of a value for stock
## - there are not births, deaths, or migration.
## We use the same method for no-region and
## with-region cases.

## HAS_TESTS
draw_values_init <- function() {
  self$draw_stock_end_init()
  self$update_counts_init()
}

PFilter$set(
  which = "public",
  name = "draw_values_init",
  value = draw_values_init
)


## make_diagnostics -----------------------------------------------------------

## Forecast population, deaths, immigration, emigration, and (if dominant)
## births for a single cohort

forecast_cohort <- function() {
    n_particle <- self$n_particle
    n_interval <- self$n_interval
    for (i_interval in seq_len(n_interval)) {
        self$stock_start <- self$counts_stock[, i_interval]
        self$draw_rates(i_interval = i_interval)
        self$draw_counts_forecast()
        self$update_counts(i_interval = i_interval)
    }
    self$resampled[] <- FALSE
    self$ess <- n_particle
    self$draw_index_output()
}

PFilter$set(
  which = "public",
  name = "forecast_cohort",
  value = forecast_cohort
)


## make_diagnostics -----------------------------------------------------------

## Create data frame with diagnostics
## for one cohort. The data frame has variables
## cohort, sexgender, time, age, ess,
## resampled, and n_unique, and has n_interval+1
## rows. The same function is used in no-region
## and with-region models.

## HAS_TESTS
make_diagnostics <- function() {
    n_interval <- self$n_interval
    cohort <- rep(self$cohort, times = n_interval + 1L)
    sexgender <- rep(self$sexgender, times = n_interval + 1L)
    time <- self$time_levels_stock
    age_events <- self$age_levels_events
    log
    ess <- self$ess
    sum_loglik <- self$sum_loglik
    sum_logtrans <- c(0, self$sum_logtrans)
    sum_logimp  <- self$sum_logimp
    sum_logwt_unnorm <- self$sum_logwt_unnorm
    resampled <- self$resampled
    n_particle <- self$n_particle
    n_thin <- self$n_thin
    n_unique <- self$n_unique
    is_popn <- self$is_popn
    is_new_cohort <- (age_events[[1L]] == 0L) && !is_popn[[1L]]
    if (is_new_cohort) {
        age <- c(-1L, age_events)
        triangle <- c("<none>",
                      rep(c("Lower", "Upper"), length.out = n_interval))
    } else {
        age <- c(age_events[[1L]], age_events)
        triangle <- c("<none>",
                      rep(c("Upper", "Lower"), length.out = n_interval))
    }
    data.frame(
        cohort = cohort,
        sexgender = sexgender,
        time = time,
        age = age,
        triangle = triangle,
        sum_loglik = sum_loglik,
        sum_logtrans = sum_logtrans,
        sum_logimp = sum_logimp,
        sum_logwt_unnorm = sum_logwt_unnorm,
        ess = ess,
        resampled = resampled,
        n_particle = n_particle,
        n_thin = n_thin,
        n_unique = n_unique
    )
}

PFilter$set(
  which = "public",
  name = "make_diagnostics",
  value = make_diagnostics
)


## make_output_events --------------------------------------------------------

## Create list with metadata and (thinned) results for
## deaths, migration, and (if 'is_forecast' is TRUE,
## and this is the dominant sex/gender, births.
## In with-region models, where the sex/gender 'is_dominant'
## is TRUE and 'is_forecast' is TRUE,
## the return value is a list with elements
## "births", "deaths", "internal_in",
## "internal_out", "immigration1", "emigration1",
## "immigration2", and "emigration2". In no-region models
## elements "internal_in" and "internal_out" are omitted.
## In models where 'is_dominant' or 'is_forecast' is FALSE,
## data frame "births" is omitted.
## In with-region models, the metadata is a data frame
## with variables cohort, sexgender, time, age, region.
## In no-region models, the region column is omitted.
## The metadata is serialized, ie turned into a 'raw' vector.
## The counts data are double vectors. Each element is
## begins with integers giving the length of the meta
## and count vectors.

## HAS_TESTS
make_output_events_noreg <- function(is_forecast) {
    n_interval <- self$n_interval
    cohort <- rep(self$cohort, times = n_interval)
    sexgender <- rep(self$sexgender, times = n_interval)
    time <- self$time_levels_events
    age <- self$age_levels_events
    meta <- data.frame(
        cohort = cohort,
        sexgender = sexgender,
        time = time,
        age = age
    )
    meta <- serialize(meta, connection = NULL)
    n_meta <- length(meta)
    index_output <- self$index_output
    is_dominant <- self$is_dominant
    count_deaths <- as.double(self$counts_deaths[index_output, ])
    count_immigration1 <- as.double(self$counts_immigration1[index_output, ])
    count_emigration1 <- as.double(self$counts_emigration1[index_output, ])
    count_immigration2 <- as.double(self$counts_immigration2[index_output, ])
    count_emigration2 <- as.double(self$counts_emigration2[index_output, ])
    n_count <- length(count_deaths)
    ans <- list(
        deaths = list(n_meta, n_count, meta, count_deaths),
        immigration1 = list(n_meta, n_count, meta, count_immigration1),
        emigration1 = list(n_meta, n_count, meta, count_emigration1),
        immigration2 = list(n_meta, n_count, meta, count_immigration2),
        emigration2 = list(n_meta, n_count, meta, count_emigration2)
    )
    if (is_forecast && is_dominant) {
        count_births <- as.double(self$counts_births[index_output, ])
        births <- list(births = list(n_meta, n_count, meta, count_births))
        ans <- c(births, ans)
    }
    ans
}

PFilterNoReg$set(
  which = "public",
  name = "make_output_events",
  value = make_output_events_noreg
  )


## HAS_TESTS
make_output_events_withreg <- function(is_forecast) {
    n_interval <- self$n_interval
    n_region <- self$n_region
    cohort <- self$cohort
    sexgender <- self$sexgender
    time_levels <- self$time_levels_events
    age_levels <- self$age_levels_events
    region_levels <- self$region_levels
    cohort <- rep(cohort, times = n_region * n_interval)
    sexgender <- rep(sexgender, times = n_region * n_interval)
    time <- rep(time_levels, each = n_region)
    age <- rep(age_levels, each = n_region)
    region <- rep(region_levels, times = n_interval)
    meta <- data.frame(
        cohort = cohort,
        sexgender = sexgender,
        time = time,
        age = age,
        region = region
    )
    meta <- serialize(meta, connection = NULL)
    n_meta <- length(meta)
    index_output <- self$index_output
    is_dominant <- self$is_dominant
    count_deaths <- as.double(self$counts_deaths[index_output, , ])
    count_internal_in <- as.double(self$counts_internal_in[index_output, , ])
    count_internal_out <- as.double(self$counts_internal_out[index_output, , ])
    count_immigration1 <- as.double(self$counts_immigration1[index_output, , ])
    count_emigration1 <- as.double(self$counts_emigration1[index_output, , ])
    count_immigration2 <- as.double(self$counts_immigration2[index_output, , ])
    count_emigration2 <- as.double(self$counts_emigration2[index_output, , ])
    n_count <- length(count_deaths)
    ans <- list(
        deaths = list(n_meta, n_count, meta, count_deaths),
        internal_in = list(n_meta, n_count, meta, count_internal_in),
        internal_out = list(n_meta, n_count, meta, count_internal_out),
        immigration1 = list(n_meta, n_count, meta, count_immigration1),
        emigration1 = list(n_meta, n_count, meta, count_emigration1),
        immigration2 = list(n_meta, n_count, meta, count_immigration2),
        emigration2 = list(n_meta, n_count, meta, count_emigration2)
    )
    if (is_forecast && is_dominant) {
        count_births <- as.double(self$counts_births[index_output, , ])
        births <- list(births = list(n_meta, n_count, meta, count_births))
        ans <- c(births, ans)
    }
    ans
}

PFilterWithReg$set(
                   which = "public",
                   name = "make_output_events",
                   value = make_output_events_withreg
               )


## make_output_population -----------------------------------------------------

## Create list with metadata and (thinned) results for
## population. The results for 'stock' include
## population and accession, so we need to drop the
## accession (using the 'in_popn' variable.)
## In no-region models, the metadata is a
## data frame with variables cohort, sexgender,
## time, and age, and in with-region models, the
## metadata also has a region variable. In both cases
## the data frame is serialised, ie turned into a 'raw' vector.
## In both models, the population data is a double vector.
## The output starts with two integers, the first giving the length
## of the metadata object, and the second giving the length
## of the count vector.

## HAS_TESTS
make_output_population_noreg <- function() {
    is_popn <- self$is_popn
    n_popn <- self$n_popn
    cohort <- rep(self$cohort, times = n_popn)
    sexgender <- rep(self$sexgender, times = n_popn)
    time <- self$time_levels_stock[is_popn]
    age <- self$age_levels_stock[is_popn]
    meta <- data.frame(
        cohort = cohort,
        sexgender = sexgender,
        time = time,
        age = age
    )
    meta <- serialize(meta, connection = NULL)
    n_meta <- length(meta)
    index_output <- self$index_output
    count <- self$counts_stock[index_output, is_popn]
    count <- as.double(count)
    n_count <- length(count)
    list(n_meta, n_count, meta, count)
}

PFilterNoReg$set(
  which = "public",
  name = "make_output_population",
  value = make_output_population_noreg
)


## HAS_TESTS
make_output_population_withreg <- function() {
    n_region <- self$n_region
    cohort <- self$cohort
    sexgender <- self$sexgender
    is_popn <- self$is_popn
    n_popn <- self$n_popn
    time_levels <- self$time_levels_stock[is_popn]
    age_levels <- self$age_levels_stock[is_popn]
    region_levels <- self$region_levels
    cohort <- rep(cohort, times = n_region * n_popn)
    sexgender <- rep(sexgender, times = n_region * n_popn)
    time <- rep(time_levels, each = n_region)
    age <- rep(age_levels, each = n_region)
    region <- rep(region_levels, times = n_popn)
    meta <- data.frame(
        cohort = cohort,
        sexgender = sexgender,
        time = time,
        age = age,
        region = region
    )
    meta <- serialize(meta, connection = NULL)
    n_meta <- length(meta)
    index_output <- self$index_output
    count <- self$counts_stock[index_output, , is_popn]
    count <- as.double(count)
    n_count <- length(count)
    list(n_meta, n_count, meta, count)
}

PFilterWithReg$set(
  which = "public",
  name = "make_output_population",
  value = make_output_population_withreg
)


## resample -------------------------------------------------------------------

## If effective sample size for interval 'i_interval',
## falls below a limit set by 'threshold',
## then resample, ie sample (with replacement)
## a subset of particles, with probability proportional
## to their weights. For efficiency, we record the
## selection via values for 'index_parent'. The elements
## of 'index_parent' show which elements of
## 'counts_stock' should be used for 'stock_start'
## in the first interval. The sampling
## procedure is identical for no-region and with-region
## models so this is a method for the PFilter superclass.

## HAS_TESTS
resample <- function(i_interval) {
    n_particle <- self$n_particle
    is_parallel_second <- (i_interval >= 1L) && self$is_parallelogram_second[[i_interval]]
    if (is_parallel_second) {
        do_resampling <- FALSE
        ess <- self$ess[i_interval]
    }
    else {
        n_interval <- self$n_interval
        threshold <- self$threshold
        logwt_unnorm <- self$logwt_unnorm[, i_interval + 1L]
        logwt_unnorm[is.infinite(logwt_unnorm)] <- -Inf ## hacky way of avoiding edge-case problems
        wt <- softmax(logwt_unnorm)
        ess <- 1 / sum(wt^2)
        ## use <= so that 'threshold = 1'
        ## implies reampling every interval
        do_resampling <- (ess / n_particle) <= threshold
        if (is.na(do_resampling)) {
            loglik <- self$calc_loglik(i_interval - 1L)
            logtrans <- self$calc_logtrans(i_interval = i_interval)
            logimp <- self$calc_logimp(is_forecast = FALSE)
            stop(sprintf("%d elements of loglik are NA\n", sum(is.na(loglik))),
                 sprintf("%d elements of loglik are infinite\n", sum(is.infinite(loglik))),
                 sprintf("%d elements of logtrans are NA\n", sum(is.na(logtrans))),
                 sprintf("%d elements of logtrans are infinite\n", sum(is.infinite(logtrans))),
                 sprintf("%d elements of logimp are NA\n", sum(is.na(logimp))),
                 sprintf("%d elements of logimp are infinite\n", sum(is.infinite(logimp))))
        }
    }
    if (do_resampling) { 
        index_parent <- draw_index_parent(wt)
    } else {
        index_parent <- seq_len(n_particle)
    }
    self$index_parent[, i_interval + 1L] <- index_parent
    self$resampled[i_interval + 1L] <- (is_parallel_second || do_resampling)
    self$ess[i_interval + 1L] <- ess
}


PFilter$set(
  which = "public",
  name = "resample",
  value = resample
)


## resample_init --------------------------------------------------------------

## If effective sample size for the initial population
## falls below a limit set by
## 'threshold', then resample, ie sample (with replacement)
## a subset of particles, with probability proportional
## to their weights. For efficiency, we record the
## selection via values for 'index_parent'. The elements
## of 'index_parent' show which elements of
## 'counts_stock' should be used for 'stock_start'
## in the first interval. The sampling
## procedure is identical for no-region and with-region
## models so this is a method for the PFilter superclass.

## HAS_TESTS
resample_init <- function() {
  self$resample(i_interval = 0L) ## '-1' in C
}

PFilter$set(
  which = "public",
  name = "resample_init",
  value = resample_init
)


## run ------------------------------------------------------------------------

## Generate values for migration, population, and
## (in forecasts only) births and deaths, plus
## associated weights, for a cohort.

## HAS_TESTS
run <- function(is_forecast) {
  n_interval <- self$n_interval
  has_stock_init <- self$has_stock_init
  if (has_stock_init) {
    self$skip_sampling_init()
  } else {
    self$draw_values_init()
    self$calc_logwt_unnorm_init()
    self$resample_init()
  }
  for (i_interval in seq_len(n_interval)) {
    self$draw_values(
      i_interval = i_interval,
      is_forecast = is_forecast
    )
    self$calc_logwt_unnorm(
      i_interval = i_interval,
      is_forecast = is_forecast
    )
    self$resample(i_interval = i_interval)
  }
  self$calc_index_ancestor()
  self$calc_trajectories()
  self$draw_index_output()
  self$calc_n_unique()
}

PFilter$set(
  which = "public",
  name = "run",
  value = run
)


## skip_sampling_init ---------------------------------------------------------

## Set values for 'logwt_norm', 'index_parent',
## 'resampled' and 'ess' to suitable
## values for the case where 'have_stock_init'
## is TRUE, so that no weighting or resampling is done.

## HAS_TESTS
skip_sampling_init <- function() {
  n_particle <- self$n_particle
  self$logwt_unnorm[, 1L] <- rep(-1 * log(n_particle), times = n_particle)
  self$index_parent[, 1L] <- seq_len(n_particle)
  self$resampled[1L] <- FALSE
  self$ess[1L] <- as.numeric(n_particle)
}


PFilter$set(
  which = "public",
  name = "skip_sampling_init",
  value = skip_sampling_init
)


## update_counts --------------------------------------------------------------

## HAS_TESTS
update_counts_noreg <- function(i_interval) {
    update_counts_noreg_inner(self, i_interval - 1L); ## expects C-style index
}

PFilterNoReg$set(
  which = "public",
  name = "update_counts",
  value = update_counts_noreg
)


## HAS_TESTS
update_counts_withreg <- function(i_interval) {
    update_counts_withreg_inner(self, i_interval - 1L); ## expects C-style index
}


PFilterWithReg$set(
  which = "public",
  name = "update_counts",
  value = update_counts_withreg
)


## update_counts_init ---------------------------------------------------------

## HAS_TESTS
update_counts_init_noreg <- function() {
  self$counts_stock[, 1L] <- self$stock_end ## when translating to C, change 1 to 0
}

PFilterNoReg$set(
  which = "public",
  name = "update_counts_init",
  value = update_counts_init_noreg
)


## HAS_TESTS
update_counts_init_withreg <- function() {
  self$counts_stock[, , 1L] <- self$stock_end ## when translating to C, change 1 to 0
}

PFilterWithReg$set(
  which = "public",
  name = "update_counts_init",
  value = update_counts_init_withreg
)


## write_results --------------------------------------------------------------

## Extract results and write to 'work_dir'

## HAS_TESTS
write_results <- function(is_forecast, work_dir) {
    cohort <- self$cohort
    sexgender <- self$sexgender
    ## outputs
    output_population <- self$make_output_population()
    output_events <- self$make_output_events(is_forecast = is_forecast)
    outputs <- c(list(population = output_population),
                 output_events)
    names(outputs) <- paste(names(outputs), cohort, sexgender, sep = "-")
    write_outputs(
        outputs = outputs,
        work_dir = work_dir
    )
    ## diagnostics
    diagnostics <- self$make_diagnostics()
    write_diagnostics(
        diagnostics = diagnostics,
        cohort = cohort,
        sexgender = sexgender,
        work_dir = work_dir
    )
}

PFilter$set(
  which = "public",
  name = "write_results",
  value = write_results
)

