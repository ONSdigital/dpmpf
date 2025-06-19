
## Functions to create simulated PFilter objects, for use in testing

## HAS_TESTS
#' Make arguments required by simulated_pf_existing_* functions
#' @noRd
make_args_simulated_existing <- function() {
  list(
    n_interval = 4L,
    time_levels_stock = c(2000L, 2000L, 2001L, 2001L, 2002L),
    time_levels_events = c(2000L, 2000L, 2001L, 2001L),
    age_levels_stock = c(20L, 21L, 21L, 22L, 22L),
    age_levels_events = c(20L, 21L, 21L, 22L),
    is_popn = c(TRUE, FALSE, TRUE, FALSE, TRUE),
    n_popn = 3L
  )
}


## HAS_TESTS
#' Make arguments required by simulated_pf_expiring_* functions
#' @noRd
make_args_simulated_expiring <- function() {
  list(
    n_interval = 1L,
    time_levels_stock = c(2000L, 2001L),
    time_levels_events = 2001L,
    age_levels_stock = c(99L, 100L),
    age_levels_events = 99L,
    is_popn = c(TRUE, FALSE),
    n_popn = 1L
  )
}


## HAS_TESTS
#' Make arguments required by simulated_pf_new_* functions
#' @noRd
make_args_simulated_new <- function() {
  list(
    n_interval = 3L,
    time_levels_stock = c(2000L, 2000L, 2001L, 201L),
    time_levels_events = c(2000L, 2001L, 2001L),
    age_levels_stock = c(0L, 0L, 1L, 1L),
    age_levels_events = c(0L, 0L, 1L),
    is_popn = c(FALSE, TRUE, FALSE, TRUE),
    n_popn = 2L
  )
}


## HAS_TESTS
#' Make a simulated inputs data frame for a particle filter with no regions
#'
#' @param has_stock_init Whether the cohort has a fixed initial value for stock
#' @param is_dominant Whether the cohort belongs to the 'dominant' sex/gender
#'   (ie the sex/gender whose population is used to calculate exposures
#'   in models of births)
#' @param has_one_imem Whether there are two immigration and two emigration
#'   flows, or just one of each.
#' @param n_particle Number of particles
#' @param n_thin Final sample consists of n_particle / n_thin particles
#' @param is_forecast Whether the particle filter is for a forecast
#' @param parallelogram Whether to use the parallelogram importance function.
#' @param cohort_type Whether the cohort (i) already existed at the start of the
#'    estimation/forecasting period (and survives past the end), (ii)
#'    is born during the estimation/forecasting period, or (iii)
#'    already existing at the start of the estimation/forecasting
#'    period and reaches the maximum age during the period
#'
#' @return A data frame
#'
#' @noRd
simulated_df_noreg <- function(has_stock_init = FALSE,
                               is_dominant = TRUE,
                               has_one_imem = FALSE,
                               n_particle = 5L,
                               n_thin = 1L,
                               is_forecast = FALSE,
                               parallelogram = FALSE,
                               cohort_type = c("existing", "new", "expiring")) {
    cohort_type <- match.arg(cohort_type)
    args_cohort <- switch(cohort_type,
                          existing = make_args_simulated_existing(),
                          new = make_args_simulated_new(),
                          expiring = make_args_simulated_expiring()
                          )
    n_interval <- args_cohort$n_interval
    if (has_stock_init) {
        counts_stock <- as.double(stats::rpois(n = n_particle, lambda = 10))
    } else {
        counts_stock <- NULL
    }
    if (is_forecast) {
        counts_births <- NULL
        counts_deaths <- NULL
    } else {
        counts_births <- as.double(stats::rpois(n_interval, lambda = 5))
        counts_deaths <- as.double(stats::rpois(n_interval, lambda = 6))
    }
    if (is_dominant) {
        rates_births <- stats::rgamma(n = n_interval, shape = 1, rate = 3)
    } else {
        rates_births <- rep(0, times = n_interval)
    }
    rates_deaths <- stats::rgamma(n = n_interval, shape = 1, rate = 3)
    rates_immigration1 <- stats::rgamma(n = n_interval, shape = 1, rate = 0.2)
    rates_emigration1 <- stats::rgamma(n = n_interval, shape = 1, rate = 3)
    if (has_one_imem) {
        rates_immigration2 <- rep(0, times = n_interval)
        rates_emigration2 <- rep(0, times = n_interval)
    } else {
        rates_immigration2 <- stats::rgamma(n = n_interval, shape = 1, rate = 0.2)
        rates_emigration2 <- stats::rgamma(n = n_interval, shape = 1, rate = 3)
    }
    counts_data_stock_1 <- as.double(stats::rpois(n = n_interval + 1L, lambda = 20))
    counts_data_stock_2 <- as.double(stats::rpois(n = n_interval + 1L, lambda = 20))
    ratio_stock_2 <- rep(1, times = n_interval + 1L)
    disp_stock_2 <- rep(1.2, times = n_interval + 1L)
    cdm_stock_1 <- new_CdmNoregPoibin(
        counts_data = counts_data_stock_1,
        prob = 0.95
    )
    cdm_stock_2 <- new_CdmNoregNbinom(
        counts_data = counts_data_stock_2,
        ratio = ratio_stock_2,
        disp = disp_stock_2
    )
    cdms_stock <- new_CdmsNoreg(list(cdm_stock_1, cdm_stock_2))
    cdms_immigration1 <- new_CdmsNoreg()
    cdms_emigration1 <- new_CdmsNoreg()
    cdms_immigration2 <- new_CdmsNoreg()
    cdms_emigration2 <- new_CdmsNoreg()
    df <- list(
        cohort = 2000L,
        sexgender = "Female",
        has_stock_init = has_stock_init,
        is_dominant = is_dominant,
        n_particle = n_particle,
        n_thin = n_thin,
        n_interval = args_cohort$n_interval
    )
    df$has_one_imem <- has_one_imem
    df$obs_zero <- 0.5
    df$is_popn <- args_cohort[["is_popn"]]
    df$n_popn <- args_cohort$n_popn
    make_parallelogram_first <- function(is_popn) {
        n <- length(is_popn)
        if (n == 2L)
            FALSE
        else {
            is_last <- seq_len(n - 1L) == n - 1L
            is_popn[-n] & !is_last
        }
    }
    make_parallelogram_second <- function(is_popn) {
        n <- length(is_popn)
        if (n == 2L)
            FALSE
        else {
            is_first <- seq_len(n - 1L) == 1L
            !is_popn[-n] & !is_first
        }
    }
    if (parallelogram) {
        df$is_parallelogram_first <- make_parallelogram_first(df$is_popn)
        df$is_parallelogram_second <- make_parallelogram_second(df$is_popn)
    }
    else {
        df$is_parallelogram_first <- rep(FALSE, times = length(df$is_popn) - 1L)
        df$is_parallelogram_second <- rep(FALSE, times = length(df$is_popn) - 1L)
    }      
    df$time_levels_stock <- args_cohort[["time_levels_stock"]]
    df$time_levels_events <- args_cohort[["time_levels_events"]]
    df$age_levels_stock <- args_cohort[["age_levels_stock"]]
    df$age_levels_events <- args_cohort[["age_levels_events"]]
    df$counts_stock <- counts_stock
    df$counts_births <- counts_births
    df$counts_deaths <- counts_deaths
    df$rates_births <- rates_births
    df$rates_deaths <- rates_deaths
    df$rates_immigration1 <- rates_immigration1
    df$rates_emigration1 <- rates_emigration1
    df$rates_immigration2 <- rates_immigration2
    df$rates_emigration2 <- rates_emigration2
    df$disp_births <- 0
    df$disp_deaths <- 0
    df$disp_immigration1 <- 0
    df$disp_emigration1 <- 0
    df$disp_immigration2 <- 0
    df$disp_emigration2 <- 0
    df$cdms_stock <- cdms_stock
    df$cdms_immigration1 <- cdms_immigration1
    df$cdms_emigration1 <- cdms_emigration1
    df$cdms_immigration2 <- cdms_immigration2
    df$cdms_emigration2 <- cdms_emigration2
    df
}


## HAS_TESTS
#' Make a simulated input dataset for a
#' particle filter for a model with regions
#'
#' @inheritParams simulated_df_noreg
#' @param n_region Number of regions
#' @param is_split_mig Whether to split migration into separate flows
#'
#' @return An object of class "PFilterWithReg"
#'
#' @noRd
simulated_df_withreg <- function(has_stock_init = FALSE,
                                 is_dominant = TRUE,
                                 n_region = 3L,
                                 has_one_imem = FALSE,
                                 is_split_mig = TRUE,
                                 n_particle = 5L,
                                 n_thin = 1L,
                                 is_forecast = FALSE,
                                 cohort_type = c("existing", "new", "expiring")) {
    cohort_type <- match.arg(cohort_type)
    args_cohort <- switch(cohort_type,
                          existing = make_args_simulated_existing(),
                          new = make_args_simulated_new(),
                          expiring = make_args_simulated_expiring()
                          )
    n_interval <- args_cohort$n_interval
    if (has_stock_init) {
        counts_stock <- matrix(as.double(stats::rpois(n = n_particle * n_region, lambda = 10)),
                               nrow = n_particle, ncol = n_region
                               )
    } else {
        counts_stock <- NULL
    }
    if (is_forecast) {
        counts_births <- NULL
        counts_deaths <- NULL
    } else {
        counts_births <- matrix(as.double(stats::rpois(n = n_region * n_interval, lambda = 5)),
                                nrow = n_region, ncol = n_interval
                                )
        counts_deaths <- matrix(as.double(stats::rpois(n = n_region * n_interval, lambda = 6)),
                                nrow = n_region, ncol = n_interval
                                )
    }
    if (is_dominant) {
        rates_births <- stats::rgamma(
                                   n = n_region * n_interval,
                                   shape = 1,
                                   rate = 3
                               )
    } else {
        rates_births <- rep(0, times = n_region * n_interval)
    }
    rates_deaths <- stats::rgamma(
                               n = n_region * n_interval,
                               shape = 1,
                               rate = 3
                           )
    rates_internal_in <- stats::rgamma(
                                    n = n_region * n_interval,
                                    shape = 1,
                                    rate = 0.2
                                )
    rates_internal_out <- stats::rgamma(
                                     n = n_region * n_interval,
                                     shape = 1,
                                     rate = 3
                                 )
    rates_immigration1 <- stats::rgamma(
                                     n = n_region * n_interval,
                                     shape = 1,
                                     rate = 0.2
                                 )
    rates_emigration1 <- stats::rgamma(
                                    n = n_region * n_interval,
                                    shape = 1,
                                    rate = 3
                                )
    rates_immigration2 <- stats::rgamma(
                                     n = n_region * n_interval,
                                     shape = 1,
                                     rate = 0.2
                                 )
    rates_emigration2 <- stats::rgamma(
                                    n = n_region * n_interval,
                                    shape = 1,
                                    rate = 3
                                )
    rates_births <- matrix(rates_births, nrow = n_region, ncol = n_interval)
    rates_deaths <- matrix(rates_deaths, nrow = n_region, ncol = n_interval)
    rates_internal_in <- matrix(rates_internal_in, nrow = n_region, ncol = n_interval)
    rates_internal_out <- matrix(rates_internal_out, nrow = n_region, ncol = n_interval)
    rates_immigration1 <- matrix(rates_immigration1, nrow = n_region, ncol = n_interval)
    rates_emigration1 <- matrix(rates_emigration1, nrow = n_region, ncol = n_interval)
    if (has_one_imem) {
        rates_immigration2 <- matrix(0, nrow = n_region, ncol = n_interval)
        rates_emigration2 <- matrix(0, nrow = n_region, ncol = n_interval)
    } else {
        rates_immigration2 <- matrix(rates_immigration2, nrow = n_region, ncol = n_interval)
        rates_emigration2 <- matrix(rates_emigration2, nrow = n_region, ncol = n_interval)
    }
    counts_data_stock_1 <- as.double(stats::rpois(
                                                 n = n_region * (n_interval + 1L),
                                                 lambda = 20
                                             ))
    counts_data_stock_2 <- as.double(stats::rpois(
                                                 n = n_region * (n_interval + 1L),
                                                 lambda = 20
                                             ))
    ratio_stock_2 <- rep(1, times = n_region * (n_interval + 1L))
    disp_stock_2 <- rep(1.2, times = n_region * (n_interval + 1L))
    counts_data_stock_1 <- matrix(counts_data_stock_1,
                                  nrow = n_region
                                  )
    counts_data_stock_2 <- matrix(counts_data_stock_2,
                                  nrow = n_region
                                  )
    ratio_stock_2 <- matrix(ratio_stock_2,
                            nrow = n_region
                            )
    disp_stock_2 <- matrix(disp_stock_2,
                           nrow = n_region
                           )
    cdm_stock_1 <- new_CdmWithregPoibin(
        counts_data = counts_data_stock_1,
        prob = 0.95
    )
    cdm_stock_2 <- new_CdmWithregNbinom(
        counts_data = counts_data_stock_2,
        ratio = ratio_stock_2,
        disp = disp_stock_2
    )
    cdms_stock <- new_CdmsWithreg(list(cdm_stock_1, cdm_stock_2))
    cdms_internal_in <- new_CdmsWithreg()
    cdms_internal_out <- new_CdmsWithreg()
    cdms_immigration1 <- new_CdmsWithreg()
    cdms_emigration1 <- new_CdmsWithreg()
    cdms_immigration2 <- new_CdmsWithreg()
    cdms_emigration2 <- new_CdmsWithreg()
    df <- list(
        cohort = 2000L,
        sexgender = "Female",
        has_stock_init = has_stock_init,
        is_dominant = is_dominant,
        n_particle = n_particle,
        n_thin = n_thin,
        n_region = n_region,
        n_interval = args_cohort$n_interval
    )
    df$obs_zero <- 0.5
    df$has_one_imem <- has_one_imem
    df$is_split_mig <- is_split_mig
    df$is_popn <- args_cohort[["is_popn"]]
    df$n_popn <- args_cohort$n_popn
    df$time_levels_stock <- args_cohort[["time_levels_stock"]]
    df$time_levels_events <- args_cohort[["time_levels_events"]]
    df$age_levels_stock <- args_cohort[["age_levels_stock"]]
    df$age_levels_events <- args_cohort[["age_levels_events"]]
    df$region_levels <- seq_len(n_region)
    df$counts_stock <- counts_stock
    df$counts_births <- counts_births
    df$counts_deaths <- counts_deaths
    df$rates_births <- rates_births
    df$rates_deaths <- rates_deaths
    df$rates_internal_in <- rates_internal_in
    df$rates_internal_out <- rates_internal_out
    df$rates_immigration1 <- rates_immigration1
    df$rates_emigration1 <- rates_emigration1
    df$rates_immigration2 <- rates_immigration2
    df$rates_emigration2 <- rates_emigration2
    df$disp_births <- 0
    df$disp_deaths <- 0
    df$disp_immigration1 <- 0
    df$disp_emigration1 <- 0
    df$disp_immigration2 <- 0
    df$disp_emigration2 <- 0
    df$cdms_stock <- cdms_stock
    df$cdms_internal_in <- cdms_internal_in
    df$cdms_internal_out <- cdms_internal_out
    df$cdms_immigration1 <- cdms_immigration1
    df$cdms_emigration1 <- cdms_emigration1
    df$cdms_immigration2 <- cdms_immigration2
    df$cdms_emigration2 <- cdms_emigration2
    df
}


## HAS_TESTS
#' Make a simulated particle filter
#' for a model with no regions
#'
#' @inheritParams simulate_df_noreg
#'
#' @return An object of class "PFilterNoReg"
#'
#' @noRd
simulated_pf_noreg <- function(has_stock_init = FALSE,
                               is_dominant = TRUE,
                               has_one_imem = FALSE,
                               n_particle = 5L,
                               n_thin = 1L,
                               threshold = 0.5,
                               is_forecast = FALSE,
                               parallelogram = FALSE,
                               cohort_type = c("existing", "new", "expiring")) {
  df <- simulated_df_noreg(
    has_stock_init = has_stock_init,
    is_dominant = is_dominant,
    has_one_imem = has_one_imem,
    n_particle = n_particle,
    n_thin = n_thin,
    is_forecast = is_forecast,
    parallelogram = parallelogram,
    cohort_type = cohort_type
  )
  make_pf(
    df = df,
    threshold = threshold,
    is_forecast = is_forecast
  )
}


## HAS_TESTS
#' Make a simulated particle filter for a model with regions
#'
#' @inheritParams simulated_df_noreg
#' @param n_region Number of regions
#' @param is_split_mig Whether to separate migration flows.
#'
#' @return An object of class "PFilterWithReg"
#'
#' @noRd
simulated_pf_withreg <- function(has_stock_init = FALSE,
                                 is_dominant = TRUE,
                                 has_one_imem = FALSE,
                                 is_split_mig = TRUE,
                                 n_region = 3L,
                                 n_particle = 5L,
                                 n_thin = 1L,
                                 threshold = 0.5,
                                 is_forecast = FALSE,
                                 cohort_type = c("existing", "new", "expiring")) {
  df <- simulated_df_withreg(
    has_stock_init = has_stock_init,
    is_dominant = is_dominant,
    has_one_imem = has_one_imem,
    is_split_mig = is_split_mig,
    n_region = n_region,
    n_particle = n_particle,
    n_thin = n_thin,
    is_forecast = is_forecast,
    cohort_type = cohort_type
  )
  make_pf(
    df = df,
    threshold = threshold,
    is_forecast = is_forecast
  )
}
