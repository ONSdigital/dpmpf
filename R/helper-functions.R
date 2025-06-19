## HAS_TESTS
draw_gross_mig_inner <- function(rates_im, rates_em, exposure, net_mig) {
    n_particle <- length(rates_im)
    lambda <- rates_im + rates_em * exposure
    lower <- abs(net_mig)
    gross_mig <- rpoistr(
        n = n_particle,
        lambda = lambda,
        lower = lower
    )
    is_diff_parity <- (gross_mig %% 2L) != (net_mig %% 2L)
    gross_mig[is_diff_parity] <- gross_mig[is_diff_parity] + 1
    logprob_high <- dpoistr(gross_mig,
                            lambda = lambda,
                            lower = lower,
                            use_log = TRUE
                            )
    logprob_low <- dpoistr(gross_mig - 1,
                           lambda = lambda,
                           lower = lower,
                           use_log = TRUE
                           )
    logprob_low[gross_mig == lower] <- -Inf
    logimp_gross_mig = log_sum_exp_2(x = logprob_high,
                                     y = logprob_low)
    list(gross_mig,
         logimp_gross_mig)
}



## HAS_TESTS
## Given a set of normalised weights, randomly select
## a new set of particles, using the "Systematic Resampling"
## algorithm, as described on page 13 of Arnaud Doucet and
## Adam M. Johansen. 2008. A Tutorial on Particle Filtering
## and Smoothing: Fifteen years later. The return value is
## an integer vector giving the index of the
## sampled particles.
draw_index_parent <- function(wt) {
  n <- length(wt)
  vec <- c(0, cumsum(wt))
  u <- stats::runif(
    n = 1L,
    min = 0,
    max = 1 / n
  )
  x <- seq(
    from = u,
    by = 1 / n,
    length.out = n
  )
  ## http://adv-r.had.co.nz/Rcpp.html discusses how to write findInterval using Rcpp,
  ## or there is a function defined somewhere in the C source code for R
  findInterval(
    x = x,
    vec = vec
  )
}


## HAS_TESTS
#' Generate demograhic rates
#'
#' Generate demographic rates by drawing from a
#' gamma distribution with a mean-dispersion
#' parameterisation.
#'
#' @param n Number of draws.
#' @param mean, disp Mean and dispersion of
#' gamma distribution. Both scalars.
#'
#' @return A vector of length n.
#'
#' @noRd
draw_rates_inner <- function(n, mean, disp) {
    eps <- 1e-10
    if (disp < eps) {
        ans <- rep(mean, times = n)
    }
    else {
        if (mean < eps) {
            ans <- rep(0, times = n)
        }
        else {
      shape <- 1 / disp # mean invariant parameterisation
      rate <- 1 / (disp*mean)
      ans <- stats::rgamma(n = n,
                           shape = shape,
                           rate = rate)
        }
    }
    ans
}


## HAS_TESTS
## Find out if inputs data frame 'df' describes
## a cohort with regions
## 'df' - A data frame
## return value - TRUE or FALSE
has_region <- function(df) {
  "n_region" %in% names(df)
}


## HAS_TESTS
## Create a particle filter. If 'df' contains regional data,
## then the particle filter has C++ class "PfWithreg"; if it
## does not, then the particle filter has C++ class "PfNoreg".
## 'df' - A data frame with input data for the cohort being modelled
## 'threshold' - Real value between 0 and 1.
## return value - An object of class "PfWithreg" or "PfNoreg"
make_pf <- function(df, n_particle, threshold, is_forecast) {
    has_region <- has_region(df)
    if (has_region)
        PFilterWithReg$new(
                           df = df,
                           threshold = threshold,
                           is_forecast = is_forecast
                       )
    else
        PFilterNoReg$new(
                         df = df,
                         threshold = threshold,
                         is_forecast = is_forecast
                     )
}


## HAS_TESTS
logprob_trans <- function(stock_start,
                          counts_bth,
                          counts_dth,
                          counts_im,
                          counts_em,
                          rates_bth,
                          rates_dth,
                          rates_im,
                          rates_em,
                          exposure,
                          is_dominant) {
    n_particle <- length(stock_start)
    ans_im <- stats::dpois(x = counts_im,
                           lambda = rates_im,
                           log = TRUE)
    ans_dth_em <- numeric(length = n_particle)
    prob_exit <- 1 - exp(-0.5 * (rates_dth + rates_em))
    prob_dth_exit <- rates_dth / (rates_dth + rates_em)
    is_im_even <- (counts_im %% 2L) == 0L
    is_im_odd <- !is_im_even
    size_mid <- stock_start + 0.5 * counts_im
    size_low <- size_mid - 0.5
    size_high <- size_mid + 0.5
    counts_exit <- counts_dth + counts_em
    is_even_valid <- is_im_even & (counts_exit <= size_mid)
    is_odd_valid <- is_im_odd & (counts_exit <= size_low)
    is_valid <- is_even_valid | is_odd_valid
    logprob_exit <- numeric(length = n_particle)
    logprob_dth_exit <- numeric(length = n_particle)
    logprob_exit[!is_valid] <- -1e300
    logprob_dth_exit[!is_valid] <- -1e300
    logprob_exit[is_even_valid] <- stats::dbinom(counts_exit[is_even_valid],
                                                 prob = prob_exit[is_even_valid],
                                                 size = size_mid[is_even_valid],
                                                 log = TRUE)
    logprob_exit[is_odd_valid] <- log_sum_exp_2(stats::dbinom(counts_exit[is_odd_valid],
                                                              prob = prob_exit[is_odd_valid],
                                                              size = size_low[is_odd_valid],
                                                              log = TRUE),
                                                stats::dbinom(counts_exit[is_odd_valid],
                                                              prob = prob_exit[is_odd_valid],
                                                              size = size_high[is_odd_valid],
                                                              log = TRUE)) - log(2)
    logprob_dth_exit[is_valid] <- stats::dbinom(counts_dth[is_valid],
                                                prob = prob_dth_exit[is_valid],
                                                size = counts_exit[is_valid],
                                                log = TRUE)
    ans_dth_em <- logprob_exit + logprob_dth_exit
    ans <- ans_im + ans_dth_em
    if (is_dominant) {
        lambda <- rates_bth * exposure
        ans_bth <- stats::dpois(x = counts_bth,
                              lambda = lambda,
                              log = TRUE)
        ans_bth[(counts_bth > 0) & (lambda == 0)] <- -1e300
        ans <- ans + ans_bth
    }
    ans[is.infinite(ans) | is.na(ans)] <- -1e300
    ans
}



## HAS_TESTS
## Message to display at start of 'create_run_destroy_pfilter'
show_start_message <- function(df) {
    str <- sprintf("  %s  starting cohort :  %4s %7s",
                   strftime(Sys.time(), format = "%H:%M:%S"),
                   df$cohort,
                   df$sexgender)
    message(str)
}


## HAS_TESTS
## Calculate softmax. Assume x is valid
## numeric vector with no NAs and
## length of at least 1. Subtracting
## max(x) helps avoid overflow.
softmax <- function(x) {
  x <- x - max(x)
  x <- exp(x)
  x / sum(x)
}


## HAS_TESTS
## Sum 'x' within groups defined by 'g'.
## 'x' - numeric vector, with length greater
##    than or equal to 1
## 'g' - integer vector, same length as 'x'
## return value - numeric vector, with length
##   equal to the number of unique elements
##   in 'g'
sum_by <- function(x, g) {
  n_x <- length(x)
  g_unique <- unique(g)
  n_g_unique <- length(g_unique)
  ans <- rep.int(0, times = n_g_unique)
  for (i in seq_len(n_x)) {
    j <- match(g[i], g_unique)
    ans[j] <- ans[j] + x[i]
  }
  ans
}


#' HAS_TESTS
#' Write a data frame of diagnostics into the 'work_dir'
#' as a RDS file.
#'
#' @param diagnostics A of data frame.
#' @param work_dir A path to a directory.
#'
#' @return TRUE if the function completes successfully.
#'
#' @noRd
write_diagnostics <- function(diagnostics, cohort, sexgender, work_dir) {
    basename <- paste0("tmp-diagnostics-", cohort, "-", sexgender, ".rds")
    file <- file.path(work_dir, basename)
    saveRDS(diagnostics, file = file, compress = FALSE) ## 'compress = FALSE' for speed
    invisible(TRUE)
}


#' HAS_TESTS
#' Given a set of output list, and the name
#' of a directory, write the data frames into
#' the directory as binary files.
#'
#' @param outputs A named list of outputs. The names
#' are used to create the file names.
#' @param work_dir A path to a directory.
#'
#' @return TRUE if the function completes successfully.
#'
#' @noRd
write_outputs <- function(outputs, work_dir) {
  nms_outputs <- names(outputs)
  for (i_output in seq_along(outputs)) {
    output <- outputs[[i_output]]
    nm_output <- nms_outputs[[i_output]]
    basename <- paste0("tmp-", nm_output, ".bin")
    file <- file.path(work_dir, basename)
    con <- file(file, "wb")
    for (i_item in seq_along(output))
        writeBin(output[[i_item]], con = con)
    close(con)
  }
  invisible(TRUE)
}
