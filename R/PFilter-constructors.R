
## globals.R needs to be run to avoid warnings about
## global variables 'self' and 'super'
## We need functions 'new_CdmsNoreg' and 'new_CdmsWithreg',
## which are in file 'cdms-constructors.R',
## for the code in this file to run. @include
## makes sure that the other files are run first.
#' @include globals.R cdms-constructors.R
NULL

## In this file we create the class objects first, and
## then add initialization methods later, using 'set'.
## In file 'PFilter-methods' we add further methods
## using 'set'.

## Incantation needed to prevent note
## "Namespace in Imports field not imported from: ‘R6’":
#' @import R6
NULL


## create class objects -------------------------------------------------------

## Superclass for PFilterNoReg and PFilterWithReg
PFilter <-
  R6::R6Class("PFilter",
    public = list(
      cohort = character(1L), # int scalar - year cohort born
      sexgender = character(1L), # string - sex/gender of cohort
      has_stock_init = logical(1L), # lgl scalar - whether has values for initial stock
      has_one_imem = logical(1L), # lgl scalar - whether there is only one set of immig/emig
      is_dominant = logical(1L), # lgl scalar - whether births attributed to this sex/gender
      obs_zero = numeric(1L), # dbl scalar - parameter for data models
      n_interval = integer(1L), # int scalar - number of intervals being estimated
      n_particle = integer(1L), # int scalar - number of sample paths
      n_thin = integer(1L), # int scalar - keep 1/thin particles
      is_popn = logical(), # lgl vaector n_interval+1 - whether stock refers to popn (not accession)
      n_popn = integer(1L), # int scalar - sum(is_popn)
      time_levels_stock = integer(), # int vector n__interval+1 - times for stock counts
      time_levels_events = integer(), # int vector n_interval - times for event counts
      age_levels_stock = integer(), # int vector n_interval+1 - ages for stock counts
      age_levels_events = integer(), # int vector n_interval - ages for event counts
      disp_births = numeric(1L), # dbl scalar - dispersion for birth rates
      disp_deaths = numeric(1L), # dbl scalar - dispersion for death rates
      disp_immigration1 = numeric(1L), # dbl scalar - dispersion for immigration rates
      disp_immigration2 = numeric(1L), # dbl scalar - dispersion for immigration rates
      disp_emigration1 = numeric(1L), # dbl scalar - dispersion for emigration rates
      disp_emigration2 = numeric(1L), # dbl scalar - dispersion for emigration rates
      logwt_unnorm = matrix(), # dbl matrix n_particle * (n_interval + 1) - log weights, unnormalised
      sum_loglik = numeric(), # dbl vector n_interval + 1 - log likelihood
      sum_logtrans = numeric(), # dbl vector n_interval - log transition probability
      sum_logimp = numeric(), # dbl vector n_interval + 1 - log importance probability
      sum_logwt_unnorm = numeric(), # dbl vector n_interval + 1 - log weights, unnormalised
      threshold = numeric(1L), # dbl scalar - threshold for resampling, between 0 and 1
      ess = numeric(), # dbl vector n_interval+1 - effective sample size, before resampling
      n_unique = integer(), # dbl vector n_interval+1 - number of unique particles
      resampled = logical(), # lgl vector n_interval+1 - whether resampling done
      index_parent = matrix(), # int matrix n_particle * (n_interval+1) - index of particle at prev interval
      index_ancestor = matrix(), # int matrix n_particle * (n_interval+1) - index of particle back to start
      index_output = integer(), # int vector n_particle / n_thin - indices of values included in output
      logimp_stock_end_init = numeric(), # dbl vector n_particle - log of importance prob of initial popn
      logimp_stock_end_net_mig = numeric(), # dbl vector n_particle - log of importance prob of population and net mig
      logimp_counts_im_em = numeric(), # dbl vector n_particle - log of importance prob of immigration, emigration
      logimp_gross_mig = numeric(), # dbl vector n_particle - log of importance prob of gross migration
      logimp_counts_dth = numeric() # dbl vector n_particle - log of importance prob of deaths
    )
  )

PFilterNoReg <-
  R6::R6Class("PFilterNoReg",
    public = list(
      ## fields - counts of people and events
      counts_stock = matrix(), # dbl matrix n_particle x (n_interval+1) - stock incl births of cohort,
      counts_births = matrix(), # dbl matrix n_particle x n_interval - births to cohort,
      counts_deaths = matrix(), # dbl matrix n_particle x n_interval
      counts_immigration1 = matrix(), # dbl matrix n_particle x n_interval
      counts_emigration1 = matrix(), # dbl matrix n_particle x n_interval
      counts_immigration2 = matrix(), # dbl matrix n_particle x n_interval
      counts_emigration2 = matrix(), # dbl matrix n_particle x n_interval
      ## fields - rates for events (where 'rates' defined to included expected counts)
      rates_births = numeric(), # dbl vector n_interval - births to cohort
      rates_deaths = numeric(), # dbl vector n_interval
      rates_immigration1 = numeric(), # dbl vector n_interval
      rates_emigration1 = numeric(), # dbl vector n_interval
      rates_immigration2 = numeric(), # dbl vector n_interval
      rates_emigration2 = numeric(), # dbl vector n_interval
      ## cohort data models
      cdms_stock = NULL, # object of class "CdmsNoreg"
      cdms_immigration1 = NULL, # object of class "CdmsNoreg"
      cdms_emigration1 = NULL, # object of class "CdmsNoreg"
      cdms_immigration2 = NULL, # object of class "CdmsNoreg"
      cdms_emigration2 = NULL, # object of class "CdmsNoreg"
      ## fields - intermediate calculations
      stock_start = numeric(), # dbl vector n_particle - stock at start of current interval
      stock_end = numeric(), # dbl vector n_particle - stock at end of current interval
      counts_bth = numeric(), # dbl vector n_particle - births during current interval
      counts_dth = numeric(), # dbl vector n_particle - deaths during current interval
      counts_im1 = numeric(), # dbl vector n_particle - immigration1 during current interval
      counts_em1 = numeric(), # dbl vector n_particle - emigration1 during current interval
      counts_im2 = numeric(), # dbl vector n_particle - immigration2 during current interval
      counts_em2 = numeric(), # dbl vector n_particle - emigration2 during current interval
      net_mig = numeric(), # dbl vector n_particle - net migration (including internal)
      gross_mig = numeric(), # dbl vector n_particle - gross migration (including internal)
      rates_bth = numeric(), # dbl vector n_particle - birth rates during current interval
      rates_dth = numeric(), # dbl vector n_particle - death rates during current interval
      rates_im1 = numeric(), # dbl vector n_particle - immigration1 rates during current interval
      rates_em1 = numeric(), # dbl vector n_particle - emigration1 rates during current interval
      rates_im2 = numeric(), # dbl vector n_particle - immigration2 rates during current interval
      rates_em2 = numeric(), # dbl vector n_particle - emigration2 rates during current interval
      rates_cim = numeric(), # dbl vector n_particle - combined immigration rate for current interval
      rates_cem = numeric(), # dbl vector n_particle - combined emigration rate for current interval
      exposure = numeric(), # dbl vector n_particle - exact exposure for current interval
      exposure_approx1 = numeric(), # dbl vector n_particle - approx exposure for current interval
      exposure_approx2 = numeric(), # dbl vector n_particle - approx exposure for current interval
      ## parallelogram importance function
      is_parallelogram_first = logical(), # indicator for importance function
      is_parallelogram_second = logical(), # indicator for importance function
      stock_end_second = numeric(), # dbl vector n_particle - stock at end of next interval
      counts_bth_second = numeric(), # dbl vector n_particle - births during second interval
      counts_dth_second = numeric(), # dbl vector n_particle - deaths during second interval
      counts_im1_second = numeric(), # dbl vector n_particle - immigration1 during second interval
      counts_em1_second = numeric(), # dbl vector n_particle - emigration1 during second interval
      net_mig_second = numeric(), # dbl vector n_particle - net migration (including internal) in second interval
      gross_mig_second = numeric(), # dbl vector n_particle - gross migration (including internal) in second interval
      net_mig_combined = numeric(), # dbl vector n_particle - net migration (including internal) in both intervals
      gross_mig_combined = numeric(), # dbl vector n_particle - gross migration (including internal) in both intervals
      rates_bth_second = numeric(), # dbl vector n_particle - birth rates during second interval
      rates_dth_second = numeric(), # dbl vector n_particle - death rates during second interval
      rates_im1_second = numeric(), # dbl vector n_particle - immigration1 rates during second interval
      rates_em1_second = numeric(), # dbl vector n_particle - emigration1 rates during second interval
      exposure_approx2_second = numeric(), # dbl vector n_particle - approx exposure for current interval
      exposure_second = numeric(), # dbl vector n_particle - exact exposure for second interval
      exposure_approx_first = numeric(), # dbl vector n_particle - approx exposure for first interval for parallelogram algorithm
      exposure_approx_second = numeric() # dbl vector n_particle - approx exposure for second interval for parallelogram algorithm
    ),
    inherit = PFilter
  )

PFilterWithReg <-
  R6::R6Class("PFilterWithReg",
              public = list(
      ## fields - estimation strategy
      is_split_mig = logical(1L), # whether to split migration into separate streams
      ## fields - description of cohort
      n_region = integer(1L), # int number of regions
      region_levels = character(), # chr vector n_region
      ## fields - counts of people and events
      counts_stock = array(), # dbl array n_particle x n_region x (n_interval+1) - stock incl births of cohort
      counts_births = array(), # dbl array n_particle x n_region x n_interval - births to cohort
      counts_deaths = array(), # dbl array n_particle x n_region x n_interval
      counts_internal_in = array(), # dbl array n_particle x n_region x n_interval
      counts_internal_out = array(), # dbl array n_particle x n_region x n_interval
      counts_immigration1 = array(), # dbl array n_particle x n_region x n_interval
      counts_emigration1 = array(), # dbl array n_particle x n_region x n_interval
      counts_immigration2 = array(), # dbl array n_particle x n_region x n_interval
      counts_emigration2 = array(), # dbl array n_particle x n_region x n_interval
      ## fields - rates for events (where 'rates' defined to included expected counts)
      rates_births = matrix(), # dbl matrix n_region x n_interval - births to cohort
      rates_deaths = matrix(), # dbl matrix n_region x n_interval
      rates_internal_in = matrix(), # dbl matrix n_region x n_interval
      rates_internal_out = matrix(), # dbl matrix n_region x n_interval
      rates_immigration1 = matrix(), # dbl matrix n_region x n_interval
      rates_emigration1 = matrix(), # dbl matrix n_region x n_interval
      rates_immigration2 = matrix(), # dbl matrix n_region x n_interval
      rates_emigration2 = matrix(), # dbl matrix n_region x n_interval
      ## fields - cohort data models
      cdms_stock = NULL, # object of class "CdmsWithreg"
      cdms_internal_in = NULL, # object of class "CdmsWithreg"
      cdms_internal_out = NULL, # object of class "CdmsWithreg"
      cdms_immigration1 = NULL, # object of class "CdmsWithreg"
      cdms_emigration1 = NULL, # object of class "CdmsWithreg"
      cdms_immigration2 = NULL, # object of class "CdmsWithreg"
      cdms_emigration2 = NULL, # object of class "CdmsWithreg"
      ## fields - intermediate calculations
      stock_start = matrix(), # dbl matrix n_particle x n_region - stock at start of interval
      stock_end = matrix(), # dbl matrix n_particle x n_region - stock at end of current interval
      counts_bth = matrix(), # dbl matrix n_particle x n_region - births during current interval
      counts_dth = matrix(), # dbl matrix n_particle x n_region - deaths during current interval
      counts_in = matrix(), # dbl matrix n_particle x n_region - internal in-mig during current interval
      counts_out = matrix(), # dbl matrix n_particle x n_region - internal out-mig during current interval
      counts_im1 = matrix(), # dbl matrix n_particle x n_region - immigration1 during current interval
      counts_em1 = matrix(), # dbl matrix n_particle x n_region - emigration1 during current interval
      counts_im2 = matrix(), # dbl matrix n_particle x n_region - immigration2 during current interval
      counts_em2 = matrix(), # dbl matrix n_particle x n_region - emigration2 during current interval
      net_mig = matrix(), # dbl matrix n_particle x n_region - net migration (including internal)
      gross_mig = matrix(), # dbl matrix n_particle x n_region - gross migration (including internal)
      rates_bth = matrix(), # dbl matrix n_particle x n_region - birth rates during current interval
      rates_dth = matrix(), # dbl matrix n_particle x n_region - death rates during current interval
      rates_in = matrix(), # dbl matrix n_particle x n_region - internal in-mig rates during interval
      rates_out = matrix(), # dbl matrix n_particle x n_region - internal out-mig rates during interval
      rates_im1 = matrix(), # dbl matrix n_particle x n_region - immigration1 rates during current interval
      rates_em1 = matrix(), # dbl matrix n_particle x n_region - emigration1 rates during current interval
      rates_im2 = matrix(), # dbl matrix n_particle x n_region - immigration2 rates during current interval
      rates_em2 = matrix(), # dbl matrix n_particle x n_region - emigration2 rates during current interval
      rates_cim = matrix(), # dbl matrix n_particle x n_region - combined immigration rate for interval
      rates_cem = matrix(), # dbl matrix n_particle x n_region - combined emigration rate for interval
      rates_cin = matrix(), # dbl matrix n_particle x n_region - combined in-migration rate for interval
      rates_cout = matrix(), # dbl matrix n_particle x n_region - combined out-migration rate for interval
      exposure = matrix(), # dbl matrix n_particle x n_region - exact exposure for current interval
      exposure_approx1 = matrix(), # dbl matrix n_particle x n_region - approx exposure for current interval
      exposure_approx2 = matrix(), # dbl matrix n_particle x n_region - approx exposure for current interval
      logimp_counts_in_out = numeric() # dbl vector n_particle log of importance prob of internal migration
    ),
    inherit = PFilter
  )




## initialisation functions ---------------------------------------------------

## These functions are called internally, and are never seen by
## users. We assume inputs have all been checked and are valid.
## The dimensions of 'counts' and 'rates' values,
## and the class of cdmsdepend on
## whether we have a no-region or with-region model

## Function 'initialize_pfilter' is called by the specific
## initialization methods for PFilterNoReg and PFilterWithReg,
## but, as noted, some of the inputs have different dimensions/classes
initialize_pfilter <- function(df, threshold) {
    n_particle <- df$n_particle
    n_thin <- df$n_thin
    n_interval <- df$n_interval
    self$cohort <- df$cohort
    self$sexgender <- df$sexgender
    self$has_stock_init <- df$has_stock_init
    self$has_one_imem <- df$has_one_imem
    self$is_dominant <- df$is_dominant
    self$obs_zero <- df$obs_zero
    self$n_interval <- n_interval
    self$n_particle <- n_particle
    self$n_thin <- n_thin
    self$threshold <- threshold
    self$is_popn <- df$is_popn
    self$n_popn <- df$n_popn
    self$time_levels_stock <- df$time_levels_stock
    self$time_levels_events <- df$time_levels_events
    self$age_levels_stock <- df$age_levels_stock
    self$age_levels_events <- df$age_levels_events
    self$logwt_unnorm <- matrix(nrow = n_particle, ncol = n_interval + 1L)
    self$sum_loglik <- numeric(n_interval + 1L)
    self$sum_logtrans <- numeric(n_interval)
    self$sum_logimp <- numeric(n_interval + 1L)
    self$sum_logwt_unnorm <- numeric(n_interval + 1L)
    self$threshold <- threshold
    self$ess <- numeric(length = n_interval + 1L)
    self$n_unique <- integer(length = n_interval + 1L)
    self$resampled <- logical(length = n_interval + 1L)
    self$index_parent <- matrix(nrow = n_particle, ncol = n_interval + 1L)
    self$index_ancestor <- matrix(nrow = n_particle, ncol = n_interval + 1)
    self$index_output <- integer(length = n_particle %/% n_thin) # integer division
    self$rates_births <- df$rates_births
    self$rates_deaths <- df$rates_deaths
    self$rates_immigration1 <- df$rates_immigration1
    self$rates_emigration1 <- df$rates_emigration1
    self$rates_immigration2 <- df$rates_immigration2
    self$rates_emigration2 <- df$rates_emigration2
    self$disp_births <- df$disp_births
    self$disp_deaths <- df$disp_deaths
    self$disp_immigration1 <- df$disp_immigration1
    self$disp_emigration1 <- df$disp_emigration1
    self$disp_immigration2 <- df$disp_immigration2
    self$disp_emigration2 <- df$disp_emigration2
    self$cdms_stock <- df$cdms_stock
    self$cdms_immigration1 <- df$cdms_immigration1
    self$cdms_emigration1 <- df$cdms_emigration1
    self$cdms_immigration2 <- df$cdms_immigration2
    self$cdms_emigration2 <- df$cdms_emigration2
    self$logimp_stock_end_init <- numeric(length = n_particle)
    self$logimp_stock_end_net_mig <- numeric(length = n_particle)
    self$logimp_counts_im_em <- numeric(length = n_particle)
    self$logimp_gross_mig <- numeric(length = n_particle)
    self$logimp_counts_dth <- numeric(length = n_particle)
}

PFilter$set(
  which = "public",
  name = "initialize",
  value = initialize_pfilter
)

initialize_pfilter_noreg <- function(df,
                                     threshold,
                                     is_forecast) {
  super$initialize(
    df = df,
    threshold = threshold
    )
  n_particle <- self$n_particle
  n_interval <- self$n_interval
  has_stock_init <- self$has_stock_init
  counts_stock <- matrix(NA_real_, nrow = n_particle, ncol = n_interval + 1L)
  if (has_stock_init) {
    counts_stock[, 1L] <- df$counts_stock
  }
  self$counts_stock <- counts_stock
  if (is_forecast) {
    counts_births <- rep(NA_real_, times = n_particle * n_interval)
    counts_deaths <- rep(NA_real_, times = n_particle * n_interval)
  } else {
    counts_births <- rep(df$counts_births, each = n_particle)
    counts_deaths <- rep(df$counts_deaths, each = n_particle)
  }
  self$counts_births <- matrix(counts_births, nrow = n_particle, ncol = n_interval)
  self$counts_deaths <- matrix(counts_deaths, nrow = n_particle, ncol = n_interval)
  self$counts_immigration1 <- matrix(NA_real_, nrow = n_particle, ncol = n_interval)
  self$counts_emigration1 <- matrix(NA_real_, nrow = n_particle, ncol = n_interval)
  self$counts_immigration2 <- matrix(NA_real_, nrow = n_particle, ncol = n_interval)
  self$counts_emigration2 <- matrix(NA_real_, nrow = n_particle, ncol = n_interval)
  self$rates_bth <- numeric(length = n_particle)
  self$rates_dth <- numeric(length = n_particle)
  self$rates_im1 <- numeric(length = n_particle)
  self$rates_em1 <- numeric(length = n_particle)
  self$rates_im2 <- numeric(length = n_particle)
  self$rates_em2 <- numeric(length = n_particle)
  self$rates_cim <- numeric(length = n_particle)
  self$rates_cem <- numeric(length = n_particle)
  self$is_parallelogram_first <- df$is_parallelogram_first
  self$is_parallelogram_second <- df$is_parallelogram_second
}

PFilterNoReg$set(
  which = "public",
  name = "initialize",
  value = initialize_pfilter_noreg
)

initialize_pfilter_withreg <- function(df,
                                       threshold,
                                       is_forecast) {
  super$initialize(
    df = df,
    threshold = threshold
    )
  n_particle <- self$n_particle
  n_interval <- self$n_interval
  n_region <- df$n_region
  self$is_split_mig <- df$is_split_mig
  self$region_levels <- df$region_levels
  self$n_region <- n_region
  has_stock_init <- self$has_stock_init
  counts_stock <- array(NA_real_, dim = c(n_particle, n_region, n_interval + 1L))
  if (has_stock_init) {
    counts_stock[, , 1L] <- df$counts_stock
  }
  self$counts_stock <- counts_stock
  if (is_forecast) {
    counts_births <- rep(NA_real_, times = n_particle * n_region * n_interval)
    counts_deaths <- rep(NA_real_, times = n_particle * n_region * n_interval)
  } else {
    counts_births <- rep(df$counts_births, each = n_particle)
    counts_deaths <- rep(df$counts_deaths, each = n_particle)
  }
  self$counts_births <- array(counts_births, dim = c(n_particle, n_region, n_interval))
  self$counts_deaths <- array(counts_deaths, dim = c(n_particle, n_region, n_interval))
  self$counts_internal_in <- array(NA_real_, dim = c(n_particle, n_region, n_interval))
  self$counts_internal_out <- array(NA_real_, dim = c(n_particle, n_region, n_interval))
  self$counts_immigration1 <- array(NA_real_, dim = c(n_particle, n_region, n_interval))
  self$counts_emigration1 <- array(NA_real_, dim = c(n_particle, n_region, n_interval))
  self$counts_immigration2 <- array(NA_real_, dim = c(n_particle, n_region, n_interval))
  self$counts_emigration2 <- array(NA_real_, dim = c(n_particle, n_region, n_interval))
  self$rates_internal_in <- df$rates_internal_in
  self$rates_internal_out <- df$rates_internal_out
  self$cdms_internal_in <- df$cdms_internal_in
  self$cdms_internal_out <- df$cdms_internal_out
  self$rates_bth <- matrix(nrow = n_particle, ncol = n_region)
  self$rates_dth <- matrix(nrow = n_particle, ncol = n_region)
  self$rates_in <- matrix(nrow = n_particle, ncol = n_region)
  self$rates_out <- matrix(nrow = n_particle, ncol = n_region)
  self$rates_im1 <- matrix(nrow = n_particle, ncol = n_region)
  self$rates_em1 <- matrix(nrow = n_particle, ncol = n_region)
  self$rates_im2 <- matrix(nrow = n_particle, ncol = n_region)
  self$rates_em2 <- matrix(nrow = n_particle, ncol = n_region)
  self$rates_cim <- matrix(nrow = n_particle, ncol = n_region)
  self$rates_cem <- matrix(nrow = n_particle, ncol = n_region)
  self$rates_cin <- matrix(nrow = n_particle, ncol = n_region)
  self$rates_cout <- matrix(nrow = n_particle, ncol = n_region)
  self$logimp_counts_in_out <- numeric(length = n_particle)
}


PFilterWithReg$set(
  which = "public",
  name = "initialize",
  value = initialize_pfilter_withreg
)




