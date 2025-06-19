
test_that("'simulated_df_noreg' creates data frame", {
  ans <- simulated_df_noreg()
  expect_true(is.list(ans))
})

test_that("'simulated_df_withreg' creates data frame", {
  ans <- simulated_df_withreg()
  expect_true(is.list(ans))
})

test_that("'simulated_pf_noreg' creates object of class PFilterNoReg", {
  ans <- simulated_pf_noreg(cohort_type = "existing")
  expect_s3_class(ans, "PFilterNoReg")
  ans <- simulated_pf_noreg(cohort_type = "new")
  expect_s3_class(ans, "PFilterNoReg")
  ans <- simulated_pf_noreg(cohort_type = "expiring")
  expect_s3_class(ans, "PFilterNoReg")
})

test_that("'simulated_pf_withreg' creates object of class PFilterWithReg", {
  ans <- dpmpf:::simulated_pf_withreg(cohort_type = "existing")
  expect_s3_class(ans, "PFilterWithReg")
  ans <- simulated_pf_withreg(cohort_type = "new")
  expect_s3_class(ans, "PFilterWithReg")
  ans <- simulated_pf_withreg(cohort_type = "expiring")
  expect_s3_class(ans, "PFilterWithReg")
})

test_that("return value for 'make_args_simulated_existing' has expected names", {
  expect_setequal(
    names(make_args_simulated_existing()),
    c(
      "n_interval",
      "time_levels_stock",
      "time_levels_events",
      "age_levels_stock",
      "age_levels_events",
      "is_popn",
      "n_popn"
    )
  )
})

test_that("return value for 'make_args_simulated_expiring' has expected names", {
  expect_setequal(
    names(make_args_simulated_expiring()),
    c(
      "n_interval",
      "time_levels_stock",
      "time_levels_events",
      "age_levels_stock",
      "age_levels_events",
      "is_popn",
      "n_popn"
    )
  )
})

test_that("return value for 'make_args_simulated_new' has expected names", {
  expect_setequal(
    names(make_args_simulated_new()),
    c(
      "n_interval",
      "time_levels_stock",
      "time_levels_events",
      "age_levels_stock",
      "age_levels_events",
      "is_popn",
      "n_popn"
    )
  )
})
