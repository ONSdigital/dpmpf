
## dpoibin1 -------------------------------------------------------------------

test_that("dpoibin1 works when x < threshold", {
  set.seed(0)
  for (i in 1:10) {
    x <- rpois(n = 1, lambda = 10)
    prob <- runif(1, min = 0.5, max = 1)
    size <- rpois(n = 1, lambda = 4)
    ans_obtained <- dpoibin1(x = x, size = size, prob = prob, use_log = FALSE)
    ans_expected <- sum(dbinom(seq(from = 0L, to = x), size = size, prob = prob) *
      dpois(seq(from = x, to = 0L), lambda = size * (1 - prob)))
    expect_equal(ans_obtained, ans_expected)
  }
})

test_that("dpoibin1 works when x > threshold", {
  set.seed(100)
  for (i in 1:10) {
    x <- rpois(n = 1, lambda = 10) + 50L
    prob <- runif(1, min = 0.5, max = 1)
    size <- as.integer(abs(rnorm(1, mean = 12, sd = 2))) + 50L
    ans_obtained <- dpoibin1(x = x, size = size, prob = prob, use_log = FALSE)
    mean <- prob * size + (1 - prob) * size
    sd <- sqrt(prob * (1 - prob) * floor(size) + (1 - prob) * size)
    ans_expected <- dnorm(x, mean = mean, sd = sd)
    expect_equal(ans_obtained, ans_expected)
  }
})

test_that("log argument for dpoibin1 works", {
  expect_equal(
    dpoibin1(x = 10L, size = 11L, prob = 0.98, use_log = TRUE),
    log(dpoibin1(x = 10L, size = 11L, prob = 0.98, use_log = FALSE))
  )
  expect_equal(
    dpoibin1(x = 1L, size = 1L, prob = 0.98, use_log = TRUE),
    log(dpoibin1(x = 1L, size = 1L, prob = 0.98, use_log = FALSE))
  )
  expect_equal(
    dpoibin1(x = 3000L, size = 3000L, prob = 0.98, use_log = TRUE),
    log(dpoibin1(x = 3000L, size = 3000L, prob = 0.98, use_log = FALSE))
  )
})

test_that("R and C versions of dpoibin1 give the same value", {
  dpoibin1_r <- function(x, size, prob, use_log = FALSE) {
    kThreshold <- 50
    lambda <- (1 - prob) * size
    if (x > kThreshold) {
      sd <- sqrt((1 - prob^2) * size)
      ans <- stats::dnorm(x, mean = size, sd = sd, log = use_log)
    } else {
      limit <- min(x, size)
      seq_binom <- seq.int(from = 0L, to = limit)
      seq_pois <- seq.int(from = x, to = x - limit)
      logimp_binom <- stats::dbinom(seq_binom, size = size, prob = prob, log = TRUE)
      logimp_pois <- stats::dpois(seq_pois, lambda = lambda, log = TRUE)
      logimp <- logimp_binom + logimp_pois
      if (use_log) {
        max_val <- max(logimp)
        ans <- max_val + log(sum(exp(logimp - max_val)))
      } else {
        ans <- sum(exp(logimp))
      }
    }
    ans
  }
  x <- round(runif(n = 100, max = 100))
  size <- pmax(x + round(runif(n = 100, min = -5, max = 5)), 0)
  prob <- runif(n = 100, min = 0.8)
  use_log <- runif(n = 100) < 0.5
  ans_r <- mapply(dpoibin1_r, x = x, size = size, prob = prob, use_log = use_log)
  ans_c <- mapply(dpoibin1, x = x, size = size, prob = prob, use_log = use_log)
  expect_equal(ans_r, ans_c)
})


## rpoibin1 -------------------------------------------------------------------

test_that("random variates generated by rpoibin1 have the expected mean and variance", {
  set.seed(100)
  prob <- 0.98
  size <- 500
  x <- replicate(n = 10000, rpoibin1(size = size, prob = prob))
  expect_equal(mean(x), size, tolerance = 0.01)
  expect_equal(var(x), (1 - prob^2) * size, tolerance = 0.01)
})

test_that("R and C versions of rpoibin1 give the same value", {
  rpoibin1_r <- function(size, prob) {
    val_binom <- stats::rbinom(n = 1L, size = size, prob = prob)
    val_pois <- stats::rpois(n = 1L, lambda = (1 - prob) * size)
    val_binom + val_pois
  }
  for (seed in 1:100) {
    set.seed(seed)
    size <- as.integer(runif(n = 1, max = 100))
    prob <- runif(1)
    set.seed(seed)
    ans_r <- rpoibin1_r(size = size, prob = prob)
    set.seed(seed)
    ans_c <- rpoibin1(size = size, prob = prob)
    expect_equal(ans_r, ans_c)
  }
})

test_that("'rpoibin1' returns doubles", {
  expect_true(is.double(rpoibin1(size = 10, prob = 0.9)))
})


## rpoistr1 -------------------------------------------------------------------

test_that("rpoistr1 passes basic sanity checks", {
  set.seed(0)
  lambda <- runif(n = 100, min = 0, max = 20)
  lower <- as.integer(lambda + runif(n = 100, min = 0, max = 5))
  ans <- mapply(rpoistr1, lambda = lambda, lower = lower)
  expect_false(anyNA(ans))
  expect_true(all(ans >= lower))
  expect_true(is.double(ans))
})

test_that("rpoistr1 gives the same distribution as using brute force when lambda is small", {
  set.seed(0)
  ans_rpoistr1 <- replicate(n = 10000, rpoistr1(lambda = 5, lower = 10L))
  ans_brute <- rpois(n = 100000, lambda = 5)
  ans_brute <- ans_brute[ans_brute >= 10L]
  expect_equal(mean(ans_rpoistr1), mean(ans_brute), tolerance = 0.1)
  expect_equal(median(ans_rpoistr1), median(ans_brute), tolerance = 0.1)
  expect_equal(var(ans_rpoistr1), var(ans_brute), tolerance = 0.1)
  expect_equal(mean(ans_rpoistr1 == 10L), mean(ans_brute == 10L), tolerance = 0.05)
})

test_that("rpoistr1 gives the same distribution as using brute force when lambda is large", {
  set.seed(0)
  ans_rpoistr1 <- replicate(n = 10000, rpoistr1(lambda = 50000, lower = 50100L))
  ans_brute <- rpois(n = 100000, lambda = 50000)
  ans_brute <- ans_brute[ans_brute >= 50100L]
  expect_equal(mean(ans_rpoistr1), mean(ans_brute), tolerance = 0.1)
  expect_equal(median(ans_rpoistr1), median(ans_brute), tolerance = 0.1)
  expect_equal(var(ans_rpoistr1), var(ans_brute), tolerance = 0.1)
  expect_equal(mean(ans_rpoistr1 == 10L), mean(ans_brute == 10L), tolerance = 0.05)
})

test_that("rpoistr1 returns 'lower' when 'p_lower' is very high", {
  set.seed(0)
  ans_obtained <- rpoistr1(lambda = 5, lower = 10000)
  ans_expected <- 10000
  expect_identical(ans_obtained, ans_expected)
})

test_that("R and C versions of rpoistr1 give the same answer", {
  rpoistr1_r <- function(lambda, lower) {
    max_attempt <- 10L
    n_attempt <- 0L
    found <- FALSE
    while ((n_attempt < max_attempt) && !found) {
      n_attempt <- n_attempt + 1L
      prop_value <- stats::rpois(n = 1L, lambda = lambda)
      found <- prop_value >= lower
    }
    if (found) {
      return(prop_value)
    }
    p_lower <- stats::ppois(
      q = lower - 1,
      lambda = lambda,
      lower.tail = TRUE,
      log.p = FALSE
    )
    p_lower_max <- 1 - 1e-10
    if (p_lower > p_lower_max) {
      return(lower)
    }
    U <- stats::runif(
      n = 1L,
      min = p_lower,
      max = 1
    )
    ans <- stats::qpois(
      p = U,
      lambda = lambda,
      lower.tail = TRUE,
      log.p = FALSE
    )
    if (ans < lower) {
      ans <- lower
    }
    ans
  }
  set.seed(0)
  lambda <- runif(n = 100, min = 0, max = 20)
  lower <- as.integer(lambda + runif(n = 100, min = 0, max = 5))
  set.seed(0)
  ans_c <- mapply(rpoistr1, lambda = lambda, lower = lower)
  set.seed(0)
  ans_r <- mapply(rpoistr1_r, lambda = lambda, lower = lower)
  expect_equal(ans_r, ans_c)
})

test_that("'rpoistr1' returns doubles", {
  expect_true(is.double(rpoistr1(lambda = 10, lower = 9)))
})



## rpoistr --------------------------------------------------------------------

test_that("R and C versions of rpoistr give the same answer", {
  rpoistr_r <- function(n, lambda, lower) {
    ans <- integer(length = n)
    for (i in seq_len(n)) {
      ans[[i]] <- rpoistr1(
        lambda = lambda[[i]],
        lower = lower[[i]]
      )
    }
    ans
  }
  set.seed(0)
  lambda <- runif(n = 100, min = 0, max = 20)
  lower <- as.integer(lambda + runif(n = 100, min = 0, max = 5))
  set.seed(0)
  ans_c <- rpoistr(n = 100, lambda = lambda, lower = lower)
  set.seed(0)
  ans_r <- rpoistr_r(n = 100, lambda = lambda, lower = lower)
  expect_equal(ans_r, ans_c)
})

test_that("'rpoistr' returns doubles", {
  expect_true(is.double(rpoistr(n = 3, lambda = 10, lower = 9)))
})


## pskel1 ---------------------------------------------------------------------

test_that("'pskel1' gives same answers as 'pskellam' function in 'skellam' package", {
  ans_obtained <- pskel1(
    q = 5,
    mu1 = 2.5,
    mu2 = 1.5,
    lower_tail = TRUE,
    log_p = FALSE
  )
  ans_expected <- 0.9846523
  expect_equal(ans_obtained, ans_expected, tolerance = 0.0001)
  ans_obtained <- pskel1(
    q = 5,
    mu1 = 2.5,
    mu2 = 1.5,
    lower_tail = FALSE,
    log_p = FALSE
  )
  ans_expected <- 0.01534768
  expect_equal(ans_obtained, ans_expected, tolerance = 0.0001)
  ans_obtained <- pskel1(
    q = -5,
    mu1 = 2.5,
    mu2 = 1.5,
    lower_tail = TRUE,
    log_p = FALSE
  )
  ans_expected <- 0.002717002
  expect_equal(ans_obtained, ans_expected, tolerance = 0.0001)
  ans_obtained <- pskel1(
    q = -5,
    mu1 = 2.5,
    mu2 = 1.5,
    lower_tail = FALSE,
    log_p = FALSE
  )
  ans_expected <- 0.997283
  expect_equal(ans_obtained, ans_expected, tolerance = 0.0001)
  ans_obtained <- pskel1(
    q = -1,
    mu1 = 2.5,
    mu2 = 4.5,
    lower_tail = TRUE,
    log_p = FALSE
  )
  ans_expected <- 0.7144911
  expect_equal(ans_obtained, ans_expected, tolerance = 0.0001)
  ans_obtained <- pskel1(
    q = 20,
    mu1 = 12.5,
    mu2 = 4.5,
    lower_tail = TRUE,
    log_p = TRUE
  )
  ans_expected <- -0.001996993
  expect_equal(ans_obtained, ans_expected, tolerance = 0.01)
})

test_that("R and C versions of pskel1 give the same answer", {
    pskel1_r <- function(q, mu1, mu2, lower_tail, log_p) {
        if (lower_tail) {
            if (q < 0) {
                ans <- pnchisq_switch(
                    q = 2 * mu2,
                    df = -2 * q,
                    ncp = 2 * mu1,
                    lower = TRUE,
                    use_log = log_p
                )
            } else {
                ans <- pnchisq_switch(
                    q = 2 * mu1,
                    df = 2 * (q + 1),
                    ncp = 2 * mu2,
                    lower = FALSE,
                    use_log = log_p
                )
            }
        } else {
            if (q < 0) {
                ans <- pnchisq_switch(
                    q = 2 * mu2,
                    df = -2 * q,
                    ncp = 2 * mu1,
                    lower = FALSE,
                    use_log = log_p
                )
            } else {
                ans <- pnchisq_switch(
                    q = 2 * mu1,
                    df = 2 * (q + 1),
                    ncp = 2 * mu2,
                    lower = TRUE,
                    use_log = log_p
                )
            }
        }
        ans
    }
    set.seed(0)
    q <- as.numeric(round(runif(n = 100, min = -10, max = 10)))
    mu1 <- runif(n = 100, min = 0, max = 20)
    mu2 <- runif(n = 100, min = 0, max = 20)
    lower_tail <- sample(c(TRUE, FALSE), size = 100, replace = TRUE)
    log_p <- sample(c(TRUE, FALSE), size = 100, replace = TRUE)
    ans_c <- mapply(pskel1_r, q = q, mu1 = mu1, mu2 = mu2, lower_tail = lower_tail, log_p = log_p)
    ans_r <- mapply(pskel1, q = q, mu1 = mu1, mu2 = mu2, lower_tail = lower_tail, log_p = log_p)
    expect_equal(ans_r, ans_c)
})


## dpoistr --------------------------------------------------------------------

test_that("'dpoistr' works when lower = 0", {
    x <- as.double(0:10)
    lambda <- rep(5, 11)
    lower <- rep(0, 11)
    ## log is TRUE
    ans_obtained <- dpoistr(x = x, lambda = lambda, lower = lower, use_log = TRUE)
    ans_expected <- dpois(x = x, lambda = lambda, log = TRUE)
    expect_equal(ans_obtained, ans_expected)
    ## log is FALSE
    ans_obtained <- dpoistr(x = x, lambda = lambda, lower = lower, use_log = FALSE)
    ans_expected <- dpois(x = x, lambda = lambda, log = FALSE)
    expect_equal(ans_obtained, ans_expected)
})

test_that("'dpoistr' works when lower = 1", {
    x <- as.double(0:10)
    lambda <- rep(5, 11)
    lower <- rep(1, 11)
    ## log is TRUE
    ans_obtained <- dpoistr(x = x, lambda = lambda, lower = lower, use_log = TRUE)
    ans_expected <- dpois(x = x, lambda = lambda, log = TRUE) -
        ppois(lower - 1, lambda = lambda, log = TRUE, lower.tail = FALSE)
    ans_expected[x < lower] <- -Inf
    expect_equal(ans_obtained, ans_expected)
    ## log is FALSE
    ans_obtained <- dpoistr(x = x, lambda = lambda, lower = lower, use_log = FALSE)
    ans_expected <- dpois(x = x, lambda = lambda, log = FALSE) /
        ppois(lower - 1, lambda = lambda, log = FALSE, lower.tail = FALSE)
    ans_expected[x < lower] <- 0
    expect_equal(ans_obtained, ans_expected)
})

test_that("'dpoistr' works when lower = 2", {
    x <- as.double(0:10)
    lambda <- rep(5, 11)
    lower <- rep(1, 11)
    ## log is TRUE
    ans_obtained <- dpoistr(x = x, lambda = lambda, lower = lower, use_log = TRUE)
    ans_expected <- dpois(x = x, lambda = lambda, log = TRUE) -
        ppois(lower - 1, lambda = lambda, log = TRUE, lower.tail = FALSE)
    ans_expected[x < lower] <- -Inf
    expect_equal(ans_obtained, ans_expected)
    ## log is FALSE
    ans_obtained <- dpoistr(x = x, lambda = lambda, lower = lower, use_log = FALSE)
    ans_expected <- dpois(x = x, lambda = lambda, log = FALSE) /
        ppois(lower - 1, lambda = lambda, log = FALSE, lower.tail = FALSE)
    ans_expected[x < lower] <- 0
    expect_equal(ans_obtained, ans_expected)
})

test_that("C and R versions of 'dpoistr' give same answer", {
    dpoistr_r <- function(x, lambda, lower, use_log) {
        log_dens_untrunc <- stats::dpois(
                                       x = x,
                                       lambda = lambda,
                                       log = TRUE
                                   )
        log_const <- stats::ppois(
                                q = lower - 1,
                                lambda = lambda,
                                lower.tail = FALSE,
                                log.p = TRUE
                            )
        log_dens_trunc <- log_dens_untrunc - log_const
        less_than_lower <- x < lower
        if (use_log) {
            ans <- log_dens_trunc
            ans[less_than_lower] <- -Inf
        } else {
            ans <- exp(log_dens_trunc)
            ans[less_than_lower] <- 0
        }
        ans
    }
    x <- rpois(n = 100, lambda = 10)
    lambda <- runif(n = 100, min = 0.01, max = 20)
    lower <- rpois(n = 100, lambda = 5)
    use_log <- sample(c(TRUE, FALSE), size = 100, replace = TRUE)
    ans_c <- mapply(dpoistr, x = x, lambda = lambda, lower = lower, use_log = use_log)
    ans_r <- mapply(dpoistr_r, x = x, lambda = lambda, lower = lower, use_log = use_log)
    expect_equal(ans_c, ans_r)    
})

test_that("'dpoistr' works when p(x) = 0", {
    expect_identical(dpoistr(x = 10, lambda = 0, lower = 8, use_log = TRUE),
                     -Inf)
    expect_identical(dpoistr(x = 10, lambda = 0, lower = 8, use_log = FALSE),
                     0)
})


## pnchisq_approx -------------------------------------------------------------

test_that("pnchisq_approx gives answers close to pchisq", {
    set.seed(0)
    n <- 100L
    q <- runif(n, 0.001, 10)
    df <- runif(n, 0.001, 10)
    ncp  <- runif(n, 0.001, 10)
    lower <- sample(c(TRUE, FALSE), size = n, replace = TRUE)
    ans_exact <- mapply(pchisq, q = q, df = df, ncp = ncp, lower = lower, log = FALSE)
    ans_approx <- mapply(pnchisq_approx, q = q, df = df, ncp = ncp,
                         lower = lower, use_log = FALSE)
    expect_equal(ans_exact, ans_approx, tolerance = 0.005)
})


## pnchisq_approx -------------------------------------------------------------

test_that("pnchisq_switch gives answers close to pchisq or to pnchisq_qpprox", {
    set.seed(0)
    n <- 100L
    q <- runif(n, 0.001, 10)
    df <- runif(n, 0.001, 10)
    ncp  <- runif(n, 0.001, 10)
    lower <- sample(c(TRUE, FALSE), size = n, replace = TRUE)
    use_approx <- (q > 20) | (df > 20) | (ncp > 20)
    ans_exact <- mapply(pchisq, q = q, df = df, ncp = ncp, lower = lower, log = FALSE)
    ans_approx <- mapply(pnchisq_approx, q = q, df = df, ncp = ncp,
                         lower = lower, use_log = FALSE)
    ans_switch <- mapply(pnchisq_switch, q = q, df = df, ncp = ncp,
                         lower = lower, use_log = FALSE)
    expect_equal(ans_switch[use_approx], ans_approx[use_approx])
    expect_equal(ans_switch[!use_approx], ans_exact[!use_approx])
})


## pskel ----------------------------------------------------------------------

test_that("R and C versions of pskel give the same answer", {
  pskel_r <- function(q, mu1, mu2, lower_tail, log_p) {
    n <- length(q)
    ans <- numeric(length = n)
    for (i in seq_len(n)) {
      ans[i] <- pskel1(
        q = q[i],
        mu1 = mu1[i],
        mu2 = mu2[i],
        lower_tail = lower_tail,
        log_p = log_p
      )
    }
    ans
  }
  set.seed(0)
  q <- as.numeric(round(runif(n = 10, min = -10, max = 10)))
  mu1 <- runif(n = 10, min = 0, max = 20)
  mu2 <- runif(n = 10, min = 0, max = 20)
  ans_c <- pskel_r(q = q, mu1 = mu1, mu2 = mu2, lower_tail = FALSE, log_p = TRUE)
  ans_r <- pskel(q = q, mu1 = mu1, mu2 = mu2, lower_tail = FALSE, log_p = TRUE)
  expect_equal(ans_r, ans_c)
})


## qskel1 ---------------------------------------------------------------------

test_that("'qskel1' gives same answers as 'qskellam' function in 'skellam' package", {
  ans_obtained <- qskel1(
    p = 0.9,
    mu1 = 12.5,
    mu2 = 10.5
  )
  ans_expected <- 8
  expect_equal(ans_obtained, ans_expected)
  ans_obtained <- qskel1(
    p = 0.9,
    mu1 = 120.5,
    mu2 = 100.5
  )
  ans_expected <- 39
  expect_equal(ans_obtained, ans_expected)
  ans_obtained <- qskel1(
    p = 0.1,
    mu1 = 120.5,
    mu2 = 100.5
  )
  ans_expected <- 1
  expect_equal(ans_obtained, ans_expected)
  ans_obtained <- qskel1(
    p = 0.00000001,
    mu1 = 120.5,
    mu2 = 100.5
  )
  ans_expected <- -63
  expect_equal(ans_obtained, ans_expected, tolerance = 0.05)
  ans_obtained <- qskel1(
    p = 0,
    mu1 = 120.5,
    mu2 = 100.5
  )
  ans_expected <- -Inf
  expect_equal(ans_obtained, ans_expected)
  ans_obtained <- qskel1(
    p = 1,
    mu1 = 120.5,
    mu2 = 100.5
  )
  ans_expected <- Inf
  expect_equal(ans_obtained, ans_expected)
  ans_obtained <- qskel1(
    p = 0.5,
    mu1 = 0,
    mu2 = 5
  )
  ans_expected <- -5
  expect_equal(ans_obtained, ans_expected)
  ans_obtained <- qskel1(
    p = 0.5,
    mu1 = 5,
    mu2 = 0
  )
  ans_expected <- 5
  expect_equal(ans_obtained, ans_expected)
})

test_that("'qskel1' and 'pskel1' are inverses", {
  mu1 <- runif(100, 0, 20)
  mu2 <- runif(100, 0, 20)
  q0 <- round(runif(100, -10, 10))
  p <- mapply(pskel1,
    q = q0,
    mu1 = mu1,
    mu2 = mu2,
    lower_tail = TRUE,
    log_p = FALSE
  )
  q1 <- mapply(qskel1,
    p = p,
    mu1 = mu1,
    mu2 = mu2
  )
  is_valid <- (p > 0) & (p < 1)
  expect_equal(q0[is_valid], q1[is_valid])
})

test_that("R and C versions of qskel1 give the same answer", {
  qskel1_r <- function(p, mu1, mu2) {
    if (mu2 == 0) {
      ans <- stats::qpois(
        p = p,
        lambda = mu1,
        lower.tail = TRUE,
        log.p = FALSE
      )
      return(ans)
    }
    if (mu1 == 0) {
      ans <- -stats::qpois(
        p = p,
        lambda = mu2,
        lower.tail = FALSE,
        log.p = FALSE
      )
      return(ans)
    }
    mu <- mu1 - mu2
    sigma_sq <- mu1 + mu2
    sigma <- sqrt(sigma_sq)
    z <- stats::qnorm(
      p = p,
      mean = 0,
      sd = 1,
      lower.tail = TRUE,
      log.p = FALSE
    )
    if (is.infinite(z)) {
      return(z)
    }
    q_prop <- mu + sigma * z
    skew <- (z^2 - 1) * mu / sigma_sq / 6
    kurt <- -(skew * mu - 2 * mu1 * mu2 * (z^2 - 3) / sigma_sq) * z / 12 / sigma
    q_prop <- round(q_prop + skew + kurt)
    p_prop <- pskel1(
      q = q_prop,
      mu1 = mu1,
      mu2 = mu2,
      lower_tail = TRUE,
      log_p = FALSE
    )
    too_high <- p_prop > p
    while (too_high) {
      q_prop <- q_prop - 1
      p_prop <- pskel1(
        q = q_prop,
        mu1 = mu1,
        mu2 = mu2,
        lower_tail = TRUE,
        log_p = FALSE
      )
      too_high <- p_prop > p
    }
    too_low <- p_prop < p
    while (too_low) {
      q_prop <- q_prop + 1
      p_prop <- pskel1(
        q = q_prop,
        mu1 = mu1,
        mu2 = mu2,
        lower_tail = TRUE,
        log_p = FALSE
      )
      too_low <- p_prop < p
    }
    q_prop
  }
  set.seed(0)
  p <- runif(n = 100)
  mu1 <- runif(n = 100, min = 0, max = 20)
  mu2 <- runif(n = 100, min = 0, max = 20)
  ans_c <- mapply(qskel1_r, p = p, mu1 = mu1, mu2 = mu2)
  ans_r <- mapply(qskel1, p = p, mu1 = mu1, mu2 = mu2)
  expect_equal(ans_r, ans_c)
})

test_that("'qskel1' returns doubles", {
  expect_true(is.double(qskel1(p = 0.3, mu1 = 5, mu2 = 5)))
})


## rskeltr1 -------------------------------------------------------------------

test_that("rskeltr1 passes basic sanity checks", {
  set.seed(0)
  mu1 <- runif(n = 100, min = 0, max = 20)
  mu2 <- runif(n = 100, min = 0, max = 20)
  lower <- as.numeric(round(runif(n = 100, min = -20, max = 20)))
  ans <- mapply(rskeltr1, mu1 = mu1, mu2 = mu2, lower = lower)
  expect_false(anyNA(ans))
  expect_true(all(ans >= lower))
  expect_true(is.double(ans))
})

test_that("rskeltr1 gives the same distribution as using brute force when mu1, mu2 are small", {
  set.seed(0)
  ans_rskeltr1 <- replicate(n = 10000, rskeltr1(mu1 = 5, mu2 = 3, lower = 5))
  ans_brute <- rpois(n = 100000, lambda = 5) - rpois(n = 10000, lambda = 3)
  ans_brute <- ans_brute[ans_brute >= 5]
  expect_equal(mean(ans_rskeltr1), mean(ans_brute), tolerance = 0.1)
  expect_equal(median(ans_rskeltr1), median(ans_brute), tolerance = 0.1)
  expect_equal(var(ans_rskeltr1), var(ans_brute), tolerance = 0.1)
  expect_equal(mean(ans_rskeltr1 == 10L), mean(ans_brute == 10L), tolerance = 0.05)
})

test_that("rskeltr1 gives the same distribution as using brute force when lambda is large", {
  set.seed(0)
  ans_rskeltr1 <- replicate(n = 10000, rskeltr1(mu1 = 60000, mu2 = 10000, lower = 50100))
  ans_brute <- rpois(n = 100000, lambda = 60000) - rpois(n = 100000, lambda = 10000)
  ans_brute <- ans_brute[ans_brute >= 50100]
  expect_equal(mean(ans_rskeltr1), mean(ans_brute), tolerance = 0.1)
  expect_equal(median(ans_rskeltr1), median(ans_brute), tolerance = 0.1)
  expect_equal(var(ans_rskeltr1), var(ans_brute), tolerance = 0.1)
  expect_equal(mean(ans_rskeltr1 == 10L), mean(ans_brute == 10L), tolerance = 0.05)
})

test_that("rskeltr1 returns 'lower' when 'p_lower' is very high", {
  set.seed(0)
  ans_obtained <- rskeltr1(mu1 = 5, mu2 = 1, lower = 10000)
  ans_expected <- 10000
  expect_identical(ans_obtained, ans_expected)
})

test_that("R and C versions of rskeltr1 give the same answer", {
  rskeltr1_r <- function(mu1, mu2, lower) {
    max_attempt <- 10L
    n_attempt <- 0L
    found <- FALSE
    while ((n_attempt < max_attempt) && !found) {
      n_attempt <- n_attempt + 1L
      x1 <- stats::rpois(n = 1L, lambda = mu1)
      x2 <- stats::rpois(n = 1L, lambda = mu2)
      prop_value <- x1 - x2
      found <- prop_value >= lower
    }
    if (found) {
      return(prop_value)
    }
    p_lower_max <- 1 - 1e-10
    p_lower <- pskel1(
      q = lower - 1,
      mu1 = mu1,
      mu2 = mu2,
      lower_tail = TRUE,
      log_p = FALSE
    )
    if (p_lower > p_lower_max) {
      return(lower)
    }
    U <- stats::runif(
      n = 1L,
      min = p_lower,
      max = 1
    )
    ans <- qskel1(
      p = U,
      mu1 = mu1,
      mu2 = mu2
    )
    if (ans < lower) {
      ans <- lower
    }
    ans
  }
  set.seed(0)
  mu1 <- runif(n = 100, min = 0, max = 20)
  mu2 <- runif(n = 100, min = 0, max = 20)
  lower <- as.numeric(round(runif(n = 100, min = -5, max = 5)))
  set.seed(0)
  ans_c <- mapply(rskeltr1, mu1 = mu1, mu2 = mu2, lower = lower)
  set.seed(0)
  ans_r <- mapply(rskeltr1_r, mu1 = mu1, mu2 = mu2, lower = lower)
  expect_equal(ans_r, ans_c)
})

test_that("'rskeltr1' returns doubles", {
  expect_true(is.double(rskeltr1(mu1 = 5, mu2 = 5, lower = 3)))
})


## rskeltr --------------------------------------------------------------------

test_that("R and C versions of rskeltr give the same answer", {
  rskeltr_r <- function(n, mu1, mu2, lower) {
    ans <- double(length = n)
    for (i in seq_len(n)) {
      ans[[i]] <- rskeltr1(
        mu1 = mu1[[i]],
        mu2 = mu2[[i]],
        lower = lower[[i]]
      )
    }
    ans
  }
  set.seed(0)
  mu1 <- runif(n = 100, min = 0, max = 20)
  mu2 <- runif(n = 100, min = 0, max = 20)
  lower <- as.numeric(round(runif(n = 100, min = -5, max = 5)))
  set.seed(0)
  ans_c <- rskeltr(n = 100, mu1 = mu1, mu2 = mu2, lower = lower)
  set.seed(0)
  ans_r <- rskeltr_r(n = 100, mu1 = mu1, mu2 = mu2, lower = lower)
  expect_equal(ans_r, ans_c)
})

test_that("'rskeltr1' returns doubles", {
  expect_true(is.double(rskeltr(n = 3, mu1 = 5, mu2 = 5, lower = 3)))
})


## dskel1 ---------------------------------------------------------------------

test_that("'dskel1' gives same answers as 'dskellam' function in 'skellam' package", {
  ans_obtained <- dskel1(
    x = 5,
    mu1 = 2.5,
    mu2 = 1.5,
    use_log = FALSE
  )
  ans_expected <- 0.02715016
  expect_equal(ans_obtained, ans_expected, tolerance = 0.0001)
  ans_obtained <- dskel1(
    x = -5,
    mu1 = 2.5,
    mu2 = 1.5,
    use_log = FALSE
  )
  ans_expected <- 0.002111197
  expect_equal(ans_obtained, ans_expected, tolerance = 0.0001)
  ans_obtained <- dskel1(
    x = 5,
    mu1 = 2.5,
    mu2 = 1.5,
    use_log = TRUE
  )
  ans_expected <- -3.606372
  expect_equal(ans_obtained, ans_expected, tolerance = 0.0001)
  ans_obtained <- dskel1(
    x = -5,
    mu1 = 2.5,
    mu2 = 1.5,
    use_log = TRUE
  )
  ans_expected <- -6.1605
  expect_equal(ans_obtained, ans_expected, tolerance = 0.0001)
})

test_that("'dskel1' gives same answers as 'dpois' when 'mu' is 0", {
  ans_obtained <- dskel1(
    x = -5,
    mu1 = 0,
    mu2 = 1.5,
    use_log = FALSE
  )
  ans_expected <- dpois(
    x = 5,
    lambda = 1.5,
    log = FALSE
  )
  expect_equal(ans_obtained, ans_expected)
  ans_obtained <- dskel1(
    x = 5,
    mu1 = 2.1,
    mu2 = 0,
    use_log = TRUE
  )
  ans_expected <- dpois(
    x = 5,
    lambda = 2.1,
    log = TRUE
  )
  expect_equal(ans_obtained, ans_expected)
})

test_that("'dskel1' copes with large values", {
    ans <- dskel1(x = 300, mu1 = 250, mu2 = 1, use_log = TRUE)
    expect_true(is.finite(ans))
    ans <- dskel1(x = -300, mu1 = 5, mu2 = 400, use_log = TRUE)
    expect_true(is.finite(ans))
    ans <- dskel1(x = 30000, mu1 = 30000, mu2 = 1, use_log = TRUE)
    expect_true(is.finite(ans))
    ans <- dskel1(x = -30000, mu1 = 5, mu2 = 400, use_log = TRUE)
    expect_true(is.finite(ans))
})

test_that("R and C versions of dskel1 give the same answer", {
  dskel1_r <- function(x, mu1, mu2, use_log) {
    if (mu1 == 0) {
      return(stats::dpois(-x, mu2, log = use_log))
    }
    if (mu2 == 0) {
      return(stats::dpois(x, mu1, log = use_log))
    }
    I <- besselI(
      x = 2 * sqrt(mu1) * sqrt(mu2),
      nu = abs(x),
      expon.scaled = FALSE
    )
    if (use_log) {
      ans <- -(mu1 + mu2) + (x / 2) * log(mu1 / mu2) + log(I)
    } else {
      ans <- exp(-(mu1 + mu2)) * (mu1 / mu2)^(x / 2) * I
    }
    ans
  }
  set.seed(0)
  x <- as.numeric(round(runif(n = 100, min = -10, max = 10)))
  mu1 <- runif(n = 100, min = 0, max = 20)
  mu2 <- runif(n = 100, min = 0, max = 20)
  ## mu1, mu2 both non-zero; use_log = TRUE
  ans_c <- mapply(dskel1, x = x, mu1 = mu1, mu2 = mu2, use_log = TRUE)
  ans_r <- mapply(dskel1_r, x = x, mu1 = mu1, mu2 = mu2, use_log = TRUE)
  expect_equal(ans_r, ans_c)
  ## mu1, mu2 both non-zero; use_log = FALSE
  ans_c <- mapply(dskel1, x = x, mu1 = mu1, mu2 = mu2, use_log = FALSE)
  ans_r <- mapply(dskel1_r, x = x, mu1 = mu1, mu2 = mu2, use_log = FALSE)
  expect_equal(ans_r, ans_c)
  ## mu1 = 0
  ans_c <- mapply(dskel1, x = x, mu1 = 0, mu2 = mu2, use_log = TRUE)
  ans_r <- mapply(dskel1_r, x = x, mu1 = 0, mu2 = mu2, use_log = TRUE)
  expect_equal(ans_r, ans_c)
  ## mu2 = 0
  ans_c <- mapply(dskel1, x = x, mu1 = mu1, mu2 = 0, use_log = TRUE)
  ans_r <- mapply(dskel1_r, x = x, mu1 = mu1, mu2 = 0, use_log = TRUE)
  expect_equal(ans_r, ans_c)
})


## dskeltr1 -------------------------------------------------------------------

test_that("'dskeltr1' values add up to 1 - lower is negative", {
  ans <- sapply(seq(from = -10, to = 40),
    dskeltr1,
    mu1 = 3.1,
    mu2 = 1.3,
    lower = -2,
    use_log = FALSE
  )
  expect_equal(sum(ans), 1)
})

test_that("'dskeltr1' values add up to 1 - lower is positive", {
  ans <- sapply(seq(from = -10, to = 40),
    dskeltr1,
    mu1 = 3.1,
    mu2 = 1.3,
    lower = 3,
    use_log = FALSE
  )
  expect_equal(sum(ans), 1)
})

test_that("'dskeltr1' gives same answer as dskel1 when lower is very small", {
  ans_dskeltr1 <- dskeltr1(
    x = 33,
    mu1 = 40,
    mu2 = 20,
    lower = -100,
    use_log = TRUE
  )
  ans_dskel1 <- dskel1(
    x = 33,
    mu1 = 40,
    mu2 = 20,
    use_log = TRUE
  )
  expect_equal(ans_dskeltr1, ans_dskel1)
})

test_that("'dskeltr1' gives higher value than dskel1 when lower is moderate", {
  ans_dskeltr1 <- dskeltr1(
    x = 33,
    mu1 = 40,
    mu2 = 20,
    lower = 5,
    use_log = TRUE
  )
  ans_dskel1 <- dskel1(
    x = 33,
    mu1 = 40,
    mu2 = 20,
    use_log = TRUE
  )
  expect_true(ans_dskeltr1 > ans_dskel1)
})

test_that("'dskeltr1' flips from 0 to positive when use_log is FALSE and x equals 'lower'", {
  val <- sapply(1:10,
    dskeltr1,
    mu1 = 40,
    mu2 = 20,
    lower = 6,
    use_log = FALSE
  )
  ans_obtained <- val > 0
  ans_expected <- rep(c(FALSE, TRUE), each = 5)
  expect_identical(ans_obtained, ans_expected)
})

test_that("'dskeltr1' flips from infinite to finite when use_log is TRUE and x equals 'lower'", {
  val <- sapply(1:10,
    dskeltr1,
    mu1 = 40,
    mu2 = 20,
    lower = 6,
    use_log = TRUE
  )
  ans_obtained <- is.finite(val)
  ans_expected <- rep(c(FALSE, TRUE), each = 5)
  expect_identical(ans_obtained, ans_expected)
})

test_that("R and C versions of dskeltr1 give the same answer", {
  dskeltr1_r <- function(x, mu1, mu2, lower, use_log) {
    if (x < lower) {
      ans <- if (use_log) -Inf else 0
      return(ans)
    }
    log_dens_no_trunc <- dskel1(
      x = x,
      mu1 = mu1,
      mu2 = mu2,
      use_log = TRUE
    )
    log_const <- pskel1(
      q = lower - 1,
      mu1 = mu1,
      mu2 = mu2,
      lower_tail = FALSE,
      log_p = TRUE
    )
    log_dens_trunc <- log_dens_no_trunc - log_const
    ans <- if (use_log) log_dens_trunc else exp(log_dens_trunc)
    ans
  }
  set.seed(0)
  x <- as.numeric(round(runif(n = 100, min = -10, max = 10)))
  mu1 <- runif(n = 100, min = 0, max = 20)
  mu2 <- runif(n = 100, min = 0, max = 20)
  lower <- as.numeric(round(runif(n = 100, min = -10, max = 5)))
  use_log <- runif(n = 100) < 0.5
  ans_c <- mapply(dskeltr1, x = x, mu1 = mu1, mu2 = mu2, lower = lower, use_log = use_log)
  ans_r <- mapply(dskeltr1_r, x = x, mu1 = mu1, mu2 = mu2, lower = lower, use_log = use_log)
  expect_equal(ans_r, ans_c)
})

## dkseltr --------------------------------------------------------------------

test_that("R and C versions of dskeltr give the same answer", {
  dskeltr_r <- function(x, mu1, mu2, lower, use_log) {
    n <- length(x)
    ans <- numeric(length = n)
    for (i in seq_len(n)) {
      ans[i] <- dskeltr1(
        x = x[i],
        mu1 = mu1[i],
        mu2 = mu2[i],
        lower = lower[i],
        use_log = use_log
      )
    }
    ans
  }
  set.seed(0)
  x <- as.numeric(round(runif(n = 100, min = -10, max = 10)))
  mu1 <- runif(n = 100, min = 0, max = 20)
  mu2 <- runif(n = 100, min = 0, max = 20)
  lower <- as.numeric(round(runif(n = 100, min = -15, max = 5)))
  ans_c <- dskeltr(x = x, mu1 = mu1, mu2 = mu2, lower = lower, use_log = TRUE)
  ans_r <- dskeltr_r(x = x, mu1 = mu1, mu2 = mu2, lower = lower, use_log = TRUE)
  expect_equal(ans_r, ans_c)
})




    

    
