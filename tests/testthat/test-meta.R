test_that("single stratum reproduces the Wald CI on the log scale", {
  lrr <- 0.5
  se  <- 0.2
  out <- meta_rate_ratio(lrr, se, weight = 100)  # default scale = "log"

  expect_equal(out$rate_ratio, exp(lrr))
  expect_equal(out$log_rr, lrr)
  expect_equal(out$se_log_rr, se)
  zc <- qnorm(0.975)
  expect_equal(out$conf.low,  exp(lrr - zc * se))
  expect_equal(out$conf.high, exp(lrr + zc * se))
  expect_identical(out$scale, "log")
})

test_that("identity scale reproduces the Appendix A weighted-average formula", {
  lrr <- c(0.41, 0.55, 0.29)
  se  <- c(0.18, 0.26, 0.15)
  w   <- c(320, 150, 410)

  out <- meta_rate_ratio(lrr, se, w, scale = "identity")

  rr_s     <- exp(lrr)
  var_rr_s <- rr_s^2 * se^2               # delta method on each stratum
  est      <- sum(w * rr_s) / sum(w)
  var_est  <- sum(w^2 * var_rr_s) / sum(w)^2

  expect_equal(out$rate_ratio, est)
  zc <- qnorm(0.975)
  expect_equal(out$conf.low,  est - zc * sqrt(var_est))
  expect_equal(out$conf.high, est + zc * sqrt(var_est))
  expect_identical(out$scale, "identity")
})

test_that("log scale reproduces the weighted-average-of-logs formula", {
  lrr <- c(0.41, 0.55, 0.29)
  se  <- c(0.18, 0.26, 0.15)
  w   <- c(320, 150, 410)

  out <- meta_rate_ratio(lrr, se, w, scale = "log")

  est     <- sum(w * lrr) / sum(w)
  var_est <- sum(w^2 * se^2) / sum(w)^2
  expect_equal(out$log_rr, est)
  expect_equal(out$se_log_rr, sqrt(var_est))
  expect_equal(out$rate_ratio, exp(est))
})

test_that("identical strata give the same pooled estimate on both scales", {
  lrr <- rep(0.3, 4)
  se  <- rep(0.2, 4)
  w   <- c(10, 20, 30, 40)
  out_log <- meta_rate_ratio(lrr, se, w, scale = "log")
  out_id  <- meta_rate_ratio(lrr, se, w, scale = "identity")
  expect_equal(out_log$rate_ratio, exp(0.3))
  expect_equal(out_id$rate_ratio, exp(0.3))
})

test_that("meta_rate_ratio validates its inputs", {
  expect_error(meta_rate_ratio(c(0.1, 0.2), c(0.1), c(1, 1)), "same length")
  expect_error(meta_rate_ratio(0.1, 0.2, -1), "strictly positive")
  expect_error(meta_rate_ratio(0.1, -0.2, 1), "non-negative")
  expect_error(meta_rate_ratio(0.1, 0.2, 1, alpha = 1.5), "between 0 and 1")
})
