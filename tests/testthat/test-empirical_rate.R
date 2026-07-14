# These tests exercise the full pipeline and therefore require RobinCar.
# They are skipped automatically where RobinCar is unavailable.

make_toy <- function() {
  data.frame(
    arm = factor(rep(c("Control", "Treatment"), each = 8),
                 levels = c("Control", "Treatment")),
    events = c(0, 1, 2, 0, 3, 1, 0, 2,      # Control: 9 events
               1, 0, 0, 1, 0, 2, 0, 0),      # Treatment: 4 events
    exposure_years = c(1.0, 0.8, 1.2, 0.9, 1.1, 1.3, 0.7, 1.0,
                       1.2, 0.9, 1.0, 1.1, 0.8, 1.3, 1.0, 0.9),
    x = c(7.5, 8.1, 8.9, 7.2, 9.0, 8.3, 7.8, 8.0,
          8.2, 7.9, 8.4, 7.6, 8.8, 7.4, 8.1, 7.7)
  )
}

test_that("unadjusted marginal rates equal the observed aggregate rates (core identity)", {
  skip_if_not_installed("RobinCar")
  d <- make_toy()

  fit <- suppressWarnings(
    empirical_rate(d, treatment = "arm", count = "events",
                   exposure = "exposure_years")
  )

  agg <- with(d, tapply(events, arm, sum) / tapply(exposure_years, arm, sum))
  est <- setNames(fit$rates$rate, fit$rates$treatment)

  expect_equal(est[["Control"]],   agg[["Control"]],   tolerance = 1e-8)
  expect_equal(est[["Treatment"]], agg[["Treatment"]], tolerance = 1e-8)
})

test_that("unadjusted rate ratio equals the raw rate ratio exactly (core identity)", {
  skip_if_not_installed("RobinCar")
  d <- make_toy()

  fit <- suppressWarnings(
    empirical_rate(d, treatment = "arm", count = "events",
                   exposure = "exposure_years", reference = "Control")
  )
  agg    <- with(d, tapply(events, arm, sum) / tapply(exposure_years, arm, sum))
  raw_rr <- agg[["Treatment"]] / agg[["Control"]]

  cmp <- fit$comparisons
  expect_equal(nrow(cmp), 1L)
  expect_equal(cmp$rate_ratio, unname(raw_rr), tolerance = 1e-8)
})

test_that("empirical_rate returns a well-formed object", {
  skip_if_not_installed("RobinCar")
  d <- make_toy()
  fit <- suppressWarnings(
    empirical_rate(d, treatment = "arm", count = "events",
                   exposure = "exposure_years")
  )

  expect_s3_class(fit, "empirical_rate")
  expect_setequal(
    names(fit),
    c("rates", "comparisons", "vcov_rate", "log_rates", "vcov_lograte",
      "robincar", "call", "formula", "covariates", "car_scheme", "strata",
      "reference", "alpha")
  )
  expect_equal(dim(fit$vcov_rate), c(2L, 2L))
  expect_true(all(fit$rates$rate > 0))
  expect_true(all(is.finite(fit$rates$se)))
  # log rates are consistent with the rate estimates.
  rate_vec <- setNames(fit$rates$rate, fit$rates$treatment)
  expect_equal(unname(fit$log_rates),
               unname(log(rate_vec[names(fit$log_rates)])),
               tolerance = 1e-10)
})

test_that("adjusted analysis runs and yields finite positive rates", {
  skip_if_not_installed("RobinCar")
  d <- make_toy()
  fit <- suppressWarnings(
    empirical_rate(d, treatment = "arm", count = "events",
                   exposure = "exposure_years", covariates = "x",
                   reference = "Control")
  )
  expect_true(all(fit$rates$rate > 0))
  expect_true(is.finite(fit$comparisons$rate_ratio))
  expect_true(grepl("x", paste(deparse(fit$formula), collapse = " ")))
})

test_that("reference orientation makes RR = other / reference", {
  skip_if_not_installed("RobinCar")
  d <- make_toy()
  fit <- suppressWarnings(
    empirical_rate(d, treatment = "arm", count = "events",
                   exposure = "exposure_years", reference = "Control")
  )
  expect_equal(fit$comparisons$group1, "Treatment")
  expect_equal(fit$comparisons$group2, "Control")
})

test_that("print and summary methods run without error", {
  skip_if_not_installed("RobinCar")
  d <- make_toy()
  fit <- suppressWarnings(
    empirical_rate(d, treatment = "arm", count = "events",
                   exposure = "exposure_years")
  )
  expect_output(print(fit), "Empirical marginal event-rate analysis")
  s <- summary(fit)
  expect_s3_class(s, "summary.empirical_rate")
  expect_output(print(s), "Pairwise rate ratios")
})

test_that("input validation catches misuse", {
  d <- make_toy()

  expect_error(empirical_rate(as.list(d), "arm", "events", "exposure_years"),
               "must be a data frame")
  expect_error(empirical_rate(d, "nope", "events", "exposure_years"),
               "not found")
  expect_error(
    empirical_rate(transform(d, events = -events), "arm", "events",
                   "exposure_years"),
    "non-negative"
  )
  expect_error(
    empirical_rate(transform(d, exposure_years = 0 * exposure_years),
                   "arm", "events", "exposure_years"),
    "strictly positive"
  )
  expect_error(
    empirical_rate(d, "arm", "events", "exposure_years", alpha = 0),
    "between 0 and 1"
  )
  expect_error(
    empirical_rate(d, "arm", "events", "exposure_years",
                   car_scheme = "permuted-block"),
    "requires `strata`"
  )
})
