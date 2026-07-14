test_that("transform_counts reproduces the aggregate rate as the arm mean of W", {
  set.seed(1)
  y   <- c(0, 2, 1, 3, 0, 5)
  exp <- c(1.0, 0.8, 1.2, 0.9, 1.1, 1.3)
  trt <- factor(c("a", "a", "a", "b", "b", "b"))

  tr <- EmpiricalRate:::transform_counts(y, exp, trt)

  # Arm-specific mean exposure.
  expect_equal(unname(tr$Tbar[["a"]]), mean(exp[trt == "a"]))
  expect_equal(unname(tr$Tbar[["b"]]), mean(exp[trt == "b"]))

  # Key identity: mean of W within an arm == sum(y) / sum(exposure).
  Wa <- tr$W[trt == "a"]
  Wb <- tr$W[trt == "b"]
  expect_equal(mean(Wa), sum(y[trt == "a"]) / sum(exp[trt == "a"]))
  expect_equal(mean(Wb), sum(y[trt == "b"]) / sum(exp[trt == "b"]))
})

test_that("transform_counts collapses to y/exposure when exposure is constant", {
  y   <- c(1, 0, 4, 2)
  exp <- rep(2, 4)
  trt <- factor(c("a", "a", "b", "b"))
  tr  <- EmpiricalRate:::transform_counts(y, exp, trt)
  expect_equal(tr$W, y / 2)
})

test_that("delta_lograte matches the closed-form Jacobian sandwich", {
  r_hat <- c(a = 0.5, b = 0.2)
  V_r   <- matrix(c(0.04, 0.01,
                    0.01, 0.02), nrow = 2,
                  dimnames = list(c("a", "b"), c("a", "b")))
  dl <- EmpiricalRate:::delta_lograte(r_hat, V_r)

  expect_equal(dl$mu_hat, log(r_hat))
  # V_mu[i,j] = V_r[i,j] / (r_i r_j).
  expect_equal(dl$V_mu[1, 1], V_r[1, 1] / (r_hat[["a"]]^2))
  expect_equal(dl$V_mu[2, 2], V_r[2, 2] / (r_hat[["b"]]^2))
  expect_equal(dl$V_mu[1, 2], V_r[1, 2] / (r_hat[["a"]] * r_hat[["b"]]))
})

test_that("delta_lograte errors on non-positive rates", {
  expect_error(
    EmpiricalRate:::delta_lograte(c(a = 0.5, b = 0), diag(2)),
    "strictly positive"
  )
  expect_error(
    EmpiricalRate:::delta_lograte(c(a = 0.5, b = -0.1), diag(2)),
    "strictly positive"
  )
})

test_that("pairwise_rr computes RR, CI, z and p correctly for two arms", {
  mu   <- c(a = log(0.4), b = log(0.1))
  V_mu <- matrix(c(0.05, 0.00,
                   0.00, 0.09), nrow = 2,
                 dimnames = list(c("a", "b"), c("a", "b")))
  out <- EmpiricalRate:::pairwise_rr(mu, V_mu, alpha = 0.05)

  expect_equal(nrow(out), 1L)
  expect_equal(out$rate_ratio, exp(mu[["a"]] - mu[["b"]]))

  se <- sqrt(V_mu[1, 1] + V_mu[2, 2])
  z  <- (mu[["a"]] - mu[["b"]]) / se
  expect_equal(out$se_log_rr, se)
  expect_equal(out$z, z)
  expect_equal(out$p_value, 2 * pnorm(-abs(z)))

  zc <- qnorm(0.975)
  expect_equal(out$conf.low,  exp((mu[["a"]] - mu[["b"]]) - zc * se))
  expect_equal(out$conf.high, exp((mu[["a"]] - mu[["b"]]) + zc * se))
})

test_that("pairwise_rr with a reference returns other/reference comparisons only", {
  mu   <- c(ctrl = log(0.3), t1 = log(0.15), t2 = log(0.6))
  V_mu <- diag(c(0.04, 0.05, 0.06))
  dimnames(V_mu) <- list(names(mu), names(mu))

  out <- EmpiricalRate:::pairwise_rr(mu, V_mu, reference = "ctrl")
  expect_equal(nrow(out), 2L)
  expect_true(all(out$group2 == "ctrl"))
  # RR = other / ctrl.
  expect_equal(out$rate_ratio[out$group1 == "t1"], exp(mu[["t1"]] - mu[["ctrl"]]))

  # Full pairwise (no reference) gives choose(3, 2) = 3 rows.
  out_all <- EmpiricalRate:::pairwise_rr(mu, V_mu)
  expect_equal(nrow(out_all), 3L)
})

test_that("pairwise_rr errors on an unknown reference", {
  mu <- c(a = 0, b = -1)
  expect_error(
    EmpiricalRate:::pairwise_rr(mu, diag(2), reference = "zzz"),
    "not one of the treatment"
  )
})
