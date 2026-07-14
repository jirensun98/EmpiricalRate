# Internal helpers -------------------------------------------------------------
# These functions contain the statistical core of the empirical method and are
# deliberately kept free of any dependency on RobinCar so that they can be unit
# tested in isolation. They are not exported.

#' Transform event counts to subject-level rate quantities
#'
#' Computes the arm-specific mean exposure \eqn{\bar d_{i\cdot}} and the
#' transformed outcome \eqn{W_{ij} = Y_{ij} / \bar d_{i\cdot}}.
#'
#' @param y Numeric vector of event counts.
#' @param exposure Numeric vector of exposures / follow-up times.
#' @param trt Factor or vector identifying the treatment arm.
#'
#' @return A list with elements `W` (the transformed outcome) and `Tbar` (a
#'   named numeric vector of arm-specific mean exposures).
#' @keywords internal
#' @noRd
transform_counts <- function(y, exposure, trt) {
  trt <- as.factor(trt)
  Tbar <- tapply(exposure, trt, mean)
  W <- y / Tbar[as.character(trt)]
  list(W = as.numeric(W), Tbar = Tbar)
}

#' Delta-method variance of log rates
#'
#' Given marginal rate estimates \eqn{\hat r} and their variance-covariance
#' matrix \eqn{V_r}, returns the log rates and the delta-method
#' variance-covariance matrix
#' \eqn{V_\mu = \mathrm{diag}(1/\hat r)\, V_r\, \mathrm{diag}(1/\hat r)}.
#'
#' @param r_hat Named numeric vector of marginal rate estimates (all > 0).
#' @param V_r Variance-covariance matrix of `r_hat`.
#'
#' @return A list with `mu_hat` (log rates) and `V_mu` (their vcov).
#' @keywords internal
#' @noRd
delta_lograte <- function(r_hat, V_r) {
  if (any(!is.finite(r_hat)) || any(r_hat <= 0)) {
    stop("All estimated rates must be finite and strictly positive to work on ",
         "the log-rate scale. Estimated rates: ",
         paste(format(r_hat), collapse = ", "), ".", call. = FALSE)
  }
  mu_hat <- log(r_hat)
  J <- diag(1 / r_hat, nrow = length(r_hat))
  V_mu <- t(J) %*% V_r %*% J
  dimnames(V_mu) <- list(names(r_hat), names(r_hat))
  list(mu_hat = mu_hat, V_mu = V_mu)
}

#' Pairwise rate ratios with Wald inference on the log scale
#'
#' @param mu_hat Named numeric vector of log rates.
#' @param V_mu Variance-covariance matrix of `mu_hat`.
#' @param alpha Significance level for two-sided confidence intervals.
#' @param reference Optional name of the reference arm. If supplied, only
#'   comparisons of each other arm versus the reference are returned (rate ratio
#'   = other / reference). If `NULL`, all unordered pairs are returned.
#'
#' @return A data frame of comparisons.
#' @keywords internal
#' @noRd
pairwise_rr <- function(mu_hat, V_mu, alpha = 0.05, reference = NULL) {
  lev <- names(mu_hat)
  K <- length(lev)
  z_crit <- stats::qnorm(1 - alpha / 2)

  # Build the list of (i, k) index pairs, RR = rate_i / rate_k.
  if (is.null(reference)) {
    pairs <- utils::combn(K, 2, simplify = FALSE)
  } else {
    if (!reference %in% lev) {
      stop("`reference` (\"", reference, "\") is not one of the treatment ",
           "levels: ", paste(lev, collapse = ", "), ".", call. = FALSE)
    }
    kref <- match(reference, lev)
    others <- setdiff(seq_len(K), kref)
    pairs <- lapply(others, function(i) c(i, kref))
  }

  rows <- lapply(pairs, function(p) {
    i <- p[1L]; k <- p[2L]
    diff_mu <- mu_hat[i] - mu_hat[k]
    var_diff <- V_mu[i, i] + V_mu[k, k] - 2 * V_mu[i, k]
    se_diff <- sqrt(var_diff)
    z <- diff_mu / se_diff
    data.frame(
      group1    = lev[i],
      group2    = lev[k],
      rate_ratio = exp(diff_mu),
      conf.low  = exp(diff_mu - z_crit * se_diff),
      conf.high = exp(diff_mu + z_crit * se_diff),
      log_rr    = as.numeric(diff_mu),
      se_log_rr = as.numeric(se_diff),
      z         = as.numeric(z),
      p_value   = 2 * stats::pnorm(-abs(z)),
      stringsAsFactors = FALSE,
      row.names = NULL
    )
  })

  out <- do.call(rbind, rows)
  out$comparison <- paste(out$group1, "vs", out$group2)
  out <- out[, c("comparison", "group1", "group2", "rate_ratio",
                 "conf.low", "conf.high", "log_rr", "se_log_rr",
                 "z", "p_value")]
  rownames(out) <- NULL
  out
}
