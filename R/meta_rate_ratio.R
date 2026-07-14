#' Pool rate ratios across strata or studies
#'
#' Combines stratum-specific (or study-specific) rate ratios into a single
#' pooled rate ratio using exposure weights, as described in Appendix A of the
#' companion paper. This is useful when event rates differ across strata (for
#' example ordinal disease-severity groups or separate studies) while the rate
#' ratio is assumed common, a situation in which a stratified linear-model
#' analysis need not be the most efficient estimator.
#'
#' @param log_rr Numeric vector of stratum-specific log rate ratios (for
#'   example the `log_rr` column returned by [empirical_rate()] for one
#'   comparison, one entry per stratum).
#' @param se_log_rr Numeric vector of the corresponding standard errors on the
#'   log scale (the `se_log_rr` column from [empirical_rate()]).
#' @param weight Numeric vector of stratum weights. Per Appendix A these are the
#'   total exposures (total follow-up time) within each stratum. Must be
#'   positive.
#' @param scale Character scalar giving the scale on which strata are pooled.
#'   `"log"` (the default) pools the log rate ratios --- the practice
#'   recommended in the paper's Discussion and the more standard choice.
#'   `"identity"` pools the rate ratios on their natural scale and reproduces
#'   the formula displayed in Appendix A exactly (natural-scale variances are
#'   obtained from `se_log_rr` by the delta method).
#' @param alpha Numeric in `(0, 1)`; significance level for the two-sided
#'   confidence interval (default `0.05`).
#'
#' @details
#' With weights \eqn{w_s} and stratum estimates \eqn{\hat\lambda_s}, the pooled
#' estimator on the identity scale is
#' \deqn{\hat\lambda = \frac{\sum_s w_s \hat\lambda_s}{\sum_s w_s}, \qquad
#'   \widehat{\mathrm{Var}}(\hat\lambda) =
#'   \frac{\sum_s w_s^2\,\widehat{\mathrm{Var}}(\hat\lambda_s)}{(\sum_s w_s)^2}.}
#' The log scale (`scale = "log"`) applies the same weighted-average and
#' variance formulas to \eqn{\log\hat\lambda_s} and then exponentiates.
#'
#' @return A one-row data frame with columns `rate_ratio`, `conf.low`,
#'   `conf.high`, `log_rr`, `se_log_rr`, and `scale`.
#'
#' @seealso [empirical_rate()].
#'
#' @references
#' Sun J, Amoafo L, Qu Y (2026). An Empirical Method for Analyzing Count Data,
#' Appendix A.
#'
#' @examples
#' # Three strata, each with a log rate ratio, its SE, and a total-exposure
#' # weight (e.g. patient-years in the stratum).
#' log_rr    <- c(0.41, 0.55, 0.29)
#' se_log_rr <- c(0.18, 0.26, 0.15)
#' weight    <- c(320, 150, 410)
#'
#' meta_rate_ratio(log_rr, se_log_rr, weight)                    # log scale
#' meta_rate_ratio(log_rr, se_log_rr, weight, scale = "identity")
#'
#' @export
meta_rate_ratio <- function(log_rr, se_log_rr, weight,
                            scale = c("log", "identity"),
                            alpha = 0.05) {
  scale <- match.arg(scale)

  if (!is.numeric(log_rr) || !is.numeric(se_log_rr) || !is.numeric(weight)) {
    stop("`log_rr`, `se_log_rr`, and `weight` must all be numeric.",
         call. = FALSE)
  }
  n <- length(log_rr)
  if (length(se_log_rr) != n || length(weight) != n) {
    stop("`log_rr`, `se_log_rr`, and `weight` must have the same length ",
         "(one entry per stratum).", call. = FALSE)
  }
  if (n < 1L) stop("At least one stratum is required.", call. = FALSE)
  if (any(!is.finite(log_rr)) || any(!is.finite(se_log_rr)) ||
      any(!is.finite(weight))) {
    stop("`log_rr`, `se_log_rr`, and `weight` must be finite.", call. = FALSE)
  }
  if (any(se_log_rr < 0)) {
    stop("`se_log_rr` values must be non-negative.", call. = FALSE)
  }
  if (any(weight <= 0)) {
    stop("`weight` values must be strictly positive.", call. = FALSE)
  }
  if (!is.numeric(alpha) || length(alpha) != 1L || alpha <= 0 || alpha >= 1) {
    stop("`alpha` must be a single number strictly between 0 and 1.",
         call. = FALSE)
  }

  z_crit <- stats::qnorm(1 - alpha / 2)
  W <- sum(weight)

  if (scale == "log") {
    est <- sum(weight * log_rr) / W
    var_est <- sum(weight^2 * se_log_rr^2) / W^2
    se_est <- sqrt(var_est)
    out <- data.frame(
      rate_ratio = exp(est),
      conf.low   = exp(est - z_crit * se_est),
      conf.high  = exp(est + z_crit * se_est),
      log_rr     = est,
      se_log_rr  = se_est,
      scale      = "log",
      stringsAsFactors = FALSE
    )
  } else {
    rr_s <- exp(log_rr)
    # delta method: Var(RR_s) = RR_s^2 * Var(log RR_s)
    var_rr_s <- rr_s^2 * se_log_rr^2
    est <- sum(weight * rr_s) / W
    var_est <- sum(weight^2 * var_rr_s) / W^2
    se_est <- sqrt(var_est)
    out <- data.frame(
      rate_ratio = est,
      conf.low   = est - z_crit * se_est,
      conf.high  = est + z_crit * se_est,
      log_rr     = log(est),
      se_log_rr  = se_est / est,   # approx SE of log(pooled RR), delta method
      scale      = "identity",
      stringsAsFactors = FALSE
    )
  }
  rownames(out) <- NULL
  out
}
