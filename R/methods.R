# S3 methods for objects of class "empirical_rate" --------------------------

#' @param x An object of class `"empirical_rate"` (or its summary).
#' @param object An object of class `"empirical_rate"`.
#' @param digits Number of significant digits for printing.
#' @param ... Currently ignored.
#' @rdname empirical_rate
#' @export
print.empirical_rate <- function(x, digits = 4, ...) {
  cat("Empirical marginal event-rate analysis\n")
  cat("Call: ", paste(deparse(x$call), collapse = " "), "\n\n", sep = "")

  adj <- if (is.null(x$covariates)) {
    "unadjusted"
  } else {
    paste0("adjusted for ", paste(x$covariates, collapse = ", "))
  }
  cat("Analysis: ", adj, "  |  working model: ",
      paste(deparse(x$formula), collapse = " "), "\n", sep = "")
  cat("Variance: car_scheme = \"", x$car_scheme, "\"", sep = "")
  if (!is.null(x$strata)) {
    cat("  (strata: ", paste(x$strata, collapse = ", "), ")", sep = "")
  }
  cat("\n\n")

  rates <- x$rates
  rates$rate <- signif(rates$rate, digits)
  rates$se   <- signif(rates$se, digits)
  cat("Marginal event rates:\n")
  print(rates, row.names = FALSE)

  cat("\nRate ratios (", format(100 * (1 - x$alpha)), "% CI):\n", sep = "")
  print(format_comparisons(x$comparisons, digits), row.names = FALSE)

  invisible(x)
}

#' @rdname empirical_rate
#' @export
summary.empirical_rate <- function(object, ...) {
  structure(
    list(
      call        = object$call,
      formula     = object$formula,
      covariates  = object$covariates,
      car_scheme  = object$car_scheme,
      strata      = object$strata,
      reference   = object$reference,
      alpha       = object$alpha,
      rates       = object$rates,
      comparisons = object$comparisons
    ),
    class = "summary.empirical_rate"
  )
}

#' @rdname empirical_rate
#' @export
print.summary.empirical_rate <- function(x, digits = 4, ...) {
  cat("Empirical marginal event-rate analysis\n")
  cat(strrep("-", 55), "\n", sep = "")
  cat("Call: ", paste(deparse(x$call), collapse = " "), "\n\n", sep = "")

  if (is.null(x$covariates)) {
    cat("Covariate adjustment : none (unadjusted)\n")
  } else {
    cat("Covariate adjustment :", paste(x$covariates, collapse = ", "), "\n")
  }
  cat("Working model        :", paste(deparse(x$formula), collapse = " "), "\n")
  cat("Randomization scheme :", x$car_scheme)
  if (!is.null(x$strata)) cat(" (strata:", paste(x$strata, collapse = ", "), ")")
  cat("\n")
  if (!is.null(x$reference)) {
    cat("Reference arm        :", x$reference, "\n")
  }
  cat("Confidence level     :", paste0(format(100 * (1 - x$alpha)), "%"), "\n\n")

  rates <- x$rates
  rates$rate <- signif(rates$rate, digits)
  rates$se   <- signif(rates$se, digits)
  cat("Marginal event rates (events per unit exposure):\n")
  print(rates, row.names = FALSE)

  cat("\nPairwise rate ratios:\n")
  print(format_comparisons(x$comparisons, digits, full = TRUE),
        row.names = FALSE)

  cat("\nEach rate ratio is rate(group1) / rate(group2); tests are two-sided",
      "Wald tests on the log-rate scale.\n")
  invisible(x)
}

# Internal: pretty-format the comparisons data frame for printing.
format_comparisons <- function(cmp, digits = 4, full = FALSE) {
  out <- data.frame(
    comparison = cmp$comparison,
    rate_ratio = signif(cmp$rate_ratio, digits),
    conf.low   = signif(cmp$conf.low, digits),
    conf.high  = signif(cmp$conf.high, digits),
    p_value    = signif(cmp$p_value, digits),
    stringsAsFactors = FALSE
  )
  if (full) {
    out$log_rr    <- signif(cmp$log_rr, digits)
    out$se_log_rr <- signif(cmp$se_log_rr, digits)
    out$z         <- signif(cmp$z, digits)
    out <- out[, c("comparison", "rate_ratio", "conf.low", "conf.high",
                   "log_rr", "se_log_rr", "z", "p_value")]
  }
  out
}
