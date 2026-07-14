#' Empirical estimation of marginal event rates and rate ratios
#'
#' Fits the empirical (distribution-free) method for comparing marginal event
#' rates between treatment groups when the endpoint is a count. Event counts are
#' transformed to subject-level rate quantities, a linear working model is fitted
#' with robust variance estimation via [RobinCar::robincar_glm()], and marginal
#' rates are back-transformed to obtain rate ratios with confidence intervals via
#' the delta method.
#'
#' @param data A data frame containing one row per subject.
#' @param treatment Character scalar. Name of the treatment-group column in
#'   `data`. It is coerced to a factor; the order of its levels determines the
#'   direction of pairwise rate ratios (see *Value*).
#' @param count Character scalar. Name of the (non-negative) event-count column.
#' @param exposure Character scalar. Name of the exposure / follow-up-time
#'   column (strictly positive, e.g. patient-years).
#' @param covariates Optional character vector of baseline-covariate column
#'   names to adjust for. `NULL` (the default) fits an unadjusted analysis.
#' @param car_scheme Character scalar passed to [RobinCar::robincar_glm()],
#'   describing the randomization scheme used in the trial. One of `"simple"`
#'   (the default), `"permuted-block"`, `"pocock-simon"`, `"biased-coin"`, or
#'   `"urn"`. Any value other than `"simple"` requires `strata`.
#' @param strata Optional character vector of randomization-stratification
#'   column names. Required when `car_scheme != "simple"`; when supplied it makes
#'   the variance estimator account for the stratified randomization scheme.
#' @param reference Optional character scalar naming the reference treatment
#'   level. If supplied, only comparisons of each other arm *versus the
#'   reference* are returned (rate ratio = other / reference). If `NULL` (the
#'   default) all unordered pairwise comparisons are returned.
#' @param alpha Numeric in `(0, 1)`. Significance level for two-sided confidence
#'   intervals (default `0.05`, i.e. 95% intervals).
#'
#' @details
#' For subject \eqn{j} in arm \eqn{i} the event count \eqn{Y_{ij}} and exposure
#' \eqn{d_{ij}} are combined into
#' \deqn{W_{ij} = Y_{ij} / \bar d_{i\cdot},}
#' where \eqn{\bar d_{i\cdot}} is the mean exposure in arm \eqn{i}. Because
#' \eqn{E[Y_{ij}] = d_{ij} r_i}, the arm-specific sample mean of \eqn{W_{ij}} is
#' exactly unbiased for the marginal rate \eqn{r_i} without any distributional
#' assumption on \eqn{Y_{ij}}. A linear working model `W ~ treatment (+
#' covariates)` is fitted and the resulting covariate-adjusted arm means
#' \eqn{\hat r_i} and their robust variance matrix are back-transformed to the
#' log-rate scale, where the rate ratio \eqn{\lambda_{ik} = r_i / r_k} is tested
#' and interval-estimated.
#'
#' In the **unadjusted** analysis (`covariates = NULL`, `car_scheme =
#' "simple"`) the estimated marginal rate for each arm equals the observed
#' aggregate rate \eqn{\sum_j Y_{ij} / \sum_j d_{ij}} exactly, and the point
#' estimate of the rate ratio equals the ratio of observed rates.
#'
#' **Working model.** Following the reference implementation, the linear
#' predictor is additive in the covariates (a common-slope ANCOVA working
#' model). Randomization-aware variance estimation --- referred to as ANHECOVA
#' in the companion paper --- is obtained by setting `car_scheme` and `strata`;
#' the working-model *form* is unchanged. Treatment-by-covariate interactions
#' (heterogeneous slopes) are not fitted.
#'
#' @return An object of class `"empirical_rate"`: a list with elements
#'   \describe{
#'     \item{`rates`}{Data frame of marginal rate estimates: `treatment`,
#'       `rate`, `se`.}
#'     \item{`comparisons`}{Data frame of rate-ratio comparisons with columns
#'       `comparison`, `group1`, `group2`, `rate_ratio`, `conf.low`,
#'       `conf.high`, `log_rr`, `se_log_rr`, `z`, `p_value`. Each rate ratio is
#'       `rate(group1) / rate(group2)`.}
#'     \item{`vcov_rate`}{Variance-covariance matrix of the rate estimates.}
#'     \item{`log_rates`, `vcov_lograte`}{Log rates and their delta-method
#'       variance-covariance matrix.}
#'     \item{`robincar`}{The full [RobinCar::robincar_glm()] fit, for advanced
#'       use.}
#'     \item{`call`, `formula`, `covariates`, `car_scheme`, `strata`,
#'       `reference`, `alpha`}{Details of the call.}
#'   }
#'
#' @seealso [meta_rate_ratio()] to pool rate ratios across strata or studies.
#'
#' @references
#' Sun J, Amoafo L, Qu Y (2026). An Empirical Method for Analyzing Count Data.
#'
#' @examplesIf requireNamespace("RobinCar", quietly = TRUE)
#' data(hypo_sim)
#'
#' # Unadjusted analysis: estimated rates equal the observed aggregate rates.
#' fit <- empirical_rate(
#'   hypo_sim,
#'   treatment = "arm",
#'   count     = "events",
#'   exposure  = "exposure_years"
#' )
#' fit
#'
#' # Adjusted analysis with the treatment arm oriented against a reference.
#' fit_adj <- empirical_rate(
#'   hypo_sim,
#'   treatment  = "arm",
#'   count      = "events",
#'   exposure   = "exposure_years",
#'   covariates = c("baseline_hba1c", "baseline_rate"),
#'   reference  = "Control"
#' )
#' summary(fit_adj)
#'
#' @export
empirical_rate <- function(data,
                           treatment,
                           count,
                           exposure,
                           covariates = NULL,
                           car_scheme = c("simple", "permuted-block",
                                          "pocock-simon", "biased-coin", "urn"),
                           strata = NULL,
                           reference = NULL,
                           alpha = 0.05) {

  call <- match.call()
  car_scheme <- match.arg(car_scheme)

  # ---- input validation ------------------------------------------------------
  if (!is.data.frame(data)) {
    stop("`data` must be a data frame.", call. = FALSE)
  }
  chk_col <- function(x, arg) {
    if (!is.character(x) || length(x) != 1L) {
      stop("`", arg, "` must be a single column name (character scalar).",
           call. = FALSE)
    }
    if (!x %in% names(data)) {
      stop("Column \"", x, "\" (given as `", arg, "`) was not found in `data`.",
           call. = FALSE)
    }
  }
  chk_col(treatment, "treatment")
  chk_col(count, "count")
  chk_col(exposure, "exposure")

  chk_cols <- function(x, arg) {
    if (is.null(x)) return(invisible())
    if (!is.character(x)) {
      stop("`", arg, "` must be a character vector of column names.",
           call. = FALSE)
    }
    miss <- setdiff(x, names(data))
    if (length(miss)) {
      stop("Column(s) named in `", arg, "` not found in `data`: ",
           paste(miss, collapse = ", "), ".", call. = FALSE)
    }
  }
  chk_cols(covariates, "covariates")
  chk_cols(strata, "strata")

  if (!is.numeric(alpha) || length(alpha) != 1L || alpha <= 0 || alpha >= 1) {
    stop("`alpha` must be a single number strictly between 0 and 1.",
         call. = FALSE)
  }

  y <- data[[count]]
  e <- data[[exposure]]
  if (!is.numeric(y) || any(is.na(y)) || any(y < 0)) {
    stop("`count` column must be non-negative and free of missing values.",
         call. = FALSE)
  }
  if (!is.numeric(e) || any(is.na(e)) || any(e <= 0)) {
    stop("`exposure` column must be strictly positive and free of missing ",
         "values.", call. = FALSE)
  }

  trt_fac <- as.factor(data[[treatment]])
  if (nlevels(trt_fac) < 2L) {
    stop("`treatment` must have at least two levels; found ", nlevels(trt_fac),
         ".", call. = FALSE)
  }
  if (!is.null(reference)) {
    if (length(reference) != 1L || !reference %in% levels(trt_fac)) {
      stop("`reference` must be one of the treatment levels: ",
           paste(levels(trt_fac), collapse = ", "), ".", call. = FALSE)
    }
  }

  if (car_scheme != "simple" && is.null(strata)) {
    stop("`car_scheme = \"", car_scheme, "\"` requires `strata` (the ",
         "randomization-stratification columns). Use `car_scheme = \"simple\"` ",
         "for an analysis that ignores the randomization scheme.",
         call. = FALSE)
  }
  if (car_scheme == "simple" && !is.null(strata)) {
    warning("`strata` is ignored when `car_scheme = \"simple\"`. Set a ",
            "`car_scheme` such as \"permuted-block\" to account for the ",
            "stratified randomization scheme.", call. = FALSE)
  }

  # ---- transform outcome -----------------------------------------------------
  tr <- transform_counts(y = y, exposure = e, trt = trt_fac)

  # Internal column names for the working model. RobinCar rejects formula terms
  # that begin with a dot (a dot is special in R formulas), so we use plain
  # names and make them unique against the user's columns to avoid clobbering.
  safe_name <- function(base, taken) {
    nm <- base
    i <- 1L
    while (nm %in% taken) { nm <- paste0(base, "_", i); i <- i + 1L }
    nm
  }
  treat_nm <- safe_name("emp_treat", names(data))
  resp_nm  <- safe_name("emp_response", c(names(data), treat_nm))

  df_emp <- data
  df_emp[[treat_nm]] <- trt_fac
  df_emp[[resp_nm]]  <- tr$W

  # ---- working-model formula (common-slope ANCOVA) ---------------------------
  if (is.null(covariates) || length(covariates) == 0L) {
    form <- stats::reformulate(treat_nm, response = resp_nm)
  } else {
    form <- stats::reformulate(c(treat_nm, covariates), response = resp_nm)
  }

  # ---- fit via RobinCar ------------------------------------------------------
  if (!requireNamespace("RobinCar", quietly = TRUE)) {
    stop("Package 'RobinCar' is required but not installed. ",
         "Install it with install.packages(\"RobinCar\").", call. = FALSE)
  }
  fit_rc <- tryCatch(
    RobinCar::robincar_glm(
      df              = df_emp,
      treat_col       = treat_nm,
      response_col    = resp_nm,
      car_strata_cols = strata,
      car_scheme      = car_scheme,
      g_family        = stats::gaussian,
      formula         = form
    ),
    error = function(cnd) {
      stop("RobinCar::robincar_glm() failed: ", conditionMessage(cnd),
           call. = FALSE)
    }
  )

  # ---- extract marginal rates and their vcov ---------------------------------
  res <- fit_rc$result
  r_hat <- as.numeric(res$estimate)
  names(r_hat) <- as.character(res$treat)

  V_r <- as.matrix(fit_rc$varcov)
  dimnames(V_r) <- list(names(r_hat), names(r_hat))

  rates <- data.frame(
    treatment = names(r_hat),
    rate      = r_hat,
    se        = as.numeric(res$se),
    row.names = NULL,
    stringsAsFactors = FALSE
  )

  # ---- delta method + pairwise rate ratios -----------------------------------
  dl <- delta_lograte(r_hat, V_r)
  comparisons <- pairwise_rr(dl$mu_hat, dl$V_mu, alpha = alpha,
                             reference = reference)

  structure(
    list(
      rates        = rates,
      comparisons  = comparisons,
      vcov_rate    = V_r,
      log_rates    = dl$mu_hat,
      vcov_lograte = dl$V_mu,
      robincar     = fit_rc,
      call         = call,
      formula      = form,
      covariates   = covariates,
      car_scheme   = car_scheme,
      strata       = strata,
      reference    = reference,
      alpha        = alpha
    ),
    class = "empirical_rate"
  )
}
