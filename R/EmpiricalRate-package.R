#' EmpiricalRate: Empirical Estimation of Event Rates and Rate Ratios
#'
#' A distribution-free method for comparing marginal event rates between
#' treatment groups when the endpoint is a count (recurrent event) such as
#' hypoglycemia. The approach transforms event counts to subject-level rate
#' quantities, fits a linear working model with robust variance estimation via
#' \pkg{RobinCar}, and back-transforms to obtain marginal rates and rate ratios
#' with confidence intervals via the delta method.
#'
#' The central object is [empirical_rate()]. For a worked introduction see the
#' package vignette with `vignette("EmpiricalRate", package = "EmpiricalRate")`
#' or `browseVignettes("EmpiricalRate")`. When installing from GitHub, build it
#' with `remotes::install_github("jirensun98/EmpiricalRate", build_vignettes =
#' TRUE)`.
#'
#' @section Method in brief:
#' For subject \eqn{j} in arm \eqn{i} with event count \eqn{Y_{ij}} and exposure
#' \eqn{d_{ij}}, define \eqn{W_{ij} = Y_{ij} / \bar d_{i\cdot}} where
#' \eqn{\bar d_{i\cdot}} is the mean exposure in arm \eqn{i}. The arm-specific
#' mean of \eqn{W_{ij}} is exactly unbiased for the marginal rate \eqn{r_i}
#' without any distributional assumption on \eqn{Y_{ij}}. A linear model is
#' fitted to \eqn{W} (optionally adjusting for baseline covariates and/or the
#' randomization scheme), and the rate ratio \eqn{\lambda_{ik} = r_i / r_k} is
#' obtained on the log scale.
#'
#' @keywords internal
#' @references
#' Sun J, Amoafo L, Qu Y (2026). An Empirical Method for Analyzing Count Data.
#'
#' Ye T, Shao J, Yi Y, Zhao Q (2023). Toward better practice of covariate
#' adjustment in analyzing randomized clinical trials. \emph{Journal of the
#' American Statistical Association}, 118(544), 2370-2382.
#'
#' Bannick MS, Shao J, Liu J, Du Y, Yi Y, Ye T (2025). A general form of
#' covariate adjustment in clinical trials under covariate-adaptive
#' randomization. \emph{Biometrika}, asaf029.
"_PACKAGE"
