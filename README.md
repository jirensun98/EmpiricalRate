# EmpiricalRate

<!-- badges: start -->
[![License: MIT](https://img.shields.io/badge/license-MIT-blue.svg)](LICENSE.md)
<!-- badges: end -->

**EmpiricalRate** implements an empirical, distribution-free method for
comparing marginal event rates between treatment groups when the endpoint is a
count (recurrent event), such as hypoglycemia episodes in a clinical trial.

The method targets the same marginal estimand as negative binomial (NB)
regression — the ratio of marginal event rates — but does not fit a parametric
count model. This avoids the numerical instability that likelihood-based
estimation can suffer when events are sparse, and it does not require the
relationship between baseline covariates and the outcome to be correctly
specified.

The package accompanies the paper:

> Sun J, Amoafo L, Qu Y (2026). *An Empirical Method for Analyzing Count Data.*

## How it works

For subject *j* in arm *i* with event count *Y*<sub>ij</sub> and exposure
*d*<sub>ij</sub>, the count is transformed to a subject-level rate quantity

*W*<sub>ij</sub> = *Y*<sub>ij</sub> / *d̄*<sub>i·</sub>,

where *d̄*<sub>i·</sub> is the **mean exposure in arm _i_**. Because
E[*Y*<sub>ij</sub>] = *d*<sub>ij</sub> *r*<sub>i</sub>, the arm-specific mean of
*W*<sub>ij</sub> is exactly unbiased for the marginal rate *r*<sub>i</sub> with
no distributional assumption on the counts. A linear working model
(`W ~ treatment (+ covariates)`) is then fitted with robust variance estimation
via the [RobinCar](https://cran.r-project.org/package=RobinCar) package, and
rate ratios with confidence intervals are obtained on the log scale by the
delta method.

A useful property: in the **unadjusted** analysis the estimated marginal rate
for each arm equals the observed aggregate rate
(&sum;*Y* / &sum;*d*) exactly, and the estimated rate ratio equals the raw
rate ratio.

## Installation

Install the development version from GitHub. Pass `build_vignettes = TRUE` so
the introductory vignette is installed alongside the package (without it,
`install_github()` skips vignettes and `vignette("EmpiricalRate")` will find
nothing):

```r
# install.packages("remotes")
remotes::install_github("jirensun98/EmpiricalRate", build_vignettes = TRUE)
```

RobinCar is the analysis engine and is installed automatically as a dependency.

## Documentation

A step-by-step introduction — covering the unadjusted, covariate-adjusted, and
stratified (ANHECOVA) analyses, and the meta-analysis helper — is provided as a
package vignette. After installing with `build_vignettes = TRUE`, open it from R
with:

```r
vignette("EmpiricalRate", package = "EmpiricalRate")   # open the vignette
browseVignettes("EmpiricalRate")                        # list vignettes in a browser
```

Prefer to read it without installing? The source renders on GitHub here:
[`vignettes/EmpiricalRate.Rmd`](vignettes/EmpiricalRate.Rmd). Function-level help
is available with `?empirical_rate`, `?meta_rate_ratio`, and
`?EmpiricalRate` (the package overview).

## Quick start

```r
library(EmpiricalRate)

data(hypo_sim)  # simulated example data shipped with the package

# Unadjusted analysis (rates equal the observed aggregate rates)
fit <- empirical_rate(
  hypo_sim,
  treatment = "arm",
  count     = "events",
  exposure  = "exposure_years",
  reference = "Control"
)
fit

# Covariate-adjusted analysis (common-slope ANCOVA working model)
fit_adj <- empirical_rate(
  hypo_sim,
  treatment  = "arm",
  count      = "events",
  exposure   = "exposure_years",
  covariates = c("baseline_hba1c", "baseline_rate"),
  reference  = "Control"
)
summary(fit_adj)
```

To account for a **stratified randomization scheme** (the ANHECOVA analysis in
the paper), pass `car_scheme` and the stratification columns via `strata`:

```r
fit_str <- empirical_rate(
  hypo_sim,
  treatment  = "arm",
  count      = "events",
  exposure   = "exposure_years",
  covariates = c("baseline_hba1c", "baseline_rate"),
  car_scheme = "permuted-block",
  strata     = "site",          # replace with your stratification factor(s)
  reference  = "Control"
)
```

## Pooling across strata or studies

For stratified analyses or meta-analyses in which the rate ratio is assumed
common across strata but event rates differ, `meta_rate_ratio()` pools the
per-stratum rate ratios with exposure weights (Appendix A of the paper). It
takes the `log_rr` and `se_log_rr` columns returned by `empirical_rate()`:

```r
meta_rate_ratio(
  log_rr    = c(0.41, 0.55, 0.29),
  se_log_rr = c(0.18, 0.26, 0.15),
  weight    = c(320, 150, 410)   # total exposure per stratum
)
```

Pooling defaults to the log scale (as recommended in the paper's Discussion);
`scale = "identity"` reproduces the natural-scale formula in Appendix A.

## Note on the example data

`hypo_sim` is **simulated** and is provided only so the examples run. It is not
the QWINT-5 trial data analysed in the paper; those data are confidential and
are not distributed with the package.

## License

MIT © Eli Lilly and Company. See [LICENSE.md](LICENSE.md).
