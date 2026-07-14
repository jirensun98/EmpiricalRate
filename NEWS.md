# EmpiricalRate 0.1.0

* Initial release.
* `empirical_rate()`: empirical (distribution-free) estimation of marginal
  event rates and pairwise rate ratios for count endpoints, with optional
  baseline-covariate adjustment (ANCOVA) and randomization-scheme-aware
  variance estimation (ANHECOVA) via RobinCar.
* `meta_rate_ratio()`: pool rate ratios across strata or studies with exposure
  weights, on the log scale (default) or the natural scale.
* `print()` and `summary()` methods for `empirical_rate` objects.
* `hypo_sim`: a simulated sparse-event example dataset.
