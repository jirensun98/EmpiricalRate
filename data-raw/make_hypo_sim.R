# Generate `hypo_sim`, a SIMULATED example dataset for the EmpiricalRate
# package. It mimics the *structure* of a sparse severe-hypoglycemia endpoint
# (two arms, heterogeneous exposure, a baseline HbA1c covariate and a baseline
# event-rate covariate) but contains no real patient data. The QWINT-5 trial
# data analysed in the companion paper are confidential and are NOT included.
#
# Run from the package root with:  Rscript data-raw/make_hypo_sim.R

set.seed(2026)

n_per_arm <- 250L
arms      <- c("Control", "Treatment")

make_arm <- function(arm, n, true_rate) {
  # Heterogeneous exposure: mixture of two uniforms, mean ~ 1 year.
  use_wide <- stats::rbinom(n, 1L, 0.5)
  exposure <- ifelse(use_wide,
                     stats::runif(n, 0.6, 1.2),
                     stats::runif(n, 0.8, 1.4))

  # Baseline HbA1c (%), mild spread around 8.
  baseline_hba1c <- round(stats::rnorm(n, mean = 8.0, sd = 0.9), 1)

  # Baseline severe-hypo rate (events/year in lead-in): sparse and skewed;
  # most participants have none, a minority have a positive baseline rate.
  had_baseline   <- stats::rbinom(n, 1L, 0.30)
  baseline_rate  <- round(had_baseline * stats::rgamma(n, shape = 1.2,
                                                       rate = 4), 3)

  # Postbaseline events: NB mean driven by the true arm rate, the exposure,
  # and a modest positive association with the baseline rate (so covariate
  # adjustment is meaningful but the relationship is not exactly log-linear).
  mu <- true_rate * exposure * exp(0.8 * baseline_rate)
  events <- stats::rnbinom(n, mu = mu, size = 0.7)

  data.frame(
    arm            = arm,
    events         = as.integer(events),
    exposure_years = round(exposure, 3),
    baseline_hba1c = baseline_hba1c,
    baseline_rate  = baseline_rate,
    stringsAsFactors = FALSE
  )
}

hypo_sim <- rbind(
  make_arm("Control",   n_per_arm, true_rate = 0.25),
  make_arm("Treatment", n_per_arm, true_rate = 0.15)
)

# Randomize row order and make `arm` a factor with Control as the first level.
hypo_sim <- hypo_sim[sample(nrow(hypo_sim)), ]
rownames(hypo_sim) <- NULL
hypo_sim$arm <- factor(hypo_sim$arm, levels = c("Control", "Treatment"))

# Quick sanity report (printed when the script is run interactively).
agg <- tapply(hypo_sim$events, hypo_sim$arm, sum) /
       tapply(hypo_sim$exposure_years, hypo_sim$arm, sum)
message("Aggregate rates by arm:")
print(round(agg, 4))
message("Observed rate ratio (Treatment / Control): ",
        round(agg[["Treatment"]] / agg[["Control"]], 4))
message("Total events by arm:")
print(tapply(hypo_sim$events, hypo_sim$arm, sum))

usethis::use_data(hypo_sim, overwrite = TRUE)
