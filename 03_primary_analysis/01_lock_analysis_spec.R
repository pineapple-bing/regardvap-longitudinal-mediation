args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L) stop("Usage: Rscript 03_primary_analysis/01_lock_analysis_spec.R config.R")
source(args[[1]]); source("R/functions_data.R")
if (!exposure_estimand %in% c("confirmed_resistant_vs_confirmed_susceptible", "resistance_detected_vs_not_detected")) {
  stop("Set exposure_estimand to one of the two prespecified values in config.R.")
}
specification <- data.frame(
  component = c("Population", "Exposure", "Mediator", "Outcome", "Effect scale", "Estimator", "Uncertainty"),
  definition = c(
    if (exposure_estimand == "confirmed_resistant_vs_confirmed_susceptible") "Eligible episodes with interpretable carbapenem susceptibility" else "Full eligible cohort; use detection language",
    exposure_estimand, "Appropriate active treatment M0, M1, M2, M3", "Death by Day 60", "Risk difference",
    "Longitudinal g-computation for interventional direct and indirect effects",
    "Patient-level nonparametric bootstrap; repeat entire fit and simulation"
  )
)
write_csv_safely(specification, "outputs/tables/03_locked_primary_specification.csv")
message("Question locked. Fit the g-formula only after timestamps and covariates are signed off.")
