# Copy this file to config.R and edit the local paths. config.R is git-ignored.

day03_path <- "data/regardvap_patient_day0_day3.xlsx"
day60_path <- "data/regardvap_patient_day0_day60.xlsx"
severity_path <- "data/regardvap_primary_severity_L0_L3.xlsx"
output_dir <- "outputs"

# This must be chosen with the clinical/statistical team before inference.
# "detected_vs_not_detected" retains culture-negative episodes in A=0.
# "resistant_vs_susceptible" requires a verified susceptible index isolate and
# excludes culture-negative or non-evaluable episodes.
exposure_estimand <- "CHOOSE_BEFORE_ANALYSIS"

# Column names in your authorised source files. Change only here when a source
# extract changes, rather than editing multiple analysis scripts.
id_var <- "subjid"
day_var <- "day_from_day0"
exposure_var <- "baseline_carba_r"
mediator_var <- "treatment_any_appropriate_active_primary"
severity_var <- "severity_lt"
outcome_var <- "death_by_day60"

# Do NOT add baseline_appropriateabx_on_symptomdate: it is essentially M0 in
# the current data and is therefore a post-exposure treatment variable.
baseline_covariates <- c("age", "sex", "country", "site", "baseline_bacteria_f")
analysis_data_path <- "outputs/analysis_dataset.rds"
