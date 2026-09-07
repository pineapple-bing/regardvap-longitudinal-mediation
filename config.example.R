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
