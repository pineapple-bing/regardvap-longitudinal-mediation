args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L) stop("Usage: Rscript 01_data_preparation/01_validate_inputs.R config.R")
source(args[[1]])
source("R/functions_data.R")
suppressPackageStartupMessages(library(readxl))

required_cfg <- c("day03_path", "day60_path", "severity_path", "exposure_estimand")
missing_cfg <- required_cfg[!vapply(required_cfg, exists, logical(1), inherits = TRUE)]
if (length(missing_cfg)) stop("Missing config values: ", paste(missing_cfg, collapse = ", "))
if (identical(exposure_estimand, "CHOOSE_BEFORE_ANALYSIS")) stop("Choose exposure_estimand in config.R before analysis.")
for (path in c(day03_path, day60_path, severity_path)) if (!file.exists(path)) stop("Missing input: ", path)

day_data <- as.data.frame(read_excel(day03_path))
stop_if_missing(day_data, c(id_var, day_var, exposure_var, mediator_var,
  "baseline_bacteria_f", "alive_at_start_of_day", "under_followup_on_day"), "day03 file")
day0 <- day_data[day_data[[day_var]] == 0, , drop = FALSE]
cat("Participants at Day 0:", length(unique(day0[[id_var]])), "\n")
cat("Exposure counts:\n"); print(table(day0[[exposure_var]], useNA = "ifany"))

culture_negative <- grepl("culture[_ ]?neg", day0$baseline_bacteria_f, ignore.case = TRUE)
n_culture_negative_a0 <- sum(culture_negative & day0[[exposure_var]] == 0, na.rm = TRUE)
if (n_culture_negative_a0) cat("FAIL/DECISION: A=0 includes", n_culture_negative_a0,
  "culture-negative episodes; do not label A=0 as susceptible.\n")

if ("baseline_appropriateabx_on_symptomdate" %in% names(day0)) {
  m0 <- as_binary(day0[[mediator_var]], mediator_var)
  baseline_tx <- as_binary(day0$baseline_appropriateabx_on_symptomdate, "baseline_appropriateabx_on_symptomdate")
  observed <- !is.na(m0) & !is.na(baseline_tx)
  agreement <- mean(m0[observed] == baseline_tx[observed])
  cat("M0 agreement with baseline appropriate-treatment field:", sprintf("%.1f%%", 100 * agreement), "\n")
  if (agreement > .95) cat("FAIL/DECISION: baseline appropriate-treatment duplicates M0; exclude it from baseline covariates.\n")
}

for (day in 0:3) {
  x <- day_data[day_data[[day_var]] == day, , drop = FALSE]
  cat("Day", day, ": n=", nrow(x), ", alive=", sum(x$alive_at_start_of_day == 1, na.rm = TRUE),
    ", followed=", sum(x$under_followup_on_day == 1, na.rm = TRUE),
    ", mediator observed=", sum(!is.na(x[[mediator_var]])), "\n", sep = "")
}
simulated_columns <- names(day_data)[grepl("_sim$|simulated", names(day_data), ignore.case = TRUE)]
if (length(simulated_columns)) cat("WARN: simulated fields detected and excluded from clinical analysis:\n- ",
  paste(simulated_columns, collapse = "\n- "), "\n", sep = "")
cat("Validation completed.\n")
