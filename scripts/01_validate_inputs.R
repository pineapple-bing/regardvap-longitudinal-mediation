args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L) stop("Usage: Rscript scripts/01_validate_inputs.R config.R")
source(args[[1]])

suppressPackageStartupMessages(library(readxl))

required <- c("day03_path", "day60_path", "severity_path", "exposure_estimand")
missing_cfg <- required[!vapply(required, exists, logical(1), inherits = TRUE)]
if (length(missing_cfg)) stop("Missing config values: ", paste(missing_cfg, collapse = ", "))
if (identical(exposure_estimand, "CHOOSE_BEFORE_ANALYSIS")) {
  stop("Choose exposure_estimand in config.R before analysis.")
}

for (p in c(day03_path, day60_path, severity_path)) {
  if (!file.exists(p)) stop("Missing input: ", p)
}

d <- as.data.frame(read_excel(day03_path))
s <- as.data.frame(read_excel(severity_path))

needed <- c(
  "subjid", "day_from_day0", "baseline_carba_r", "baseline_bacteria_f",
  "baseline_appropriateabx_on_symptomdate",
  "treatment_any_appropriate_active_primary", "alive_at_start_of_day",
  "under_followup_on_day"
)
missing_cols <- setdiff(needed, names(d))
if (length(missing_cols)) stop("Missing columns: ", paste(missing_cols, collapse = ", "))

day0 <- d[d$day_from_day0 == 0, , drop = FALSE]
cat("Participants at Day 0:", length(unique(day0$subjid)), "\n")
cat("Exposure counts:\n")
print(table(day0$baseline_carba_r, useNA = "ifany"))
cat("\nOrganism group by exposure:\n")
print(with(day0, table(baseline_bacteria_f, baseline_carba_r, useNA = "ifany")))

culture_negative <- grepl("culture[_ ]?neg", day0$baseline_bacteria_f, ignore.case = TRUE)
n_cn_a0 <- sum(culture_negative & day0$baseline_carba_r == 0, na.rm = TRUE)
if (n_cn_a0 > 0) {
  cat("\nFAIL/DECISION: A=0 contains", n_cn_a0,
      "culture-negative episodes. It cannot be labelled susceptible without a verified eligibility rule.\n")
}

m0 <- suppressWarnings(as.numeric(day0$treatment_any_appropriate_active_primary))
b0 <- suppressWarnings(as.numeric(day0$baseline_appropriateabx_on_symptomdate))
both <- !is.na(m0) & !is.na(b0)
agreement <- if (any(both)) mean(m0[both] == b0[both]) else NA_real_
cat("\nDay-0 treatment agreement with baseline appropriate-treatment variable:",
    sprintf("%.1f%%", 100 * agreement), "\n")
if (is.finite(agreement) && agreement > 0.95) {
  cat("FAIL/DECISION: baseline appropriate treatment duplicates M0 and should not be used as a baseline confounder.\n")
}

days <- 0:3
qc <- do.call(rbind, lapply(days, function(day) {
  x <- d[d$day_from_day0 == day, , drop = FALSE]
  data.frame(
    day = day,
    n = nrow(x),
    alive = sum(x$alive_at_start_of_day == 1, na.rm = TRUE),
    followed = sum(x$under_followup_on_day == 1, na.rm = TRUE),
    treatment_observed = sum(!is.na(x$treatment_any_appropriate_active_primary))
  )
}))
cat("\nAt-risk and treatment observation checks:\n")
print(qc, row.names = FALSE)

if (all(c("subjid", "day_from_day0", "severity_lt_complete") %in% names(s))) {
  sev <- do.call(rbind, lapply(days, function(day) {
    x <- s[s$day_from_day0 == day, , drop = FALSE]
    data.frame(day = day, n = nrow(x), severity_complete = sum(x$severity_lt_complete == 1, na.rm = TRUE))
  }))
  cat("\nSeverity completeness:\n")
  print(sev, row.names = FALSE)
}

sim_cols <- names(d)[grepl("_sim$", names(d))]
if (length(sim_cols)) {
  cat("\nWARN: simulated columns detected and excluded from the planned clinical pipeline:\n")
  cat(paste0("- ", sim_cols, collapse = "\n"), "\n")
}

cat("\nValidation finished. Resolve all FAIL/DECISION messages before fitting causal models.\n")
