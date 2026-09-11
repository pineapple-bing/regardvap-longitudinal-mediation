#!/usr/bin/env Rscript
# Complete verified baseline Table 1 for the REGARD-VAP longitudinal mediation study.
#
# Target population: every eligible participant with a recorded baseline exposure
# (baseline_carba_r) and 60-day mortality outcome in the current integrated
# analysis extract. The table is stratified only by the recorded binary exposure
# used by the analysis: A = 0 versus A = 1.
#
# Deliberate exclusions: 3GCR subgroup labels, WGS variables, plasmid variables,
# and any other unvalidated category are NOT reconstructed or shown. This script
# does not infer resistance class from organism name alone.
#
# Usage:
# Rscript 02_complete_baseline_table.R <day0_day3_integrated.xlsx> <output_directory>
#
# Outputs:
#   table1_complete_verified.csv      complete publication-ready summary
#   table1_patient_level_verified.csv exact participant-level source data
#   table1_cohort_audit.csv           cohort and group-count audit
#   table1_complete_verified.html     formatted table, if package gt is installed

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) {
  stop("Usage: Rscript 02_complete_baseline_table.R <day0_day3_integrated.xlsx> <output_directory>", call. = FALSE)
}
input_path <- normalizePath(args[[1]], mustWork = TRUE)
output_dir <- args[[2]]
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
if (!requireNamespace("readxl", quietly = TRUE)) stop("Package `readxl` is required.", call. = FALSE)

clean_num <- function(x) suppressWarnings(as.numeric(x))
clean_chr <- function(x) {
  out <- trimws(as.character(x))
  out[out %in% c("", "NA", "NaN")] <- NA_character_
  out
}
fmt_median_iqr <- function(x, digits = 1L) {
  x <- clean_num(x); x <- x[!is.na(x)]
  if (!length(x)) return("NA")
  q <- stats::quantile(x, c(.25, .5, .75), names = FALSE)
  sprintf(paste0("%.", digits, "f [%.", digits, "f, %.", digits, "f]"), q[2], q[1], q[3])
}
fmt_n_pct <- function(x) {
  x <- as.integer(x); denom <- sum(!is.na(x))
  if (!denom) return("NA")
  sprintf("%d (%.1f%%)", sum(x == 1L, na.rm = TRUE), 100 * sum(x == 1L, na.rm = TRUE) / denom)
}

# Absolute standardized difference between A=0 and A=1. It is a descriptive
# imbalance measure, not a p value. No SMD is manufactured if data are absent.
smd_continuous <- function(x, A) {
  x0 <- x[A == 0 & !is.na(x)]; x1 <- x[A == 1 & !is.na(x)]
  if (length(x0) < 2L || length(x1) < 2L) return(NA_character_)
  pooled_sd <- sqrt((stats::var(x0) + stats::var(x1)) / 2)
  if (!is.finite(pooled_sd) || pooled_sd == 0) return(NA_character_)
  sprintf("%.2f", abs(mean(x1) - mean(x0)) / pooled_sd)
}
smd_binary <- function(x, A) {
  x0 <- x[A == 0 & !is.na(x)]; x1 <- x[A == 1 & !is.na(x)]
  if (!length(x0) || !length(x1)) return(NA_character_)
  p0 <- mean(x0 == 1L); p1 <- mean(x1 == 1L)
  denom <- sqrt((p0 * (1 - p0) + p1 * (1 - p1)) / 2)
  if (!is.finite(denom) || denom == 0) return(NA_character_)
  sprintf("%.2f", abs(p1 - p0) / denom)
}

summary_continuous <- function(data, variable, label) {
  x <- data[[variable]]
  data.frame(
    section = NA_character_, variable = label, row_type = "continuous",
    overall = fmt_median_iqr(x),
    `Non-resistant (A=0)` = fmt_median_iqr(x[data$A == 0]),
    `Carbapenem-resistant (A=1)` = fmt_median_iqr(x[data$A == 1]),
    `Max pairwise SMD` = smd_continuous(x, data$A),
    check.names = FALSE, stringsAsFactors = FALSE
  )
}
summary_binary <- function(data, variable, label) {
  x <- as.integer(data[[variable]])
  data.frame(
    section = NA_character_, variable = label, row_type = "binary",
    overall = fmt_n_pct(x),
    `Non-resistant (A=0)` = fmt_n_pct(x[data$A == 0]),
    `Carbapenem-resistant (A=1)` = fmt_n_pct(x[data$A == 1]),
    `Max pairwise SMD` = smd_binary(x, data$A),
    check.names = FALSE, stringsAsFactors = FALSE
  )
}
summary_levels <- function(data, variable, label, levels_to_show) {
  x <- clean_chr(data[[variable]])
  do.call(rbind, lapply(levels_to_show, function(level) {
    data$indicator <- as.integer(!is.na(x) & x == level)
    summary_binary(data, "indicator", paste0(label, ": ", level))
  }))
}
section_row <- function(label) {
  data.frame(
    section = label, variable = label, row_type = "section",
    overall = NA_character_, `Non-resistant (A=0)` = NA_character_,
    `Carbapenem-resistant (A=1)` = NA_character_, `Max pairwise SMD` = NA_character_,
    check.names = FALSE, stringsAsFactors = FALSE
  )
}
organism_indicator <- function(x, pattern) as.integer(!is.na(x) & grepl(pattern, x, ignore.case = TRUE, perl = TRUE))

day0 <- readxl::read_excel(input_path)
needed <- c(
  "subjid", "day_from_day0", "baseline_carba_r", "outcome_death60_fixed",
  "baseline_age", "baseline_male", "baseline_charlson", "baseline_sofa",
  "baseline_map_min", "baseline_hr_max", "baseline_spo2fio2_min",
  "baseline_country_f", "baseline_icu_type", "baseline_intub_days",
  "baseline_intub_reason", "baseline_bacteria_f", "baseline_index_culture_organisms_all"
)
missing <- setdiff(needed, names(day0))
if (length(missing)) stop("Input is missing required columns: ", paste(missing, collapse = ", "), call. = FALSE)

base <- day0[day0$day_from_day0 == 0L, , drop = FALSE]
base <- base[!duplicated(base$subjid), , drop = FALSE]
base$A <- clean_num(base$baseline_carba_r)
base$death60 <- clean_num(base$outcome_death60_fixed)
base <- base[base$A %in% c(0, 1) & !is.na(base$death60), , drop = FALSE]
if (!nrow(base)) stop("No eligible participants after applying recorded A and 60-day outcome criteria.", call. = FALSE)

# Verified, explicitly named baseline/onset variables in the current extract.
base$age <- clean_num(base$baseline_age)
base$male <- clean_num(base$baseline_male)
base$charlson <- clean_num(base$baseline_charlson)
base$sofa <- clean_num(base$baseline_sofa)
base$map_min <- clean_num(base$baseline_map_min)
base$hr_max <- clean_num(base$baseline_hr_max)
base$spo2fio2_min <- clean_num(base$baseline_spo2fio2_min)
base$country <- clean_chr(base$baseline_country_f)
base$icu_type <- clean_chr(base$baseline_icu_type)
base$intub_days <- clean_num(base$baseline_intub_days)
base$intub_reason <- clean_chr(base$baseline_intub_reason)
base$bacteria_group <- clean_chr(base$baseline_bacteria_f)
base$index_organisms <- clean_chr(base$baseline_index_culture_organisms_all)
base$polymicrobial <- as.integer(!is.na(base$index_organisms) & grepl("\\|", base$index_organisms))
base$acinetobacter <- organism_indicator(base$index_organisms, "acito")
base$pseudomonas <- organism_indicator(base$index_organisms, "psudo")
base$enterobacterales <- organism_indicator(base$index_organisms, "ecoli|kleb|serat|prote|morg|entbc|provid")
base$staph_aureus <- organism_indicator(base$index_organisms, "staphy")

n0 <- sum(base$A == 0); n1 <- sum(base$A == 1)
n_row <- data.frame(
  section = "Header", variable = "Patients, n", row_type = "n",
  overall = as.character(nrow(base)), `Non-resistant (A=0)` = as.character(n0),
  `Carbapenem-resistant (A=1)` = as.character(n1), `Max pairwise SMD` = NA_character_,
  check.names = FALSE, stringsAsFactors = FALSE
)

table1 <- rbind(
  n_row,
  section_row("Demographics and comorbidity"),
  summary_continuous(base, "age", "Age, years, median [IQR]"),
  summary_binary(base, "male", "Male sex, n (%)"),
  summary_continuous(base, "charlson", "Charlson comorbidity index, median [IQR]"),
  section_row("Healthcare setting and pre-VAP course"),
  summary_levels(base, "country", "Country", c("Singapore", "Thailand", "Nepal")),
  summary_levels(base, "icu_type", "ICU type", c("Medical ICU", "Surgical ICU")),
  summary_continuous(base, "intub_days", "Days intubated before VAP onset, median [IQR]"),
  summary_levels(base, "intub_reason", "Indication for intubation", c(
    "Respiratory failure", "Neurological failure", "Trauma/airway obstruction",
    "Post-operative care", "Cardiovascular failure", "Metabolic acidosis"
  )),
  section_row("Clinical severity at VAP onset"),
  summary_continuous(base, "sofa", "SOFA score, median [IQR]"),
  summary_continuous(base, "map_min", "Lowest mean arterial pressure, mmHg, median [IQR]"),
  summary_continuous(base, "hr_max", "Maximum heart rate, beats/min, median [IQR]"),
  summary_continuous(base, "spo2fio2_min", "SpO2/FiO2 ratio, median [IQR]"),
  section_row("Index-culture characteristics"),
  summary_binary(base, "polymicrobial", "Polymicrobial index culture, n (%)"),
  summary_binary(base, "acinetobacter", "Any Acinetobacter spp., n (%)"),
  summary_binary(base, "pseudomonas", "Any Pseudomonas spp., n (%)"),
  summary_binary(base, "enterobacterales", "Any Enterobacterales, n (%)"),
  summary_binary(base, "staph_aureus", "Any Staphylococcus aureus, n (%)"),
  section_row("Outcome"),
  summary_binary(base, "death60", "60-day mortality, n (%)")
)

# Culture-negative is directly derived from the verified baseline bacteria field.
x_culture_negative <- as.integer(base$bacteria_group == "Culture_neg")
culture_negative_row <- data.frame(
  section = NA_character_, variable = "Culture-negative, n (%)", row_type = "binary",
  overall = fmt_n_pct(x_culture_negative),
  `Non-resistant (A=0)` = fmt_n_pct(x_culture_negative[base$A == 0]),
  `Carbapenem-resistant (A=1)` = fmt_n_pct(x_culture_negative[base$A == 1]),
  `Max pairwise SMD` = smd_binary(x_culture_negative, base$A),
  check.names = FALSE, stringsAsFactors = FALSE
)
table1 <- rbind(table1[seq_len(nrow(table1) - 2L), ], culture_negative_row, table1[(nrow(table1) - 1L):nrow(table1), ])

cohort_audit <- data.frame(
  stage = c("Day-0 person records", "Eligible: recorded A and 60-day outcome", "A=0", "A=1"),
  n = c(sum(day0$day_from_day0 == 0L), nrow(base), n0, n1),
  stringsAsFactors = FALSE
)
utils::write.csv(table1, file.path(output_dir, "table1_complete_verified.csv"), row.names = FALSE, na = "")
utils::write.csv(base, file.path(output_dir, "table1_patient_level_verified.csv"), row.names = FALSE)
utils::write.csv(cohort_audit, file.path(output_dir, "table1_cohort_audit.csv"), row.names = FALSE)

if (requireNamespace("gt", quietly = TRUE)) {
  display <- table1[, c("variable", "overall", "Non-resistant (A=0)", "Carbapenem-resistant (A=1)", "Max pairwise SMD")]
  names(display)[1] <- "Characteristic"
  section_rows <- which(table1$row_type == "section")
  gt_tab <- gt::gt(display) |>
    gt::tab_header(
      title = "Table 1. Baseline and onset characteristics by recorded carbapenem-resistance status",
      subtitle = "All eligible REGARD-VAP participants with recorded A and 60-day outcome. Values are median [IQR] or n (%)."
    ) |>
    gt::tab_style(style = gt::cell_fill(color = "#0D7182"), locations = gt::cells_body(rows = section_rows)) |>
    gt::tab_style(style = gt::cell_text(color = "white", weight = "bold"), locations = gt::cells_body(rows = section_rows)) |>
    gt::tab_style(style = gt::cell_fill(color = "#F9EEE9"), locations = gt::cells_body(rows = which(table1$variable == "60-day mortality, n (%)"))) |>
    gt::tab_source_note("SMD is the absolute standardized difference for A=0 versus A=1; it is descriptive and is not a p value.") |>
    gt::tab_options(table.font.size = gt::px(12), data_row.padding = gt::px(3))
  gt::gtsave(gt_tab, file.path(output_dir, "table1_complete_verified.html"))
}

cat("Verified complete baseline Table 1 written to: ", normalizePath(output_dir), "\n", sep = "")
