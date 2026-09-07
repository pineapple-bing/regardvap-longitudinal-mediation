args <- commandArgs(trailingOnly = TRUE)

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (length(script_arg) != 1L) stop("Run this file with Rscript.", call. = FALSE)
module_dir <- dirname(normalizePath(sub("^--file=", "", script_arg), mustWork = TRUE))
project_root <- normalizePath(file.path(module_dir, ".."), mustWork = TRUE)
source(file.path(project_root, "01_study_population_and_baseline", "R", "01_cohort_diagnostics.R"))
source(file.path(project_root, "02_longitudinal_trajectories", "R", "01_observed_trajectory_diagnostics.R"))
source(file.path(module_dir, "R", "01_analysis_contract.R"))
source(file.path(module_dir, "R", "02_diagnostics.R"))
source(file.path(module_dir, "R", "03_bootstrap.R"))

suppressPackageStartupMessages({
  library(readxl)
})

default_day03_path <- "/Users/bing/Desktop/中介分析 最新数据（最新sofa定义/副本regardvap_patient_day0_day3_integrated_longitudinal_itt460_v.xlsx"
default_day60_path <- "/Users/bing/Desktop/中介分析 最新数据（最新sofa定义/副本regardvap_patient_day0_day60_integrated_longitudinal_itt460_.xlsx"
default_severity_path <- "/Users/bing/Desktop/中介分析 最新数据（最新sofa定义/副本regardvap_primary_severity_L0_L3_itt460.xlsx"
default_itt_rds_path <- "/Users/bing/Desktop/中介分析 最新数据（最新sofa定义/itt.dataset.RDS"
default_out_dir <- "/Users/bing/Documents/中介效应/outputs/regardvap_gformula_latest_v3"

day03_path <- if (length(args) >= 1) args[1] else default_day03_path
day60_path <- if (length(args) >= 2) args[2] else default_day60_path
severity_path <- if (length(args) >= 3) args[3] else default_severity_path
itt_rds_path <- if (length(args) >= 4) args[4] else default_itt_rds_path
out_dir <- if (length(args) >= 5) args[5] else default_out_dir

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

step_dirs <- c(
  step00 = file.path(out_dir, "step00_data_prep"),
  step01 = file.path(out_dir, "step01_baseline"),
  step02 = file.path(out_dir, "step02_longitudinal_summary"),
  step03 = file.path(out_dir, "step03_main_gformula"),
  step04 = file.path(out_dir, "step04_sensitivity"),
  step05 = file.path(out_dir, "step05_hte")
)
invisible(lapply(step_dirs, dir.create, recursive = TRUE, showWarnings = FALSE))

clean_num <- function(x) {
  out <- suppressWarnings(as.numeric(x))
  out[is.nan(out)] <- NA_real_
  out
}

clean_binary01 <- function(x) {
  out <- suppressWarnings(as.numeric(x))
  out[!(out %in% c(0, 1))] <- NA_real_
  out
}

safe_char <- function(x) {
  out <- as.character(x)
  out[out %in% c("", "NA", "NaN")] <- NA_character_
  out
}

fill_missing_factor <- function(x, missing_label = "Missing") {
  out <- safe_char(x)
  out[is.na(out)] <- missing_label
  factor(out)
}

fmt_mean_sd <- function(x, digits = 1) {
  x <- x[!is.na(x)]
  if (length(x) == 0) return("NA")
  sprintf(paste0("%.", digits, "f (%.", digits, "f)"), mean(x), stats::sd(x))
}

fmt_median_iqr <- function(x, digits = 1) {
  x <- x[!is.na(x)]
  if (length(x) == 0) return("NA")
  qs <- stats::quantile(x, c(0.25, 0.5, 0.75), names = FALSE)
  sprintf(
    paste0("%.", digits, "f [%.", digits, "f, %.", digits, "f]"),
    qs[2], qs[1], qs[3]
  )
}

fmt_n_pct <- function(n, denom) {
  if (is.na(denom) || denom == 0) return("0 (0.0%)")
  sprintf("%d (%.1f%%)", n, 100 * n / denom)
}

rhs_join <- function(terms) {
  terms <- unique(terms[!is.na(terms) & nzchar(terms)])
  if (length(terms) == 0) {
    "1"
  } else {
    paste(terms, collapse = " + ")
  }
}

is_informative_variable <- function(x) {
  x <- x[!is.na(x)]
  if (length(x) == 0) return(FALSE)
  if (is.factor(x)) x <- droplevels(x)
  length(unique(as.character(x))) >= 2
}

drop_invariant_columns <- function(df, vars) {
  keep <- vars[vapply(vars, function(v) is_informative_variable(df[[v]]), logical(1))]
  drop <- setdiff(vars, keep)
  if (length(keep) > 0) {
    for (v in keep) {
      if (is.factor(df[[v]])) {
        df[[v]] <- droplevels(df[[v]])
      }
    }
  }
  list(data = df, keep = keep, drop = drop)
}

score_temp_component <- function(x) {
  out <- rep(NA_real_, length(x))
  out[!is.na(x) & x < 36] <- 1
  out[!is.na(x) & x >= 36 & x < 39] <- 0
  out[!is.na(x) & x >= 39] <- 1
  out
}

score_map_component <- function(x) {
  out <- rep(NA_real_, length(x))
  out[!is.na(x) & x >= 70] <- 0
  out[!is.na(x) & x < 70 & x >= 60] <- 1
  out[!is.na(x) & x < 60] <- 2
  out
}

score_hr_component <- function(x) {
  out <- rep(NA_real_, length(x))
  out[!is.na(x) & x < 110] <- 0
  out[!is.na(x) & x >= 110 & x < 130] <- 1
  out[!is.na(x) & x >= 130] <- 2
  out
}

score_spo2fio2_component <- function(x) {
  out <- rep(NA_real_, length(x))
  out[!is.na(x) & x >= 315] <- 0
  out[!is.na(x) & x < 315 & x >= 235] <- 1
  out[!is.na(x) & x < 235 & x >= 150] <- 2
  out[!is.na(x) & x < 150] <- 3
  out
}

score_lt_from_components <- function(temp, map, hr, inotrope, mech_vent, spo2fio2) {
  temp_score <- score_temp_component(temp)
  map_score <- score_map_component(map)
  hr_score <- score_hr_component(hr)
  spo2fio2_score <- score_spo2fio2_component(spo2fio2)
  inotrope_score <- ifelse(is.na(inotrope), NA_real_, 2 * inotrope)
  mech_vent_score <- ifelse(is.na(mech_vent), NA_real_, 1 * mech_vent)
  parts <- cbind(temp_score, map_score, hr_score, spo2fio2_score, inotrope_score, mech_vent_score)
  out <- rowSums(parts, na.rm = TRUE)
  out[rowSums(!is.na(parts)) == 0] <- NA_real_
  out
}

fit_binomial_model <- function(formula, data) {
  glm(formula, data = data, family = stats::binomial())
}

fit_gaussian_model <- function(formula, data) {
  stats::lm(formula, data = data)
}

predict_binomial_prob <- function(model, newdata) {
  p <- suppressWarnings(stats::predict(model, newdata = newdata, type = "response"))
  p <- as.numeric(p)
  p[is.na(p)] <- mean(model$y, na.rm = TRUE)
  pmin(pmax(p, 0.001), 0.999)
}

predict_gaussian_draw <- function(model, newdata) {
  mu <- suppressWarnings(stats::predict(model, newdata = newdata))
  mu <- as.numeric(mu)
  mu[is.na(mu)] <- mean(model$model[[1]], na.rm = TRUE)
  sigma <- summary(model)$sigma
  if (is.na(sigma) || sigma <= 0) sigma <- 0.25
  stats::rnorm(nrow(newdata), mean = mu, sd = sigma)
}

estimate_risk <- function(sim_df) {
  mean(sim_df$Y, na.rm = TRUE)
}

summarize_continuous <- function(df, var, label, stat = c("mean_sd", "median_iqr")) {
  stat <- match.arg(stat)
  g0 <- df[df$A == 0, var]
  g1 <- df[df$A == 1, var]
  overall <- df[[var]]
  out0 <- if (stat == "mean_sd") fmt_mean_sd(g0) else fmt_median_iqr(g0)
  out1 <- if (stat == "mean_sd") fmt_mean_sd(g1) else fmt_median_iqr(g1)
  out_all <- if (stat == "mean_sd") fmt_mean_sd(overall) else fmt_median_iqr(overall)
  data.frame(
    section = "Continuous",
    variable = label,
    level = "",
    overall = out_all,
    resistant_0 = out0,
    resistant_1 = out1,
    stringsAsFactors = FALSE
  )
}

summarize_binary <- function(df, var, label, positive_value = 1) {
  d <- df[[var]]
  denom_all <- sum(!is.na(d))
  denom0 <- sum(!is.na(d[df$A == 0]))
  denom1 <- sum(!is.na(d[df$A == 1]))
  n_all <- sum(d == positive_value, na.rm = TRUE)
  n0 <- sum(d[df$A == 0] == positive_value, na.rm = TRUE)
  n1 <- sum(d[df$A == 1] == positive_value, na.rm = TRUE)
  data.frame(
    section = "Binary",
    variable = label,
    level = "Yes",
    overall = fmt_n_pct(n_all, denom_all),
    resistant_0 = fmt_n_pct(n0, denom0),
    resistant_1 = fmt_n_pct(n1, denom1),
    stringsAsFactors = FALSE
  )
}

summarize_categorical <- function(df, var, label) {
  x <- safe_char(df[[var]])
  levs <- sort(unique(x[!is.na(x)]))
  out <- vector("list", length(levs))
  for (i in seq_along(levs)) {
    lv <- levs[i]
    denom_all <- sum(!is.na(x))
    denom0 <- sum(!is.na(x[df$A == 0]))
    denom1 <- sum(!is.na(x[df$A == 1]))
    n_all <- sum(x == lv, na.rm = TRUE)
    n0 <- sum(x[df$A == 0] == lv, na.rm = TRUE)
    n1 <- sum(x[df$A == 1] == lv, na.rm = TRUE)
    out[[i]] <- data.frame(
      section = "Categorical",
      variable = label,
      level = lv,
      overall = fmt_n_pct(n_all, denom_all),
      resistant_0 = fmt_n_pct(n0, denom0),
      resistant_1 = fmt_n_pct(n1, denom1),
      stringsAsFactors = FALSE
    )
  }
  do.call(rbind, out)
}

binary_any_micro_info <- function(prior_state, prior_respiratory_state, sampled_lag1, carb_r_lag1, baseline_index_n, day) {
  prior_state <- safe_char(prior_state)
  prior_respiratory_state <- safe_char(prior_respiratory_state)
  sampled_lag1 <- clean_binary01(sampled_lag1)
  carb_r_lag1 <- clean_binary01(carb_r_lag1)
  baseline_index_n <- clean_num(baseline_index_n)

  has_prior_state <- ifelse(
    !is.na(prior_state) | !is.na(prior_respiratory_state) | sampled_lag1 == 1 | carb_r_lag1 == 1,
    1,
    0
  )
  has_baseline_index <- ifelse(!is.na(baseline_index_n) & baseline_index_n > 0, 1, 0)
  ifelse(day == 0, has_baseline_index, has_prior_state)
}

if (!file.exists(day03_path)) stop("Missing file: ", day03_path)
if (!file.exists(day60_path)) stop("Missing file: ", day60_path)
if (!file.exists(severity_path)) stop("Missing file: ", severity_path)
if (!file.exists(itt_rds_path)) stop("Missing file: ", itt_rds_path)

day03_raw <- read_excel(day03_path)
day60_raw <- read_excel(day60_path)
severity_raw <- read_excel(severity_path)
itt_raw <- readRDS(itt_rds_path)

required_day_cols <- c(
  "subjid", "day_from_day0", "calendar_date", "alive_at_start_of_day", "under_followup_on_day",
  "baseline_age", "baseline_male", "baseline_charlson", "baseline_sofa",
  "baseline_country_f", "baseline_site_f", "baseline_icu_type", "baseline_bacteria_f",
  "baseline_carba_r", "baseline_appropriateabx_on_symptomdate",
  "treatment_any_appropriate_active_primary", "treatment_cum_appropriate_days",
  "outcome_death60_fixed", "outcome_death_on_day", "outcome_death_by_end_of_day",
  "micro_sampled_lag1", "micro_any_carbapenem_R_lag1",
  "prior_collected_micro_observation_state", "prior_collected_micro_respiratory_state",
  "baseline_index_culture_record_n"
)
missing_day_cols <- setdiff(required_day_cols, names(day03_raw))
if (length(missing_day_cols) > 0) {
  stop("Day0-3 integrated file is missing columns: ", paste(missing_day_cols, collapse = ", "))
}

required_severity_cols <- c(
  "subjid", "day_from_day0", "severity_lt_definition", "severity_l0_sofa",
  "severity_lt_temp_max_c", "severity_lt_map_min", "severity_lt_hr_max",
  "severity_lt_any_inotrope", "severity_lt_mechanical_ventilation",
  "severity_lt_spo2fio2_min", "severity_lt_complete",
  "severity_temp_max_c", "severity_map_min", "severity_hr_max",
  "severity_any_inotrope", "severity_mechanical_ventilation", "severity_spo2fio2_min"
)
missing_sev_cols <- setdiff(required_severity_cols, names(severity_raw))
if (length(missing_sev_cols) > 0) {
  stop("Primary severity file is missing columns: ", paste(missing_sev_cols, collapse = ", "))
}

day03 <- as.data.frame(day03_raw, stringsAsFactors = FALSE)
day60 <- as.data.frame(day60_raw, stringsAsFactors = FALSE)
severity <- as.data.frame(severity_raw, stringsAsFactors = FALSE)

day03$subjid <- safe_char(day03$subjid)
day03$day <- clean_num(day03$day_from_day0)
day03$A <- clean_binary01(day03$baseline_carba_r)
day03$M <- clean_binary01(day03$treatment_any_appropriate_active_primary)
day03$cum_appropriate <- clean_num(day03$treatment_cum_appropriate_days)
day03$alive_at_start_of_day <- clean_binary01(day03$alive_at_start_of_day)
day03$under_followup_on_day <- clean_binary01(day03$under_followup_on_day)
day03$Y_day <- clean_binary01(day03$outcome_death_by_end_of_day)
day03$Y_day_incident <- clean_binary01(day03$outcome_death_on_day)
day03$O_info <- binary_any_micro_info(
  prior_state = day03$prior_collected_micro_observation_state,
  prior_respiratory_state = day03$prior_collected_micro_respiratory_state,
  sampled_lag1 = day03$micro_sampled_lag1,
  carb_r_lag1 = day03$micro_any_carbapenem_R_lag1,
  baseline_index_n = day03$baseline_index_culture_record_n,
  day = day03$day
)
day03$O_micro_carbR <- ifelse(day03$day == 0, day03$A, clean_binary01(day03$micro_any_carbapenem_R_lag1))

day03$C <- ifelse(
  is.na(day03$alive_at_start_of_day) | is.na(day03$under_followup_on_day),
  NA_real_,
  ifelse(day03$alive_at_start_of_day == 1 & day03$under_followup_on_day == 1, 0, 1)
)

day03 <- day03[day03$day %in% 0:3, , drop = FALSE]

severity$subjid <- safe_char(severity$subjid)
severity$day <- clean_num(severity$day_from_day0)
severity <- severity[severity$day %in% 0:3, , drop = FALSE]

severity$L_main <- ifelse(
  severity$day == 0,
  clean_num(severity$severity_l0_sofa),
  score_lt_from_components(
    clean_num(severity$severity_lt_temp_max_c),
    clean_num(severity$severity_lt_map_min),
    clean_num(severity$severity_lt_hr_max),
    clean_binary01(severity$severity_lt_any_inotrope),
    clean_binary01(severity$severity_lt_mechanical_ventilation),
    clean_num(severity$severity_lt_spo2fio2_min)
  )
)

severity$L_alt <- ifelse(
  severity$day == 0,
  clean_num(severity$severity_l0_sofa),
  score_lt_from_components(
    clean_num(severity$severity_temp_max_c),
    clean_num(severity$severity_map_min),
    clean_num(severity$severity_hr_max),
    clean_binary01(severity$severity_any_inotrope),
    clean_binary01(severity$severity_mechanical_ventilation),
    clean_num(severity$severity_spo2fio2_min)
  )
)

severity$lt_complete <- clean_binary01(severity$severity_lt_complete)
severity$lt_definition <- safe_char(severity$severity_lt_definition)

panel <- merge(
  day03,
  severity[, c("subjid", "day", "L_main", "L_alt", "lt_complete", "lt_definition")],
  by = c("subjid", "day"),
  all.x = TRUE,
  sort = FALSE
)

outcome_day60 <- day60[day60$day_from_day0 %in% c(60, "60"), c("subjid", "outcome_death60_fixed"), drop = FALSE]
if (nrow(outcome_day60) == 0) {
  outcome_day60 <- day60[!duplicated(day60$subjid), c("subjid", "outcome_death60_fixed"), drop = FALSE]
}
outcome_day60$subjid <- safe_char(outcome_day60$subjid)
outcome_day60$Y <- clean_binary01(outcome_day60$outcome_death60_fixed)
outcome_day60 <- outcome_day60[!duplicated(outcome_day60$subjid), c("subjid", "Y"), drop = FALSE]

panel <- merge(panel, outcome_day60, by = "subjid", all.x = TRUE, sort = FALSE)

panel$age <- clean_num(panel$baseline_age)
panel$male <- clean_binary01(panel$baseline_male)
panel$charlson <- clean_num(panel$baseline_charlson)
panel$country <- fill_missing_factor(panel$baseline_country_f)
panel$site <- fill_missing_factor(panel$baseline_site_f)
panel$icu_type <- fill_missing_factor(panel$baseline_icu_type)
panel$bacteria <- fill_missing_factor(panel$baseline_bacteria_f)
panel$baseline_sofa_raw <- clean_num(panel$baseline_sofa)
panel$severity_group <- factor(ifelse(panel$baseline_sofa_raw >= 6, "SOFA>=6", "SOFA<6"))
panel$country_simple <- fill_missing_factor(panel$baseline_country_f)
panel$bacteria_simple <- fill_missing_factor(ifelse(
  safe_char(panel$baseline_bacteria_f) %in% c("Klebsiella", "Acinetobacter", "Pseudomonas"),
  safe_char(panel$baseline_bacteria_f),
  "Other_or_Culture_neg"
))

panel <- panel[!is.na(panel$subjid) & !is.na(panel$A), , drop = FALSE]
panel <- panel[order(panel$subjid, panel$day), , drop = FALSE]
validate_day_index(panel)

panel$M[panel$alive_at_start_of_day == 0 | panel$under_followup_on_day == 0] <- NA_real_
panel$L_main[panel$alive_at_start_of_day == 0 | panel$under_followup_on_day == 0] <- NA_real_
panel$L_alt[panel$alive_at_start_of_day == 0 | panel$under_followup_on_day == 0] <- NA_real_
panel$O_info[panel$alive_at_start_of_day == 0 | panel$under_followup_on_day == 0] <- NA_real_
panel$O_micro_carbR[panel$alive_at_start_of_day == 0 | panel$under_followup_on_day == 0] <- NA_real_
panel$Y_day[panel$alive_at_start_of_day == 0 | panel$under_followup_on_day == 0] <- NA_real_
panel$Y_day_incident[panel$alive_at_start_of_day == 0 | panel$under_followup_on_day == 0] <- NA_real_

patient_level <- panel[panel$day == 0, c("subjid", "Y", "age", "male", "charlson", "country", "site", "icu_type", "bacteria"), drop = FALSE]
patient_level <- patient_level[!duplicated(patient_level$subjid), , drop = FALSE]

wide_keep <- c("subjid", "day", "A", "C", "M", "L_main", "L_alt", "O_info", "O_micro_carbR", "Y_day", "Y_day_incident")
wide <- reshape(panel[, wide_keep], idvar = "subjid", timevar = "day", direction = "wide")
names(wide) <- gsub("\\.", "", names(wide))
analysis <- merge(wide, patient_level, by = "subjid", all.x = TRUE, sort = FALSE)

required_wide <- c(
  "A0", "M0", "M1", "M2", "M3", "L_main0", "L_main1", "L_main2", "L_main3",
  "O_info0", "O_info1", "O_info2", "O_info3", "Y"
)
missing_wide <- setdiff(required_wide, names(analysis))
if (length(missing_wide) > 0) {
  stop("Wide panel is missing required variables: ", paste(missing_wide, collapse = ", "))
}

analysis$A <- analysis$A0
analysis$C0 <- analysis$C0
analysis$C1 <- analysis$C1
analysis$C2 <- analysis$C2
analysis$C3 <- analysis$C3
analysis$L0 <- analysis$L_main0
analysis$L1 <- analysis$L_main1
analysis$L2 <- analysis$L_main2
analysis$L3 <- analysis$L_main3
analysis$L0_alt <- analysis$L_alt0
analysis$L1_alt <- analysis$L_alt1
analysis$L2_alt <- analysis$L_alt2
analysis$L3_alt <- analysis$L_alt3
analysis$O0 <- analysis$O_info0
analysis$O1 <- analysis$O_info1
analysis$O2 <- analysis$O_info2
analysis$O3 <- analysis$O_info3

analysis <- analysis[!is.na(analysis$A) & !is.na(analysis$Y), , drop = FALSE]
assert_clinical_analysis_only(names(analysis))
analysis$baseline_sofa_raw <- clean_num(analysis$L0)
analysis$severity_group <- factor(ifelse(analysis$baseline_sofa_raw >= 6, "SOFA>=6", "SOFA<6"))
analysis$country_simple <- fill_missing_factor(analysis$country)
analysis$bacteria_simple <- fill_missing_factor(ifelse(
  safe_char(analysis$bacteria) %in% c("Klebsiella", "Acinetobacter", "Pseudomonas"),
  safe_char(analysis$bacteria),
  "Other_or_Culture_neg"
))
make_table1 <- function(panel_df) {
  base <- panel_df[panel_df$day == 0, , drop = FALSE]
  base <- base[!duplicated(base$subjid), , drop = FALSE]
  rows <- list(
    data.frame(
      section = "Header",
      variable = "N",
      level = "",
      overall = as.character(nrow(base)),
      resistant_0 = as.character(sum(base$A == 0, na.rm = TRUE)),
      resistant_1 = as.character(sum(base$A == 1, na.rm = TRUE)),
      stringsAsFactors = FALSE
    ),
    summarize_continuous(base, "age", "Age, years", stat = "mean_sd"),
    summarize_binary(base, "male", "Male sex"),
    summarize_continuous(base, "charlson", "Charlson comorbidity index", stat = "median_iqr"),
    summarize_continuous(base, "L_main", "Baseline severity (L0 / SOFA)", stat = "median_iqr"),
    summarize_categorical(base, "country", "Country"),
    summarize_categorical(base, "site", "Site"),
    summarize_categorical(base, "icu_type", "ICU type"),
    summarize_categorical(base, "bacteria", "Baseline bacteria group"),
    summarize_binary(base, "baseline_appropriate", "Appropriate antibiotics on symptom day"),
    summarize_binary(base, "Y", "60-day mortality")
  )
  table1 <- do.call(rbind, rows)
  rownames(table1) <- NULL
  list(base = base, table = table1)
}

make_table2 <- function(panel_df) {
  use <- panel_df[panel_df$day %in% 0:3, , drop = FALSE]
  mk_row <- function(df, a_value, day_value) {
    d <- df[df$A == a_value & df$day == day_value, , drop = FALSE]
    data.frame(
      group = ifelse(a_value == 1, "Resistant", "Non-resistant"),
      day = day_value,
      n_rows = nrow(d),
      n_treatment_observed = sum(!is.na(d$M)),
      n_severity_observed = sum(!is.na(d$L_main)),
      n_micro_info_observed = sum(!is.na(d$O_info)),
      appropriate_treatment = fmt_n_pct(sum(d$M == 1, na.rm = TRUE), sum(!is.na(d$M))),
      daily_severity_main = fmt_mean_sd(d$L_main, digits = 2),
      daily_severity_alt = fmt_mean_sd(d$L_alt, digits = 2),
      microbiology_info_available = fmt_n_pct(sum(d$O_info == 1, na.rm = TRUE), sum(!is.na(d$O_info))),
      prior_carbapenem_r_info = fmt_n_pct(sum(d$O_micro_carbR == 1, na.rm = TRUE), sum(!is.na(d$O_micro_carbR))),
      cumulative_appropriate_days = fmt_median_iqr(d$cum_appropriate, digits = 2),
      stringsAsFactors = FALSE
    )
  }
  rows <- list()
  for (a_value in c(0, 1)) {
    for (day_value in 0:3) {
      rows[[length(rows) + 1]] <- mk_row(use, a_value, day_value)
    }
  }
  do.call(rbind, rows)
}

make_trajectory_table <- function(panel_df) {
  use <- panel_df[panel_df$day %in% 0:3, c("subjid", "A", "day", "M", "O_info"), drop = FALSE]
  wide_m <- reshape(use, idvar = c("subjid", "A"), timevar = "day", direction = "wide")
  names(wide_m) <- sub("^M\\.", "M", names(wide_m))
  names(wide_m) <- sub("^O_info\\.", "O", names(wide_m))
  wide_m$treatment_trajectory <- paste0(wide_m$M0, wide_m$M1, wide_m$M2, wide_m$M3)
  wide_m$micro_info_trajectory <- paste0(wide_m$O0, wide_m$O1, wide_m$O2, wide_m$O3)
  out <- as.data.frame(
    table(
      group = ifelse(wide_m$A == 1, "Resistant", "Non-resistant"),
      treatment_trajectory = wide_m$treatment_trajectory,
      micro_info_trajectory = wide_m$micro_info_trajectory
    ),
    stringsAsFactors = FALSE
  )
  names(out) <- c("group", "treatment_trajectory", "micro_info_trajectory", "n")
  out[out$n > 0, , drop = FALSE]
}

draw_figure3_panels <- function(sum_df) {
  cols <- c("Non-resistant" = "#1f77b4", "Resistant" = "#d62728")
  tmp0 <- sum_df[sum_df$group == "Non-resistant", ]
  tmp1 <- sum_df[sum_df$group == "Resistant", ]
  ylim_l <- range(sum_df$L_main, na.rm = TRUE)

  graphics::plot(tmp0$day, tmp0$M, type = "o", pch = 16, col = cols["Non-resistant"], ylim = c(0, 1),
                 xlab = "Day from baseline", ylab = "Observed appropriate-treatment rate",
                 main = "A. Daily appropriate-treatment rate", xaxt = "n", lwd = 2)
  graphics::axis(1, at = 0:3, labels = paste("Day", 0:3))
  graphics::lines(tmp1$day, tmp1$M, type = "o", pch = 17, col = cols["Resistant"], lwd = 2)
  graphics::legend("bottomright", legend = c("Non-resistant", "Resistant"), col = cols,
                   pch = c(16, 17), lwd = 2, bty = "n")

  graphics::plot(tmp0$day, tmp0$L_main, type = "o", pch = 16, col = cols["Non-resistant"], ylim = ylim_l,
                 xlab = "Day from baseline", ylab = "Observed daily severity score",
                 main = "B. Daily severity trajectory", xaxt = "n", lwd = 2)
  graphics::axis(1, at = 0:3, labels = paste("Day", 0:3))
  graphics::lines(tmp1$day, tmp1$L_main, type = "o", pch = 17, col = cols["Resistant"], lwd = 2)
  graphics::legend("topright", legend = c("Non-resistant", "Resistant"), col = cols,
                   pch = c(16, 17), lwd = 2, bty = "n")

  graphics::plot(tmp0$day, tmp0$O_info, type = "o", pch = 16, col = cols["Non-resistant"], ylim = c(0, 1),
                 xlab = "Day from baseline", ylab = "Prior microbiology info available",
                 main = "C. O_t information availability", xaxt = "n", lwd = 2)
  graphics::axis(1, at = 0:3, labels = paste("Day", 0:3))
  graphics::lines(tmp1$day, tmp1$O_info, type = "o", pch = 17, col = cols["Resistant"], lwd = 2)
  graphics::legend("bottomright", legend = c("Non-resistant", "Resistant"), col = cols,
                   pch = c(16, 17), lwd = 2, bty = "n")
}

make_figure3 <- function(panel_df, png_path, svg_path = NULL) {
  use <- panel_df[panel_df$day %in% 0:3, , drop = FALSE]
  sum_df <- stats::aggregate(
    cbind(M = use$M, L_main = use$L_main, O_info = use$O_info) ~ A + day,
    data = use,
    FUN = function(x) mean(x, na.rm = TRUE)
  )
  sum_df$group <- ifelse(sum_df$A == 1, "Resistant", "Non-resistant")

  grDevices::png(png_path, width = 1800, height = 600, res = 150)
  op <- graphics::par(mfrow = c(1, 3), mar = c(5, 5, 4, 2) + 0.1, family = "serif")
  draw_figure3_panels(sum_df)
  graphics::par(op)
  grDevices::dev.off()

  if (!is.null(svg_path)) {
    tryCatch({
      suppressWarnings(grDevices::svg(svg_path, width = 12, height = 4.2, pointsize = 12, family = "serif"))
      op <- graphics::par(mfrow = c(1, 3), mar = c(5, 5, 4, 2) + 0.1, family = "serif")
      draw_figure3_panels(sum_df)
      graphics::par(op)
      grDevices::dev.off()
    }, error = function(e) {
      message("Skipping SVG export for Figure 3: ", conditionMessage(e))
      if (names(grDevices::dev.cur()) != "null device") {
        grDevices::dev.off()
      }
    })
  }

  sum_df
}

run_gformula_legacy <- function(analysis_df, mediator_history = c("lagged", "none"), lt_variant = c("main", "alt"), include_ot = FALSE, nsim = 4000, seed = 20260903) {
  mediator_history <- match.arg(mediator_history)
  lt_variant <- match.arg(lt_variant)
  set.seed(seed)

  df <- analysis_df
  if (lt_variant == "alt") {
    df$L0 <- df$L0_alt
    df$L1 <- df$L1_alt
    df$L2 <- df$L2_alt
    df$L3 <- df$L3_alt
  }

  # Do not adjust for baseline_appropriateabx_on_symptomdate: in this extract it
  # is essentially M0 and would adjust for the mediator itself.
  covars <- c("age", "male", "charlson", "country", "site", "icu_type", "bacteria")
  ot_terms <- if (include_ot) c("O0", "O1", "O2", "O3") else character(0)
  fit_vars <- unique(c("A", "Y", "M0", "M1", "M2", "M3", "L0", "L1", "L2", "L3", covars, ot_terms))
  fit_df <- df[, fit_vars, drop = FALSE]
  fit_df <- fit_df[stats::complete.cases(fit_df), , drop = FALSE]
  if (nrow(fit_df) < 100) stop("Too few complete cases for g-formula fitting.")
  if (!is_informative_variable(fit_df$A)) stop("Exposure A has fewer than 2 observed levels in this subset.")

  pruned <- drop_invariant_columns(fit_df, c(covars, ot_terms))
  fit_df <- pruned$data
  covars <- intersect(covars, pruned$keep)
  ot_terms <- intersect(ot_terms, pruned$keep)

  ot_term_for_day <- function(day) {
    term <- paste0("O", day)
    if (term %in% ot_terms) term else character(0)
  }

  mediator_formula <- function(t) {
    ot_term <- ot_term_for_day(t)
    if (t == 0) {
      stats::as.formula(paste("M0 ~", rhs_join(c("A", "L0", ot_term, covars))))
    } else {
      prev_m <- paste0("M", t - 1)
      if (mediator_history == "lagged") {
        stats::as.formula(paste0("M", t, " ~ ", rhs_join(c("A", paste0("L", t), prev_m, ot_term, "L0", covars))))
      } else {
        stats::as.formula(paste0("M", t, " ~ ", rhs_join(c("A", paste0("L", t), ot_term, "L0", covars))))
      }
    }
  }

  l_formula <- function(t) {
    prev_l <- paste0("L", t - 1)
    prev_m <- paste0("M", t - 1)
    ot_term <- ot_term_for_day(t - 1)
    stats::as.formula(paste0("L", t, " ~ ", rhs_join(c("A", prev_l, prev_m, ot_term, covars))))
  }

  y_terms <- c("A", "L0", "L1", "L2", "L3", covars)
  if (mediator_history == "lagged") {
    y_terms <- c(y_terms, "M0", "M1", "M2", "M3")
  } else {
    y_terms <- c(y_terms, "M3")
  }
  y_terms <- c(y_terms, ot_terms)
  y_formula <- stats::as.formula(paste("Y ~", rhs_join(y_terms)))

  m0_fit <- fit_binomial_model(mediator_formula(0), fit_df)
  m1_fit <- fit_binomial_model(mediator_formula(1), fit_df)
  m2_fit <- fit_binomial_model(mediator_formula(2), fit_df)
  m3_fit <- fit_binomial_model(mediator_formula(3), fit_df)
  l1_fit <- fit_gaussian_model(l_formula(1), fit_df)
  l2_fit <- fit_gaussian_model(l_formula(2), fit_df)
  l3_fit <- fit_gaussian_model(l_formula(3), fit_df)
  y_fit <- fit_binomial_model(y_formula, fit_df)

  base_cols <- unique(c(covars, "L0", ot_terms))
  sim_base <- fit_df[rep(seq_len(nrow(fit_df)), length.out = nsim), base_cols, drop = FALSE]

  simulate_regime <- function(a_for_y, a_for_m) {
    sim <- sim_base
    sim$A <- a_for_y

    m_dat0 <- sim
    m_dat0$A <- a_for_m
    p0 <- predict_binomial_prob(m0_fit, m_dat0)
    sim$M0 <- stats::rbinom(nrow(sim), 1, p0)

    sim$L1 <- predict_gaussian_draw(l1_fit, sim)
    m_dat1 <- sim
    m_dat1$A <- a_for_m
    p1 <- predict_binomial_prob(m1_fit, m_dat1)
    sim$M1 <- stats::rbinom(nrow(sim), 1, p1)

    sim$L2 <- predict_gaussian_draw(l2_fit, sim)
    m_dat2 <- sim
    m_dat2$A <- a_for_m
    p2 <- predict_binomial_prob(m2_fit, m_dat2)
    sim$M2 <- stats::rbinom(nrow(sim), 1, p2)

    sim$L3 <- predict_gaussian_draw(l3_fit, sim)
    m_dat3 <- sim
    m_dat3$A <- a_for_m
    p3 <- predict_binomial_prob(m3_fit, m_dat3)
    sim$M3 <- stats::rbinom(nrow(sim), 1, p3)

    y_dat <- sim
    y_dat$A <- a_for_y
    py <- predict_binomial_prob(y_fit, y_dat)
    sim$Y <- stats::rbinom(nrow(sim), 1, py)
    sim
  }

  r11 <- estimate_risk(simulate_regime(1, 1))
  r10 <- estimate_risk(simulate_regime(1, 0))
  r00 <- estimate_risk(simulate_regime(0, 0))

  data.frame(
    lt_variant = lt_variant,
    mediator_history = mediator_history,
    include_ot = include_ot,
    nsim = nsim,
    n_complete = nrow(fit_df),
    active_covars = paste(covars, collapse = ";"),
    dropped_covars = paste(pruned$drop, collapse = ";"),
    R_1_G1 = r11,
    R_1_G0 = r10,
    R_0_G0 = r00,
    TE = r11 - r00,
    IDE = r10 - r00,
    IIE = r11 - r10,
    stringsAsFactors = FALSE
  )
}

run_hte_gformula <- function(analysis_df, modifier_var, modifier_label, levels_keep = NULL, min_n = 40, min_complete = 60) {
  df <- analysis_df
  x <- df[[modifier_var]]
  if (is.null(levels_keep)) {
    levels_keep <- sort(unique(as.character(x[!is.na(x)])))
  }

  rows <- list()
  for (lv in levels_keep) {
    sub_df <- df[as.character(x) == lv, , drop = FALSE]
    a_levels <- sort(unique(sub_df$A[!is.na(sub_df$A)]))
    if (nrow(sub_df) < min_n || length(a_levels) < 2) {
      rows[[length(rows) + 1]] <- data.frame(
        modifier = modifier_label,
        level = lv,
        n = nrow(sub_df),
        n_A0 = sum(sub_df$A == 0, na.rm = TRUE),
        n_A1 = sum(sub_df$A == 1, na.rm = TRUE),
        status = "skipped_insufficient_variation",
        n_complete = NA_integer_,
        R_1_G1 = NA_real_,
        R_1_G0 = NA_real_,
        R_0_G0 = NA_real_,
        TE = NA_real_,
        IDE = NA_real_,
        IIE = NA_real_,
        stringsAsFactors = FALSE
      )
      next
    }

    fit_check <- sub_df[, c("A", "Y", "M0", "M1", "M2", "M3", "L0", "L1", "L2", "L3", "age", "male", "charlson", "country", "site", "icu_type", "bacteria"), drop = FALSE]
    n_complete <- sum(stats::complete.cases(fit_check))
    if (n_complete < min_complete) {
      rows[[length(rows) + 1]] <- data.frame(
        modifier = modifier_label,
        level = lv,
        n = nrow(sub_df),
        n_A0 = sum(sub_df$A == 0, na.rm = TRUE),
        n_A1 = sum(sub_df$A == 1, na.rm = TRUE),
        status = "skipped_low_complete_cases",
        n_complete = n_complete,
        R_1_G1 = NA_real_,
        R_1_G0 = NA_real_,
        R_0_G0 = NA_real_,
        TE = NA_real_,
        IDE = NA_real_,
        IIE = NA_real_,
        stringsAsFactors = FALSE
      )
      next
    }

    fit <- tryCatch(
      run_gformula(sub_df, mediator_history = "lagged", lt_variant = "main", include_ot = FALSE, nsim = 3000, seed = 20260903),
      error = function(e) e
    )

    if (inherits(fit, "error")) {
      rows[[length(rows) + 1]] <- data.frame(
        modifier = modifier_label,
        level = lv,
        n = nrow(sub_df),
        n_A0 = sum(sub_df$A == 0, na.rm = TRUE),
        n_A1 = sum(sub_df$A == 1, na.rm = TRUE),
        status = paste0("failed: ", fit$message),
        n_complete = n_complete,
        R_1_G1 = NA_real_,
        R_1_G0 = NA_real_,
        R_0_G0 = NA_real_,
        TE = NA_real_,
        IDE = NA_real_,
        IIE = NA_real_,
        stringsAsFactors = FALSE
      )
    } else {
      rows[[length(rows) + 1]] <- data.frame(
        modifier = modifier_label,
        level = lv,
        n = nrow(sub_df),
        n_A0 = sum(sub_df$A == 0, na.rm = TRUE),
        n_A1 = sum(sub_df$A == 1, na.rm = TRUE),
        status = "ok",
        n_complete = fit$n_complete,
        R_1_G1 = fit$R_1_G1,
        R_1_G0 = fit$R_1_G0,
        R_0_G0 = fit$R_0_G0,
        TE = fit$TE,
        IDE = fit$IDE,
        IIE = fit$IIE,
        stringsAsFactors = FALSE
      )
    }
  }

  do.call(rbind, rows)
}

table1_obj <- make_table1(panel)
table2 <- make_table2(panel)
trajectory_tab <- make_trajectory_table(panel)
figure3_sum <- make_figure3(
  panel,
  file.path(step_dirs[["step02"]], "figure3_observed_processes.png"),
  file.path(step_dirs[["step02"]], "figure3_observed_processes.svg")
)

source_audit <- data.frame(
  dataset = c("day0_3_integrated", "day0_60_integrated", "primary_severity", "itt_rds_screen", "itt_rds_randomisation", "itt_rds_vitalsigns"),
  source_path = c(day03_path, day60_path, severity_path, itt_rds_path, itt_rds_path, itt_rds_path),
  n_rows = c(
    nrow(day03),
    nrow(day60),
    nrow(severity_raw),
    if ("screen" %in% names(itt_raw)) nrow(itt_raw$screen) else NA_integer_,
    if ("randomisation" %in% names(itt_raw)) nrow(itt_raw$randomisation) else NA_integer_,
    if ("vitalsigns" %in% names(itt_raw)) nrow(itt_raw$vitalsigns) else NA_integer_
  ),
  note = c(
    "Main longitudinal panel for Day 0-3 and mediator M_t",
    "60-day follow-up file used for outcome Y",
    "Primary source for L_t definition in the mediation window",
    "Original CRF module bundle",
    "Original treatment allocation module bundle",
    "Original daily physiology source behind severity construction"
  ),
  stringsAsFactors = FALSE
)

censoring_summary <- data.frame(
  day = 0:3,
  n_censored_at_start = sapply(0:3, function(d) sum(panel$C[panel$day == d] == 1, na.rm = TRUE)),
  n_death_by_end_of_day = sapply(0:3, function(d) sum(panel$Y_day[panel$day == d] == 1, na.rm = TRUE)),
  n_micro_info_available = sapply(0:3, function(d) sum(panel$O_info[panel$day == d] == 1, na.rm = TRUE)),
  stringsAsFactors = FALSE
)

lt_qc <- stats::aggregate(
  cbind(L_main = panel$L_main, L_alt = panel$L_alt, O_info = panel$O_info) ~ day,
  data = panel,
  FUN = function(x) sum(!is.na(x))
)
lt_complete_tab <- as.data.frame(with(severity, table(day, lt_complete, useNA = "ifany")))

source(file.path(module_dir, "R", "04_gformula_engine.R"))
source(file.path(module_dir, "R", "05_sensitivity_runner.R"))

main_gf <- run_gformula(analysis, mediator_history = "lagged", lt_variant = "main", include_ot = FALSE, nsim = 4000)
table3 <- main_gf[, c("n_complete", "R_1_G1", "R_1_G0", "R_0_G0", "TE", "IDE", "IIE")]
table4 <- run_prespecified_sensitivities(analysis, nsim = 4000)
write_cohort_diagnostics(panel, analysis, step_dirs[["step00"]])
write_observed_trajectory_diagnostics(panel, step_dirs[["step02"]])
write_positivity_diagnostics(analysis, step_dirs[["step03"]])

bootstrap_n <- suppressWarnings(as.integer(Sys.getenv("REGARDVAP_N_BOOT", unset = "0")))
if (!is.na(bootstrap_n) && bootstrap_n > 0L) {
  bootstrap_nsim <- suppressWarnings(as.integer(Sys.getenv("REGARDVAP_BOOT_NSIM", unset = "4000")))
  if (is.na(bootstrap_nsim) || bootstrap_nsim < 100L) bootstrap_nsim <- 4000L
  bootstrap_results <- bootstrap_gformula(analysis, B = bootstrap_n, nsim = bootstrap_nsim)
  utils::write.csv(bootstrap_results, file.path(step_dirs[["step03"]], "bootstrap_effect_estimates.csv"), row.names = FALSE)
  utils::write.csv(bootstrap_percentile_ci(bootstrap_results), file.path(step_dirs[["step03"]], "bootstrap_percentile_ci.csv"), row.names = FALSE)
}

utils::write.csv(panel, file.path(step_dirs[["step00"]], "analysis_panel_long_day0_day3.csv"), row.names = FALSE)
utils::write.csv(analysis, file.path(step_dirs[["step00"]], "analysis_panel_wide_for_gformula.csv"), row.names = FALSE)
utils::write.csv(source_audit, file.path(step_dirs[["step00"]], "data_sources_used.csv"), row.names = FALSE)
utils::write.csv(lt_qc, file.path(step_dirs[["step00"]], "lt_nonmissing_counts_by_day.csv"), row.names = FALSE)
utils::write.csv(lt_complete_tab, file.path(step_dirs[["step00"]], "lt_complete_table.csv"), row.names = FALSE)
utils::write.csv(censoring_summary, file.path(step_dirs[["step00"]], "censoring_and_daily_event_summary.csv"), row.names = FALSE)

utils::write.csv(table1_obj$base, file.path(step_dirs[["step01"]], "table1_patient_level_dataset.csv"), row.names = FALSE)
utils::write.csv(table1_obj$table, file.path(step_dirs[["step01"]], "table1_baseline_summary.csv"), row.names = FALSE)
utils::write.csv(table2, file.path(step_dirs[["step02"]], "table2_longitudinal_day0_day3_summary.csv"), row.names = FALSE)
utils::write.csv(trajectory_tab, file.path(step_dirs[["step02"]], "table2_observed_trajectory_distribution.csv"), row.names = FALSE)
utils::write.csv(figure3_sum, file.path(step_dirs[["step02"]], "figure3_daywise_summary.csv"), row.names = FALSE)

utils::write.csv(table3, file.path(step_dirs[["step03"]], "table3_main_gformula_estimates.csv"), row.names = FALSE)
utils::write.csv(table4, file.path(step_dirs[["step04"]], "table4_sensitivity_analyses.csv"), row.names = FALSE)

notes <- c(
  "# REGARD-VAP longitudinal mediation pipeline v3",
  "",
  "Pipeline structure:",
  "- step00_data_prep: analysis-ready long and wide panels, source audit, missingness checks",
  "- step01_baseline: baseline Table 1 style summary by baseline carbapenem resistance",
  "- step02_longitudinal_summary: Day 0-3 observed treatment, severity, and O_t information summary",
  "- step03_main_gformula: primary interventional direct and indirect effect estimates",
  "- step04_sensitivity: alternative L_t definition, no-lag mediator history, and optional O_t-adjusted model",
  "- step05_hte: reserved for future heterogeneity analyses; no results are generated",
  "",
  "Variable system used in this v3 draft:",
  "- V = baseline covariates: age, sex, Charlson, country, site, ICU type, and bacteria group",
  "- A = baseline carbapenem resistance",
  "- L_t = daily clinical state; primary model uses Day 0 SOFA plus Day 1-3 pragmatic LT score",
  "- O_t = binary proxy for microbiology information available before each daily treatment decision",
  "- M_t = daily appropriate active treatment",
  "- Y = 60-day mortality",
  "",
  "Interpretation note:",
  "- The main model still matches your current binary A = 0/1 workflow.",
  "- The O_t piece is included as a draft data-structure extension inspired by the newer slide deck, but it should be reviewed before becoming the final manuscript model.",
  "- Heterogeneity analyses are not run in this version.",
  "",
  "Inputs:",
  paste0("- Day 0-3 integrated longitudinal file: ", day03_path),
  paste0("- Day 0-60 integrated longitudinal file: ", day60_path),
  paste0("- Primary severity file: ", severity_path),
  paste0("- ITT original RDS bundle: ", itt_rds_path)
)
writeLines(notes, file.path(out_dir, "README_regardvap_gformula_latest_v3.md"))

cat("Analysis finished.\n")
cat("Output directory:\n")
cat(out_dir, "\n")
