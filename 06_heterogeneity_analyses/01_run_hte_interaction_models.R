args <- commandArgs(trailingOnly = TRUE)

suppressPackageStartupMessages({
  library(readxl)
  library(ggplot2)
  library(splines)
})

default_day03_path <- "/Users/bing/Desktop/中介分析 最新数据（最新sofa定义/副本regardvap_patient_day0_day3_integrated_longitudinal_itt460_v.xlsx"
default_day60_path <- "/Users/bing/Desktop/中介分析 最新数据（最新sofa定义/副本regardvap_patient_day0_day60_integrated_longitudinal_itt460_.xlsx"
default_out_dir <- "/Users/bing/Documents/中介效应/outputs/hte_interaction_model_v1"

day03_path <- if (length(args) >= 1) args[1] else default_day03_path
day60_path <- if (length(args) >= 2) args[2] else default_day60_path
out_dir <- if (length(args) >= 3) args[3] else default_out_dir

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

safe_char <- function(x) {
  out <- as.character(x)
  out[out %in% c("", "NA", "NaN")] <- NA_character_
  out
}

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

trim_levels <- function(x) {
  factor(safe_char(x))
}

predict_prob <- function(model, newdata) {
  out <- suppressWarnings(stats::predict(model, newdata = newdata, type = "response"))
  as.numeric(out)
}

safe_lrt_p <- function(full_model, reduced_model) {
  cmp <- tryCatch(stats::anova(reduced_model, full_model, test = "LRT"), error = function(e) NULL)
  if (is.null(cmp) || nrow(cmp) < 2 || !"Pr(>Chi)" %in% names(cmp)) return(NA_real_)
  as.numeric(cmp[2, "Pr(>Chi)"])
}

bootstrap_effect_curve <- function(data, fit_fn, effect_fn, n_boot = 120, seed = 20260903) {
  set.seed(seed)
  boot_list <- vector("list", n_boot)
  n <- nrow(data)
  for (b in seq_len(n_boot)) {
    idx <- sample.int(n, size = n, replace = TRUE)
    boot_data <- data[idx, , drop = FALSE]
    fit <- tryCatch(fit_fn(boot_data), error = function(e) NULL)
    if (is.null(fit)) next
    out <- tryCatch(effect_fn(fit, boot_data), error = function(e) NULL)
    if (is.null(out)) next
    boot_list[[b]] <- out
  }
  ok <- boot_list[!vapply(boot_list, is.null, logical(1))]
  if (length(ok) == 0) return(NULL)
  do.call(rbind, Map(function(df, i) transform(df, boot_id = i), ok, seq_along(ok)))
}

if (!file.exists(day03_path)) stop("Missing file: ", day03_path)
if (!file.exists(day60_path)) stop("Missing file: ", day60_path)

day03_raw <- as.data.frame(read_excel(day03_path), stringsAsFactors = FALSE)
day60_raw <- as.data.frame(read_excel(day60_path), stringsAsFactors = FALSE)

base <- day03_raw[day03_raw$day_from_day0 %in% c(0, "0"), , drop = FALSE]
base <- base[!duplicated(base$subjid), , drop = FALSE]

outcome_day60 <- day60_raw[day60_raw$day_from_day0 %in% c(60, "60"), c("subjid", "outcome_death60_fixed"), drop = FALSE]
if (nrow(outcome_day60) == 0) {
  outcome_day60 <- day60_raw[!duplicated(day60_raw$subjid), c("subjid", "outcome_death60_fixed"), drop = FALSE]
}
outcome_day60$subjid <- safe_char(outcome_day60$subjid)
outcome_day60$Y <- clean_binary01(outcome_day60$outcome_death60_fixed)
outcome_day60 <- outcome_day60[!duplicated(outcome_day60$subjid), c("subjid", "Y"), drop = FALSE]

base$subjid <- safe_char(base$subjid)
base <- merge(base, outcome_day60, by = "subjid", all.x = TRUE, sort = FALSE)

base$A <- clean_binary01(base$baseline_carba_r)
base$baseline_sofa <- clean_num(base$baseline_sofa)
base$baseline_age <- clean_num(base$baseline_age)
base$baseline_male <- clean_binary01(base$baseline_male)
base$baseline_charlson <- clean_num(base$baseline_charlson)
base$baseline_appropriate <- clean_binary01(base$baseline_appropriateabx_on_symptomdate)
base$baseline_bacteria_f <- trim_levels(base$baseline_bacteria_f)
base$baseline_country_f <- trim_levels(base$baseline_country_f)

# Keep the HTE dataset close to the mediator analysis target population:
# baseline-defined, non-missing, culture-positive patients.
model_df <- base[
  !is.na(base$A) &
    !is.na(base$Y) &
    !is.na(base$baseline_sofa) &
    !is.na(base$baseline_age) &
    !is.na(base$baseline_male) &
    !is.na(base$baseline_charlson) &
    !is.na(base$baseline_appropriate) &
    !is.na(base$baseline_bacteria_f) &
    !is.na(base$baseline_country_f) &
    as.character(base$baseline_bacteria_f) != "Culture_neg",
  c("subjid", "Y", "A", "baseline_sofa", "baseline_age", "baseline_male",
    "baseline_charlson", "baseline_appropriate", "baseline_bacteria_f", "baseline_country_f"),
  drop = FALSE
]

model_df$baseline_bacteria_f <- droplevels(model_df$baseline_bacteria_f)
model_df$baseline_country_f <- droplevels(model_df$baseline_country_f)

fit_sofa_model <- function(df) {
  glm(
    Y ~ A * ns(baseline_sofa, df = 3) +
      baseline_age + baseline_male + baseline_charlson +
      baseline_appropriate + baseline_bacteria_f + baseline_country_f,
    data = df,
    family = stats::binomial()
  )
}

fit_sofa_reduced_model <- function(df) {
  glm(
    Y ~ A + ns(baseline_sofa, df = 3) +
      baseline_age + baseline_male + baseline_charlson +
      baseline_appropriate + baseline_bacteria_f + baseline_country_f,
    data = df,
    family = stats::binomial()
  )
}

fit_bacteria_model <- function(df) {
  glm(
    Y ~ A * baseline_bacteria_f +
      baseline_sofa + baseline_age + baseline_male + baseline_charlson +
      baseline_appropriate + baseline_country_f,
    data = df,
    family = stats::binomial()
  )
}

fit_bacteria_reduced_model <- function(df) {
  glm(
    Y ~ A + baseline_bacteria_f +
      baseline_sofa + baseline_age + baseline_male + baseline_charlson +
      baseline_appropriate + baseline_country_f,
    data = df,
    family = stats::binomial()
  )
}

fit_country_model <- function(df) {
  glm(
    Y ~ A * baseline_country_f +
      baseline_sofa + baseline_age + baseline_male + baseline_charlson +
      baseline_appropriate + baseline_bacteria_f,
    data = df,
    family = stats::binomial()
  )
}

fit_country_reduced_model <- function(df) {
  glm(
    Y ~ A + baseline_country_f +
      baseline_sofa + baseline_age + baseline_male + baseline_charlson +
      baseline_appropriate + baseline_bacteria_f,
    data = df,
    family = stats::binomial()
  )
}

make_sofa_effect_curve <- function(model, data, n_points = 80) {
  sofa_grid <- seq(
    stats::quantile(data$baseline_sofa, 0.05, na.rm = TRUE),
    stats::quantile(data$baseline_sofa, 0.95, na.rm = TRUE),
    length.out = n_points
  )
  out <- vector("list", length(sofa_grid))
  for (i in seq_along(sofa_grid)) {
    sofa_value <- sofa_grid[i]
    nd1 <- data
    nd0 <- data
    nd1$A <- 1
    nd0$A <- 0
    nd1$baseline_sofa <- sofa_value
    nd0$baseline_sofa <- sofa_value
    p1 <- predict_prob(model, nd1)
    p0 <- predict_prob(model, nd0)
    out[[i]] <- data.frame(
      baseline_sofa = sofa_value,
      risk_A1 = mean(p1, na.rm = TRUE),
      risk_A0 = mean(p0, na.rm = TRUE),
      risk_difference = mean(p1 - p0, na.rm = TRUE)
    )
  }
  do.call(rbind, out)
}

make_factor_effect_table <- function(model, data, factor_var) {
  levs <- levels(data[[factor_var]])
  out <- vector("list", length(levs))
  for (i in seq_along(levs)) {
    lv <- levs[i]
    nd1 <- data
    nd0 <- data
    nd1$A <- 1
    nd0$A <- 0
    nd1[[factor_var]] <- factor(lv, levels = levs)
    nd0[[factor_var]] <- factor(lv, levels = levs)
    p1 <- predict_prob(model, nd1)
    p0 <- predict_prob(model, nd0)
    out[[i]] <- data.frame(
      modifier = factor_var,
      level = lv,
      n_observed = sum(as.character(data[[factor_var]]) == lv, na.rm = TRUE),
      risk_A1 = mean(p1, na.rm = TRUE),
      risk_A0 = mean(p0, na.rm = TRUE),
      risk_difference = mean(p1 - p0, na.rm = TRUE),
      stringsAsFactors = FALSE
    )
  }
  do.call(rbind, out)
}

sofa_fit <- fit_sofa_model(model_df)
bacteria_fit <- fit_bacteria_model(model_df)
country_fit <- fit_country_model(model_df)
sofa_fit_reduced <- fit_sofa_reduced_model(model_df)
bacteria_fit_reduced <- fit_bacteria_reduced_model(model_df)
country_fit_reduced <- fit_country_reduced_model(model_df)

sofa_curve <- make_sofa_effect_curve(sofa_fit, model_df)
bacteria_effects <- make_factor_effect_table(bacteria_fit, model_df, "baseline_bacteria_f")
country_effects <- make_factor_effect_table(country_fit, model_df, "baseline_country_f")

interaction_test_table <- data.frame(
  modifier = c("baseline_sofa", "baseline_bacteria_f", "baseline_country_f"),
  model_type = c("continuous spline interaction", "categorical interaction", "categorical interaction"),
  n = nrow(model_df),
  global_interaction_p = c(
    safe_lrt_p(sofa_fit, sofa_fit_reduced),
    safe_lrt_p(bacteria_fit, bacteria_fit_reduced),
    safe_lrt_p(country_fit, country_fit_reduced)
  ),
  stringsAsFactors = FALSE
)

bacteria_balance <- as.data.frame.matrix(with(model_df, table(baseline_bacteria_f, A)))
bacteria_balance$level <- rownames(bacteria_balance)
rownames(bacteria_balance) <- NULL
names(bacteria_balance)[names(bacteria_balance) == "0"] <- "n_A0"
names(bacteria_balance)[names(bacteria_balance) == "1"] <- "n_A1"
bacteria_balance$modifier <- "baseline_bacteria_f"
bacteria_balance <- bacteria_balance[, c("modifier", "level", "n_A0", "n_A1")]

country_balance <- as.data.frame.matrix(with(model_df, table(baseline_country_f, A)))
country_balance$level <- rownames(country_balance)
rownames(country_balance) <- NULL
names(country_balance)[names(country_balance) == "0"] <- "n_A0"
names(country_balance)[names(country_balance) == "1"] <- "n_A1"
country_balance$modifier <- "baseline_country_f"
country_balance <- country_balance[, c("modifier", "level", "n_A0", "n_A1")]

balance_table <- rbind(bacteria_balance, country_balance)

sofa_boot <- bootstrap_effect_curve(
  data = model_df,
  fit_fn = fit_sofa_model,
  effect_fn = function(model, data) make_sofa_effect_curve(model, data, n_points = 80)[, c("baseline_sofa", "risk_difference")],
  n_boot = 120,
  seed = 20260903
)

if (!is.null(sofa_boot)) {
  sofa_split <- split(sofa_boot$risk_difference, sofa_boot$baseline_sofa)
  sofa_curve$rd_low <- vapply(
    as.character(sofa_curve$baseline_sofa),
    function(k) stats::quantile(sofa_split[[k]], 0.025, na.rm = TRUE, names = FALSE),
    numeric(1)
  )
  sofa_curve$rd_high <- vapply(
    as.character(sofa_curve$baseline_sofa),
    function(k) stats::quantile(sofa_split[[k]], 0.975, na.rm = TRUE, names = FALSE),
    numeric(1)
  )
} else {
  sofa_curve$rd_low <- NA_real_
  sofa_curve$rd_high <- NA_real_
}

factor_boot_ci <- function(data, fit_fn, factor_var, seed_offset = 0) {
  boot <- bootstrap_effect_curve(
    data = data,
    fit_fn = fit_fn,
    effect_fn = function(model, data) make_factor_effect_table(model, data, factor_var)[, c("level", "risk_difference")],
    n_boot = 120,
    seed = 20260903 + seed_offset
  )
  if (is.null(boot)) return(NULL)
  stats::aggregate(risk_difference ~ level, data = boot, FUN = function(x) stats::quantile(x, c(0.025, 0.975), na.rm = TRUE, names = FALSE))
}

bacteria_ci <- factor_boot_ci(model_df, fit_bacteria_model, "baseline_bacteria_f", seed_offset = 10)
country_ci <- factor_boot_ci(model_df, fit_country_model, "baseline_country_f", seed_offset = 20)

if (!is.null(bacteria_ci)) {
  bacteria_map <- setNames(bacteria_ci$risk_difference, bacteria_ci$level)
  bacteria_effects$rd_low <- vapply(as.character(bacteria_effects$level), function(k) bacteria_map[[k]][1], numeric(1))
  bacteria_effects$rd_high <- vapply(as.character(bacteria_effects$level), function(k) bacteria_map[[k]][2], numeric(1))
} else {
  bacteria_effects$rd_low <- NA_real_
  bacteria_effects$rd_high <- NA_real_
}

if (!is.null(country_ci)) {
  country_map <- setNames(country_ci$risk_difference, country_ci$level)
  country_effects$rd_low <- vapply(as.character(country_effects$level), function(k) country_map[[k]][1], numeric(1))
  country_effects$rd_high <- vapply(as.character(country_effects$level), function(k) country_map[[k]][2], numeric(1))
} else {
  country_effects$rd_low <- NA_real_
  country_effects$rd_high <- NA_real_
}

sofa_plot <- ggplot(sofa_curve, aes(x = baseline_sofa, y = risk_difference)) +
  geom_ribbon(aes(ymin = rd_low, ymax = rd_high), fill = "#9ecae1", alpha = 0.35, na.rm = TRUE) +
  geom_line(color = "#203A57", linewidth = 1.1) +
  geom_hline(yintercept = 0, linetype = 2, color = "gray50") +
  labs(
    title = "Interaction Model HTE: baseline SOFA",
    subtitle = "Predicted AMR effect on 60-day mortality across baseline SOFA",
    x = "Baseline SOFA",
    y = "Predicted risk difference: A=1 minus A=0"
  ) +
  theme_minimal(base_size = 15)

bacteria_plot <- ggplot(bacteria_effects, aes(x = risk_difference, y = reorder(level, risk_difference))) +
  geom_vline(xintercept = 0, linetype = 2, color = "gray50") +
  geom_errorbarh(aes(xmin = rd_low, xmax = rd_high), height = 0.18, color = "#6c757d", na.rm = TRUE) +
  geom_point(size = 3, color = "#0E7C86") +
  labs(
    title = "Interaction Model HTE: baseline bacteria group",
    subtitle = "Predicted AMR effect on 60-day mortality by bacteria group",
    x = "Predicted risk difference: A=1 minus A=0",
    y = NULL
  ) +
  theme_minimal(base_size = 15)

country_plot <- ggplot(country_effects, aes(x = risk_difference, y = reorder(level, risk_difference))) +
  geom_vline(xintercept = 0, linetype = 2, color = "gray50") +
  geom_errorbarh(aes(xmin = rd_low, xmax = rd_high), height = 0.18, color = "#6c757d", na.rm = TRUE) +
  geom_point(size = 3, color = "#D98E04") +
  labs(
    title = "Interaction Model HTE: baseline country",
    subtitle = "Predicted AMR effect on 60-day mortality by country",
    x = "Predicted risk difference: A=1 minus A=0",
    y = NULL
  ) +
  theme_minimal(base_size = 15)

ggsave(file.path(out_dir, "figure_hte_interaction_baseline_sofa.png"), sofa_plot, width = 8.5, height = 5.5, dpi = 220)
ggsave(file.path(out_dir, "figure_hte_interaction_bacteria.png"), bacteria_plot, width = 8.5, height = 5.5, dpi = 220)
ggsave(file.path(out_dir, "figure_hte_interaction_country.png"), country_plot, width = 8.5, height = 4.8, dpi = 220)

utils::write.csv(model_df, file.path(out_dir, "analysis_dataset_for_interaction_hte.csv"), row.names = FALSE)
utils::write.csv(sofa_curve, file.path(out_dir, "table_hte_interaction_sofa_curve.csv"), row.names = FALSE)
utils::write.csv(bacteria_effects, file.path(out_dir, "table_hte_interaction_bacteria.csv"), row.names = FALSE)
utils::write.csv(country_effects, file.path(out_dir, "table_hte_interaction_country.csv"), row.names = FALSE)
utils::write.csv(interaction_test_table, file.path(out_dir, "table_hte_interaction_global_tests.csv"), row.names = FALSE)
utils::write.csv(balance_table, file.path(out_dir, "table_hte_interaction_balance_counts.csv"), row.names = FALSE)

model_notes <- c(
  "# Interaction-model HTE draft",
  "",
  "Date: 2026-09-03",
  "",
  "Population:",
  "- Baseline-defined culture-positive patients with non-missing A, Y, baseline SOFA, age, sex, Charlson, baseline appropriate treatment, bacteria group, and country.",
  "",
  "Outcome and exposure:",
  "- Y = 60-day mortality",
  "- A = baseline carbapenem resistance (1 resistant, 0 non-resistant)",
  "",
  "Interaction models:",
  "- SOFA model: Y ~ A * ns(baseline_sofa, 3) + age + sex + Charlson + baseline appropriate treatment + bacteria + country",
  "- Bacteria model: Y ~ A * baseline_bacteria_f + baseline_sofa + age + sex + Charlson + baseline appropriate treatment + country",
  "- Country model: Y ~ A * baseline_country_f + baseline_sofa + age + sex + Charlson + baseline appropriate treatment + bacteria",
  "",
  "Plot interpretation:",
  "- The y-axis is the model-based predicted risk difference for AMR (A=1 minus A=0).",
  "- For baseline SOFA, the line shows how the predicted AMR effect changes across the SOFA range.",
  "- For bacteria and country, the points show the predicted AMR effect within each category.",
  "",
  "Companion outputs:",
  "- table_hte_interaction_global_tests.csv gives the overall interaction p-value for each modifier.",
  "- table_hte_interaction_balance_counts.csv shows the A=0 and A=1 counts within each bacteria/country level to assess sparsity.",
  "",
  "Important note:",
  "- This is an interaction-model HTE analysis for the final binary outcome, not the longitudinal mediational g-formula decomposition of TE / IDE / IIE."
)
writeLines(model_notes, file.path(out_dir, "README_hte_interaction_model_v1.md"))

cat("Wrote interaction-model HTE outputs to:\n")
cat(out_dir, "\n")
