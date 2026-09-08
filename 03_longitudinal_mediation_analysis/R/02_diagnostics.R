# Diagnostics separated from model fitting so they are visible and auditable.

write_positivity_diagnostics <- function(
    analysis_df, output_dir,
    exposure_covariates = c("age", "male", "charlson", "site", "icu_type", "bacteria")) {
  required <- c("A", "M0", "M1", "M2", "M3")
  if (!all(required %in% names(analysis_df))) {
    stop("Cannot calculate positivity diagnostics: missing A or daily mediator variables.", call. = FALSE)
  }

  day_rows <- lapply(0:3, function(day) {
    mediator <- paste0("M", day)
    tab <- as.data.frame(table(A = analysis_df$A, M = analysis_df[[mediator]], useNA = "no"))
    names(tab)[names(tab) == "Freq"] <- "n"
    tab$day <- day
    tab$proportion_within_A <- ave(tab$n, tab$A, FUN = function(x) x / sum(x))
    tab
  })

  complete_history <- stats::complete.cases(analysis_df[c("A", "M0", "M1", "M2", "M3")])
  history <- analysis_df[complete_history, c("A", "M0", "M1", "M2", "M3"), drop = FALSE]
  history$history <- paste0(history$M0, history$M1, history$M2, history$M3)
  history_tab <- as.data.frame(table(A = history$A, history = history$history, useNA = "no"))
  names(history_tab)[names(history_tab) == "Freq"] <- "n"
  history_tab$proportion_within_A <- ave(history_tab$n, history_tab$A, FUN = function(x) x / sum(x))

  exposure_support <- function(stratum, label) {
    stratum <- as.character(stratum)
    stratum[is.na(stratum) | !nzchar(stratum)] <- "Missing"
    tab <- table(stratum = stratum, A = factor(analysis_df$A, levels = c(0, 1)))
    out <- data.frame(
      variable = label,
      stratum = rownames(tab),
      n_A0 = as.integer(tab[, "0"]),
      n_A1 = as.integer(tab[, "1"]),
      stringsAsFactors = FALSE
    )
    out$n_total <- out$n_A0 + out$n_A1
    out$both_exposure_levels_observed <- out$n_A0 > 0 & out$n_A1 > 0
    out$minimum_exposure_cell <- pmin(out$n_A0, out$n_A1)
    out
  }

  support_vars <- intersect(c("country", "site", "icu_type", "bacteria"), names(analysis_df))
  exposure_support_marginal <- do.call(rbind, lapply(support_vars, function(v) {
    exposure_support(analysis_df[[v]], v)
  }))
  if (all(c("site", "bacteria") %in% names(analysis_df))) {
    joint_stratum <- interaction(
      as.character(analysis_df$site), as.character(analysis_df$bacteria),
      sep = " | ", drop = TRUE
    )
    exposure_support_joint <- exposure_support(joint_stratum, "site_x_bacteria")
  } else {
    exposure_support_joint <- data.frame()
  }

  missing_exposure_covariates <- setdiff(exposure_covariates, names(analysis_df))
  if (length(missing_exposure_covariates)) {
    stop(
      "Cannot calculate exposure propensity: missing covariates ",
      paste(missing_exposure_covariates, collapse = ", "), call. = FALSE
    )
  }
  exposure_df <- analysis_df[, unique(c("subjid", "A", exposure_covariates)), drop = FALSE]
  exposure_df <- exposure_df[stats::complete.cases(exposure_df), , drop = FALSE]
  pruned <- drop_invariant_columns(exposure_df, exposure_covariates)
  exposure_df <- pruned$data
  exposure_covariates_used <- intersect(exposure_covariates, pruned$keep)
  exposure_formula <- stats::as.formula(paste("A ~", rhs_join(exposure_covariates_used)))
  exposure_fit <- fit_binomial_model(exposure_formula, exposure_df)
  exposure_df$propensity_A1 <- predict_binomial_prob(exposure_fit, exposure_df, lower = 0, upper = 1)
  exposure_df$outside_01_99 <- exposure_df$propensity_A1 < 0.01 | exposure_df$propensity_A1 > 0.99
  exposure_df$outside_025_975 <- exposure_df$propensity_A1 < 0.025 | exposure_df$propensity_A1 > 0.975
  exposure_df$outside_05_95 <- exposure_df$propensity_A1 < 0.05 | exposure_df$propensity_A1 > 0.95

  propensity_summary <- do.call(rbind, lapply(c("overall", "A0", "A1"), function(group) {
    use <- if (group == "overall") exposure_df else exposure_df[exposure_df$A == as.numeric(sub("A", "", group)), , drop = FALSE]
    p <- use$propensity_A1
    data.frame(
      group = group, n = length(p), min = min(p), p01 = stats::quantile(p, 0.01, names = FALSE),
      p05 = stats::quantile(p, 0.05, names = FALSE), median = stats::median(p),
      p95 = stats::quantile(p, 0.95, names = FALSE), p99 = stats::quantile(p, 0.99, names = FALSE),
      max = max(p), n_outside_01_99 = sum(use$outside_01_99),
      n_outside_025_975 = sum(use$outside_025_975), n_outside_05_95 = sum(use$outside_05_95),
      stringsAsFactors = FALSE
    )
  }))
  exposure_fit_type <- if (inherits(exposure_fit, "ridge_binomial_model")) {
    exposure_fit$fit_type
  } else {
    attr(exposure_fit, "fit_type")
  }
  exposure_model <- data.frame(
    formula = paste(deparse(exposure_formula), collapse = ""),
    model_type = exposure_fit_type,
    n_fit = nrow(exposure_df),
    active_covariates = paste(exposure_covariates_used, collapse = ";"),
    dropped_covariates = paste(pruned$drop, collapse = ";"),
    stringsAsFactors = FALSE
  )

  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  utils::write.csv(do.call(rbind, day_rows), file.path(output_dir, "diagnostic_positivity_by_day.csv"), row.names = FALSE)
  utils::write.csv(history_tab, file.path(output_dir, "diagnostic_positivity_by_history.csv"), row.names = FALSE)
  utils::write.csv(exposure_support_marginal, file.path(output_dir, "diagnostic_exposure_support_by_stratum.csv"), row.names = FALSE)
  utils::write.csv(exposure_support_joint, file.path(output_dir, "diagnostic_exposure_support_site_by_bacteria.csv"), row.names = FALSE)
  utils::write.csv(exposure_df[, c("subjid", "A", "propensity_A1", "outside_01_99", "outside_025_975", "outside_05_95")],
    file.path(output_dir, "diagnostic_exposure_propensity_individual.csv"), row.names = FALSE)
  utils::write.csv(propensity_summary, file.path(output_dir, "diagnostic_exposure_propensity_summary.csv"), row.names = FALSE)
  utils::write.csv(exposure_model, file.path(output_dir, "diagnostic_exposure_propensity_model.csv"), row.names = FALSE)
  invisible(list(
    by_day = do.call(rbind, day_rows), by_history = history_tab,
    exposure_support_marginal = exposure_support_marginal,
    exposure_support_joint = exposure_support_joint,
    exposure_propensity = exposure_df,
    exposure_propensity_summary = propensity_summary,
    exposure_model = exposure_model
  ))
}
