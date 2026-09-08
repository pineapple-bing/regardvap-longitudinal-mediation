# Prespecified robustness analyses for the longitudinal mediation model.
#
# A sensitivity analysis should target a specific source of uncertainty. This
# runner records the analysis family, estimand status, population restriction,
# covariate set, support diagnostics, and failures rather than returning only
# successful point estimates.

primary_covariates <- function() {
  # Site already determines country in this six-centre dataset. Including both
  # creates an exactly rank-deficient design matrix.
  c("age", "male", "charlson", "site", "icu_type", "bacteria")
}

country_only_covariates <- function() {
  c("age", "male", "charlson", "country", "icu_type", "bacteria")
}

no_centre_covariates <- function() {
  c("age", "male", "charlson", "icu_type", "bacteria")
}

make_sensitivity_specifications <- function() {
  list(
    list(
      id = "S1_alternative_severity", family = "measurement", same_estimand = TRUE,
      description = "Alternative daily severity definition", subset = "all",
      mediator_history = "lagged", mediator_severity_timing = "same_day", outcome_mediator_summary = "full_history",
      lt_variant = "alt", horizon = 3L, covariates = primary_covariates(),
      probability_bounds = c(0.001, 0.999)
    ),
    list(
      id = "S2_transition_without_lagged_M", family = "treatment_process_model", same_estimand = TRUE,
      description = "Mediator transition models omit prior-day treatment; outcome retains full M0-M3 history",
      subset = "all", mediator_history = "none", mediator_severity_timing = "same_day", outcome_mediator_summary = "full_history",
      lt_variant = "main", horizon = 3L, covariates = primary_covariates(),
      probability_bounds = c(0.001, 0.999)
    ),
    list(
      id = "S3_cumulative_treatment_outcome", family = "outcome_model", same_estimand = TRUE,
      description = "Outcome model summarises M0-M3 as cumulative appropriate-treatment days",
      subset = "all", mediator_history = "lagged", mediator_severity_timing = "same_day", outcome_mediator_summary = "cumulative_days",
      lt_variant = "main", horizon = 3L, covariates = primary_covariates(),
      probability_bounds = c(0.001, 0.999)
    ),
    list(
      id = "S4_ever_treated_outcome", family = "outcome_model", same_estimand = TRUE,
      description = "Outcome model summarises M0-M3 as any appropriate treatment",
      subset = "all", mediator_history = "lagged", mediator_severity_timing = "same_day", outcome_mediator_summary = "ever_treated",
      lt_variant = "main", horizon = 3L, covariates = primary_covariates(),
      probability_bounds = c(0.001, 0.999)
    ),
    list(
      id = "S5_early_treatment_window", family = "alternative_estimand", same_estimand = FALSE,
      description = "Alternative mediator window through Day 1 only", subset = "all",
      mediator_history = "lagged", mediator_severity_timing = "same_day", outcome_mediator_summary = "full_history",
      lt_variant = "main", horizon = 1L, covariates = primary_covariates(),
      probability_bounds = c(0.001, 0.999)
    ),
    list(
      id = "S6_culture_positive", family = "target_population", same_estimand = FALSE,
      description = "Restrict to culture-positive participants", subset = "culture_positive",
      mediator_history = "lagged", mediator_severity_timing = "same_day", outcome_mediator_summary = "full_history",
      lt_variant = "main", horizon = 3L, covariates = primary_covariates(),
      probability_bounds = c(0.001, 0.999)
    ),
    list(
      id = "S7_major_pathogens", family = "target_population", same_estimand = FALSE,
      description = "Restrict to Klebsiella, Acinetobacter, or Pseudomonas", subset = "major_pathogens",
      mediator_history = "lagged", mediator_severity_timing = "same_day", outcome_mediator_summary = "full_history",
      lt_variant = "main", horizon = 3L, covariates = primary_covariates(),
      probability_bounds = c(0.001, 0.999)
    ),
    list(
      id = "S8_country_without_site", family = "centre_structure", same_estimand = TRUE,
      description = "Adjustment includes country but omits nested site indicators", subset = "all",
      mediator_history = "lagged", mediator_severity_timing = "same_day", outcome_mediator_summary = "full_history",
      lt_variant = "main", horizon = 3L, covariates = country_only_covariates(),
      probability_bounds = c(0.001, 0.999)
    ),
    list(
      id = "S9_without_centre_indicators", family = "centre_structure", same_estimand = TRUE,
      description = "Adjustment omits both site and country indicators", subset = "all",
      mediator_history = "lagged", mediator_severity_timing = "same_day", outcome_mediator_summary = "full_history",
      lt_variant = "main", horizon = 3L, covariates = no_centre_covariates(),
      probability_bounds = c(0.001, 0.999)
    ),
    list(
      id = "S10_probability_bounds_01", family = "positivity_numerical_stability", same_estimand = TRUE,
      description = "Prediction probabilities truncated to 0.01-0.99", subset = "all",
      mediator_history = "lagged", mediator_severity_timing = "same_day", outcome_mediator_summary = "full_history",
      lt_variant = "main", horizon = 3L, covariates = primary_covariates(),
      probability_bounds = c(0.01, 0.99)
    ),
    list(
      id = "S11_probability_bounds_025", family = "positivity_numerical_stability", same_estimand = TRUE,
      description = "Prediction probabilities truncated to 0.025-0.975", subset = "all",
      mediator_history = "lagged", mediator_severity_timing = "same_day", outcome_mediator_summary = "full_history",
      lt_variant = "main", horizon = 3L, covariates = primary_covariates(),
      probability_bounds = c(0.025, 0.975)
    ),
    list(
      id = "S12_marginal_exposure_support", family = "exposure_positivity", same_estimand = FALSE,
      description = "Restrict to site and bacteria levels with at least five participants in each exposure group", subset = "marginal_exposure_support",
      mediator_history = "lagged", mediator_severity_timing = "same_day", outcome_mediator_summary = "full_history",
      lt_variant = "main", horizon = 3L, covariates = primary_covariates(),
      probability_bounds = c(0.001, 0.999)
    ),
    list(
      id = "S13_joint_exposure_support", family = "exposure_positivity", same_estimand = FALSE,
      description = "Restrict to site-by-bacteria cells with at least five participants in each exposure group", subset = "joint_exposure_support",
      mediator_history = "lagged", mediator_severity_timing = "same_day", outcome_mediator_summary = "full_history",
      lt_variant = "main", horizon = 3L, covariates = primary_covariates(),
      probability_bounds = c(0.001, 0.999)
    ),
    list(
      id = "S14_prior_day_severity_for_treatment", family = "temporal_ordering", same_estimand = TRUE,
      description = "Mediator models use prior-day rather than same-day severity", subset = "all",
      mediator_history = "lagged", mediator_severity_timing = "prior_day", outcome_mediator_summary = "full_history",
      lt_variant = "main", horizon = 3L, covariates = primary_covariates(),
      probability_bounds = c(0.001, 0.999)
    )
  )
}

apply_sensitivity_subset <- function(analysis_df, subset_name) {
  bacteria <- as.character(analysis_df$bacteria)
  site <- as.character(analysis_df$site)
  levels_with_exposure_support <- function(x, minimum_per_group = 5L) {
    tab <- table(x = x, A = factor(analysis_df$A, levels = c(0, 1)))
    rownames(tab)[tab[, "0"] >= minimum_per_group & tab[, "1"] >= minimum_per_group]
  }
  supported_sites <- levels_with_exposure_support(site)
  supported_bacteria <- levels_with_exposure_support(bacteria)
  joint <- interaction(site, bacteria, sep = " | ", drop = TRUE)
  supported_joint <- levels_with_exposure_support(joint)
  keep <- switch(
    subset_name,
    all = rep(TRUE, nrow(analysis_df)),
    culture_positive = !is.na(bacteria) & bacteria != "Culture_neg",
    major_pathogens = !is.na(bacteria) & bacteria %in% c("Klebsiella", "Acinetobacter", "Pseudomonas"),
    marginal_exposure_support = site %in% supported_sites & bacteria %in% supported_bacteria,
    joint_exposure_support = as.character(joint) %in% supported_joint,
    stop("Unknown sensitivity-analysis subset: ", subset_name, call. = FALSE)
  )
  droplevels(analysis_df[keep, , drop = FALSE])
}

sensitivity_support_row <- function(spec, dat) {
  a0 <- sum(dat$A == 0, na.rm = TRUE)
  a1 <- sum(dat$A == 1, na.rm = TRUE)
  data.frame(
    sensitivity_id = spec$id, analysis_family = spec$family, subset = spec$subset,
    n_input = nrow(dat), n_A0 = a0, n_A1 = a1,
    n_deaths = sum(dat$Y == 1, na.rm = TRUE),
    smallest_exposure_group = min(a0, a1), stringsAsFactors = FALSE
  )
}

failed_sensitivity_row <- function(spec, dat, message) {
  data.frame(
    sensitivity_id = spec$id, analysis_family = spec$family,
    same_estimand_as_primary = spec$same_estimand, description = spec$description,
    subset = spec$subset, status = "failed", error_message = message,
    lt_variant = spec$lt_variant, mediator_history = spec$mediator_history,
    mediator_severity_timing = spec$mediator_severity_timing,
    outcome_mediator_summary = spec$outcome_mediator_summary, include_ot = FALSE,
    horizon = spec$horizon, probability_lower = spec$probability_bounds[1],
    probability_upper = spec$probability_bounds[2], n_input = nrow(dat),
    nsim = NA_integer_, n_complete = NA_integer_,
    active_covars = paste(spec$covariates, collapse = ";"), dropped_covars = "",
    R_1_G1 = NA_real_, R_1_G0 = NA_real_, R_0_G0 = NA_real_,
    TE = NA_real_, IDE = NA_real_, IIE = NA_real_, stringsAsFactors = FALSE
  )
}

run_one_sensitivity <- function(analysis_df, spec, nsim, seed) {
  dat <- apply_sensitivity_subset(analysis_df, spec$subset)
  if (nrow(dat) < 100L || length(unique(dat$A[!is.na(dat$A)])) < 2L) {
    return(failed_sensitivity_row(spec, dat, "Insufficient sample size or exposure variation after restriction."))
  }
  fit <- tryCatch(
    run_gformula(
      dat, mediator_history = spec$mediator_history,
      mediator_severity_timing = spec$mediator_severity_timing,
      outcome_mediator_summary = spec$outcome_mediator_summary,
      lt_variant = spec$lt_variant, include_ot = FALSE, horizon = spec$horizon,
      covariates = spec$covariates, probability_bounds = spec$probability_bounds,
      nsim = nsim, seed = seed
    ),
    error = function(e) e
  )
  if (inherits(fit, "error")) {
    return(failed_sensitivity_row(spec, dat, conditionMessage(fit)))
  }
  cbind(
    sensitivity_id = spec$id, analysis_family = spec$family,
    same_estimand_as_primary = spec$same_estimand, description = spec$description,
    subset = spec$subset, status = "ok", error_message = "", fit,
    stringsAsFactors = FALSE
  )
}

run_mi_severity_sensitivity <- function(analysis_df, m = 10L, maxit = 10L, nsim = 10000L, seed = 20260903L) {
  spec <- list(
    id = "S15_multiple_imputation_severity", family = "missing_data", same_estimand = TRUE,
    description = paste0("Predictive mean matching for missing L1-L3 (m=", m, ")"),
    subset = "all", mediator_history = "lagged", mediator_severity_timing = "same_day", outcome_mediator_summary = "full_history",
    lt_variant = "main", horizon = 3L, covariates = primary_covariates(),
    probability_bounds = c(0.001, 0.999)
  )
  if (!requireNamespace("mice", quietly = TRUE)) {
    return(list(summary = failed_sensitivity_row(spec, analysis_df, "Package 'mice' is not installed."),
      imputations = data.frame(), logged_events = data.frame()))
  }

  mi_vars <- unique(c("A", "Y", paste0("M", 0:3), paste0("L", 0:3), primary_covariates()))
  mi_data <- analysis_df[, mi_vars, drop = FALSE]
  methods <- rep("", ncol(mi_data)); names(methods) <- names(mi_data)
  targets <- intersect(c("L1", "L2", "L3"), names(mi_data))
  methods[targets] <- "pmm"
  predictor_matrix <- mice::make.predictorMatrix(mi_data)
  predictor_matrix[,] <- 0
  # Site is nested within country in this dataset. Use site, not both site and
  # country, as an imputation predictor to avoid deterministic collinearity.
  mi_predictors <- unique(c(
    "A", "Y", paste0("M", 0:3), paste0("L", 0:3),
    "age", "male", "charlson", "site", "icu_type", "bacteria"
  ))
  for (target in targets) predictor_matrix[target, setdiff(mi_predictors, target)] <- 1

  imp <- tryCatch(
    suppressWarnings(mice::mice(
      mi_data, m = as.integer(m), maxit = as.integer(maxit), method = methods,
      predictorMatrix = predictor_matrix, printFlag = FALSE, seed = seed
    )),
    error = function(e) e
  )
  if (inherits(imp, "error")) {
    return(list(summary = failed_sensitivity_row(spec, analysis_df,
      paste("Multiple imputation failed:", conditionMessage(imp))),
      imputations = data.frame(), logged_events = data.frame()))
  }

  fits <- lapply(seq_len(m), function(i) {
    completed <- mice::complete(imp, action = i)
    fit <- tryCatch(
      run_gformula(
        completed, mediator_history = "lagged", mediator_severity_timing = "same_day", outcome_mediator_summary = "full_history",
        lt_variant = "main", include_ot = FALSE, horizon = 3L,
        covariates = primary_covariates(), probability_bounds = c(0.001, 0.999),
        nsim = nsim, seed = seed + i
      ), error = function(e) e
    )
    if (inherits(fit, "error")) {
      data.frame(imputation = i, status = "failed", error_message = conditionMessage(fit),
        R_1_G1 = NA_real_, R_1_G0 = NA_real_, R_0_G0 = NA_real_, TE = NA_real_, IDE = NA_real_, IIE = NA_real_)
    } else {
      data.frame(imputation = i, status = "ok", error_message = "",
        R_1_G1 = fit$R_1_G1, R_1_G0 = fit$R_1_G0, R_0_G0 = fit$R_0_G0,
        TE = fit$TE, IDE = fit$IDE, IIE = fit$IIE)
    }
  })
  fits <- do.call(rbind, fits)
  ok <- fits[fits$status == "ok", , drop = FALSE]
  if (!nrow(ok)) {
    return(list(summary = failed_sensitivity_row(spec, analysis_df,
      "Every imputed-data g-formula fit failed."), imputations = fits,
      logged_events = if (is.null(imp$loggedEvents)) data.frame() else imp$loggedEvents))
  }

  estimates <- c("R_1_G1", "R_1_G0", "R_0_G0", "TE", "IDE", "IIE")
  averaged <- vapply(estimates, function(v) mean(ok[[v]], na.rm = TRUE), numeric(1))
  summary <- data.frame(
    sensitivity_id = spec$id, analysis_family = spec$family,
    same_estimand_as_primary = spec$same_estimand, description = spec$description,
    subset = spec$subset, status = if (nrow(ok) == m) "ok" else "partial_success",
    error_message = if (nrow(ok) == m) "" else paste(m - nrow(ok), "imputed fits failed"),
    lt_variant = "main", mediator_history = "lagged", mediator_severity_timing = "same_day",
    outcome_mediator_summary = "full_history", include_ot = FALSE, horizon = 3L,
    probability_lower = 0.001, probability_upper = 0.999,
    n_input = nrow(analysis_df), nsim = nsim, n_complete = nrow(analysis_df),
    active_covars = paste(primary_covariates(), collapse = ";"), dropped_covars = "",
    R_1_G1 = averaged["R_1_G1"], R_1_G0 = averaged["R_1_G0"],
    R_0_G0 = averaged["R_0_G0"], TE = averaged["TE"],
    IDE = averaged["IDE"], IIE = averaged["IIE"], stringsAsFactors = FALSE
  )
  list(
    summary = summary,
    imputations = fits,
    logged_events = if (is.null(imp$loggedEvents)) data.frame() else imp$loggedEvents
  )
}

run_prespecified_sensitivities <- function(
    analysis_df, nsim = 4000, seed = 20260903,
    run_multiple_imputation = TRUE, mi_m = 10L, mi_maxit = 10L,
    mi_nsim = min(nsim, 10000L)) {
  specifications <- make_sensitivity_specifications()
  rows <- lapply(seq_along(specifications), function(i) {
    run_one_sensitivity(analysis_df, specifications[[i]], nsim = nsim, seed = seed + i)
  })
  support <- do.call(rbind, lapply(specifications, function(spec) {
    sensitivity_support_row(spec, apply_sensitivity_subset(analysis_df, spec$subset))
  }))

  mi_details <- data.frame()
  mi_logged_events <- data.frame()
  if (isTRUE(run_multiple_imputation)) {
    mi <- run_mi_severity_sensitivity(analysis_df, m = mi_m, maxit = mi_maxit,
      nsim = mi_nsim, seed = seed + 1000L)
    rows[[length(rows) + 1L]] <- mi$summary
    mi_details <- mi$imputations
    mi_logged_events <- mi$logged_events
    support <- rbind(support, sensitivity_support_row(
      list(id = "S15_multiple_imputation_severity", family = "missing_data", subset = "all"),
      analysis_df
    ))
  }

  result <- do.call(rbind, rows)
  rownames(result) <- NULL
  attr(result, "support_diagnostics") <- support
  attr(result, "mi_imputation_estimates") <- mi_details
  attr(result, "mi_logged_events") <- mi_logged_events
  result
}
