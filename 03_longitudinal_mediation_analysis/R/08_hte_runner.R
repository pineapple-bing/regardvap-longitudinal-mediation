# Exploratory heterogeneity-of-treatment-effect (HTE) analysis.
#
# This module deliberately does not label candidate modifiers as confirmatory.
# Validity of the AMR exposure and the primary target population must be
# confirmed before any HTE result can be interpreted as a clinical finding.

make_hte_specifications <- function() {
  list(
    list(id = "icu_type", modifier_var = "icu_type", modifier_label = "ICU type",
         levels_keep = c("Medical ICU", "Surgical ICU")),
    list(id = "severity_group", modifier_var = "severity_group",
         modifier_label = "Baseline severity (SOFA>=6 vs <6)",
         levels_keep = c("SOFA<6", "SOFA>=6"))
  )
}

hte_fit_variables <- function(covariates = primary_covariates()) {
  unique(c("A", "Y", paste0("M", 0:3), paste0("L", 0:3), covariates))
}

hte_gformula_arguments <- function(nsim, seed) {
  list(
    mediator_history = "lagged", mediator_severity_timing = "same_day",
    outcome_mediator_summary = "full_history", lt_variant = "main",
    include_ot = FALSE, horizon = 3L, covariates = primary_covariates(),
    probability_bounds = c(0.001, 0.999), nsim = nsim, seed = seed
  )
}

# Exposure support is reported separately in each modifier level. The screen is
# descriptive: having both observed A levels is necessary, but not sufficient,
# for causal identification.
hte_exposure_support <- function(analysis_df, specifications,
                                 covariates = primary_covariates(), seed = 20260911L) {
  rows <- list()
  cell_index <- 0L
  for (spec in specifications) {
    for (lv in spec$levels_keep) {
      cell_index <- cell_index + 1L
      dat <- analysis_df[as.character(analysis_df[[spec$modifier_var]]) == lv, , drop = FALSE]
      vars <- hte_fit_variables(covariates)
      complete <- dat[stats::complete.cases(dat[, vars, drop = FALSE]), , drop = FALSE]
      base <- data.frame(
        modifier = spec$modifier_label, level = lv, n_input = nrow(dat),
        n_complete = nrow(complete), n_complete_A0 = sum(complete$A == 0),
        n_complete_A1 = sum(complete$A == 1), stringsAsFactors = FALSE
      )
      if (nrow(complete) == 0L || length(unique(complete$A)) < 2L) {
        rows[[length(rows) + 1L]] <- cbind(base, support_status = "failed_no_complete_exposure_overlap",
          propensity_min = NA_real_, propensity_p05 = NA_real_, propensity_median = NA_real_,
          propensity_p95 = NA_real_, propensity_max = NA_real_, n_outside_05_95 = NA_integer_,
          exposure_model_type = NA_character_)
        next
      }
      x <- complete[, unique(c("A", covariates)), drop = FALSE]
      pruned <- drop_invariant_columns(x, covariates)
      active <- intersect(covariates, pruned$keep)
      set.seed(seed + cell_index)
      fit <- tryCatch(fit_binomial_model(stats::as.formula(paste("A ~", rhs_join(active))), pruned$data), error = function(e) e)
      if (inherits(fit, "error")) {
        rows[[length(rows) + 1L]] <- cbind(base, support_status = paste0("failed_propensity_model: ", fit$message),
          propensity_min = NA_real_, propensity_p05 = NA_real_, propensity_median = NA_real_,
          propensity_p95 = NA_real_, propensity_max = NA_real_, n_outside_05_95 = NA_integer_,
          exposure_model_type = NA_character_)
        next
      }
      p <- predict_binomial_prob(fit, pruned$data, lower = 0, upper = 1)
      fit_type <- if (inherits(fit, "ridge_binomial_model")) fit$fit_type else attr(fit, "fit_type")
      rows[[length(rows) + 1L]] <- cbind(base, support_status = "review_required",
        propensity_min = min(p), propensity_p05 = unname(stats::quantile(p, .05)),
        propensity_median = stats::median(p), propensity_p95 = unname(stats::quantile(p, .95)),
        propensity_max = max(p), n_outside_05_95 = sum(p < .05 | p > .95),
        exposure_model_type = fit_type)
    }
  }
  do.call(rbind, rows)
}

run_hte_gformula <- function(analysis_df, modifier_var, modifier_label, levels_keep,
                             min_complete_per_exposure = 20L, nsim = 4000L,
                             seed = 20260903L) {
  rows <- list()
  vars <- hte_fit_variables()
  for (i in seq_along(levels_keep)) {
    lv <- levels_keep[[i]]
    sub <- analysis_df[as.character(analysis_df[[modifier_var]]) == lv, , drop = FALSE]
    complete <- sub[stats::complete.cases(sub[, vars, drop = FALSE]), , drop = FALSE]
    n_a0 <- sum(complete$A == 0); n_a1 <- sum(complete$A == 1)
    base <- list(modifier = modifier_label, level = lv, n_input = nrow(sub), n_complete = nrow(complete),
                 n_complete_A0 = n_a0, n_complete_A1 = n_a1)
    blank <- list(R_1_G1 = NA_real_, R_1_G0 = NA_real_, R_0_G0 = NA_real_,
                  TE = NA_real_, IDE = NA_real_, IIE = NA_real_)
    if (min(n_a0, n_a1) < min_complete_per_exposure) {
      rows[[i]] <- data.frame(c(base, list(status = "skipped_low_complete_exposure_cell"), blank), stringsAsFactors = FALSE)
      next
    }
    fit <- tryCatch(do.call(run_gformula, c(list(analysis_df = sub), hte_gformula_arguments(nsim, seed + i))), error = function(e) e)
    if (inherits(fit, "error")) {
      rows[[i]] <- data.frame(c(base, list(status = paste0("failed: ", fit$message)), blank), stringsAsFactors = FALSE)
    } else {
      rows[[i]] <- data.frame(c(base, list(status = "ok"), as.list(fit[1, c("R_1_G1", "R_1_G0", "R_0_G0", "TE", "IDE", "IIE")])), stringsAsFactors = FALSE)
    }
  }
  do.call(rbind, rows)
}

run_all_hte <- function(analysis_df, specifications, nsim = 4000L, seed = 20260903L) {
  do.call(rbind, lapply(specifications, function(spec) run_hte_gformula(
    analysis_df, spec$modifier_var, spec$modifier_label, spec$levels_keep, nsim = nsim, seed = seed
  )))
}

# Uses the same participant bootstrap sample for both levels, then contrasts
# their estimates. This is the appropriate resampling unit for an interaction
# contrast; comparing overlap of two separate subgroup CIs is not a test.
bootstrap_hte_difference <- function(analysis_df, specifications, B = 500L,
                                     nsim = 10000L, seed = 20260910L,
                                     checkpoint_path = NULL, progress_every = 25L) {
  set.seed(seed)
  indices <- lapply(seq_len(B), function(x) sample.int(nrow(analysis_df), nrow(analysis_df), replace = TRUE))
  required <- c("bootstrap_id", "status", "message")
  for (spec in specifications) for (tag in c("level1", "level2")) {
    required <- c(required, paste0(spec$id, "_", tag, "_", c("TE", "IDE", "IIE")))
  }
  rows <- vector("list", B)
  if (!is.null(checkpoint_path) && file.exists(checkpoint_path)) {
    prior <- utils::read.csv(checkpoint_path, stringsAsFactors = FALSE)
    if (!all(required %in% names(prior))) stop("Invalid HTE-bootstrap checkpoint.", call. = FALSE)
    for (b in unique(prior$bootstrap_id[prior$bootstrap_id <= B])) rows[[b]] <- prior[prior$bootstrap_id == b, required, drop = FALSE][1, ]
  }
  for (b in seq_len(B)) {
    if (!is.null(rows[[b]])) next
    sampled <- analysis_df[indices[[b]], , drop = FALSE]
    values <- list(bootstrap_id = b, message = "")
    messages <- character(0); ok <- TRUE
    for (spec in specifications) for (j in seq_along(spec$levels_keep)) {
      tag <- c("level1", "level2")[[j]]
      sub <- sampled[as.character(sampled[[spec$modifier_var]]) == spec$levels_keep[[j]], , drop = FALSE]
      fit <- tryCatch(do.call(run_gformula, c(list(analysis_df = sub), hte_gformula_arguments(nsim, seed + b * 100L + j))), error = function(e) e)
      if (inherits(fit, "error")) { ok <- FALSE; messages <- c(messages, paste(spec$id, tag, fit$message, sep = ":")) }
      for (est in c("TE", "IDE", "IIE")) values[[paste0(spec$id, "_", tag, "_", est)]] <- if (inherits(fit, "error")) NA_real_ else fit[[est]]
    }
    values$status <- if (ok) "ok" else "failed"; values$message <- paste(messages, collapse = " | ")
    rows[[b]] <- as.data.frame(values, stringsAsFactors = FALSE)
    if (!is.null(checkpoint_path) && (b %% progress_every == 0L || b == B)) utils::write.csv(do.call(rbind, rows[!vapply(rows, is.null, logical(1))]), checkpoint_path, row.names = FALSE)
  }
  do.call(rbind, rows)
}

hte_difference_percentile_ci <- function(draws, specifications) {
  ok <- draws[draws$status == "ok", , drop = FALSE]
  ci <- function(x) if (sum(is.finite(x))) stats::quantile(x[is.finite(x)], c(.025, .975), names = FALSE) else c(NA_real_, NA_real_)
  do.call(rbind, lapply(specifications, function(spec) do.call(rbind, lapply(c("TE", "IDE", "IIE"), function(est) {
    x <- ok[[paste0(spec$id, "_level1_", est)]] - ok[[paste0(spec$id, "_level2_", est)]]
    z <- ci(x)
    data.frame(modifier = spec$modifier_label, estimand = est, level1 = spec$levels_keep[[1]], level2 = spec$levels_keep[[2]],
      delta_bootstrap_mean = mean(x, na.rm = TRUE), delta_lower_95 = z[1], delta_upper_95 = z[2],
      n_success = sum(is.finite(x)), n_requested = nrow(draws), stringsAsFactors = FALSE)
  }))))
}
