# Core longitudinal g-formula engine.
#
# This module is intentionally data-source agnostic: it expects the analysis
# dataset assembled by the master runner and uses the modelling helpers defined
# there. Keeping the engine here makes the estimand and simulation code easier
# to audit independently of data cleaning.

run_gformula <- function(analysis_df, mediator_history = c("lagged", "none"), lt_variant = c("main", "alt"), include_ot = FALSE, horizon = 3L, nsim = 4000, seed = 20260903) {
  mediator_history <- match.arg(mediator_history)
  lt_variant <- match.arg(lt_variant)
  horizon <- as.integer(horizon)
  if (is.na(horizon) || !(horizon %in% 1:3)) {
    stop("horizon must be 1, 2, or 3 (the final mediator day).", call. = FALSE)
  }
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
  days <- 0:horizon
  mediator_terms <- paste0("M", days)
  severity_terms <- paste0("L", days)
  ot_terms <- if (include_ot) paste0("O", days) else character(0)
  fit_vars <- unique(c("A", "Y", mediator_terms, severity_terms, covars, ot_terms))
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

  y_terms <- c("A", severity_terms, covars)
  if (mediator_history == "lagged") {
    y_terms <- c(y_terms, mediator_terms)
  } else {
    y_terms <- c(y_terms, paste0("M", horizon))
  }
  y_formula <- stats::as.formula(paste("Y ~", rhs_join(c(y_terms, ot_terms))))

  m_fits <- lapply(days, function(t) fit_binomial_model(mediator_formula(t), fit_df))
  names(m_fits) <- as.character(days)
  l_fits <- lapply(seq_len(horizon), function(t) fit_gaussian_model(l_formula(t), fit_df))
  names(l_fits) <- as.character(seq_len(horizon))
  y_fit <- fit_binomial_model(y_formula, fit_df)

  fit_type <- function(x) {
    if (inherits(x, "ridge_binomial_model")) x$fit_type else attr(x, "fit_type")
  }
  lambda_1se <- function(x) {
    if (inherits(x, "ridge_binomial_model")) x$lambda_1se else NA_real_
  }
  binomial_models <- c(m_fits, list(Y = y_fit))
  nuisance_model_diagnostics <- data.frame(
    model = names(binomial_models),
    model_type = vapply(binomial_models, fit_type, character(1)),
    lambda_1se = vapply(binomial_models, lambda_1se, numeric(1)),
    n_fit = nrow(fit_df),
    stringsAsFactors = FALSE
  )
  conditional_positivity <- do.call(rbind, lapply(days, function(t) {
    p <- predict_binomial_prob(m_fits[[as.character(t)]], fit_df)
    do.call(rbind, lapply(sort(unique(fit_df$A)), function(a) {
      pa <- p[fit_df$A == a]
      data.frame(day = t, A = a, n = length(pa), min_probability = min(pa), max_probability = max(pa),
        n_below_0_01 = sum(pa < 0.01), n_above_0_99 = sum(pa > 0.99), stringsAsFactors = FALSE)
    }))
  }))

  base_cols <- unique(c(covars, "L0", ot_terms))
  sim_base <- fit_df[rep(seq_len(nrow(fit_df)), length.out = nsim), base_cols, drop = FALSE]

  prediction_log <- list()
  record_prediction <- function(x, model, regime, diagnostic_type) {
    d <- attr(x, "prediction_diagnostics")
    if (is.null(d)) return(invisible(NULL))
    prediction_log[[length(prediction_log) + 1L]] <<- data.frame(
      regime = regime, model = model, diagnostic_type = diagnostic_type,
      n_na_prediction = unname(d["n_na_prediction"]),
      n_truncated = unname(d[grep("truncated", names(d))[1]]),
      stringsAsFactors = FALSE
    )
    invisible(NULL)
  }

  simulate_regime <- function(a_for_y, a_for_m) {
    regime <- paste0("R(", a_for_y, ",G", a_for_m, ")")
    sim <- sim_base
    sim$A <- a_for_y
    for (t in days) {
      m_dat <- sim
      m_dat$A <- a_for_m
      m_fit <- m_fits[[as.character(t)]]
      p_m <- predict_binomial_prob(m_fit, m_dat)
      record_prediction(p_m, paste0("M", t), regime, "probability")
      sim[[paste0("M", t)]] <- stats::rbinom(nrow(sim), 1, p_m)
      if (t < horizon) {
        l_draw <- predict_gaussian_draw(l_fits[[as.character(t + 1)]], sim, lower = 0, upper = 11)
        record_prediction(l_draw, paste0("L", t + 1), regime, "severity_support")
        sim[[paste0("L", t + 1)]] <- l_draw
      }
    }
    y_dat <- sim
    y_dat$A <- a_for_y
    p_y <- predict_binomial_prob(y_fit, y_dat)
    record_prediction(p_y, "Y", regime, "probability")
    sim$Y <- stats::rbinom(nrow(sim), 1, p_y)
    sim
  }

  r11 <- estimate_risk(simulate_regime(1, 1))
  r10 <- estimate_risk(simulate_regime(1, 0))
  r00 <- estimate_risk(simulate_regime(0, 0))
  result <- data.frame(lt_variant = lt_variant, mediator_history = mediator_history,
    include_ot = include_ot, horizon = horizon, nsim = nsim, n_complete = nrow(fit_df),
    active_covars = paste(covars, collapse = ";"), dropped_covars = paste(pruned$drop, collapse = ";"),
    R_1_G1 = r11, R_1_G0 = r10, R_0_G0 = r00, TE = r11 - r00,
    IDE = r10 - r00, IIE = r11 - r10, stringsAsFactors = FALSE)
  attr(result, "nuisance_model_diagnostics") <- nuisance_model_diagnostics
  attr(result, "conditional_positivity") <- conditional_positivity
  attr(result, "prediction_diagnostics") <- do.call(rbind, prediction_log)
  result
}
