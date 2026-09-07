# Core longitudinal g-formula engine.
#
# This module is intentionally data-source agnostic: it expects the analysis
# dataset assembled by the master runner and uses the modelling helpers defined
# there. Keeping the engine here makes the estimand and simulation code easier
# to audit independently of data cleaning.

run_gformula <- function(analysis_df, mediator_history = c("lagged", "none"), lt_variant = c("main", "alt"), include_ot = FALSE, nsim = 4000, seed = 20260903) {
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
  y_formula <- stats::as.formula(paste("Y ~", rhs_join(c(y_terms, ot_terms))))

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
    for (t in 0:3) {
      m_dat <- sim
      m_dat$A <- a_for_m
      m_fit <- list(m0_fit, m1_fit, m2_fit, m3_fit)[[t + 1]]
      sim[[paste0("M", t)]] <- stats::rbinom(nrow(sim), 1, predict_binomial_prob(m_fit, m_dat))
      if (t < 3) sim[[paste0("L", t + 1)]] <- predict_gaussian_draw(list(l1_fit, l2_fit, l3_fit)[[t + 1]], sim)
    }
    y_dat <- sim
    y_dat$A <- a_for_y
    sim$Y <- stats::rbinom(nrow(sim), 1, predict_binomial_prob(y_fit, y_dat))
    sim
  }

  r11 <- estimate_risk(simulate_regime(1, 1))
  r10 <- estimate_risk(simulate_regime(1, 0))
  r00 <- estimate_risk(simulate_regime(0, 0))
  data.frame(lt_variant = lt_variant, mediator_history = mediator_history,
    include_ot = include_ot, nsim = nsim, n_complete = nrow(fit_df),
    active_covars = paste(covars, collapse = ";"), dropped_covars = paste(pruned$drop, collapse = ";"),
    R_1_G1 = r11, R_1_G0 = r10, R_0_G0 = r00, TE = r11 - r00,
    IDE = r10 - r00, IIE = r11 - r10, stringsAsFactors = FALSE)
}
