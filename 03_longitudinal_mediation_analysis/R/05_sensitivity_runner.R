# Prespecified alternatives to the primary longitudinal mediation specification.

run_prespecified_sensitivities <- function(analysis_df, nsim = 4000, seed = 20260903) {
  specifications <- list(
    list(id = "S1_alternative_severity", description = "Alternative daily severity definition", mediator_history = "lagged", lt_variant = "alt", include_ot = FALSE, horizon = 3L),
    list(id = "S2_no_lagged_mediator_history", description = "Mediator model without lagged mediator history", mediator_history = "none", lt_variant = "main", include_ot = FALSE, horizon = 3L),
    list(id = "S3_early_treatment_window", description = "Early treatment window through Day 1 only", mediator_history = "lagged", lt_variant = "main", include_ot = FALSE, horizon = 1L)
  )
  rows <- lapply(seq_along(specifications), function(i) {
    spec <- specifications[[i]]
    fit <- run_gformula(analysis_df, mediator_history = spec$mediator_history,
      lt_variant = spec$lt_variant, include_ot = spec$include_ot,
      horizon = spec$horizon, nsim = nsim, seed = seed + i)
    cbind(sensitivity_id = spec$id, description = spec$description, fit,
      stringsAsFactors = FALSE)
  })
  do.call(rbind, rows)
}
