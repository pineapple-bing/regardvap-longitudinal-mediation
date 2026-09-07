# Prespecified alternatives to the primary longitudinal mediation specification.

run_prespecified_sensitivities <- function(analysis_df, nsim = 4000, seed = 20260903) {
  do.call(rbind, list(
    run_gformula(analysis_df, mediator_history = "lagged", lt_variant = "alt", include_ot = FALSE, nsim = nsim, seed = seed + 1L),
    run_gformula(analysis_df, mediator_history = "none", lt_variant = "main", include_ot = FALSE, nsim = nsim, seed = seed + 2L),
    run_gformula(analysis_df, mediator_history = "lagged", lt_variant = "main", include_ot = TRUE, nsim = nsim, seed = seed + 3L)
  ))
}
