# Monte Carlo diagnostics for the primary g-formula specification.

run_monte_carlo_stability <- function(analysis_df, seeds = c(20260911L, 20260912L, 20260913L), nsim = 4000L) {
  rows <- lapply(seeds, function(seed) {
    fit <- run_gformula(analysis_df, mediator_history = "lagged", lt_variant = "main",
      include_ot = FALSE, horizon = 3L, nsim = nsim, seed = seed)
    cbind(seed = seed, fit, stringsAsFactors = FALSE)
  })
  do.call(rbind, rows)
}
