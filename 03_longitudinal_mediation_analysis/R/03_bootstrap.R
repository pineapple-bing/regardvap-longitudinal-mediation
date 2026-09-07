# Patient-level bootstrap for uncertainty around longitudinal g-formula effects.

bootstrap_gformula <- function(analysis_df, B = 200L, nsim = 4000L, seed = 20260907L) {
  if (!exists("run_gformula", mode = "function", inherits = TRUE)) {
    stop("Define run_gformula before calling bootstrap_gformula.", call. = FALSE)
  }
  set.seed(seed)
  n <- nrow(analysis_df)
  if (n < 2) stop("Bootstrap requires at least two participants.", call. = FALSE)

  rows <- vector("list", B)
  for (b in seq_len(B)) {
    sampled <- analysis_df[sample.int(n, size = n, replace = TRUE), , drop = FALSE]
    fit <- tryCatch(
      run_gformula(sampled, mediator_history = "lagged", lt_variant = "main", include_ot = FALSE,
        nsim = nsim, seed = seed + b),
      error = function(e) e
    )
    if (inherits(fit, "error")) {
      rows[[b]] <- data.frame(bootstrap_id = b, status = "failed", message = fit$message,
        TE = NA_real_, IDE = NA_real_, IIE = NA_real_)
    } else {
      rows[[b]] <- data.frame(bootstrap_id = b, status = "ok", message = "",
        TE = fit$TE, IDE = fit$IDE, IIE = fit$IIE)
    }
  }
  do.call(rbind, rows)
}

bootstrap_percentile_ci <- function(bootstrap_df) {
  effects <- c("TE", "IDE", "IIE")
  do.call(rbind, lapply(effects, function(effect) {
    x <- bootstrap_df[[effect]][bootstrap_df$status == "ok"]
    data.frame(effect = effect, n_success = sum(!is.na(x)),
      lower_95 = stats::quantile(x, 0.025, na.rm = TRUE, names = FALSE),
      upper_95 = stats::quantile(x, 0.975, na.rm = TRUE, names = FALSE))
  }))
}
