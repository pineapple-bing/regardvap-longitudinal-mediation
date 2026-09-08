# Patient-level bootstrap for uncertainty around longitudinal g-formula effects.

bootstrap_gformula <- function(
    analysis_df, B = 500L, nsim = 10000L, seed = 20260907L,
    checkpoint_path = NULL, progress_every = 25L) {
  if (!exists("run_gformula", mode = "function", inherits = TRUE)) {
    stop("Define run_gformula before calling bootstrap_gformula.", call. = FALSE)
  }
  set.seed(seed)
  n <- nrow(analysis_df)
  if (n < 2) stop("Bootstrap requires at least two participants.", call. = FALSE)

  # Generate all resample indices before any g-formula simulation changes the
  # global RNG state. This keeps the patient resamples identical when nsim is
  # changed, allowing a clean Monte Carlo convergence comparison.
  bootstrap_indices <- lapply(seq_len(B), function(b) {
    sample.int(n, size = n, replace = TRUE)
  })
  rows <- vector("list", B)
  if (!is.null(checkpoint_path) && file.exists(checkpoint_path)) {
    prior <- utils::read.csv(checkpoint_path, stringsAsFactors = FALSE)
    required <- c("bootstrap_id", "status", "message", "R_1_G1", "R_1_G0", "R_0_G0", "TE", "IDE", "IIE")
    if (!all(required %in% names(prior))) stop("Invalid primary-bootstrap checkpoint.", call. = FALSE)
    prior <- prior[prior$bootstrap_id %in% seq_len(B), required, drop = FALSE]
    for (b in unique(prior$bootstrap_id)) rows[[b]] <- prior[prior$bootstrap_id == b, , drop = FALSE][1, ]
    message("Resuming primary bootstrap with ", sum(!vapply(rows, is.null, logical(1))), " of ", B, " resamples complete.")
  }
  for (b in seq_len(B)) {
    if (!is.null(rows[[b]])) next
    sampled <- analysis_df[bootstrap_indices[[b]], , drop = FALSE]
    fit <- tryCatch(
      run_gformula(sampled, mediator_history = "lagged", lt_variant = "main", include_ot = FALSE,
        nsim = nsim, seed = seed + b),
      error = function(e) e
    )
    if (inherits(fit, "error")) {
      rows[[b]] <- data.frame(bootstrap_id = b, status = "failed", message = fit$message,
        R_1_G1 = NA_real_, R_1_G0 = NA_real_, R_0_G0 = NA_real_,
        TE = NA_real_, IDE = NA_real_, IIE = NA_real_)
    } else {
      rows[[b]] <- data.frame(bootstrap_id = b, status = "ok", message = "",
        R_1_G1 = fit$R_1_G1, R_1_G0 = fit$R_1_G0, R_0_G0 = fit$R_0_G0,
        TE = fit$TE, IDE = fit$IDE, IIE = fit$IIE)
    }
    if (!is.null(checkpoint_path) && (b %% progress_every == 0L || b == B)) {
      completed <- do.call(rbind, rows[!vapply(rows, is.null, logical(1))])
      utils::write.csv(completed, checkpoint_path, row.names = FALSE)
      message("Primary bootstrap: ", nrow(completed), "/", B, " resamples complete.")
    }
  }
  do.call(rbind, rows)
}

bootstrap_percentile_ci <- function(bootstrap_df) {
  estimands <- c("R_1_G1", "R_1_G0", "R_0_G0", "TE", "IDE", "IIE")
  do.call(rbind, lapply(estimands, function(estimand) {
    x <- bootstrap_df[[estimand]][bootstrap_df$status == "ok"]
    x <- x[!is.na(x)]
    interval <- if (length(x)) {
      stats::quantile(x, c(0.025, 0.975), na.rm = TRUE, names = FALSE)
    } else {
      c(NA_real_, NA_real_)
    }
    data.frame(
      estimand = estimand, n_success = length(x),
      lower_95 = interval[1], upper_95 = interval[2]
    )
  }))
}
