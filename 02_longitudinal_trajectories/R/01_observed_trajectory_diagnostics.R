# Descriptive checks that make the observed longitudinal process transparent.

write_observed_trajectory_diagnostics <- function(panel, output_dir) {
  required <- c("day", "A", "M", "L_main", "O_info")
  if (!all(required %in% names(panel))) {
    stop("Panel is missing variables for observed trajectory diagnostics.", call. = FALSE)
  }
  groups <- split(panel, interaction(panel$day, panel$A, drop = TRUE))
  out <- do.call(rbind, lapply(groups, function(d) {
    data.frame(
      day = d$day[1], A = d$A[1], n_records = nrow(d),
      mediator_observed = sum(!is.na(d$M)),
      appropriate_active_n = sum(d$M == 1, na.rm = TRUE),
      appropriate_active_proportion = mean(d$M == 1, na.rm = TRUE),
      severity_observed = sum(!is.na(d$L_main)),
      severity_mean = mean(d$L_main, na.rm = TRUE),
      micro_info_available_n = sum(d$O_info == 1, na.rm = TRUE),
      stringsAsFactors = FALSE
    )
  }))
  out$appropriate_active_proportion[is.nan(out$appropriate_active_proportion)] <- NA_real_
  out$severity_mean[is.nan(out$severity_mean)] <- NA_real_
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  utils::write.csv(out, file.path(output_dir, "diagnostic_observed_trajectory_by_day.csv"), row.names = FALSE)
  invisible(out)
}
