# Diagnostics separated from model fitting so they are visible and auditable.

write_positivity_diagnostics <- function(analysis_df, output_dir) {
  required <- c("A", "M0", "M1", "M2", "M3")
  if (!all(required %in% names(analysis_df))) {
    stop("Cannot calculate positivity diagnostics: missing A or daily mediator variables.", call. = FALSE)
  }

  day_rows <- lapply(0:3, function(day) {
    mediator <- paste0("M", day)
    tab <- as.data.frame(table(A = analysis_df$A, M = analysis_df[[mediator]], useNA = "no"))
    names(tab)[names(tab) == "Freq"] <- "n"
    tab$day <- day
    tab$proportion_within_A <- ave(tab$n, tab$A, FUN = function(x) x / sum(x))
    tab
  })

  complete_history <- stats::complete.cases(analysis_df[c("A", "M0", "M1", "M2", "M3")])
  history <- analysis_df[complete_history, c("A", "M0", "M1", "M2", "M3"), drop = FALSE]
  history$history <- paste0(history$M0, history$M1, history$M2, history$M3)
  history_tab <- as.data.frame(table(A = history$A, history = history$history, useNA = "no"))
  names(history_tab)[names(history_tab) == "Freq"] <- "n"
  history_tab$proportion_within_A <- ave(history_tab$n, history_tab$A, FUN = function(x) x / sum(x))

  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  utils::write.csv(do.call(rbind, day_rows), file.path(output_dir, "diagnostic_positivity_by_day.csv"), row.names = FALSE)
  utils::write.csv(history_tab, file.path(output_dir, "diagnostic_positivity_by_history.csv"), row.names = FALSE)
  invisible(list(by_day = do.call(rbind, day_rows), by_history = history_tab))
}
