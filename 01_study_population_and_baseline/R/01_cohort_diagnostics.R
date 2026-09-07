# Cohort-definition diagnostics used before causal-model fitting.

write_cohort_diagnostics <- function(panel, analysis, output_dir) {
  day0 <- panel[panel$day == 0, , drop = FALSE]
  flow <- data.frame(
    stage = c(
      "Day-0 records with non-missing participant ID",
      "Day-0 records with observed exposure A",
      "Person-level records with observed Day-60 outcome Y",
      "Complete daily mediator history M0-M3",
      "Complete primary severity history L0-L3"
    ),
    n = c(
      sum(!is.na(day0$subjid)),
      sum(!is.na(day0$subjid) & !is.na(day0$A)),
      sum(!is.na(analysis$Y)),
      sum(stats::complete.cases(analysis[c("M0", "M1", "M2", "M3")])),
      sum(stats::complete.cases(analysis[c("L0", "L1", "L2", "L3")]))
    ),
    stringsAsFactors = FALSE
  )
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  utils::write.csv(flow, file.path(output_dir, "diagnostic_cohort_flow.csv"), row.names = FALSE)
  invisible(flow)
}
