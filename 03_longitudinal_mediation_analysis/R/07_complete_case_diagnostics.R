# Audit the explicitly defined primary complete-case cohort.

write_complete_case_diagnostics <- function(analysis_df, output_dir, horizon = 3L) {
  covars <- c("age", "male", "charlson", "country", "site", "icu_type", "bacteria")
  required <- unique(c("A", "Y", paste0("M", 0:horizon), paste0("L", 0:horizon), covars))
  baseline_eligible <- !is.na(analysis_df$A) & !is.na(analysis_df$Y)
  complete_primary <- baseline_eligible & stats::complete.cases(analysis_df[, required, drop = FALSE])
  missing_counts <- data.frame(
    variable = required,
    n_missing_among_baseline_eligible = vapply(required, function(v) sum(is.na(analysis_df[[v]]) & baseline_eligible), integer(1)),
    stringsAsFactors = FALSE
  )
  flow <- data.frame(
    stage = c("All rows in participant-level analysis data", "Eligible: observed A and Day-60 Y", "Primary complete-case cohort"),
    n = c(nrow(analysis_df), sum(baseline_eligible), sum(complete_primary)),
    stringsAsFactors = FALSE
  )
  inclusion_by_ay <- as.data.frame(with(analysis_df[baseline_eligible, , drop = FALSE], table(A, Y, complete_case = complete_primary[baseline_eligible], useNA = "ifany")))
  model_n <- data.frame(
    model = c(paste0("M", 0:horizon), paste0("L", seq_len(horizon)), "Y"),
    n_primary_complete_case = sum(complete_primary),
    stringsAsFactors = FALSE
  )
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  utils::write.csv(flow, file.path(output_dir, "diagnostic_primary_complete_case_flow.csv"), row.names = FALSE)
  utils::write.csv(missing_counts, file.path(output_dir, "diagnostic_primary_missingness_by_variable.csv"), row.names = FALSE)
  utils::write.csv(inclusion_by_ay, file.path(output_dir, "diagnostic_complete_case_inclusion_by_A_Y.csv"), row.names = FALSE)
  utils::write.csv(model_n, file.path(output_dir, "diagnostic_nuisance_model_sample_sizes.csv"), row.names = FALSE)
  invisible(list(required = required, complete_primary = complete_primary, flow = flow))
}
