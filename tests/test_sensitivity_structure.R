source(file.path("03_longitudinal_mediation_analysis", "R", "05_sensitivity_runner.R"))

specs <- make_sensitivity_specifications()
stopifnot(length(specs) == 14L)
stopifnot(sum(vapply(specs, `[[`, logical(1), "same_estimand")) == 9L)
stopifnot(sum(!vapply(specs, `[[`, logical(1), "same_estimand")) == 5L)
stopifnot(identical(make_primary_specification()$id, "PRIMARY"))

set.seed(1)
test_df <- data.frame(
  A = rep(c(0, 1), 20), Y = rep(c(0, 0, 1, 1), 10),
  bacteria = factor(rep(c("Klebsiella", "Culture_neg"), 20)),
  site = factor(rep(c("S1", "S2"), each = 20))
)
population <- make_population_description(
  test_df,
  list(make_primary_specification(), specs[[6]])
)
stopifnot(nrow(population) == 2L)
stopifnot(population$n[population$sensitivity_id == "PRIMARY"] == 40L)
stopifnot(population$n[population$sensitivity_id == "S6_culture_positive"] == 20L)

draws <- data.frame(
  bootstrap_id = rep(1:10, each = 2),
  sensitivity_id = rep(c("PRIMARY", "S1_alternative_severity"), 10),
  status = "ok", error_message = "",
  R_1_G1 = seq_len(20) / 100, R_1_G0 = seq_len(20) / 110,
  R_0_G0 = seq_len(20) / 120, TE = seq_len(20) / 130,
  IDE = seq_len(20) / 140, IIE = seq_len(20) / 150
)
ci <- sensitivity_bootstrap_percentile_ci(draws)
stopifnot(nrow(ci) == 12L)
stopifnot(all(ci$successful_resamples == 10L))
stopifnot(all(ci$lower_95 <= ci$upper_95))

cat("Sensitivity structure tests passed.\n")
