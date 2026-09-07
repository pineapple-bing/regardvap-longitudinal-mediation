args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L) stop("Usage: Rscript 02_descriptive_analysis/01_cohort_flow.R config.R")
source(args[[1]]); source("R/functions_data.R")
data <- readRDS(analysis_data_path)
flow <- data.frame(
  stage = c("All Day-0 records", "Known exposure", "Observed 60-day outcome", "Complete M0–M3"),
  n = c(nrow(data), sum(!is.na(data$A)), sum(!is.na(data$Y)), sum(stats::complete.cases(data[c("M0", "M1", "M2", "M3")])))
)
write_csv_safely(flow, "outputs/tables/01_cohort_flow.csv")
