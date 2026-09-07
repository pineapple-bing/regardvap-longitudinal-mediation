args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L) stop("Usage: Rscript 02_descriptive_analysis/02_treatment_by_exposure.R config.R")
source(args[[1]]); source("R/functions_data.R")
data <- readRDS(analysis_data_path)
summarise_binary <- function(x, label) data.frame(
  endpoint = label, A = c(0, 1),
  n = c(sum(data$A == 0 & !is.na(x)), sum(data$A == 1 & !is.na(x))),
  proportion = c(mean(x[data$A == 0], na.rm = TRUE), mean(x[data$A == 1], na.rm = TRUE))
)
table_out <- do.call(rbind, list(
  summarise_binary(data$appropriate_by_end_day1, "Appropriate by end of Day 1"),
  summarise_binary(data$appropriate_by_end_day2, "Appropriate by end of Day 2"),
  summarise_binary(data$appropriate_by_end_day3, "Appropriate by end of Day 3")
))
write_csv_safely(table_out, "outputs/tables/02_early_treatment_by_exposure.csv")
