args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L) stop("Usage: Rscript 01_data_preparation/02_build_analysis_dataset.R config.R")
source(args[[1]])
source("R/functions_data.R")
suppressPackageStartupMessages({ library(readxl); library(dplyr); library(tidyr) })

day_data <- as.data.frame(read_excel(day03_path))
severity_data <- as.data.frame(read_excel(severity_path))
outcome_data <- as.data.frame(read_excel(day60_path))
stop_if_missing(day_data, c(id_var, day_var, exposure_var, mediator_var), "day03 file")
stop_if_missing(severity_data, c(id_var, day_var, severity_var), "severity file")
stop_if_missing(outcome_data, c(id_var, outcome_var), "day60 file")

# One row per person. Do not add future microbiology variables here.
panel <- day_data |>
  filter(.data[[day_var]] %in% 0:3) |>
  transmute(!!id_var := .data[[id_var]], day = .data[[day_var]],
    A = as_binary(.data[[exposure_var]], exposure_var), M = as_binary(.data[[mediator_var]], mediator_var),
    alive = as_binary(.data[["alive_at_start_of_day"]], "alive_at_start_of_day"),
    followed = as_binary(.data[["under_followup_on_day"]], "under_followup_on_day"),
    across(any_of(baseline_covariates)))
severity_panel <- severity_data |>
  filter(.data[[day_var]] %in% 0:3) |>
  transmute(!!id_var := .data[[id_var]], day = .data[[day_var]], L = .data[[severity_var]])
wide <- panel |>
  left_join(severity_panel, by = c(id_var, "day")) |>
  select(all_of(id_var), day, A, M, L, alive, followed, any_of(baseline_covariates)) |>
  pivot_wider(names_from = day, values_from = c(M, L, alive, followed), names_glue = "{.value}{day}") |>
  left_join(outcome_data |> select(all_of(id_var), Y = all_of(outcome_var)), by = id_var) |>
  mutate(appropriate_by_end_day1 = any_observed_one(M0, M1),
    appropriate_by_end_day2 = any_observed_one(M0, M1, M2),
    appropriate_by_end_day3 = any_observed_one(M0, M1, M2, M3),
    appropriate_days_day0_3 = rowSums(across(c(M0, M1, M2, M3)), na.rm = TRUE))
dir.create(dirname(analysis_data_path), recursive = TRUE, showWarnings = FALSE)
saveRDS(wide, analysis_data_path)
message("Wrote person-level analysis dataset: ", analysis_data_path)
