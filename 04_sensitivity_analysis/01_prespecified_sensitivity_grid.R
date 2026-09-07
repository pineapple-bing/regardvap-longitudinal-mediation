args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L) stop("Usage: Rscript 04_sensitivity_analysis/01_prespecified_sensitivity_grid.R config.R")
source(args[[1]]); source("R/functions_data.R")
grid <- data.frame(
  analysis_id = sprintf("S%02d", 1:8),
  domain = c("Population", "Mediator window", "Treatment definition", "Severity", "Missingness", "Site", "Organism", "Model fit"),
  alternative = c(
    "Resistance detected vs not detected; label contrast correctly",
    "Day 0–1 and Day 0–2 alongside Day 0–3 primary sequence",
    "Newly started appropriate treatment; prescription-end handling",
    "Alternative pre-decision severity construction",
    "Longitudinal multiple imputation versus complete case",
    "Country fixed effects plus site-cluster bootstrap / hierarchical sensitivity",
    "Monomicrobial and major organism-group exploratory analyses",
    "Flexible terms/interactions and observed-versus-simulated trajectory checks"
  )
)
write_csv_safely(grid, "outputs/tables/04_prespecified_sensitivity_grid.csv")
