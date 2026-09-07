# README: MEDIATION Analysis Guide

## Introduction

This repository contains the R workflow for the REGARD-VAP clinical longitudinal mediation analysis. The analysis examines baseline carbapenem resistance, appropriate active antibiotic treatment during Day 0--3, and 60-day mortality.

This README provides a step-by-step guide for running the available clinical workflow. It covers data preparation, study-population checks, observed longitudinal summaries, longitudinal g-formula estimation, and the currently implemented sensitivity analyses.

---

## Step-by-Step Guide

### Preparation

#### Step 1: Install R and RStudio

Make sure R and RStudio are installed on your computer:

- Download R from: <https://cran.r-project.org/>
- Download RStudio from: <https://posit.co/download/rstudio-desktop/>

Install the R packages used by this workflow:

```r
install.packages(c("readxl", "ggplot2", "gt"))
```

#### Step 2: Obtain the approved clinical data files

Use the approved study-data process to obtain:

```text
regardvap_patient_day0_day3_integrated_longitudinal_itt460_v.xlsx
regardvap_patient_day0_day60_integrated_longitudinal_itt460_.xlsx
regardvap_primary_severity_L0_L3_itt460.xlsx
itt.dataset.RDS
```

**Important:** Do not modify the source files or change their names before running the workflow. Do not upload participant-level data to this repository.

#### Step 3: Download the code

Go to the `<> Code` menu in this repository and download the ZIP file, or clone the repository with Git. Keep the source data in an approved local directory outside the repository.

---

### Preparing the analysis data

#### Step 1: Select an output directory

Choose an approved local output directory. The workflow creates:

```text
step00_data_prep/
step01_baseline/
step02_longitudinal_summary/
step03_main_gformula/
step04_sensitivity/
step05_hte/
```

#### Step 2: Run the clinical workflow

From the repository root, run:

```bash
Rscript scripts/00_run_regardvap_manuscript_pipeline_no_wgs.R \
  /approved/path/day0_day3.xlsx \
  /approved/path/day0_day60.xlsx \
  /approved/path/severity.xlsx \
  /approved/path/itt.dataset.RDS \
  /approved/path/outputs/regardvap_manuscript_run
```

The runner reads the source files, constructs Day 0--3 and participant-level analysis datasets, and then runs the available analysis stages.

#### Step 3: Data-preparation outputs

The following files are written to `step00_data_prep/`:

- `data_sources_used.csv`: source-file paths and row counts.
- `analysis_panel_long_day0_day3.csv`: Day 0--3 participant-day panel.
- `analysis_panel_wide_for_gformula.csv`: participant-level analysis dataset.
- `censoring_and_daily_event_summary.csv`: daily censoring, death and information summaries.

---

### Study population and baseline characteristics

The study-population diagnostic code is located in:

```text
01_study_population_and_baseline/R/01_cohort_diagnostics.R
```

It checks the cohort and writes baseline outputs to `step01_baseline/`.

#### Baseline characteristics

The workflow writes:

- `table1_baseline_summary.csv`
- `table1_patient_level_dataset.csv`
- `cohort_diagnostic_summary.csv`

Review the source-data audit, participant/day structure, missingness and censoring summaries before interpreting subsequent outputs.

---

### Longitudinal treatment and clinical summaries

The observed longitudinal diagnostic code is located in:

```text
02_longitudinal_trajectories/R/01_observed_trajectory_diagnostics.R
```

It describes the Day 0--3 observed process and writes results to `step02_longitudinal_summary/`.

#### Outputs

- `table2_longitudinal_day0_day3_summary.csv`
- `table2_observed_trajectory_distribution.csv`
- `figure3_daywise_summary.csv`
- `figure3_observed_processes.svg`

Review these outputs to confirm observed daily severity, treatment and microbiology-information patterns before running or reporting the mediation analysis.

---

### Longitudinal mediation analysis

The master analysis script is:

```text
03_longitudinal_mediation_analysis/01_run_regardvap_longitudinal_mediation_v3_no_wgs.R
```

The analysis uses these supporting files:

```text
03_longitudinal_mediation_analysis/R/01_analysis_contract.R
03_longitudinal_mediation_analysis/R/02_diagnostics.R
03_longitudinal_mediation_analysis/R/03_bootstrap.R
03_longitudinal_mediation_analysis/R/04_gformula_engine.R
03_longitudinal_mediation_analysis/R/05_sensitivity_runner.R
```

#### Main analysis

The main analysis writes g-formula estimates and diagnostic tables to `step03_main_gformula/`, including:

- `table3_main_gformula_estimates.csv`
- `diagnostic_positivity_by_day.csv`
- `diagnostic_positivity_by_history.csv`

#### Sensitivity analyses

The workflow writes the implemented sensitivity-analysis results to:

```text
step04_sensitivity/table4_sensitivity_analyses.csv
```

The current script runs three alternatives: the alternative severity variable; a model without lagged mediator history; and a model including the daily microbiology-information proxy.

#### Bootstrap

Bootstrap is not run unless requested. To run it, set the number of participant-level resamples and Monte Carlo simulations before starting the main script:

```bash
export REGARDVAP_N_BOOT=500
export REGARDVAP_BOOT_NSIM=4000
```

The bootstrap outputs are written to `step03_main_gformula/`.

---

### Heterogeneity of treatment effect

The directory `04_heterogeneity_of_treatment_effect/` is reserved for future work. No HTE analysis is run by the current workflow.

---

### Troubleshooting

If the workflow stops, first check that the four supplied input paths exist, that the source-file names have not changed, and that the required R packages are installed. Review `data_sources_used.csv` and the R error message before changing code or source data.
