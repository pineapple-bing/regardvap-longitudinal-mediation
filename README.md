# README: MEDIATION Analysis Guide

## Introduction

This repository contains the clinical-analysis workflow for REGARD-VAP. The study asks whether the association between baseline carbapenem resistance and 60-day mortality is mediated by the sequence of receiving appropriate active antibiotic treatment during Day 0--3.

This is the complete project guide. It follows one sequence: verify the source data, build the cohort, describe the observed longitudinal process, and estimate the longitudinal g-formula. The repository contains **code and documentation only**; individual-level data and generated results stay in an approved local location and are never committed to GitHub.

## Study question and causal ordering

Among eligible VAP episodes, estimate the total effect of baseline carbapenem resistance (`A`) on Day-60 mortality (`Y`) and decompose it into an interventional direct effect (IDE) and interventional indirect effect (IIE) through appropriate active treatment on Day 0--3 (`M0`--`M3`).

```text
Baseline V, A -> L0, O0 -> M0 -> L1, O1 -> M1 -> L2, O2 -> M2 -> L3, O3 -> M3 -> Y
```

`L_t` is clinical severity measured before the treatment decision on day `t`; `O_t` is microbiology information available before that decision; and `V` is the baseline covariate set. This is a longitudinal interventional mediation analysis, not a one-time-point mediation model.

## Repository structure

```text
01_study_population_and_baseline/       # cohort and baseline diagnostics
02_longitudinal_trajectories/           # observed Day 0--3 process diagnostics
03_longitudinal_mediation_analysis/     # g-formula, sensitivity and bootstrap modules
04_heterogeneity_of_treatment_effect/   # reserved; HTE is not run
data/                                   # no participant data
scripts/                                # end-to-end clinical runner
```

The workflow entry point is `scripts/00_run_regardvap_manuscript_pipeline_no_wgs.R`.

## Software and local data

Install the R packages used by the scripts:

```r
install.packages(c("readxl", "ggplot2", "gt"))
```

Keep these approved clinical sources outside the repository:

```text
regardvap_patient_day0_day3_integrated_longitudinal_itt460_v.xlsx
regardvap_patient_day0_day60_integrated_longitudinal_itt460_.xlsx
regardvap_primary_severity_L0_L3_itt460.xlsx
itt.dataset.RDS
```

Do not upload these files, identifiers, raw WGS files or generated outputs to GitHub. Do not rename a source file without recording the change in the source-data audit.

Run the full clinical workflow from the repository root:

```bash
Rscript scripts/00_run_regardvap_manuscript_pipeline_no_wgs.R \
  /approved/path/day0_day3.xlsx \
  /approved/path/day0_day60.xlsx \
  /approved/path/severity.xlsx \
  /approved/path/itt.dataset.RDS \
  /approved/path/outputs/regardvap_manuscript_run
```

The runner creates local numbered output folders.

---

## Step 1. Study population and baseline characteristics

**Code:** `01_study_population_and_baseline/R/01_cohort_diagnostics.R`

### Purpose

Establish the eligible analytic cohort before fitting a causal model.

### Workflow

1. Read the Day 0--3, Day 0--60, severity and original CRF sources.
2. Record source-file provenance and row counts.
3. Build a person-day `panel` for days 0--3.
4. Build a person-level `analysis` dataset for the g-formula.
5. Check participant/day uniqueness, missingness, death and censoring.
6. Generate baseline characteristics by resistance status `A`.

### Outputs to review

```text
step00_data_prep/data_sources_used.csv
step00_data_prep/analysis_panel_long_day0_day3.csv
step00_data_prep/analysis_panel_wide_for_gformula.csv
step00_data_prep/censoring_and_daily_event_summary.csv
step01_baseline/table1_baseline_summary.csv
step01_baseline/cohort_diagnostic_summary.csv
```

### Interpretation boundary

In the current extract, `A=0` includes culture-negative or non-evaluable episodes. Do not label it “carbapenem susceptible” unless the population is restricted to confirmed susceptible index isolates.

---

## Step 2. Observed longitudinal trajectories

**Code:** `02_longitudinal_trajectories/R/01_observed_trajectory_diagnostics.R`

### Purpose

Describe the observed Day 0--3 clinical and treatment process before interpreting causal estimates.

### Workflow

1. Summarise daily severity, appropriate active treatment and microbiology-information availability.
2. Tabulate observed treatment histories.
3. Create the observed treatment/severity trajectory figure.
4. Export distributions used to assess model plausibility.

### Outputs to review

```text
step02_longitudinal_summary/table2_longitudinal_day0_day3_summary.csv
step02_longitudinal_summary/table2_observed_trajectory_distribution.csv
step02_longitudinal_summary/figure3_observed_processes.svg
step02_longitudinal_summary/figure3_daywise_summary.csv
```

### Required clinical check

Confirm that `L_t` and `O_t` are known before the same-day `M_t` decision. A day after death/loss to follow-up is not an ordinary treatment day; values remain missing rather than being set to zero.

---

## Step 3. Longitudinal mediation analysis

**Master code:** `03_longitudinal_mediation_analysis/01_run_regardvap_longitudinal_mediation_v3_no_wgs.R`

### Modules

```text
R/01_analysis_contract.R     blocks WGS/simulated fields and checks the time index
R/02_diagnostics.R           exposure and treatment-history positivity tables
R/03_bootstrap.R             patient-level bootstrap and percentile intervals
R/04_gformula_engine.R       nuisance models and counterfactual simulation
R/05_sensitivity_runner.R    prespecified alternative specifications
```

### Primary nuisance models

```text
P(M0 | A, L0, V)
P(L1 | A, L0, M0, V)
P(M1 | A, L1, M0, L0, V)
P(L2 | A, L1, M1, V)
P(M2 | A, L2, M1, L0, V)
P(L3 | A, L2, M2, V)
P(M3 | A, L3, M2, L0, V)
P(Y  | A, L0:L3, M0:M3, V)
```

`V` includes age, sex, Charlson score, country, site, ICU type and bacteria group. `baseline_appropriateabx_on_symptomdate` is excluded because it nearly duplicates `M0`; using it as a baseline covariate would adjust for the mediator.

### Counterfactual regimes and effects

| Regime | Outcome exposure | Mediator process | Interpretation |
| --- | --- | --- | --- |
| `R(1,G1)` | `A=1` | drawn under `A=1` | resistant world |
| `R(1,G0)` | `A=1` | drawn under `A=0` | resistant world with mediator distribution from `A=0` |
| `R(0,G0)` | `A=0` | drawn under `A=0` | reference world |

```text
TE  = R(1,G1) - R(0,G0)
IDE = R(1,G0) - R(0,G0)
IIE = R(1,G1) - R(1,G0)
```

All effects are risk differences.

### Primary outputs

```text
step03_main_gformula/table3_main_gformula_estimates.csv
step03_main_gformula/diagnostic_positivity_by_day.csv
step03_main_gformula/diagnostic_positivity_by_history.csv
step04_sensitivity/table4_sensitivity_analyses.csv
```

### Prespecified sensitivity analyses

1. Alternative daily severity definition (`L_alt`).
2. A mediator model without lagged mediator history.
3. Inclusion of the time-varying microbiology-information proxy (`O_t`).

These are reproducibility checks. Report the analytic sample and active covariates with every sensitivity result.

### Bootstrap confidence intervals

Bootstrap is off by default. After checking a successful main run:

```bash
export REGARDVAP_N_BOOT=500
export REGARDVAP_BOOT_NSIM=4000
```

Rerun the master workflow. It resamples participants, repeats fitting/simulation and writes `bootstrap_effect_estimates.csv` and `bootstrap_percentile_ci.csv`.

### Required checks before reporting causal estimates

1. Review positivity by exposure, day and treatment history.
2. Compare observed and simulated severity/treatment distributions.
3. Check nuisance-model calibration and separation.
4. Constrain simulated severity to its feasible support before final modelling.
5. Use patient-level bootstrap inference.
6. Assess sensitivity to unmeasured mediator--outcome confounding.

---

## Step 4. Heterogeneity of treatment effect

This folder is a placeholder. HTE is **not** implemented or reported yet. Start it only after the primary cohort, time ordering, positivity, outcome model and bootstrap inference are finalised.

## Scope and data-sharing rules

### WGS and plasmid data

This is a clinical mediation repository. Simulated WGS fields, plasmid-cluster variables and carbapenemase variables are excluded from the executable code. A real *Pseudomonas aeruginosa* WGS substudy, if pursued, must be a separate analysis with its own question, data dictionary and methods.

### Sharing a package with Danny

Only share a governance-approved, de-identified clinical analysis extract with a data dictionary, timeline/DAG and focused review questions. Do not share simulated WGS data, raw WGS, identifiers, or preliminary estimates presented as final results. Use the study-approved transfer route.

## Reproducibility record

For every shared result, record date, source-file versions, R/package versions, exact command, random seed, Monte Carlo draws, successful bootstrap resamples, and any deviation from this guide.
