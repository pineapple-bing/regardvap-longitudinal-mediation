# Step 1. Study population and baseline characteristics

## Purpose

Define the analytic cohort before fitting a model. This stage is where the
clinical source files are audited, the Day 0--3 panel is checked, and the
meaning of the exposure comparison is made explicit.

## Inputs and checks

| Check | Why it matters |
| --- | --- |
| one Day-0 record per participant | establishes time zero and baseline `A` |
| four records per participant where at risk | verifies the daily observation process |
| `A` and Day-60 `Y` observed | defines the current analysis population |
| `M_t`, `L_t`, `O_t` availability by day | determines modelled histories and complete-case selection |
| death/loss-to-follow-up status | prevents post-event variables being treated as observed zeros |

## Exposure decision

There are two scientifically different choices. Select one before interpreting
effect estimates.

1. **Etiologic contrast:** confirmed carbapenem-resistant versus confirmed
   carbapenem-susceptible index isolates. Culture-negative/non-evaluable
   episodes are excluded.
2. **Surveillance contrast:** resistance detected versus not detected in the
   full eligible cohort. Here `A=0` is not “susceptible”.

The current v3 code retains the latter binary structure; the first requires a
verified eligibility rule.

## Expected local outputs

```text
_scratch/longitudinal_core/step00_data_prep/
├── data_sources_used.csv
├── analysis_panel_long_day0_day3.csv
├── analysis_panel_wide_for_gformula.csv
├── censoring_and_daily_event_summary.csv
├── lt_nonmissing_counts_by_day.csv
└── lt_complete_table.csv

3.1_study_population_and_baseline/
└── table1_baseline_characteristics.csv
```

## Reporting boundary

Table 1 is descriptive. Do not use p-values in Table 1 as a criterion for
confounder selection, and do not include baseline appropriate treatment because
it is almost the same variable as `M0`.

`R/01_cohort_diagnostics.R` is executed by the master workflow and writes an
explicit cohort-flow diagnostic alongside the local analysis panel.
