# Step 2. Longitudinal trajectories

## Purpose

Make the observed clinical pathway and temporal assumptions inspectable before
the longitudinal g-formula is run. This is a descriptive section: it does not
claim a causal effect.

## Variables by day

| Component | Day 0 | Day 1--3 | Role |
| --- | --- | --- | --- |
| `A` | baseline carbapenem resistance | fixed | exposure |
| `L_t` | baseline SOFA | pragmatic daily severity score | time-varying clinical state |
| `O_t` | index-culture information | prior microbiology information | decision information |
| `M_t` | appropriate active treatment | appropriate active treatment | mediator sequence |
| observation status | alive/followed at start of day | alive/followed at start of day | at-risk process |
| `Y` | not assessed | death by Day 60 | outcome |

## Required figures and tables

```text
Figure 1  Measurement timeline
Figure 2  Longitudinal DAG
Table 2   Day 0--3 severity, treatment, and information summaries
Figure 3  Observed daily trajectories
Support   Distribution of observed treatment histories
```

## Clinical checks before Step 3

1. Confirm Day 0 is the clinically defensible symptom/diagnosis time zero.
2. Confirm the severity measure used as `L_t` is available before the same-day
   treatment decision.
3. Confirm whether `O_t` contains only information known before treatment, not
   a later final microbiology result.
4. Describe date-level treatment endpoints as “by end of Day 1”, rather than
   “within 48 hours”, unless timestamps support the latter.

## Expected local outputs

```text
3.2_longitudinal_trajectories/
├── figure1_timeline.svg
├── figure2_dag.svg
├── table2_day0_day3_treatment_severity_summary.csv
├── table2_day0_day3_treatment_severity_summary_pretty.html
├── table2_observed_trajectory_distribution.csv
├── figure3_observed_trajectories.png
└── figure3_daywise_summary.csv
```

`R/01_observed_trajectory_diagnostics.R` is executed by the master workflow and
writes a transparent day-by-exposure diagnostic table used to verify Table 2
and Figure 3 inputs.
