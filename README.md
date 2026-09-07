# REGARD-VAP longitudinal mediation analysis

Reproducible code scaffold for studying whether early appropriate antibiotic
treatment mediates the association between baseline carbapenem resistance and
60-day mortality.

## Current status

Method-development repository. The previous v3 estimates are preliminary and
must not be reported as final causal estimates. The main unresolved decision is
whether the exposure contrast is:

1. carbapenem-resistant infection versus confirmed carbapenem-susceptible
   infection; or
2. detection of carbapenem resistance versus no detection in the full cohort.

These are different target populations and estimands. In the current integrated
file, `baseline_carba_r = 0` includes culture-negative episodes, so it must not
be labelled "carbapenem susceptible" without additional eligibility rules.

## Proposed primary timeline

- Day 0: first VAP symptom/diagnosis date.
- `M_t`: appropriate active treatment on each calendar day, Day 0 through Day 3.
- `L_t`: clinical severity measured before that day's treatment decision.
- Outcome: death by Day 60.
- Descriptive early-treatment endpoint: appropriate active treatment by the end
  of Day 1. With date-only data, call this "by end of Day 1", not an exact
  "within 48 hours" endpoint.

The ACORN-HAI empirical window (Day -1 to Day +2) is useful as a secondary
descriptive definition. It is not a direct substitute for the daily mediator
sequence in the longitudinal g-formula.

## Repository safety

Participant-level data, raw WGS files, simulated WGS variables, generated
outputs, PowerPoint files, and Word documents are intentionally excluded.

## First checks

1. Copy `config.example.R` to `config.R` and set authorised local paths.
2. Run the scripts in numeric order, starting with
   `Rscript 01_data_preparation/01_validate_inputs.R config.R`.
3. Resolve every `FAIL` and review every `WARN` before model fitting.

See `docs/analysis_plan.md` for the proposed primary and sensitivity analyses.

## Project layout

```
01_data_preparation/     # validate inputs and build a person-level Day 0–3 panel
02_descriptive_analysis/ # cohort flow and observed treatment summaries
03_primary_analysis/     # locked analysis specification and modelling entry point
04_sensitivity_analysis/ # prespecified alternatives, not post-hoc fishing
05_figures/              # report tables and figures
R/                        # small reusable helper functions
data/                     # documentation only; no participant-level files
outputs/                  # local derived data/tables/figures, ignored by git
```

The numbered folders mirror the transparent ACORN-HAI structure. Unlike its
baseline-treatment analyses, this project preserves `L_t -> M_t -> L_(t+1) -> Y`
for longitudinal mediation.
