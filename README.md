# REGARD-VAP: longitudinal mediation of carbapenem resistance and mortality

This is the codebase for the REGARD-VAP clinical longitudinal-mediation study,
not a generic mediation template. It reconstructs the v3 pipeline into a
manuscript-facing workflow with its figures and tables preserved.

## Clinical question

Among eligible VAP episodes, what is the total effect of baseline carbapenem
resistance (`A`) on 60-day mortality (`Y`), and what component is mediated by
the Day 0--3 sequence of appropriate active treatment (`M0`--`M3`)?

The model uses daily clinical severity (`L_t`) and pre-decision microbiology
information (`O_t`) as time-varying history. The primary estimand is an
interventional direct/indirect-effect decomposition from a longitudinal
g-formula.

## Repository map

```text
data/                              # documentation only; no participant data
scripts/                           # one reproducible manuscript runner
01_study_population_and_baseline/  # cohort definition, source audit, Table 1
02_longitudinal_trajectories/      # timeline, DAG, Table 2, trajectory figure
03_longitudinal_mediation_analysis/# v3 clinical g-formula; Tables 3--4, Figure 4
04_heterogeneity_of_treatment_effect/ # reserved; no results yet
docs/                              # analysis plan, code map, decisions
outputs/                           # local only and git-ignored
```

The central executable is
`scripts/00_run_regardvap_manuscript_pipeline_no_wgs.R`. It calls the
complete v3 clinical longitudinal-mediation code in
`03_longitudinal_mediation_analysis/01_run_regardvap_longitudinal_mediation_v3_no_wgs.R`
and creates manuscript-facing output sections:

1. Study population and baseline characteristics
2. Longitudinal trajectories (timeline, DAG, Table 2, Figure 3)
3. Longitudinal mediation analysis (Table 3, Table 4, Figure 4)
4. Heterogeneity of treatment effect (future work only)

## WGS policy

All simulated WGS/plasmid/carbapenemase fields have been removed from the
clinical pipeline and from HTE. The current repository does not execute a WGS
analysis. A future real *Pseudomonas aeruginosa* WGS substudy must be rebuilt
from the verified real dataset in a separate, explicitly labelled module.

## Safety

No participant-level data, raw WGS files, simulation fields, generated tables,
figures, PowerPoint files, or Word documents are versioned. Copy
`config.example.R` to a local ignored `config.R` before running.

## Important current methods notes

- `baseline_appropriateabx_on_symptomdate` is excluded from adjustment because
  it nearly duplicates `M0` in the current extract.
- The `A=0` group includes culture-negative episodes. It must not be called
  “carbapenem susceptible” unless eligibility is restricted to confirmed
  susceptible index isolates.
- The legacy v3 severity simulator is retained for reproducibility only; the
  final analysis should constrain simulated severity to its feasible support.
