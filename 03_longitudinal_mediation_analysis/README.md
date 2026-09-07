# Step 3. Longitudinal mediation analysis

## Purpose

This module implements the clinical v3 longitudinal g-formula after WGS,
plasmid, carbapenemase, and all simulated WGS variables have been removed. It
estimates the interventional analogue of direct and indirect effects of
baseline carbapenem resistance on 60-day mortality through Day 0--3 appropriate
active treatment.

The executable is preserved as a complete, auditable v3 source file while this
guide makes its internal stages explicit.

```text
01_run_regardvap_longitudinal_mediation_v3_no_wgs.R
├── A. Read, validate and harmonise approved clinical sources
├── B. Build Day 0--3 long and person-level wide panels
├── C. Describe observed baseline and longitudinal data
├── D. Fit daily treatment, severity and outcome nuisance models
├── E. Simulate three counterfactual regimes
├── F. Calculate TE, IDE and IIE
├── G. Run prespecified sensitivity specifications
└── H. Export local tables, figures and diagnostics
```

## Inputs

| Source | Role | Required timing |
| --- | --- | --- |
| Day 0--3 integrated clinical file | `A`, daily `M_t`, observation status, candidate `O_t`, baseline covariates | Day 0--3 |
| Day 0--60 integrated file | Day-60 mortality `Y` | follow-up |
| Primary severity file | `L0` and daily `L1`--`L3` construction | before `M_t` |
| ITT RDS bundle | provenance/source audit only | baseline and daily CRF source |

## A. Data construction

The code creates two local objects:

- `panel`: one row per participant-day (`t=0,1,2,3`) for observed trajectories;
- `analysis`: one row per participant with `M0`--`M3`, `L0`--`L3`, `O0`--`O3`,
  `A`, `Y`, and baseline covariates for the g-formula.

Participants without observed `A` or `Y` are excluded from the current wide
analysis. Days after death or loss to follow-up are missing, not ordinary
treatment zeros.

## B. Nuisance models

For the primary specification, the code fits:

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

`V` includes age, sex, Charlson score, country, site, ICU type and baseline
bacteria group. It excludes `baseline_appropriateabx_on_symptomdate`, because
that field nearly duplicates `M0`.

## C. Counterfactual simulation

| Regime | Exposure in outcome model | Mediator process | Interpretation |
| --- | --- | --- | --- |
| `R(1,G1)` | `A=1` | drawn under `A=1` | resistant world |
| `R(1,G0)` | `A=1` | drawn under `A=0` | resistant world with non-resistant mediator distribution |
| `R(0,G0)` | `A=0` | drawn under `A=0` | non-resistant world |

```text
TE  = R(1,G1) − R(0,G0)
IDE = R(1,G0) − R(0,G0)
IIE = R(1,G1) − R(1,G0)
```

Effects are risk differences. Monte Carlo draws approximate each regime risk;
Monte Carlo stability must be checked before interpretation.

## D. Prespecified sensitivity analyses

The legacy v3 runner currently exports:

1. an alternative daily severity construction (`L_alt`);
2. an altered mediator-history specification; and
3. inclusion of the time-varying microbiology-information proxy (`O_t`).

These are reproducibility outputs, not a claim that every variation is a clean
one-factor sensitivity analysis. Report the analytic sample and active terms
with every result.

## E. Diagnostics required before reporting

The code already exports completeness, censoring and observed-trajectory
summaries. Before final reporting, add and review:

- positivity tables by exposure and daily treatment history;
- observed-versus-simulated severity and treatment distributions;
- nuisance-model calibration;
- constrained modelling of bounded/discrete severity;
- patient-level nonparametric bootstrap confidence intervals; and
- sensitivity analysis for unmeasured mediator--outcome confounding.

## Outputs

```text
3.3_longitudinal_mediation_analysis/
├── table3_counterfactual_risks_te_ide_iie.csv
├── table3_counterfactual_risks_te_ide_iie_pretty.html
├── table4_sensitivity_analyses.csv
├── table4_sensitivity_analyses_pretty.html
├── figure4_te_decomposition.png
└── figure4_te_decomposition_data.csv
```

## What this module does not do

- It does not analyse WGS, plasmid clusters, or carbapenemase genes.
- It does not imply `A=0` is confirmed susceptibility.
- It does not provide final confidence intervals until bootstrap inference is
  implemented.
- It does not replace clinical confirmation of same-day ordering.
