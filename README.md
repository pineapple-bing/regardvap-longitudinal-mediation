<p align="center">
  <strong>REGARD-VAP</strong><br>
  <em>Longitudinal mediation analysis of carbapenem resistance, appropriate treatment, and 60-day mortality</em>
</p>

---

# REGARD-VAP analysis guide

## Introduction

This repository contains the reproducible clinical-analysis workflow for the
REGARD-VAP longitudinal mediation study. The study asks whether the association
between baseline carbapenem resistance and 60-day mortality is partly explained
by the daily sequence of receiving appropriate active antibiotic treatment.

The repository is deliberately organised as a step-by-step guide. Every
analysis section states its input data, its analytic purpose, the expected
outputs, and the interpretation boundary. The structure is inspired by the
transparent organisation of the public ACORN-HAI analysis repository, while the
causal model, variables, and estimand are specific to REGARD-VAP.

## Study question and estimand

Among eligible VAP episodes, estimate the total effect of baseline carbapenem
resistance (`A`) on death by Day 60 (`Y`) and decompose it into an
interventional direct effect (IDE) and an interventional indirect effect (IIE)
through appropriate active treatment on Day 0--3 (`M0`--`M3`).

The daily causal ordering is:

```text
Baseline V, A ──► L0, O0 ──► M0 ──► L1, O1 ──► M1 ──► L2, O2 ──► M2 ──► L3, O3 ──► M3 ──► Y
```

`L_t` is clinical severity measured before the day-t treatment decision and
`O_t` is microbiology information available before that decision. This is an
interventional longitudinal-mediation analysis, not a conventional single-time
point mediation model.

## Repository structure

```text
data/                                  # data dictionary and local-data instructions only
scripts/                               # one reproducible manuscript runner
01_study_population_and_baseline/      # source audit, eligibility, cohort flow, Table 1
02_longitudinal_trajectories/          # timeline, DAG, Table 2, observed trajectories
03_longitudinal_mediation_analysis/    # longitudinal g-formula, Tables 3--4, Figure 4
04_heterogeneity_of_treatment_effect/  # reserved; no analysis is run yet
docs/                                  # analysis decisions, code map, collaboration handoff
outputs/                               # generated locally; never committed
```

## Before you run anything

### Software

Install a current version of R and the packages used by the scripts:

```r
install.packages(c("readxl", "ggplot2", "gt"))
```

### Local data layout

Participant-level files are **not** stored in this repository. In an approved
local location, keep the following sources unchanged:

```text
regardvap_patient_day0_day3_integrated_longitudinal_itt460_v.xlsx
regardvap_patient_day0_day60_integrated_longitudinal_itt460_.xlsx
regardvap_primary_severity_L0_L3_itt460.xlsx
itt.dataset.RDS
```

Do not upload any of these files to GitHub. Do not rename a source file without
updating the run command and recording the change in the data-source audit.

### Run the full clinical workflow

From the repository root, run:

```bash
Rscript scripts/00_run_regardvap_manuscript_pipeline_no_wgs.R \
  /approved/path/day0_day3.xlsx \
  /approved/path/day0_day60.xlsx \
  /approved/path/severity.xlsx \
  /approved/path/itt.dataset.RDS \
  /approved/path/outputs/regardvap_manuscript_run
```

The runner creates four numbered output directories. Sections 1--3 contain
analysis outputs; Section 4 only records that HTE is not yet run.

---

## Step-by-step analysis guide

### Step 1. Study population and baseline characteristics

**Folder:** `01_study_population_and_baseline/`

**Purpose:** Establish the eligible analytic cohort before any causal model is
fitted. Audit the four source files, confirm unique participant identifiers,
check the Day 0--3 observation structure, and document missingness and
censoring.

**Key outputs:**

- source-data audit;
- long and wide analysis-panel files (local only);
- cohort-flow and censoring summaries; and
- Table 1 baseline characteristics by `A`.

**Interpretation rule:** In the current extract `A=0` includes
culture-negative/non-evaluable episodes. It must not be labelled “carbapenem
susceptible” unless the population is restricted to confirmed susceptible
index isolates.

See [Step 1 guide](01_study_population_and_baseline/README.md).

### Step 2. Longitudinal trajectories

**Folder:** `02_longitudinal_trajectories/`

**Purpose:** Describe, rather than causally estimate, the observed Day 0--3
clinical pathway. This section makes the temporal assumptions visible before
the g-formula is interpreted.

**Key outputs:**

- Figure 1: measurement timeline;
- Figure 2: longitudinal DAG;
- Table 2: daily severity, treatment and microbiology-information summaries;
- Figure 3: observed treatment/severity trajectories; and
- a trajectory-distribution table for model checking.

See [Step 2 guide](02_longitudinal_trajectories/README.md).

### Step 3. Longitudinal mediation analysis

**Folder:** `03_longitudinal_mediation_analysis/`

**Purpose:** Estimate counterfactual risks under three regimes:

```text
R(1, G1): A = resistant; mediator process drawn under A = resistant
R(1, G0): A = resistant; mediator process drawn under A = non-resistant
R(0, G0): A = non-resistant; mediator process drawn under A = non-resistant
```

The effect decomposition is:

```text
TE  = R(1, G1) - R(0, G0)
IDE = R(1, G0) - R(0, G0)
IIE = R(1, G1) - R(1, G0)
```

**Core executable:**
`03_longitudinal_mediation_analysis/01_run_regardvap_longitudinal_mediation_v3_no_wgs.R`

**Key outputs:**

- Table 3: counterfactual risks, total effect, IDE and IIE;
- Table 4: prespecified sensitivity analyses; and
- Figure 4: risk and effect-decomposition display.

The Step 3 guide separates the g-formula into data construction, nuisance
models, simulation, diagnostics, uncertainty, and sensitivity analyses so that
the large v3 script can be audited rather than treated as a black box.

See [Step 3 guide](03_longitudinal_mediation_analysis/README.md).

### Step 4. Heterogeneity of treatment effect

**Folder:** `04_heterogeneity_of_treatment_effect/`

This is intentionally a placeholder. HTE will not be run until the primary
population, time ordering, positivity diagnostics, outcome model, and bootstrap
uncertainty have been finalised.

---

## Methods safeguards and current decisions

1. **Do not adjust for `baseline_appropriateabx_on_symptomdate`.** It agrees
   with `M0` for almost all observed episodes and is therefore a treatment/
   mediator variable, not a baseline confounder.
2. **No WGS in the clinical pipeline.** All simulated WGS, plasmid-cluster and
   carbapenemase variables have been removed from the executable clinical code.
   A real *Pseudomonas aeruginosa* WGS substudy, if pursued, must be separate
   and labelled as such.
3. **Do not call an exact 48-hour window from date-only data.** The current
   descriptive endpoint is “appropriate active treatment by the end of Day 1”.
4. **The legacy v3 severity simulation is retained for reproducibility, not as
   the final model.** The final version should constrain simulated severity to
   its feasible support and report observed-versus-simulated checks.
5. **Uncertainty remains incomplete without a patient-level bootstrap.** The
   current v3 core is a reproducible starting point, not a final reportable
   causal analysis.

## Working with Danny

Before sending any data or results, use the
[Danny handoff checklist](docs/danny_data_handoff.md). The appropriate package
is a governance-approved clinical analysis extract plus data dictionary and
specific review questions. Do **not** include simulated WGS fields, raw WGS,
unapproved identifiers, or claim that the mediation estimates are final.

## Reproducibility and support

The repository contains code and documentation only. Generated outputs are
local and ignored by Git. Record software versions, source-file dates, and any
analytic deviation in `docs/` before a result is shared externally.
