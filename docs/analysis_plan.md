# Analysis plan and code audit

## Primary question

Among eligible REGARD-VAP participants, estimate the total effect of baseline
carbapenem resistance on 60-day mortality and the interventional indirect effect
operating through the Day 0-3 sequence of appropriate active treatment.

The effect is an interventional analogue, not a natural direct/indirect effect.
The exact stochastic intervention for the mediator must be written explicitly
before the final model is run.

## Decisions required before inference

### Exposure and target population

Choose one and do not switch after viewing results:

- **Preferred etiologic contrast:** resistant versus confirmed susceptible
  infection among episodes with an eligible index organism and interpretable
  carbapenem susceptibility. Culture-negative/non-evaluable episodes are
  excluded.
- **Alternative surveillance contrast:** resistance detected versus not detected
  in the full VAP cohort. Here A=0 is not "susceptible" and conclusions must use
  detection language.

### Time zero and within-day order

Day 0 is the first VAP symptom/diagnosis date in the current files. For every
day, construct variables so that severity and microbiology information available
before the treatment decision precede `M_t`. Same-day summaries that mix values
recorded after treatment initiation can introduce temporal ambiguity.

### Treatment windows

- Primary mediator: daily appropriate active treatment, `M0`-`M3`.
- Secondary summary: any appropriate active treatment by end of Day 1.
- Sensitivity: any appropriate active treatment by end of Day 2.
- ACORN-aligned descriptive window: any empirical treatment active from Day -1
  through Day +2, only if Day -1 data are complete enough.
- Definitive treatment on Day 3-5 is a separate treatment-transition analysis,
  not the primary mediator.

## Problems in the previous v3 code

1. Simulated WGS fields were relabelled as observed WGS and included in summaries
   and heterogeneity analyses. They must be removed from the clinical pipeline.
2. `baseline_appropriateabx_on_symptomdate` was included as a baseline covariate,
   although it agrees with `M0` for nearly all participants. This duplicates or
   adjusts for the mediator and risks separation/overadjustment.
3. `country` and `site` were both entered as fixed categorical effects even
   though site is nested within country. Use country fixed effects plus a
   prespecified site strategy, or a hierarchical model when feasible.
4. Complete-case fitting retained 404 of 460 participants because daily severity
   was missing. All 460 were alive and under follow-up through Day 3 in the
   current extract, so this is missing-covariate selection rather than early-death
   selection. Missingness still requires a justified strategy.
5. The `include_ot=TRUE` run retained only 153 participants and then dropped all
   O variables as invariant. Its difference from the main result is therefore a
   selected-sample result, not an O-adjusted sensitivity analysis.
6. The `mediator_history="none"` run simultaneously removed lagged mediator
   terms from mediator models and replaced the full mediator history in the
   outcome model with `M3`. It changes two features and is not a clean one-factor
   sensitivity analysis.
7. Severity was simulated with an unbounded Gaussian model even though the score
   is bounded and discrete, permitting impossible simulated values.
8. Only Monte Carlo variability was used. There are no bootstrap confidence
   intervals for total, direct, or indirect effects.
9. Positivity, model calibration, simulated-versus-observed trajectory checks,
   and Monte Carlo error were not reported.

## Recommended analysis set

### Primary analysis

- Prespecified exposure contrast and eligible population.
- Daily `L_t -> M_t` ordering confirmed from timestamps or data construction.
- Flexible models for treatment and severity; severity predictions constrained
  to its support.
- Interventional direct and indirect effects on the risk-difference scale.
- Patient-level nonparametric bootstrap confidence intervals, with the entire
  model-fitting and simulation procedure repeated inside each bootstrap sample.
- Monte Carlo sample large enough that Monte Carlo error is negligible relative
  to bootstrap uncertainty.

### Core sensitivity analyses

1. **Exposure population:** confirmed resistant versus confirmed susceptible;
   separately report the detection contrast if scientifically relevant.
2. **Mediator window:** Day 0-1, Day 0-2, and Day 0-3 definitions.
3. **Treatment definition:** active appropriate treatment versus newly started
   appropriate treatment; known-end versus open-ended prescription handling.
4. **Severity definition:** primary daily score versus alternative construction.
5. **Missing severity:** multiple imputation compatible with longitudinal
   structure versus complete case, with missingness tables by A, day, site, and Y.
6. **Site heterogeneity:** country fixed effects; site-cluster bootstrap or a
   hierarchical sensitivity model, subject to only six sites and sparse sites.
7. **Organism eligibility:** monomicrobial episodes only and major organism
   groups separately; these are exploratory if numbers are small.
8. **Microbiology information:** model `O_t` as a time-varying process or define
   it clearly as pre-decision information. Do not condition on observed future
   `O_t` values during counterfactual simulation.
9. **Model specification:** nonlinear continuous terms and selected interactions;
   compare simulated and observed daily treatment/severity distributions.
10. **Unmeasured mediator-outcome confounding:** present a quantitative bias or
    tipping-point analysis rather than claiming it is solved by g-computation.

### Secondary benchmark

Use an ACORN-HAI-style baseline propensity-score/IPW analysis to estimate the
total effect of resistance on mortality. This is a useful benchmark, but it does
not estimate longitudinal mediation and cannot replace the primary g-formula.

## Current v3 sensitivity analyses, decoded

- `lt_variant="alt"`: substitutes an alternative daily severity source. In the
  saved results it produced exactly the same estimates as the primary analysis,
  so the two inputs or retained values should be verified.
- `mediator_history="none"`: removes prior-M terms and uses only `M3` in the
  outcome model; this is not a pure lag-history sensitivity.
- `include_ot=TRUE`: intended to add daily microbiology-information variables,
  but used only 153 complete cases and dropped O0-O3 as invariant. It does not
  test the intended question.

These analyses should be replaced rather than described as the final sensitivity
analysis package.
