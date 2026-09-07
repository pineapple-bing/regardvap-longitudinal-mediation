# v3 code and output map

The previous v3 code was a 995-line monolithic runner. Its complete clinical
logic is retained in `04_longitudinal_mediation/01_run_regardvap_longitudinal_mediation_v3_no_wgs.R`, with the simulated WGS code removed.

| Manuscript section | v3 outputs retained | Code owner |
| --- | --- | --- |
| 1. Study population and baseline | input audit, inclusion/exclusion checks, Table 1 | data preparation + baseline modules |
| 2. Longitudinal trajectories | timeline, DAG, Table 2, observed trajectory figure | trajectory module |
| 3. Longitudinal mediation | counterfactual risks, TE/IDE/IIE, sensitivity table, effect-decomposition figure | g-formula module |
| 4. Heterogeneity | severity, organism group, and country analyses | HTE module |

Excluded from this clinical repository: all `*_sim` WGS variables, plasmid-cluster analyses, carbapenemase-gene analyses, and WGS-derived figures/tables.
