# Danny data handoff checklist

## Is it appropriate to send data now?

**Yes, for data and method review; no, as a final causal-results package.**
The current priority is to have Danny verify the clinical data structure,
treatment definition, and timing. The g-formula estimates are still under
methods review and should not be presented as final.

## Before transfer

- Confirm Danny is covered by the study's ethics approval, data-sharing
  agreement, and access permissions.
- Use an approved institutional storage/transfer route, not personal email or a
  public GitHub repository.
- Remove direct identifiers and any field unnecessary for the agreed review.
- Exclude all `*_sim` fields, simulated WGS/plasmid variables, raw WGS files,
  and unvalidated derived results.
- Record the extract date, row count, variable list, and a data version label.

## Recommended review package

1. **Clinical analysis extract**: one row per participant-day for Day 0--3,
   plus a patient-level Day-60 outcome file, if this level of access is
   approved.
2. **Data dictionary**: variable name, definition, source, units/coding,
   measurement time, and whether it is available before `M_t`.
3. **Study timeline and DAG**: the current editable figures, labelled as draft.
4. **Specific questions for Danny**:
   - Is Day 0 clinically defensible as time zero?
   - Does each daily severity field precede the treatment decision?
   - Is the “appropriate active treatment” definition clinically correct?
   - Which microbiology information is actually available before each decision?
   - Which episodes should define the susceptible comparison group?

## Suggested message

> I am sharing a de-identified clinical analysis extract for review of the
> longitudinal data structure and treatment timing. The current mediation
> estimates are preliminary and should not be interpreted as final. Simulated
> WGS fields and WGS analyses are excluded. I would especially value your view
> on time zero, pre-treatment daily severity, the definition of appropriate
> active treatment, and the eligible comparison population.
