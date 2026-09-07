# Data-contract checks shared by the clinical longitudinal-mediation workflow.

assert_clinical_analysis_only <- function(variable_names) {
  forbidden <- variable_names[grepl("wgs|plasmid|carbapenemase|_simulated|_sim$", variable_names, ignore.case = TRUE)]
  if (length(forbidden)) {
    stop(
      "Clinical analysis object contains prohibited WGS/simulated fields: ",
      paste(forbidden, collapse = ", "),
      call. = FALSE
    )
  }
  invisible(TRUE)
}

validate_day_index <- function(panel, id_var = "subjid", day_var = "day") {
  required <- c(id_var, day_var)
  if (!all(required %in% names(panel))) {
    stop("Panel is missing day-index variables.", call. = FALSE)
  }
  if (anyDuplicated(panel[c(id_var, day_var)])) {
    stop("More than one record exists for a participant-day.", call. = FALSE)
  }
  if (any(!panel[[day_var]] %in% 0:3)) {
    stop("Clinical mediation panel must contain only Day 0--3.", call. = FALSE)
  }
  invisible(TRUE)
}

assert_v3_complete_risk_set <- function(panel, risk_var = "R") {
  if (!(risk_var %in% names(panel))) {
    stop("Panel is missing R_t (alive and under follow-up at the beginning of day t).", call. = FALSE)
  }
  r <- panel[[risk_var]]
  if (any(!is.na(r) & !(r %in% c(0, 1)))) {
    stop("R_t must be coded 0/1.", call. = FALSE)
  }
  if (any(is.na(r))) {
    stop("R_t is missing for at least one person-day; resolve risk-set status before running V3.", call. = FALSE)
  }
  if (any(r == 0)) {
    stop(
      "This V3 implementation does not yet simulate death/censoring within Day 0-3. ",
      "It must not silently complete-case exclude non-risk person-days; extend the survival/censoring model before using data with R_t = 0.",
      call. = FALSE
    )
  }
  invisible(TRUE)
}
