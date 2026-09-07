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
