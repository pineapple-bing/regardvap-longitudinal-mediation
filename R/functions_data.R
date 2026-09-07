# Shared, deliberately small helpers for the REGARD-VAP analysis scripts.

stop_if_missing <- function(data, variables, object_name = deparse(substitute(data))) {
  missing <- setdiff(variables, names(data))
  if (length(missing)) stop(object_name, " is missing: ", paste(missing, collapse = ", "), call. = FALSE)
}

as_binary <- function(x, name) {
  value <- suppressWarnings(as.numeric(x))
  invalid <- !is.na(value) & !(value %in% c(0, 1))
  if (any(invalid)) stop(name, " must contain only 0, 1, or NA.", call. = FALSE)
  value
}

any_observed_one <- function(...) {
  values <- cbind(...)
  ifelse(rowSums(!is.na(values)) == 0, NA_real_, as.numeric(rowSums(values == 1, na.rm = TRUE) > 0))
}

write_csv_safely <- function(x, file) {
  dir.create(dirname(file), recursive = TRUE, showWarnings = FALSE)
  utils::write.csv(x, file, row.names = FALSE, na = "")
  message("Wrote: ", file)
}
