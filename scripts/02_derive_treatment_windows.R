# Reusable helpers for date-only daily treatment data.
# These functions deliberately avoid calling Day 0-1 an exact 48-hour window.

validate_binary <- function(x, name) {
  bad <- !is.na(x) & !(x %in% c(0, 1))
  if (any(bad)) stop(name, " contains values outside 0/1/NA")
  as.numeric(x)
}

any_observed_one <- function(...) {
  x <- cbind(...)
  ifelse(rowSums(!is.na(x)) == 0, NA_real_, as.numeric(rowSums(x == 1, na.rm = TRUE) > 0))
}

derive_treatment_windows <- function(wide) {
  required <- paste0("M", 0:3)
  missing <- setdiff(required, names(wide))
  if (length(missing)) stop("Missing daily mediator columns: ", paste(missing, collapse = ", "))
  for (v in required) wide[[v]] <- validate_binary(wide[[v]], v)

  wide$appropriate_by_end_day1 <- any_observed_one(wide$M0, wide$M1)
  wide$appropriate_by_end_day2 <- any_observed_one(wide$M0, wide$M1, wide$M2)
  wide$appropriate_by_end_day3 <- any_observed_one(wide$M0, wide$M1, wide$M2, wide$M3)
  wide$appropriate_days_day0_3 <- rowSums(wide[required], na.rm = TRUE)
  wide$appropriate_days_day0_3[rowSums(!is.na(wide[required])) == 0] <- NA_real_
  wide
}
