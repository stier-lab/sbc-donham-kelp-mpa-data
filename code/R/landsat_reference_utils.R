# Remove documented empty spatial slots, never zero biomass in occupied slots.
remove_empty_landsat_references <- function(data, corrections) {
  required <- c("MPA", "status", "rep")
  stopifnot(all(required %in% names(data)), all(required %in% names(corrections)))
  years <- grep("^X?[0-9]{4}$", names(data), value = TRUE)
  if (!length(years)) stop("Landsat annual columns are missing.")
  key <- function(x) paste(x$MPA, x$status, x$rep, sep = "|")
  data_keys <- key(data)
  correction_keys <- key(corrections)
  if (anyDuplicated(data_keys) || anyDuplicated(correction_keys) ||
      nrow(corrections) != 6L || any(corrections$status != "reference") ||
      any(corrections$rep != 2L)) stop("Unexpected empty-reference correction registry.")
  idx <- match(correction_keys, data_keys)
  if (anyNA(idx)) stop("A documented empty-reference slot is missing; inspect the source version.")
  values <- as.matrix(data[idx, years, drop = FALSE])
  if (any(!is.finite(values)) || any(values != 0)) {
    stop("A documented empty slot is no longer all zero; do not remove it automatically.")
  }
  data[-idx, , drop = FALSE]
}
