# Internal exposure-value handling ------------------------------------------

# Convert the caller-selected exposure and reference values to the package's
# internal 1/0 representation.  This is deliberately unconditional: a data
# column containing 0/1 must still be flipped when exp_value = 0 and
# ref_value = 1.
.canonicalize_exposure <- function(x, exp_value, ref_value,
                                   exposure_var = "exposure",
                                   require_both = TRUE,
                                   warn_unmatched = TRUE) {
  if (length(exp_value) != 1L || is.na(exp_value) ||
      length(ref_value) != 1L || is.na(ref_value)) {
    stop("exp_value and ref_value must each be a single non-missing value")
  }

  exp_key <- as.character(exp_value)
  ref_key <- as.character(ref_value)
  if (identical(exp_key, ref_key)) {
    stop("exp_value and ref_value must differ")
  }

  x_key <- as.character(x)
  is_exp <- !is.na(x) & x_key == exp_key
  is_ref <- !is.na(x) & x_key == ref_key

  if (require_both && !any(is_exp)) {
    stop(paste("exp_value", exp_value, "not found in", exposure_var))
  }
  if (require_both && !any(is_ref)) {
    stop(paste("ref_value", ref_value, "not found in", exposure_var))
  }

  unmatched <- !is.na(x) & !is_exp & !is_ref
  n_unmatched <- sum(unmatched)
  if (warn_unmatched && n_unmatched > 0L) {
    warning(sprintf(
      "%d observations had values other than exp_value or ref_value and were set to NA",
      n_unmatched
    ))
  }

  value <- rep(NA_integer_, length(x))
  value[is_ref] <- 0L
  value[is_exp] <- 1L

  list(value = value, n_unmatched = n_unmatched)
}
