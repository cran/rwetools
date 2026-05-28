# TBL1 HELPERS- build_table1 & related functions ###########
#### table1 Package & Helper Functions ########

#' Format a numeric value as a fixed-point string
#'
#' @param x Numeric scalar to format.
#' @param digits Integer. Number of decimal places (default 3).
#' @param use_abs Logical. If TRUE, format the absolute value (default FALSE).
#' @return Character string, or "NA" when \code{x} is \code{NA}.
#' @keywords internal
fmt_num <- function(x, digits = 3, use_abs = FALSE) {
  if (is.na(x)) return("NA")
  if (use_abs) x <- abs(x)
  formatC(x, format = "f", digits = digits)
}

#' Format mean and standard deviation as "mean (SD)"
#'
#' @param m Numeric scalar. The mean value.
#' @param s Numeric scalar. The standard deviation.
#' @param mean_digits Integer. Decimal places for the mean (default 2).
#' @param sd_digits Integer. Decimal places for the SD (default 2).
#' @param use_abs Logical. If TRUE, use absolute values (default FALSE).
#' @return Character string in the form "mean (sd)", or "NA" when \code{m} is \code{NA}.
#' @keywords internal
fmt_mean_sd <- function(m, s, mean_digits = 2, sd_digits = 2, use_abs = FALSE) {
  if (is.na(m)) return("NA")
  if (use_abs) {
    m <- abs(m)
    s <- abs(s)
  }
  paste0(formatC(m, format = "f", digits = mean_digits), 
         " (", formatC(s, format = "f", digits = sd_digits), ")")
}

#' Format count and percentage as "n (pct%)"
#'
#' @param n Numeric scalar. The count (numerator).
#' @param d Numeric scalar. The denominator.
#' @param n_digits Integer. Decimal places for the count (default 0).
#' @param pct_digits Integer. Decimal places for the percentage (default 1).
#' @param use_abs Logical. If TRUE, use absolute values (default FALSE).
#' @return Character string in the form "n (pct%)", or "NA" when inputs are
#'   missing or \code{d == 0}.
#' @keywords internal
fmt_n_pct <- function(n, d, n_digits = 0, pct_digits = 1, use_abs = FALSE) {
  if (is.na(n) | is.na(d) | d == 0) return("NA")
  if (use_abs) {
    n <- abs(n)
    d <- abs(d)
  }
  pct <- 100 * n / d
  paste0(formatC(n, format = "f", digits = n_digits), 
         " (", format(round(pct, pct_digits), nsmall = pct_digits), "%)")
}

#' Format a proportion as a percentage string
#'
#' @param p Numeric scalar. A proportion (0 to 1 scale).
#' @param digits Integer. Decimal places for the percentage (default 1).
#' @param use_abs Logical. If TRUE, use absolute value (default FALSE).
#' @return Character string in the form "xx.x%", or "--" when \code{p}
#'   is \code{NA}.
#' @keywords internal
fmt_pct <- function(p, digits = 1, use_abs = FALSE) {
  if (is.na(p)) return("--")
  if (use_abs) p <- abs(p)
  sprintf(paste0("%.", digits, "f%%"), p * 100)
}

#' Canonize levels for consistent factor conversion
#'
#' Converts a vector to a factor with canonized levels. Logical vectors
#' become \code{factor("0", "1")}, numeric vectors are coerced via
#' \code{as.character}, and character vectors are trimmed.
#'
#' @param x A vector (logical, numeric, integer, character, or factor).
#' @return A factor with canonical levels.
#' @keywords internal
canonize_levels <- function(x) {
  if (is.logical(x)) {
    f <- factor(ifelse(is.na(x), NA, ifelse(x, "1", "0")), levels = c("0","1"))
  } else if (is.numeric(x) || is.integer(x)) {
    f <- factor(as.character(x))
  } else {
    f <- factor(trimws(as.character(x)))
  }
  f
}

#' Get non-baseline factor levels
#'
#' Returns all factor levels except those commonly used as a baseline
#' (e.g., "0", "No", "None", "FALSE").
#'
#' @param f A factor.
#' @return Character vector of non-baseline levels.
#' @keywords internal
levels_excl_zero <- function(f) {
  stopifnot(is.factor(f))
  base_like <- c("0","No","no","None","none","FALSE","False")
  setdiff(levels(f), base_like)
}

#' Classify variables as continuous, binary, or categorical
#'
#' Inspects one or more data frames to determine the type of each variable
#' using a combination of pattern-based rules and data-driven heuristics.
#'
#' @param df_list A data frame or list of data frames to inspect.
#' @param max_cat_levels Integer. Numeric variables with more than this many
#'   unique values are classified as continuous (default 12).
#' @return A named list with elements:
#'   \describe{
#'     \item{cat_vars}{Sorted character vector of categorical variable names.}
#'     \item{cat_vars_w_levels}{Named list mapping each categorical variable to
#'       its sorted unique levels.}
#'     \item{binary_vars}{Sorted character vector of binary (0/1) variable names.}
#'     \item{cont_vars}{Sorted character vector of continuous variable names.}
#'   }
#' @keywords internal
get_var_types <- function(df_list, max_cat_levels = 12) {
  
  # Ensure list input
  if (is.data.frame(df_list)) df_list <- list(df_list)
  if (!is.list(df_list)) stop("df_list must be a list of dataframes")
  
  # Collect all variable names across all dataframes
  unique_vars <- unique(unlist(lapply(df_list, names)))
  unique_vars_work <- unique_vars
  
  # Initialize outputs
  cat_vars_w_levels <- list()  # named list: var -> levels
  cat_vars <- character()       # just variable names
  cont_vars <- character()      # just variable names
  binary_vars <- character()      # just variable names
  
  #' ----------------------------- #
  # PASS 1: Pattern-based rules           
  #' ----------------------------- #
  for (var in unique_vars) {
    if (!(var %in% unique_vars_work)) next
    
    base_match_cat <- grepl("female|sex|gender|busin_data_type|business_type|data_type", var)
    suffix_match_cat <- grepl("_cat$", var) | grepl("_Cat$", var)
    
    matched_cat <- base_match_cat | suffix_match_cat
    matched_cont <- grepl("count|number", var)
    matched_binary <- grepl("osteoporosis_nofx_n|colonoscopy", var)
    
    if (matched_cat) {
      cat_vars <- c(cat_vars, var)
      cat_vars_w_levels[[var]] <- NA
      unique_vars_work <- setdiff(unique_vars_work, var)
      next
    }
    
    if (matched_cont) {
      cont_vars <- c(cont_vars, var)
      unique_vars_work <- setdiff(unique_vars_work, var)
      next
    }
    
    if (matched_binary) {
      binary_vars <- c(binary_vars, var)
      unique_vars_work <- setdiff(unique_vars_work, var)
      next
    }
    
    # Character / factor check
    is_char_factor <- FALSE
    for (df in df_list) {
      if (var %in% names(df)) {
        x <- df[[var]]
        if (is.character(x) || is.factor(x)) {
          is_char_factor <- TRUE
          break
        }
      }
    }
    
    if (is_char_factor) {
      cat_vars <- c(cat_vars, var)
      cat_vars_w_levels[[var]] <- NA
      unique_vars_work <- setdiff(unique_vars_work, var)
    }
  }
  
  #' ----------------------------- #
  # PASS 2: Continuous early-stop          
  #' ----------------------------- #
  
  seen_values <- vector("list", length(unique_vars_work))
  names(seen_values) <- unique_vars_work
  
  for (df in df_list) {
    vars_now <- intersect(unique_vars_work, names(df))
    
    for (var in vars_now) {
      x <- df[[var]]
      if (!(is.numeric(x) || is.integer(x))) next
      
      ux <- unique(stats::na.omit(x))
      seen_values[[var]] <- unique(c(seen_values[[var]], ux))
      
      if (length(seen_values[[var]]) > max_cat_levels) {
        cont_vars <- c(cont_vars, var)
        unique_vars_work <- setdiff(unique_vars_work, var)
      }
    }
  }
  
  #' ----------------------------- #
  # PASS 3: Remaining vars -> binary or categorical            
  #' ----------------------------- #
  for (var in unique_vars_work) {
    
    # reuse value in pass 2 if it exists
    if (!is.null(seen_values[[var]])) {
      uniq <- seen_values[[var]]
    } else {
      # else calculate it here
      all_values <- c()
      for (df in df_list) {
        if (var %in% names(df)) {
          x <- df[[var]]
          ux <- unique(stats::na.omit(x))
          all_values <- c(all_values, ux)
        }
      }
      uniq <- unique(all_values)
    }
    
    # skip variable with no data (NA) 
    if (length(uniq) == 0) {
      next 
    }
    # binary criteria
    else if (all(uniq %in% c(0,1))) {
      binary_vars <- c(binary_vars, var)
    } 
    # others are Categorical
    else {
      cat_vars <- c(cat_vars, var)
      cat_vars_w_levels[[var]] <- NA
    }
  }
  
  #' ----------------------------- #
  # PASS 4 (optional): fill categorical levels           
  #' ----------------------------- #
  for (var in names(cat_vars_w_levels)) {
    levels_all <- c()
    for (df in df_list) {
      if (var %in% names(df)) {
        x <- df[[var]]
        if (is.factor(x)) {
          levels_all <- c(levels_all, levels(x))
        } else {
          levels_all <- c(levels_all, unique(stats::na.omit(x)))
        }
      }
    }
    cat_vars_w_levels[[var]] <- sort(unique(levels_all))
  }
  
  #' ----------------------------- #
  # Final output          
  #' ----------------------------- #
  list(
    cat_vars = sort(unique(cat_vars)),
    cat_vars_w_levels = cat_vars_w_levels,
    binary_vars = sort(unique(binary_vars)),
    cont_vars = sort(unique(cont_vars))
  )
}

#' Combine named lists, keeping unique values
#'
#' Merges multiple named lists by name, concatenating and de-duplicating the
#' values for each shared key. Useful for combining \code{cat_vars_w_levels}
#' outputs from multiple calls to \code{\link{get_var_types}}.
#'
#' @param ... Named lists to combine.
#' @return A single named list with unique, combined values per key.
#' @keywords internal
combine_named_list <- function(...) {
  lsts <- list(...)
  all_names <- unique(unlist(lapply(lsts, names)))
  out <- stats::setNames(vector("list", length(all_names)), all_names)
  for (nm in all_names) {
    vals <- lapply(lsts, function(x) x[[nm]])
    vals <- vals[!sapply(vals, is.null)]
    out[[nm]] <- unique(unlist(vals))
  }
  out
}

#' Compute weighted mean, SD, and total weight from a survey design
#'
#' @param des A \code{survey.design} object (from \pkg{survey}).
#' @param var_name Character string. Name of the variable in \code{des$variables}.
#' @return A named numeric vector with elements \code{mean}, \code{sd}, and
#'   \code{sum_w}. Returns \code{c(mean = NA, sd = NA, sum_w = 0)} when all
#'   values are \code{NA}.
#' @keywords internal
get_weighted_stats <- function(des, var_name) {
  var_formula <- stats::as.formula(paste0("~", var_name))
  x <- des$variables[[var_name]]
  if (all(is.na(x))) return(c(mean = NA, sd = NA, sum_w = 0))
  
  m <- suppressWarnings(as.numeric(survey::svymean(var_formula, design = des, na.rm = TRUE)))
  v <- suppressWarnings(as.numeric(survey::svyvar(var_formula, design = des, na.rm = TRUE)))
  s <- if (is.na(v)) NA else sqrt(v)
  w_all <- suppressWarnings(stats::weights(des, type = "sampling"))
  sw <- sum(w_all[!is.na(x)], na.rm = TRUE)
  c(mean = m, sd = s, sum_w = sw)
}

#' Compute weighted proportion for a single factor level
#'
#' @param design A \code{survey.design} object (from \pkg{survey}).
#' @param var Character string. Name of the categorical variable.
#' @param level Character string. The factor level whose proportion is needed.
#' @return Numeric scalar: the weighted proportion, 0 if the level is absent,
#'   or \code{NA_real_} on error.
#' @keywords internal
get_weighted_prop <- function(design, var, level) {
  tryCatch({
    tab <- survey::svytable(stats::as.formula(paste0("~", var)), design)
    if (level %in% names(tab)) {
      prop <- tab[[level]] / sum(tab)
      return(as.numeric(prop))
    } else {
      return(0)
    }
  }, error = function(e) {
    return(NA_real_)
  })
}

#' Process a single categorical or binary variable for Table 1
#'
#' Computes weighted counts, percentages, crude differences, and standardized
#' mean differences for each level of a categorical or binary variable.
#' Uses \code{stats::xtabs} to compute the full weighted cross-tabulation
#' once per variable, then extracts per-level results efficiently.
#'
#' @param v Character string. Variable name.
#' @param dat Data frame containing \code{v}, a \code{.grp} column (0/1),
#'   and optionally a \code{..w} weights column.
#' @param levels_to_show Character. Either \code{"all"} (show every level) or
#'   \code{"non_zero"} (show only non-baseline levels).
#' @param predefined_levels Character vector of levels to display, or \code{NULL}
#'   to use levels found in the data.
#' @param N_all Integer. Total number of patients.
#' @param N_ref Integer. Number of reference-group patients.
#' @param N_tx Integer. Number of exposed-group patients.
#' @param idx_ref Integer vector. Row indices for the reference group.
#' @param idx_tx Integer vector. Row indices for the exposed group.
#' @param w_vec Numeric vector. Pre-computed weight vector for all rows.
#' @param grp_vec Integer vector. Pre-computed group vector (0/1) for all rows.
#' @param using_w Logical. Whether weights are active (any weight > 1).
#' @param n_decimal Integer. Decimal places for counts (default 0).
#' @param pct_decimal Integer. Decimal places for percentages (default 1).
#' @param use_abs Logical. Use absolute values for differences (default FALSE).
#' @param calculate_existing Logical. If TRUE, compute non-missing counts
#'   (default FALSE).
#' @param verbose Logical. Print progress messages (default TRUE).
#' @return A matrix of list rows (one per level), or \code{NULL} if the
#'   variable is not in \code{dat} or has no displayable levels.
#' @keywords internal
process_catbin_direct <- function(v, dat, levels_to_show = c("all", "non_zero"), 
                                  predefined_levels = NULL,
                                  N_all, N_ref, N_tx,
                                  idx_ref, idx_tx,
                                  w_vec, grp_vec,
                                  using_w,
                                  n_decimal = 0, pct_decimal = 1, 
                                  use_abs = FALSE,
                                  calculate_existing = FALSE,
                                  verbose = TRUE) {
  levels_to_show <- match.arg(levels_to_show)
  if (!v %in% names(dat)) return(NULL)
  
  if (verbose) message(sprintf("  Processing %s variable: %s", 
                               ifelse(levels_to_show == "non_zero", "binary", "categorical"), v))
  
  # ---- Per-variable computation (ONCE) ----
  
  # Canonize levels
  x_all_full <- canonize_levels(dat[[v]])
  lv_all <- levels(x_all_full)
  
  # Determine which levels to show
  if (levels_to_show == "all") {
    if (!is.null(predefined_levels)) {
      lv_show <- predefined_levels
      if (verbose) message(sprintf("    Using predefined levels: %s", paste(lv_show, collapse = ", ")))
    } else {
      lv_show <- lv_all
    }
  } else {
    lv_show <- levels_excl_zero(factor(x_all_full, levels = lv_all))
  }
  
  if (length(lv_show) == 0) return(NULL)
  
  # Non-missing mask
  nonmiss <- !is.na(x_all_full)
  
  # Missing counts — computed ONCE per variable (same for all levels)
  miss_all <- sum(!nonmiss)
  miss_ref <- sum(!nonmiss[idx_ref])
  miss_tx  <- sum(!nonmiss[idx_tx])
  
  # Existing counts — computed ONCE per variable if requested
  existing_all <- existing_ref <- existing_tx <- NULL
  if (calculate_existing) {
    existing_all <- sum(nonmiss)
    existing_ref <- sum(nonmiss[idx_ref])
    existing_tx  <- sum(nonmiss[idx_tx])
  }
  
  # Weighted cross-tabulation — computed ONCE per variable
  # as.character() ensures integer 0/1 matches factor levels "0"/"1"
  grp_factor <- factor(as.character(grp_vec[nonmiss]), levels = c("0", "1"))
  tab_w <- stats::xtabs(w_vec[nonmiss] ~ grp_factor + x_all_full[nonmiss])
  
  # Denominators from cross-tab (sum across all levels per group)
  den_ref <- sum(tab_w["0", ])
  den_tx  <- sum(tab_w["1", ])
  den_all <- den_ref + den_tx
  
  # ---- Per-level loop (now just table lookups) ----
  do.call(rbind, lapply(lv_show, function(L) {
    
    # Look up weighted counts from pre-computed cross-tab
    if (L %in% colnames(tab_w)) {
      num_ref <- tab_w["0", L]
      num_tx  <- tab_w["1", L]
    } else {
      num_ref <- 0
      num_tx  <- 0
    }
    num_all <- num_ref + num_tx
    
    # Calculate proportions
    p_all <- if (den_all > 0) num_all / den_all else NA_real_
    p_ref <- if (den_ref > 0) num_ref / den_ref else NA_real_
    p_tx  <- if (den_tx  > 0) num_tx  / den_tx  else NA_real_
    
    # Calculate crude difference and SMD
    crude_diff_pp <- 100 * (p_tx - p_ref)
    denom_smd <- sqrt((p_tx * (1 - p_tx) + p_ref * (1 - p_ref)) / 2)
    smd <- if (is.na(denom_smd) || denom_smd == 0) NA_real_ else (p_tx - p_ref) / denom_smd
    
    # Format output strings
    if (using_w) {
      Total_str <- paste0(formatC(round(num_all), format = "f", digits = n_decimal), 
                          " (", format(round(100 * p_all, pct_decimal), nsmall = pct_decimal), "%)")
      Tx_str    <- paste0(formatC(round(num_tx), format = "f", digits = n_decimal), 
                          " (", format(round(100 * p_tx, pct_decimal), nsmall = pct_decimal), "%)")
      Ref_str   <- paste0(formatC(round(num_ref), format = "f", digits = n_decimal), 
                          " (", format(round(100 * p_ref, pct_decimal), nsmall = pct_decimal), "%)")
    } else {
      Total_str <- fmt_n_pct(round(num_all), round(den_all), n_decimal, pct_decimal)
      Tx_str    <- fmt_n_pct(round(num_tx), round(den_tx), n_decimal, pct_decimal)
      Ref_str   <- fmt_n_pct(round(num_ref), round(den_ref), n_decimal, pct_decimal)
    }
    
    # Variable name formatting
    var_name_str <- if (levels_to_show == "all") {
      paste0(v, " (= ", L, ")")
    } else {
      v
    }
    
    var_type_str <- if (levels_to_show == "all") {
      "Categorical"
    } else {
      "Binary"
    }
    
    # Build return list
    result <- list(
      Variable = var_name_str,
      Type = var_type_str,
      Total = Total_str,
      Exp = Tx_str,
      Ref = Ref_str,
      Crude_diff = fmt_num(as.numeric(crude_diff_pp), digits = pct_decimal, use_abs = use_abs),
      Std_diff = fmt_num(as.numeric(smd), digits = 3, use_abs = use_abs),
      Missing_Total_N_Pct = fmt_n_pct(miss_all, N_all, n_decimal, pct_decimal),
      Missing_Exp_N_Pct = fmt_n_pct(miss_tx, N_tx, n_decimal, pct_decimal),
      Missing_Ref_N_Pct = fmt_n_pct(miss_ref, N_ref, n_decimal, pct_decimal)
    )
    
    # Add existing counts if calculated
    if (calculate_existing) {
      result$Existing_Total_N <- formatC(existing_all, format = "f", digits = n_decimal)
      result$Existing_Exp_N <- formatC(existing_tx, format = "f", digits = n_decimal)
      result$Existing_Ref_N <- formatC(existing_ref, format = "f", digits = n_decimal)
    }
    
    result
  }))
}
