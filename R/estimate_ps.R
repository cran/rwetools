#  PS Calculate and Add PS to Dataset ################ 
#' Calculate propensity scores and add them to the dataset
#'
#' This function calculates propensity scores using logistic regression and adds them 
#' as a new column to the dataset. It can output both an R data frame and/or CSV file.
#' Additionally, it can save odds ratio table from the PS model to Excel or data frame.
#'
#' @param in_df Data frame containing the input data (optional if in_csvpath provided)
#' @param in_csvpath Character string. Path to input CSV file (optional if in_df provided)
#' @param out_csvpath Character string. Path for output CSV file (optional, if NULL no CSV is saved)
#' @param out_xlsxpath_odds_ratio Character string. Path for output Excel file with OR table (optional)
#' @param exposure_var Character string. Name of the binary exposure/treatment variable column (default: "exposure")
#' @param exp_value Value representing the exposed/treated group (default: 1)
#' @param ref_value Value representing the reference/control group (default: 0)
#' @param class_vars Character vector. Names of categorical/factor variables to include in PS model
#' @param cont_vars Character vector. Names of continuous variables to include in PS model
#' @param ps_var Character string. Name for the calculated PS variable column (default: "ps")
#' @param interactions Character string. Interaction terms to include (e.g., "var1:var2 + var1:var3")
#' @param exclude_vars_w_extreme_distribution Logical. If TRUE, automatically excludes variables
#'   with extreme distribution (categorical: only 1 unique level in either group, or any
#'   level with count < 5 in either group; continuous: constant/single unique value
#'   (SD < 1e-6) in either group, or SMD > 1.5 between groups)
#'   instead of throwing error (default: FALSE). The categorical screen unions levels
#'   across the two exposure groups, so a level present in one group but absent in the
#'   other is counted as n = 0 and flagged by the cnt < 5 rule.
#' @param verbose Logical. Print progress messages (default TRUE).
#'
#' @return A data frame with the PS column added. Also saves to CSV if path specified.
#'
#' @section Side Effects:
#' \itemize{
#'   \item Writes a CSV file when \code{out_csvpath} is provided.
#'   \item Writes an Excel file with odds-ratio table when
#'     \code{out_xlsxpath_odds_ratio} is provided.
#' }
#'
#' @export
#'
#' @examples
#' csv_path <- system.file("extdata", "sample_data.csv", package = "rwetools")
#' df <- read.csv(csv_path)
#'
#' # Basic usage: estimate PS and add it as a new column
#' result <- estimate_ps(
#'   in_df        = df,
#'   exposure_var = "exposure",
#'   class_vars   = c("cat1", "cat2", "cat3", "cat4"),
#'   cont_vars    = c("cont1", "cont2", "cont3"),
#'   verbose      = FALSE
#' )
#' head(result$ps)
#'
#' \donttest{
#' # With CSV output and Excel OR table (requires openxlsx)
#' if (requireNamespace("openxlsx", quietly = TRUE)) {
#'   out_csv  <- tempfile(fileext = ".csv")
#'   out_xlsx <- tempfile(fileext = ".xlsx")
#'   result2 <- estimate_ps(
#'     in_df                               = df,
#'     out_csvpath                         = out_csv,
#'     out_xlsxpath_odds_ratio             = out_xlsx,
#'     exposure_var                        = "exposure",
#'     class_vars                          = c("cat1", "cat2"),
#'     cont_vars                           = c("cont1", "cont2"),
#'     exclude_vars_w_extreme_distribution = TRUE,
#'     verbose                             = FALSE
#'   )
#' }
#' }
estimate_ps <- function(
    in_df = NULL,
    in_csvpath = NULL,
    out_csvpath = NULL,
    out_xlsxpath_odds_ratio = NULL,
    exposure_var = "exposure",
    exp_value = 1,
    ref_value = 0,
    class_vars = NULL,
    cont_vars = NULL,
    ps_var = "ps",
    interactions = NULL,
    exclude_vars_w_extreme_distribution = FALSE,
    verbose = TRUE) {
  
  if (is.null(in_df) && is.null(in_csvpath)) {
    stop("Either in_df or in_csvpath must be provided")
  }
  if (!is.null(in_df) && !is.null(in_csvpath)) {
    warning("Both in_df and in_csvpath provided. Using in_df.")
  }
  if (is.null(in_df)) {
    in_df <- utils::read.csv(in_csvpath, stringsAsFactors = FALSE)
  }
  
  # Check required packages
  if (!is.null(out_xlsxpath_odds_ratio)) {
    if (!requireNamespace("openxlsx", quietly = TRUE)) {
      stop("Package 'openxlsx' is required for Excel output but not installed.")
    }
  }
  
  # Validate inputs
  if (!is.data.frame(in_df)) {
    stop("in_df must be a data frame")
  }
  
  if (!exposure_var %in% names(in_df)) {
    stop(paste("Exposure variable", exposure_var, "not found in the dataset"))
  }
  
  if (is.null(class_vars) && is.null(cont_vars)) {
    stop("At least one of class_vars or cont_vars must be specified")
  }
  if (exp_value == ref_value) {
    stop("exp_value and ref_value must be different")
  }
  
  # Check if variables exist in dataset
  all_vars <- c(class_vars, cont_vars)
  missing_vars <- setdiff(all_vars, names(in_df))
  if (length(missing_vars) > 0) {
    stop(paste("Variables not found in dataset:", paste(missing_vars, collapse = ", ")))
  }
  
  if (verbose) message("\n========================================")
  if (verbose) message("PROPENSITY SCORE CALCULATION")
  if (verbose) message("========================================")
  
  # Make a copy of the data
  data_work <- in_df
  
  # *** Step 0: Check for extreme distribution variables ***
  if (verbose) message("Step 0: Checking variables for extreme distribution")
  if (verbose) message("-------------------------------------------------------------------")
  
  # Save original input variables (before any exclusion)
  input_class_vars <- class_vars
  input_cont_vars  <- cont_vars
  
  # Track excluded variables
  excluded_class_vars <- character(0)
  excluded_cont_vars <- character(0)
  
  # Split data by exposure group for distribution checks
  exp_df <- data_work[data_work[[exposure_var]] == exp_value, , drop = FALSE]
  ref_df <- data_work[data_work[[exposure_var]] == ref_value, , drop = FALSE]
  
  # --- Check categorical variables for extreme distribution ---
  if (!is.null(class_vars) && length(class_vars) > 0) {
    extreme_class <- character(0)
    
    for (v in class_vars) {
      # Rule 0 (Single-level / constant): only 1 unique non-NA value in either group
      n_levels_exp <- length(unique(stats::na.omit(exp_df[[v]])))
      n_levels_ref <- length(unique(stats::na.omit(ref_df[[v]])))
      if (n_levels_exp <= 1 || n_levels_ref <= 1) {
        extreme_class <- c(extreme_class, v)
        next
      }
      
      # Rule 1 (Sparse level): any level with count < 5 in either group.
      # Union levels across exp/ref so zero-count cells (level present in one
      # group but absent in the other) are visible to the cnt<5 check.
      # (Without the union, table() silently drops missing levels and 0-count
      # cells can sneak past the screen -> quasi-separation in glm.)
      vals_e   <- exp_df[[v]][!is.na(exp_df[[v]])]
      vals_r   <- ref_df[[v]][!is.na(ref_df[[v]])]
      all_lvls <- sort(unique(c(vals_e, vals_r)))
      count_exp <- as.integer(table(factor(vals_e, levels = all_lvls)))
      count_ref <- as.integer(table(factor(vals_r, levels = all_lvls)))
      flag_exp <- any(count_exp < 5)
      flag_ref <- any(count_ref < 5)
      
      if (flag_exp || flag_ref) {
        extreme_class <- c(extreme_class, v)
      }
    }
    
    if (length(extreme_class) > 0) {
      if (exclude_vars_w_extreme_distribution) {
        if (verbose) message("[OK] Auto-excluding the following class variables with extreme distribution:")
        if (verbose) message("  ", paste(extreme_class, collapse = ", "))
        if (verbose) message("  (Reason: only 1 unique level in either group, or any level with count < 5 in either group)")
        excluded_class_vars <- extreme_class
        class_vars <- setdiff(class_vars, extreme_class)
      } else {
        if (verbose) message("[X] ERROR: The following class variables have extreme distribution:")
        if (verbose) message("  ", paste(extreme_class, collapse = ", "))
        if (verbose) message("  (Only 1 unique level in either group, or any level with count < 5 in either exposure group)")
        if (verbose) message("This may cause PS model instability.")
        if (verbose) message("Please remove these variables from class_vars before calling this function,")
        if (verbose) message("or set exclude_vars_w_extreme_distribution = TRUE to automatically exclude them.")
        stop("PS model cannot be fitted: class variables with extreme distribution detected")
      }
    } else {
      if (verbose) message("[OK] All class variables pass extreme distribution check")
    }
  }
  
  # --- Check continuous variables for extreme distribution ---
  if (!is.null(cont_vars) && length(cont_vars) > 0) {
    extreme_cont <- character(0)
    
    for (v in cont_vars) {
      vals_e <- exp_df[[v]][!is.na(exp_df[[v]])]
      vals_r <- ref_df[[v]][!is.na(ref_df[[v]])]
      
      sd_e  <- stats::sd(vals_e)
      sd_r  <- stats::sd(vals_r)
      
      # Rule 1 (Constant): SD < 1e-6 in either group (effectively a constant)
      # Also flag if SD is NA (fewer than 2 non-NA observations)
      rule1 <- is.na(sd_e) || is.na(sd_r) || (sd_e < 1e-6) || (sd_r < 1e-6)
      
      # Rule 2 (Separation): SMD > 1.5 (extreme mean difference between groups)
      if (rule1) {
        rule2 <- FALSE  # already caught by rule1
      } else {
        pooled_sd <- sqrt((sd_e^2 + sd_r^2) / 2)
        if (pooled_sd < 1e-6) {
          rule2 <- FALSE
        } else {
          smd <- abs(mean(vals_e) - mean(vals_r)) / pooled_sd
          rule2 <- (smd > 1.5)
        }
      }
      
      if (rule1 || rule2) {
        extreme_cont <- c(extreme_cont, v)
      }
    }
    
    if (length(extreme_cont) > 0) {
      if (exclude_vars_w_extreme_distribution) {
        if (verbose) message("[OK] Auto-excluding the following continuous variables with extreme distribution:")
        if (verbose) message("  ", paste(extreme_cont, collapse = ", "))
        if (verbose) message("  (Reason: constant/single unique value (SD < 1e-6) in either group, or SMD > 1.5 between groups)")
        excluded_cont_vars <- extreme_cont
        cont_vars <- setdiff(cont_vars, extreme_cont)
      } else {
        if (verbose) message("[X] ERROR: The following continuous variables have extreme distribution:")
        if (verbose) message("  ", paste(extreme_cont, collapse = ", "))
        if (verbose) message("  (Constant/single unique value (SD < 1e-6) in either group, or SMD > 1.5 between groups)")
        if (verbose) message("This may cause PS model instability.")
        if (verbose) message("Please remove these variables from cont_vars before calling this function,")
        if (verbose) message("or set exclude_vars_w_extreme_distribution = TRUE to automatically exclude them.")
        stop("PS model cannot be fitted: continuous variables with extreme distribution detected")
      }
    } else {
      if (verbose) message("[OK] All continuous variables pass extreme distribution check")
    }
  }
  
  # Check if at least one valid variable remains after exclusion
  total_valid_vars <- length(c(class_vars, cont_vars))
  if (total_valid_vars == 0) {
    stop("No valid variables remain for PS model after excluding extreme-distribution variables.\n",
         "At least one variable with sufficient variation is required.")
  }
  
  if (exclude_vars_w_extreme_distribution && (length(excluded_class_vars) > 0 || length(excluded_cont_vars) > 0)) {
    if (verbose) message(sprintf("\n[OK] Proceeding with %d valid variable(s) for PS model", total_valid_vars))
    if (length(excluded_class_vars) > 0) {
      if (verbose) message(sprintf("  Excluded %d class variable(s): %s", 
                                   length(excluded_class_vars), paste(excluded_class_vars, collapse = ", ")))
    }
    if (length(excluded_cont_vars) > 0) {
      if (verbose) message(sprintf("  Excluded %d continuous variable(s): %s", 
                                   length(excluded_cont_vars), paste(excluded_cont_vars, collapse = ", ")))
    }
  }
  
  if (verbose) message("")
  
  # Step 1: Recode exposure variable to 0/1 if needed          ######################  
  unique_exposure_vals <- unique(data_work[[exposure_var]])
  
  if (!(all(unique_exposure_vals %in% c(0, 1)))) {
    if (verbose) message("Step 1: Recoding exposure variable")
    if (verbose) message("----------------------------------------")
    if (verbose) message(sprintf("Original exposure values: %s", paste(unique_exposure_vals, collapse = ", ")))
    if (verbose) message(sprintf("Recoding: %s (exposed) -> 1, %s (reference) -> 0", exp_value, ref_value))
    
    # Check if exp_value and ref_value exist
    if (!(exp_value %in% unique_exposure_vals)) {
      stop(paste("exp_value", exp_value, "not found in", exposure_var))
    }
    if (!(ref_value %in% unique_exposure_vals)) {
      stop(paste("ref_value", ref_value, "not found in", exposure_var))
    }
    
    # Recode
    data_work[[exposure_var]] <- ifelse(data_work[[exposure_var]] == exp_value, 1,
                                        ifelse(data_work[[exposure_var]] == ref_value, 0, NA))
    
    # Check for NAs created
    n_na <- sum(is.na(data_work[[exposure_var]]))
    if (n_na > 0) {
      warning(paste(n_na, "observations had values other than exp_value or ref_value and were set to NA"))
    }
    
    # Show recoded distribution
    recoded_table <- table(data_work[[exposure_var]], useNA = "ifany")
    if (verbose) message("\nRecoded exposure distribution:")
    if (verbose) print(recoded_table)
    if (verbose) message("")
  }
  
  # Step 2: Build PS model formula
  if (verbose) message("Step 2: Building propensity score model")
  if (verbose) message("----------------------------------------")
  
  formula_parts <- c()
  
  # Add continuous variables
  if (!is.null(cont_vars) && length(cont_vars) > 0) {
    formula_parts <- c(formula_parts, cont_vars)
    if (verbose) message(sprintf("Continuous variables (%d): %s", 
                                 length(cont_vars), 
                                 paste(cont_vars, collapse = ", ")))
  }
  
  # Add categorical variables
  if (!is.null(class_vars) && length(class_vars) > 0) {
    # Convert to factors if not already
    for (var in class_vars) {
      if (!is.factor(data_work[[var]])) {
        data_work[[var]] <- as.factor(data_work[[var]])
      }
    }
    formula_parts <- c(formula_parts, class_vars)
    if (verbose) message(sprintf("Class (binary or categorical) variables (%d): %s", 
                                 length(class_vars), 
                                 paste(class_vars, collapse = ", ")))
  }
  
  # Add interactions
  if (!is.null(interactions) && interactions != "") {
    formula_parts <- c(formula_parts, interactions)
    if (verbose) message(sprintf("Interaction terms: %s", interactions))
  }
  
  # Create formula
  if (length(formula_parts) > 0) {
    formula_str <- paste(exposure_var, "~", paste(formula_parts, collapse = " + "))
  } else {
    stop("No variables specified for PS model")
  }
  
  if (verbose) message(sprintf("\nModel formula: %s\n", formula_str))
  formula_obj <- stats::as.formula(formula_str)
  
  # Step 3: Fit logistic regression model
  if (verbose) message("Step 3: Fitting logistic regression model")
  if (verbose) message("----------------------------------------")
  
  ps_model <- stats::glm(formula_obj, 
                         data = data_work, 
                         family = stats::binomial(link = "logit"))
  
  # Calculate propensity scores
  data_work[[ps_var]] <- stats::predict(ps_model, type = "response")
  
  # Calculate C-statistic
  c_stat_result <- calc_c_statistic(data_work[[exposure_var]], data_work[[ps_var]], label = "PS Model", verbose = verbose)
  c_stat <- if (!is.null(c_stat_result)) c_stat_result$c_stat else NA_real_
  
  if (verbose) message(sprintf("Model converged: %s", ps_model$converged))
  if (verbose) message(sprintf("Number of observations: %d", nrow(stats::model.frame(ps_model))))
  if (verbose) message(sprintf("C-statistic: %.3f", c_stat))
  
  # Show PS distribution summary
  ps_summary <- summary(data_work[[ps_var]])
  if (verbose) message("\nPropensity score distribution:")
  if (verbose) print(ps_summary)
  
  # Check for extreme PS values
  extreme_low <- sum(data_work[[ps_var]] < 0.01, na.rm = TRUE)
  extreme_high <- sum(data_work[[ps_var]] > 0.99, na.rm = TRUE)
  if (extreme_low > 0 || extreme_high > 0) {
    if (verbose) message(sprintf("\nWarning: %d observations with PS < 0.01, %d with PS > 0.99", 
                                 extreme_low, extreme_high))
  }
  
  #  Step 3.5: Create Odds Ratio Table (if requested)              ######################
  if (!is.null(out_xlsxpath_odds_ratio)) {
    if (verbose) message("\nStep 3.5: Creating Odds Ratio Table")
    if (verbose) message("----------------------------------------")
    
    # Extract coefficients and create odds ratio table
    coef_summary <- summary(ps_model)$coefficients
    
    # Calculate odds ratios and confidence intervals
    or_table <- data.frame(
      Variable = rownames(coef_summary),
      Coefficient = coef_summary[, "Estimate"],
      SE = coef_summary[, "Std. Error"],
      z_value = coef_summary[, "z value"],
      p_value = coef_summary[, "Pr(>|z|)"],
      stringsAsFactors = FALSE
    )
    or_table$Odds_Ratio    <- exp(or_table$Coefficient)
    or_table$OR_CI_Lower   <- exp(or_table$Coefficient - 1.96 * or_table$SE)
    or_table$OR_CI_Upper   <- exp(or_table$Coefficient + 1.96 * or_table$SE)
    or_table$OR_95CI       <- sprintf("%.2f (%.2f, %.2f)", or_table$Odds_Ratio, or_table$OR_CI_Lower, or_table$OR_CI_Upper)
    or_table$p_formatted   <- ifelse(or_table$p_value < 0.001, "<0.001", sprintf("%.3f", or_table$p_value))
    
    # Format for display (remove intercept for cleaner output)
    ps_model_results <- or_table[or_table$Variable != "(Intercept)",
                                 c("Variable", "Odds_Ratio", "OR_CI_Lower", "OR_CI_Upper", "OR_95CI", "p_value", "p_formatted"),
                                 drop = FALSE]
    rownames(ps_model_results) <- NULL
    
    # Reformat class variable names: "female1" -> "female (= 1)", "raceBlack" -> "race (= Black)"
    # Sort class_vars by name length (longest first) to avoid partial prefix matches
    # Only reformat if the suffix is an actual factor level of that variable
    if (!is.null(class_vars) && length(class_vars) > 0) {
      sorted_cvars <- class_vars[order(nchar(class_vars), decreasing = TRUE)]
      for (i in seq_len(nrow(ps_model_results))) {
        coef_name <- ps_model_results$Variable[i]
        for (cv in sorted_cvars) {
          if (startsWith(coef_name, cv) && nchar(coef_name) > nchar(cv)) {
            level_part <- substring(coef_name, nchar(cv) + 1)
            # Verify suffix is an actual factor level (avoids false matches like
            # cont_var "glucose_test_monitor_count" matching class_var "glucose_test_monitor")
            if (is.factor(data_work[[cv]]) && level_part %in% levels(data_work[[cv]])) {
              ps_model_results$Variable[i] <- paste0(cv, " (= ", level_part, ")")
              break
            }
          }
        }
      }
    }
    
    if (verbose) message(sprintf("  OR table includes %d variables (excluding intercept)", nrow(ps_model_results)))
    
    # Save to Excel (2 sheets: PS Model Odds Ratios + PS Model Summary)
    wb <- openxlsx::createWorkbook()
    
    # Sheet 1: PS Model Odds Ratios (clean table only) 
    openxlsx::addWorksheet(wb, "PS Model Odds Ratios")
    
    ps_table_formatted <- ps_model_results[, c("Variable", "OR_95CI", "p_formatted"), drop = FALSE]
    names(ps_table_formatted) <- c("Variables included in PS model", "Odds Ratio (95% CI)", "p-value")
    
    openxlsx::writeDataTable(wb, "PS Model Odds Ratios", ps_table_formatted, startRow = 1)
    openxlsx::setColWidths(wb, "PS Model Odds Ratios", cols = 1:3, widths = c(40, 25, 15))
    
    # Sheet 2: PS Model Summary (2-column metadata table)
    openxlsx::addWorksheet(wb, "PS Model Summary")
    
    n_total <- nrow(stats::model.frame(ps_model))
    n_exp <- sum(data_work[[exposure_var]] == 1, na.rm = TRUE)
    n_ref <- sum(data_work[[exposure_var]] == 0, na.rm = TRUE)
    
    # PS distribution quantiles
    ps_vals <- data_work[[ps_var]]
    ps_q <- stats::quantile(ps_vals, probs = c(0, 0.25, 0.5, 0.75, 1), na.rm = TRUE)
    
    # Helper to collapse variable vectors into comma-separated strings
    .collapse_vars <- function(v) {
      if (is.null(v) || length(v) == 0) "None" else paste(v, collapse = ", ")
    }
    
    summary_labels <- c(
      "Number of observations (Total)",
      "Number of observations (Exposure)",
      "Number of observations (Reference)",
      "PS Model Converged",
      "PS Model C-statistic",
      "PS Model AIC",
      "PS Min",
      "PS Q1",
      "PS Median",
      "PS Mean",
      "PS Q3",
      "PS Max",
      "Number of observations with PS < 0.01",
      "Number of observations with PS > 0.99",
      "Class variables provided as input",
      "Continuous variables provided as input",
      "Class variables included in PS model",
      "Continuous variables included in PS model",
      "Class variables excluded from PS model (due to extreme distribution)",
      "Continuous variables excluded from PS model (due to extreme distribution)"
    )
    
    summary_values <- c(
      as.character(n_total),
      as.character(n_exp),
      as.character(n_ref),
      as.character(ps_model$converged),
      sprintf("%.3f", c_stat),
      sprintf("%.1f", stats::AIC(ps_model)),
      sprintf("%.6f", ps_q[[1]]),
      sprintf("%.6f", ps_q[[2]]),
      sprintf("%.6f", ps_q[[3]]),
      sprintf("%.6f", mean(ps_vals, na.rm = TRUE)),
      sprintf("%.6f", ps_q[[4]]),
      sprintf("%.6f", ps_q[[5]]),
      as.character(extreme_low),
      as.character(extreme_high),
      .collapse_vars(input_class_vars),
      .collapse_vars(input_cont_vars),
      .collapse_vars(class_vars),
      .collapse_vars(cont_vars),
      .collapse_vars(excluded_class_vars),
      .collapse_vars(excluded_cont_vars)
    )
    
    summary_df <- data.frame(Item = summary_labels, Value = summary_values, stringsAsFactors = FALSE)
    openxlsx::writeDataTable(wb, "PS Model Summary", summary_df, startRow = 1)
    
    # Style: wrap text + wider column for Value
    wrap_style <- openxlsx::createStyle(wrapText = TRUE, valign = "top")
    openxlsx::addStyle(wb, "PS Model Summary", style = wrap_style,
                       rows = 1:(nrow(summary_df) + 1), cols = 2, gridExpand = TRUE)
    openxlsx::setColWidths(wb, "PS Model Summary", cols = 1, widths = 60)
    openxlsx::setColWidths(wb, "PS Model Summary", cols = 2, widths = 80)
    
    # Save workbook
    openxlsx::saveWorkbook(wb, out_xlsxpath_odds_ratio, overwrite = TRUE)
    if (verbose) message(sprintf("  OR table saved to Excel: %s", out_xlsxpath_odds_ratio))
  }
  
  # Step 4: Save outputs          
  if (verbose) message("\nStep 4: Saving outputs")
  if (verbose) message("----------------------------------------")
  
  # Save to CSV if path specified
  if (!is.null(out_csvpath)) {
    utils::write.csv(data_work, file = out_csvpath, row.names = FALSE)
    if (verbose) message(sprintf("Data saved to CSV: %s", out_csvpath))
  }
  
  if (verbose) message("\n========================================")
  if (verbose) message("PROPENSITY SCORE CALCULATION COMPLETE")
  if (verbose) message("========================================")
  
  # Return the data with PS added
  return(invisible(data_work))
}