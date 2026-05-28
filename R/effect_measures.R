# HR IR Incidence Rates and Hazard Ratios ############################################
#' Estimate Incidence Rates and Hazard Ratios
#'
#' This function estimates incidence rates and hazard ratios for time-to-event outcomes.
#' It can perform unweighted analysis, weighted (propensity score adjusted) analysis, or both,
#' using potentially different datasets for each analysis type.
#' Incidence rate CIs are estimated using Poisson distribution.
#'
#' @param in_df_unwt Data frame for unweighted analysis. Set to NULL to skip unweighted analysis.
#'                                  This could be the full cohort, matched data, or any specific population
#' @param in_df_wted Data frame for weighted analysis. Set to NULL to skip weighted analysis.
#'                                  This could be the same as unweighted data or a different population (e.g., trimmed)
#' @param out_xlsxpath Character string. Path for output Excel file with all results
#' @param exposure_var Character string. Name of binary exposure/treatment variable (default: "exp")
#' @param exp_value Value representing the exposed/treated group (default: 1)
#' @param ref_value Value representing the reference/control group (default: 0)
#' @param outcome_var Character string. Name of the event indicator variable (0=censored, 1=event)
#' @param weight_var Character string. Name of weight variable in the weighted dataset.
#'                   Required when in_df_wted is not NULL.
#' @param survival_time Character string. Name of survival/follow-up time variable (required)
#' @param time_unit Character string. Unit of time for survival analysis: "days" or "years".
#'                  Default: "days"
#' @param ir_per_pyears Numeric. Multiplier for incidence rate per person-years.
#'                      Allowed values: 1, 100, 1000, 10000, 100000.
#'                      Default: 1000 (i.e., IR per 1,000 person-years).
#'                      Column name is fixed as IR_per_Npy_Unwt / IR_per_Npy_Wted. The actual
#'                      multiplier value is documented in the Analysis Summary sheet.
#' @param confidence_level Numeric. Confidence level for confidence intervals (default: 0.95 for 95% CI)
#' @param readme_text Character string. Optional message to include in README sheet at beginning of Excel file.
#'                    Use this to document the analysis parameters, data sources, or any notes
#' @param stratify_by Character string. Name of the stratification variable.
#'                    When specified:
#'                    - Cox model uses strata() in the formula (separate baseline hazards per stratum)
#'                    - IR/IRD are computed via direct standardization using total cohort person-time
#'                      distribution as weights
#'                    Default: NULL (no stratification)
#' @param verbose Logical. Print progress messages (default TRUE).
#'
#' @return Invisibly returns a list containing:
#'   - incidence_rates: Data frame with event counts, person-years, and incidence rates
#'   - hazard_ratios: Data frame with unweighted and/or weighted HR estimates
#'   - unweighted_model: The unweighted Cox model object (if unweighted analysis performed)
#'   - weighted_model: The weighted Cox model object (if weighted analysis performed)
#'   - stratum_details_unwt: Data frame with stratum-level details (if stratified, unweighted)
#'   - stratum_details_wted: Data frame with stratum-level details (if stratified, weighted)
#'   Also saves results to Excel file when out_xlsxpath is provided
#'
#' @section Side Effects:
#' Creates the output directory and writes an Excel workbook when
#' \code{out_xlsxpath} is provided.
#'
#' @details
#' The function allows flexible analysis configurations:
#'   - Unweighted only: Set in_df_unwt = data, in_df_wted = NULL
#'   - Weighted only: Set in_df_unwt = NULL, in_df_wted = data
#'   - Both: Provide both datasets (can be the same or different populations)
#'
#' For unweighted analysis:
#'   - Uses standard Cox regression without weights
#'   - Can be applied to matched data, trimmed data, or full cohort
#'
#' For weighted analysis:
#'   - Uses weighted Cox regression with specified weights
#'   - Typically applied to PS-weighted data (IPTW, MW, OW, FS weights)
#'
#' Stratification (stratify_by parameter):
#'   When stratify_by is specified:
#'   - Cox model: Fits a stratified Cox proportional hazards model using strata() in the
#'     model formula. This allows separate baseline hazard functions for each stratum while
#'     estimating a common treatment effect (HR) across strata.
#'   - IR/IRD: Uses direct standardization with total cohort person-time distribution as
#'     standard weights.
#'       - IR_std = sum(Rate_k * W_k) where W_k = PT_total_k / PT_total
#'       - Variance via delta method: Var = sum(W_k^2 * d_k / PT_k^2)
#'       - IRD_std = IR_std_exposed - IR_std_unexposed
#'       - CI for standardized IR uses log-transformation to ensure positive bounds
#'
#' ir_per_pyears parameter:
#'   Controls the person-year multiplier for incidence rates and incidence rate differences.
#'   For example:
#'   - ir_per_pyears = 1000: rates expressed per 1,000 person-years (default)
#'   - ir_per_pyears = 100000: rates expressed per 100,000 person-years
#'   Output column names are fixed (e.g., IR_per_Npy_Unwt, IR_per_Npy_Wted). The multiplier
#'   value is recorded in the Analysis Summary sheet.
#'
#' @export
#'
#' @examples
#' csv_path <- system.file("extdata", "sample_data.csv", package = "rwetools")
#' df <- read.csv(csv_path)
#'
#' # Unweighted analysis only
#' result <- estimate_hr_ir(
#'   in_df_unwt    = df,
#'   in_df_wted    = NULL,
#'   exposure_var  = "exposure",
#'   outcome_var   = "outcome",
#'   survival_time = "follow_up_days",
#'   time_unit     = "days",
#'   ir_per_pyears = 1000,
#'   verbose       = FALSE
#' )
#' result$incidence_rates
#'
#' \donttest{
#' # Weighted analysis with Excel output (requires openxlsx)
#' if (requireNamespace("openxlsx", quietly = TRUE)) {
#'   df_ps <- estimate_ps(
#'     in_df        = df,
#'     exposure_var = "exposure",
#'     class_vars   = c("cat1", "cat2", "cat3", "cat4"),
#'     cont_vars    = c("cont1", "cont2", "cont3"),
#'     verbose      = FALSE
#'   )
#'   df_wt <- create_ps_weights(
#'     in_df         = df_ps,
#'     exposure_var  = "exposure",
#'     ps_var        = "ps",
#'     weight_method = "mw",
#'     weight_var    = "mw_wt",
#'     verbose       = FALSE
#'   )
#'   out_xlsx <- tempfile(fileext = ".xlsx")
#'   result2 <- estimate_hr_ir(
#'     in_df_unwt    = df,
#'     in_df_wted    = df_wt,
#'     out_xlsxpath  = out_xlsx,
#'     exposure_var  = "exposure",
#'     outcome_var   = "outcome",
#'     survival_time = "follow_up_days",
#'     weight_var    = "mw_wt",
#'     ir_per_pyears = 1000,
#'     verbose       = FALSE
#'   )
#' }
#' }
estimate_hr_ir <- function(
    in_df_unwt = NULL,
    in_df_wted = NULL,
    out_xlsxpath = NULL,
    exposure_var = "exp",
    exp_value = 1,
    ref_value = 0,
    outcome_var = NULL,
    weight_var = NULL,
    survival_time = NULL,
    time_unit = c("days", "years"),
    ir_per_pyears = 1000,
    confidence_level = 0.95,
    readme_text = NULL,
    stratify_by = NULL,
    verbose = TRUE) {



  # Check required packages
  if (!requireNamespace("survival", quietly = TRUE)) {
    stop("Package 'survival' is required but not installed.")
  }
  if (!is.null(out_xlsxpath)) {
    if (!requireNamespace("openxlsx", quietly = TRUE)) {
      stop("Package 'openxlsx' is required for Excel output but not installed.")
    }
  }

  if (verbose) message("\n========================================")
  if (verbose) message("INCIDENCE RATE & HAZARD RATIO ESTIMATION")
  if (verbose) message("========================================")

  # Check that at least one dataset is provided
  if (is.null(in_df_unwt) && is.null(in_df_wted)) {
    stop("At least one of in_df_unwt or in_df_wted must be provided")
  }

  # Validate IR_per_pyears
  allowed_multipliers <- c(1, 100, 1000, 10000, 100000)
  if (!ir_per_pyears %in% allowed_multipliers) {
    stop(paste0("ir_per_pyears must be one of: ", paste(allowed_multipliers, collapse = ", "),
                ". Got: ", ir_per_pyears))
  }

  # Fixed column name base
  ir_col_base <- "IR_per_Npy"

  # Determine which analyses to perform
  do_unweighted <- !is.null(in_df_unwt)
  do_weighted <- !is.null(in_df_wted)
  do_stratify <- !is.null(stratify_by)

  # Validate required inputs
  if (is.null(outcome_var)) {
    stop("outcome_var must be specified")
  }

  if (is.null(survival_time)) {
    stop("survival_time must be specified")
  }
  if (!is.numeric(confidence_level) || length(confidence_level) != 1 ||
      confidence_level <= 0 || confidence_level >= 1) {
    stop("confidence_level must be a single numeric value between 0 and 1 (e.g., 0.95)")
  }

  # Match time_unit argument
  time_unit <- match.arg(time_unit)

  # Check variables exist in provided datasets
  if (do_unweighted) {
    required_vars <- c(outcome_var, exposure_var, survival_time)
    missing_vars <- required_vars[!required_vars %in% names(in_df_unwt)]
    if (length(missing_vars) > 0) {
      stop(paste("Variables not found in unweighted dataset:", paste(missing_vars, collapse = ", ")))
    }

    if (do_stratify && !stratify_by %in% names(in_df_unwt)) {
      stop(paste("Stratification variable", stratify_by, "not found in unweighted dataset"))
    }
  }

  if (do_weighted) {
    required_vars <- c(outcome_var, exposure_var, survival_time)
    missing_vars <- required_vars[!required_vars %in% names(in_df_wted)]
    if (length(missing_vars) > 0) {
      stop(paste("Variables not found in weighted dataset:", paste(missing_vars, collapse = ", ")))
    }

    if (is.null(weight_var)) {
      stop("weight_var must be specified when in_df_wted is provided")
    }
    if (!weight_var %in% names(in_df_wted)) {
      stop(paste("Weight variable", weight_var, "not found in weighted dataset"))
    }

    if (do_stratify && !stratify_by %in% names(in_df_wted)) {
      stop(paste("Stratification variable", stratify_by, "not found in weighted dataset"))
    }
  }

  # Calculate z-score for confidence intervals
  z_score <- stats::qnorm(1 - (1 - confidence_level) / 2)

  # Initialize results storage
  results <- list()
  incidence_rates <- NULL
  hazard_ratios <- NULL
  stratum_details_unwt <- NULL
  stratum_details_wted <- NULL

  if (verbose) message(sprintf("Outcome: %s", outcome_var))
  if (verbose) message(sprintf("Survival Time: %s (%s)", survival_time, time_unit))
  if (verbose) message(sprintf("IR per person-years multiplier: %s (rates expressed per %s person-years)",
                               format(ir_per_pyears, big.mark = ","), format(ir_per_pyears, big.mark = ",")))
  if (do_stratify) {
    if (verbose) message(sprintf("Stratification Variable: %s", stratify_by))
    if (verbose) message("  Cox model: strata()")
    if (verbose) message("  IR/IRD: Direct Standardization")
  }
  if (verbose) message(sprintf("Confidence Level: %.0f%%", confidence_level * 100))
  if (verbose) message(sprintf("Unweighted Analysis: %s", ifelse(do_unweighted, "Yes", "No")))
  if (verbose) message(sprintf("Weighted Analysis: %s", ifelse(do_weighted, "Yes", "No")))
  if (do_weighted) {
    if (verbose) message(sprintf("Weight Variable: %s", weight_var))
  }
  if (verbose) message("")

  # STEP 1: INCIDENCE RATE CALCULATION         #####################
  if (verbose) message("Step 1: Calculating Incidence Rates")
  if (verbose) message("------------------------------------")

  ##### Process unweighted data #######
  if (do_unweighted) {
    data_unwt <- in_df_unwt

    # Recode exposure if needed
    unique_exposure_vals <- unique(data_unwt[[exposure_var]])
    if (!(all(unique_exposure_vals %in% c(0, 1)))) {
      if (verbose) message("Recoding exposure variable for unweighted data...")
      data_unwt[[exposure_var]] <- ifelse(data_unwt[[exposure_var]] == exp_value, 1,
                                          ifelse(data_unwt[[exposure_var]] == ref_value, 0, NA))
    }

    # Convert time unit if needed
    if (tolower(time_unit) == "days") {
      data_unwt$person_years <- as.numeric(data_unwt[[survival_time]]) / 365.25
      if (verbose) message("Converting survival time from days to years (unweighted)")
    } else {
      data_unwt$person_years <- as.numeric(data_unwt[[survival_time]])
      if (verbose) message("Using survival time in years (unweighted)")
    }

    # Remove missing values
    n_missing <- sum(is.na(data_unwt$person_years))
    if (n_missing > 0) {
      warning(paste(n_missing, "observations with missing survival time removed from unweighted data"))
      data_unwt <- data_unwt[!is.na(data_unwt$person_years), ]
    }

    # Calculate unweighted incidence rates
    if (do_stratify) {
      if (verbose) message("Calculating STANDARDIZED incidence rates (unweighted)...")
      unwt_ir_result <- calc_standardized_ir_ird(
        df = data_unwt, exp_var = exposure_var, out_var = outcome_var, py_var = "person_years",
        strat_var = stratify_by, conf_level = confidence_level, z = z_score, suffix = "Unwt",
        multiplier = ir_per_pyears, ir_base = ir_col_base
      )
      incidence_rates <- unwt_ir_result$summary
      stratum_details_unwt <- unwt_ir_result$stratum_detail
    } else {
      if (verbose) message("Calculating crude incidence rates (unweighted)...")
      incidence_rates <- calc_crude_ir_ird(
        df = data_unwt, exp_var = exposure_var, out_var = outcome_var, py_var = "person_years",
        conf_level = confidence_level, z = z_score, suffix = "Unwt",
        multiplier = ir_per_pyears, ir_base = ir_col_base
      )
    }
  }

  ##### Process weighted data #######
  if (do_weighted) {
    data_wted <- in_df_wted

    # Recode exposure if needed
    unique_exposure_vals <- unique(data_wted[[exposure_var]])
    if (!(all(unique_exposure_vals %in% c(0, 1)))) {
      if (verbose) message("Recoding exposure variable for weighted data...")
      data_wted[[exposure_var]] <- ifelse(data_wted[[exposure_var]] == exp_value, 1,
                                          ifelse(data_wted[[exposure_var]] == ref_value, 0, NA))
    }

    # Convert time unit if needed
    if (tolower(time_unit) == "days") {
      data_wted$person_years <- as.numeric(data_wted[[survival_time]]) / 365.25
      if (verbose) message("Converting survival time from days to years (weighted)")
    } else {
      data_wted$person_years <- as.numeric(data_wted[[survival_time]])
      if (verbose) message("Using survival time in years (weighted)")
    }

    # Remove missing values
    n_missing <- sum(is.na(data_wted$person_years))
    if (n_missing > 0) {
      warning(paste(n_missing, "observations with missing survival time removed from weighted data"))
      data_wted <- data_wted[!is.na(data_wted$person_years), ]
    }

    # Calculate weighted incidence rates
    if (do_stratify) {
      if (verbose) message("Calculating STANDARDIZED incidence rates (weighted)...")
      wted_ir_result <- calc_standardized_ir_ird_weighted(
        df = data_wted, exp_var = exposure_var, out_var = outcome_var, py_var = "person_years",
        wt_var = weight_var, strat_var = stratify_by,
        conf_level = confidence_level, z = z_score, suffix = "Wted",
        multiplier = ir_per_pyears, ir_base = ir_col_base
      )
      wted_summary <- wted_ir_result$summary
      stratum_details_wted <- wted_ir_result$stratum_detail
    } else {
      if (verbose) message("Calculating crude incidence rates (weighted)...")
      wted_summary <- calc_crude_ir_ird_weighted(
        df = data_wted, exp_var = exposure_var, out_var = outcome_var, py_var = "person_years",
        wt_var = weight_var, conf_level = confidence_level, z = z_score, suffix = "Wted",
        multiplier = ir_per_pyears, ir_base = ir_col_base
      )
    }

    # Combine unweighted and weighted results
    if (do_unweighted) {
      wted_cols <- setdiff(names(wted_summary), exposure_var)
      incidence_rates <- cbind(incidence_rates, wted_summary[, wted_cols, drop = FALSE])
    } else {
      incidence_rates <- wted_summary
    }
  }

  # Add Exposure_Group labels
  incidence_rates$Exposure_Group <- ifelse(
    incidence_rates[[exposure_var]] == 99, "Total",
    ifelse(incidence_rates[[exposure_var]] == 1, "Exposed",
           ifelse(incidence_rates[[exposure_var]] == 0, "Reference",
                  as.character(incidence_rates[[exposure_var]])))
  )

  # Reorder columns: Exposure_Group first, remove raw exposure
  col_order <- c("Exposure_Group", setdiff(names(incidence_rates), c("Exposure_Group", exposure_var)))
  incidence_rates <- incidence_rates[, col_order]

  # Print incidence summary
  if (verbose) message(sprintf("\nIncidence Rate Summary (per %s person-years):", format(ir_per_pyears, big.mark = ",")))
  if (do_stratify) {
    if (verbose) message("  IR/IRD Method: Direct Standardization")
  }
  if (verbose) print(incidence_rates, digits = 3)
  if (verbose) message("")

  # STEP 2: HAZARD RATIO ESTIMATION         #####################
  if (verbose) message("Step 2: Fitting Cox Proportional Hazards Models")
  if (do_stratify) {
    if (verbose) message(sprintf("        (Stratified by: %s)", stratify_by))
  }
  if (verbose) message("------------------------------------------------")

  if (do_unweighted) {
    if (verbose) message("Unweighted Cox model...")

    # Build formula with or without stratification
    if (do_stratify) {
      unwt_formula <- stats::as.formula(paste("survival::Surv(", survival_time, ",", outcome_var, ") ~",
                                              exposure_var, "+ survival::strata(", stratify_by, ")"))
      if (verbose) message(sprintf("  Formula: Surv(%s, %s) ~ %s + strata(%s)",
                                   survival_time, outcome_var, exposure_var, stratify_by))
    } else {
      unwt_formula <- stats::as.formula(paste("survival::Surv(", survival_time, ",", outcome_var, ") ~", exposure_var))
      if (verbose) message(sprintf("  Formula: Surv(%s, %s) ~ %s", survival_time, outcome_var, exposure_var))
    }

    unwt_model <- survival::coxph(unwt_formula, data = data_unwt)

    # Extract HR from Cox model
    unwt_s <- summary(unwt_model, conf.int = confidence_level)
    unwt_coef_df <- data.frame(
      term      = rownames(unwt_s$coefficients),
      estimate  = unwt_s$conf.int[, "exp(coef)"],
      conf.low  = unwt_s$conf.int[, 3],
      conf.high = unwt_s$conf.int[, 4],
      std.error = unwt_s$coefficients[, "se(coef)"],
      p.value   = unwt_s$coefficients[, "Pr(>|z|)"],
      stringsAsFactors = FALSE,
      row.names = NULL
    )
    unwt_row <- unwt_coef_df[unwt_coef_df$term == exposure_var, , drop = FALSE]
    unwt_results <- data.frame(
      Analysis = ifelse(do_stratify,
                        paste0("Unweighted (stratified by ", stratify_by, ")"),
                        "Unweighted"),
      HR_CI = sprintf("%.2f (%.2f, %.2f)", unwt_row$estimate, unwt_row$conf.low, unwt_row$conf.high),
      HR = unwt_row$estimate,
      HR_LCI = unwt_row$conf.low,
      HR_UCI = unwt_row$conf.high,
      lnHR = log(unwt_row$estimate),
      lnHR_SE = unwt_row$std.error,
      Pvalue = unwt_row$p.value,
      stringsAsFactors = FALSE,
      row.names = NULL
    )

    hazard_ratios <- unwt_results
    results$unweighted_model <- unwt_model
  }

  if (do_weighted) {
    if (verbose) message("Weighted Cox model...")

    # Build formula with or without stratification
    if (do_stratify) {
      wted_formula <- stats::as.formula(paste("survival::Surv(", survival_time, ",", outcome_var, ") ~",
                                              exposure_var, "+ survival::strata(", stratify_by, ")"))
      if (verbose) message(sprintf("  Formula: Surv(%s, %s) ~ %s + strata(%s)",
                                   survival_time, outcome_var, exposure_var, stratify_by))
    } else {
      wted_formula <- stats::as.formula(paste("survival::Surv(", survival_time, ",", outcome_var, ") ~", exposure_var))
      if (verbose) message(sprintf("  Formula: Surv(%s, %s) ~ %s", survival_time, outcome_var, exposure_var))
    }

    weight_values <- data_wted[[weight_var]]
    weighted_model <- survival::coxph(wted_formula, data = data_wted, weights = weight_values)

    # Extract HR from Cox model
    wted_s <- summary(weighted_model, conf.int = confidence_level)
    wted_coef_df <- data.frame(
      term      = rownames(wted_s$coefficients),
      estimate  = wted_s$conf.int[, "exp(coef)"],
      conf.low  = wted_s$conf.int[, 3],
      conf.high = wted_s$conf.int[, 4],
      std.error = wted_s$coefficients[, "se(coef)"],
      p.value   = wted_s$coefficients[, "Pr(>|z|)"],
      stringsAsFactors = FALSE,
      row.names = NULL
    )
    wted_row <- wted_coef_df[wted_coef_df$term == exposure_var, , drop = FALSE]
    weighted_results <- data.frame(
      Analysis = ifelse(do_stratify,
                        paste0("Weighted (stratified by ", stratify_by, ")"),
                        "Weighted"),
      HR_CI = sprintf("%.2f (%.2f, %.2f)", wted_row$estimate, wted_row$conf.low, wted_row$conf.high),
      HR = wted_row$estimate,
      HR_LCI = wted_row$conf.low,
      HR_UCI = wted_row$conf.high,
      lnHR = log(wted_row$estimate),
      lnHR_SE = wted_row$std.error,
      Pvalue = wted_row$p.value,
      stringsAsFactors = FALSE,
      row.names = NULL
    )

    if (do_unweighted) {
      hazard_ratios <- rbind(unwt_results, weighted_results)
    } else {
      hazard_ratios <- weighted_results
    }
    results$weighted_model <- weighted_model
  }

  if (verbose) message("\nHazard Ratio Estimates:")
  if (verbose) print(hazard_ratios, digits = 3)
  if (verbose) message("")

  # STEP 3: SAVE RESULTS         ######################
  if (verbose) message("Step 3: Saving Results")
  if (verbose) message("----------------------")

  # Save to Excel if specified
  if (!is.null(out_xlsxpath)) {
    wb <- openxlsx::createWorkbook()

    # Add README if provided
    if (!is.null(readme_text)) {
      add_readme_sheet(wb, readme_text, verbose = verbose)
    }

    # Analysis Summary Sheet
    openxlsx::addWorksheet(wb, "Analysis Summary")

    # Determine sample sizes
    unwt_sample_size <- ifelse(do_unweighted, nrow(data_unwt), "N/A")
    wted_sample_size <- ifelse(do_weighted, nrow(data_wted), "N/A")

    summary_info <- data.frame(
      Parameter = c("Effect Measure", "Outcome", "Exposure",
                    "Survival Time Variable", "Time Unit",
                    "IR per Person-Years Multiplier",
                    "Confidence Level",
                    "Unweighted Analysis", "Unweighted Sample Size",
                    "Weighted Analysis", "Weighted Sample Size", "Weight Variable"),
      Value = c("Hazard Ratio", outcome_var, exposure_var,
                survival_time, time_unit,
                format(ir_per_pyears, big.mark = ","),
                paste0(confidence_level * 100, "%"),
                ifelse(do_unweighted, "Yes", "No"),
                as.character(unwt_sample_size),
                ifelse(do_weighted, "Yes", "No"),
                as.character(wted_sample_size),
                ifelse(do_weighted, weight_var, "N/A"))
    )

    # Add stratification info if applicable
    if (do_stratify) {
      stratify_levels_unwt <- "N/A"
      stratify_levels_wted <- "N/A"

      if (do_unweighted) {
        unique_vals_unwt <- sort(unique(data_unwt[[stratify_by]]))
        stratify_levels_unwt <- paste(unique_vals_unwt, collapse = ", ")
      }
      if (do_weighted) {
        unique_vals_wted <- sort(unique(data_wted[[stratify_by]]))
        stratify_levels_wted <- paste(unique_vals_wted, collapse = ", ")
      }

      summary_info <- rbind(
        summary_info,
        data.frame(Parameter = "Stratification Variable", Value = stratify_by),
        data.frame(Parameter = "Stratification Levels (Unweighted)", Value = stratify_levels_unwt),
        data.frame(Parameter = "Stratification Levels (Weighted)", Value = stratify_levels_wted),
        data.frame(Parameter = "Cox Stratification", Value = "strata()"),
        data.frame(Parameter = "IR/IRD Stratification Method", Value = "Direct Standardization")
      )
    }

    openxlsx::writeDataTable(wb, "Analysis Summary", summary_info, startRow = 1)

    # Incidence Rates Sheet
    openxlsx::addWorksheet(wb, "Incidence Rates")
    openxlsx::writeDataTable(wb, "Incidence Rates", incidence_rates, startRow = 1)
    if (verbose) message("  Added Incidence Rates sheet")

    # Hazard Ratios Sheet
    openxlsx::addWorksheet(wb, "Hazard Ratios")
    openxlsx::writeDataTable(wb, "Hazard Ratios", hazard_ratios, startRow = 1)
    if (verbose) message("  Added Hazard Ratios sheet")

    # Stratum Details Sheets (if stratified)
    if (!is.null(stratum_details_unwt)) {
      sheet_name <- "Stratum Detail (Unwt Std)"
      openxlsx::addWorksheet(wb, sheet_name)
      openxlsx::writeDataTable(wb, sheet_name, stratum_details_unwt, startRow = 1)
      if (verbose) message(sprintf("  Added %s sheet", sheet_name))
    }

    if (!is.null(stratum_details_wted)) {
      sheet_name <- "Stratum Detail (Wted Std)"
      openxlsx::addWorksheet(wb, sheet_name)
      openxlsx::writeDataTable(wb, sheet_name, stratum_details_wted, startRow = 1)
      if (verbose) message(sprintf("  Added %s sheet", sheet_name))
    }

    # Check if output directory exists, create if it doesn't
    out_dir <- dirname(out_xlsxpath)
    if (!dir.exists(out_dir)) {
      dir.create(out_dir, recursive = TRUE)
      if (verbose) message(sprintf("  Created output directory: %s", out_dir))
    }

    # Save workbook
    openxlsx::saveWorkbook(wb, out_xlsxpath, overwrite = TRUE)
    if (verbose) message(sprintf("\u2713 Results saved to: %s", out_xlsxpath))
  }

  if (verbose) message("\n========================================")
  if (verbose) message("ANALYSIS COMPLETE")
  if (verbose) message("========================================")

  # Store results
  results$incidence_rates <- incidence_rates
  results$hazard_ratios <- hazard_ratios
  results$ir_per_pyears <- ir_per_pyears
  if (do_stratify) {
    results$stratify_by <- stratify_by
  }
  if (!is.null(stratum_details_unwt)) results$stratum_details_unwt <- stratum_details_unwt
  if (!is.null(stratum_details_wted)) results$stratum_details_wted <- stratum_details_wted

  return(invisible(results))
}


# RR Risk Ratios and Risk Differences ##############
#' Estimate Risk Ratios and Risk Differences using Cumulative Incidence
#'
#' This function estimates risk ratios (RR) and risk differences (RD) for time-to-event outcomes
#' using cumulative incidence at a specified timepoint.
#' Two estimators are supported:
#' \describe{
#'   \item{Kaplan-Meier (KM)}{Standard 1 - S(t) cumulative incidence. Appropriate when there
#'     are no competing risks or competing events are treated as censored.}
#'   \item{Aalen-Johansen (AJ)}{Cumulative incidence function accounting for competing risks.
#'     Produces unbiased risk estimates in the presence of competing events by treating them
#'     as a separate absorbing state rather than censoring. Requires a competing event indicator
#'     via \code{competing_event_var}. Uses multi-state \code{\link[survival]{survfit}} internally
#'     (\code{survival} >= 3.1).}
#' }
#' It can perform unweighted analysis, weighted (propensity score adjusted) analysis, or both,
#' using potentially different datasets for each analysis type.
#' Confidence intervals can be estimated using bootstrapping, analytical (delta method / Greenwood),
#' or both.
#'
#' When \code{stratify_by} is specified, the function performs direct standardization:
#' stratum-specific risks are combined using the total population distribution as the standard.
#' Bootstrap resampling is stratified to preserve stratum structure.
#'
#' @param in_df_unwt Data frame for unweighted analysis. Set to NULL to skip unweighted analysis.
#'                   This could be the full cohort, matched data, or any specific population.
#' @param in_df_wted Data frame for weighted analysis. Set to NULL to skip weighted analysis.
#'                   This could be the same as unweighted data or a different population (e.g., trimmed).
#' @param out_xlsxpath Character string. Path for output Excel file with all results.
#' @param exposure_var Character string. Name of binary exposure/treatment variable (default: "exp").
#' @param exp_value Value representing the exposed/treated group (default: 1).
#' @param ref_value Value representing the reference/control group (default: 0).
#' @param outcome_var Character string. Name of the event indicator variable (0=censored, 1=event).
#' @param weight_var Character string. Name of weight variable in the weighted dataset.
#'                   Required when \code{in_df_wted} is not NULL.
#' @param survival_time Character string. Name of survival/follow-up time variable (required).
#' @param time_unit Character string. Unit of time for survival analysis: "days", "months", or "years".
#'                  Default: "days".
#' @param rr_rd_at_timepoint Numeric. Timepoint at which to calculate cumulative incidence.
#'                           Units depend on \code{time_unit} parameter. Default: 365.
#' @param rr_rd_per_individuals Numeric. Denominator for expressing risk and risk difference.
#'                              Default: 1000.
#' @param confidence_level Numeric. Confidence level for confidence intervals (default: 0.95).
#' @param conf_int_method Character string. Method for confidence interval estimation:
#'   \describe{
#'     \item{"bootstrap"}{Percentile bootstrap CIs only.}
#'     \item{"analytical"}{Analytical CIs only (delta method for RR, Greenwood-based for RD).
#'       \code{bootstrap_count} is ignored.}
#'     \item{"both"}{Compute both bootstrap and analytical CIs.}
#'   }
#'   Default: "bootstrap".
#' @param risk_estimator Character string. Cumulative incidence estimator:
#'   \describe{
#'     \item{"KM"}{Kaplan-Meier (1 - survival). Aliases: "Kaplan-Meier", "kaplan-meier".}
#'     \item{"AJ"}{Aalen-Johansen cumulative incidence function for competing risks.
#'       Requires \code{competing_event_var}. Aliases: "Aalen-Johansen", "aalen-johansen".}
#'   }
#'   Default: "KM".
#' @param competing_event_var Character string. Name of the competing event indicator variable
#'   (1 = competing event occurred, e.g., death; 0 = no competing event).
#'   Required when \code{risk_estimator = "AJ"}. Ignored when \code{risk_estimator = "KM"}.
#'   The multi-state status is constructed internally:
#'   \code{outcome_var == 1} -> event of interest (1),
#'   \code{outcome_var == 0 & competing_event_var == 1} -> competing event (2),
#'   \code{outcome_var == 0 & competing_event_var == 0} -> censored (0).
#' @param bootstrap_count Integer. Number of bootstrap iterations for CI estimation (default: 500).
#'                        Ignored when \code{conf_int_method = "analytical"}.
#' @param stratify_by Character string. Name of stratification variable for direct standardization
#'                    (default: NULL). When specified, the function:
#'                    (1) calculates stratum weights (w_k = N_k / N_total) from the total population,
#'                    (2) computes risk within each stratum using the selected estimator,
#'                    (3) standardizes: Risk_std = Sum(Risk_k * w_k),
#'                    (4) derives RR and RD from standardized risks.
#'                    Bootstrap CIs use stratified resampling (within each stratum).
#'                    Analytical SEs use Greenwood formula: Var(Risk_std) = Sum(w_k^2 * Var(R_k)).
#' @param n_cores Integer. Number of CPU cores for parallel bootstrap. Default uses all
#'                available cores minus 1. Set to 1 to disable parallelization.
#'                Ignored when \code{conf_int_method = "analytical"}.
#' @param seed Integer or NULL. Optional seed for the bootstrap random number
#'             generator. When NULL (default), the caller's current RNG state
#'             is used and not modified. When an integer is supplied, it is
#'             passed to \code{\link[parallel]{clusterSetRNGStream}} (parallel)
#'             or \code{\link{set.seed}} (sequential) to make the bootstrap
#'             CIs reproducible. In the sequential case the caller's
#'             \code{.Random.seed} is saved and restored on exit.
#'             Ignored when \code{conf_int_method = "analytical"}.
#' @param readme_text Character string. Optional message to include in README sheet of Excel output.
#' @param verbose Logical. Print progress messages (default TRUE).
#'
#' @return Invisibly returns a list containing:
#'   \describe{
#'     \item{estimates}{Data frame with RR, RD, confidence intervals, and
#'       \code{Risk_Estimator} (\code{"KM"} or \code{"AJ"}).}
#'     \item{cumulative_incidence}{Data frame with cumulative incidence by group.
#'       Columns: \code{Risk} (probability), \code{Risk_LCI} / \code{Risk_UCI}
#'       (log-transformed CI), \code{Risk_SE} (SE on 0-1 scale),
#'       \code{Risk_Var} (variance), \code{RiskperN} (scaled risk),
#'       \code{RiskperN_LCI} / \code{RiskperN_UCI} (scaled CI),
#'       and \code{Risk_Estimator}.}
#'     \item{stratum_details}{(if \code{stratify_by} specified) Data frame with
#'       stratum-level results and \code{Risk_Estimator}.}
#'     \item{unweighted_bootstrap}{(if applicable) Bootstrap results matrix.}
#'     \item{weighted_bootstrap}{(if applicable) Bootstrap results matrix.}
#'   }
#'   Also saves results to Excel file when \code{out_xlsxpath} is provided.
#'
#' @section Side Effects:
#' Creates the output directory and writes an Excel workbook when
#' \code{out_xlsxpath} is provided.
#'
#' @details
#' The function proceeds as follows:
#'   \enumerate{
#'     \item Fits cumulative incidence curves by exposure group using the selected estimator
#'     \item Extracts cumulative incidence (risk) at the specified timepoint
#'     \item Computes Risk Ratio (RR) = Risk_exposed / Risk_reference
#'     \item Computes Risk Difference (RD) = Risk_exposed - Risk_reference
#'     \item Estimates confidence intervals using the selected method
#'   }
#'
#' \strong{Kaplan-Meier (KM) estimator:}
#' Risk = 1 - S(t), where S(t) is the KM survival estimate. Competing events, if present,
#' are treated as censored, which can overestimate cumulative incidence.
#'
#' \strong{Aalen-Johansen (AJ) estimator:}
#' Cumulative incidence function from a multi-state model via
#' \code{\link[survival]{survfit}} with \code{Surv(time, factor(status))}.
#' Competing events are modelled as a separate absorbing state, producing
#' unbiased risk estimates. Variance is based on the counting-process
#' variance estimator (Aalen, 1978).
#'
#' \strong{Analytical CI method:}
#' \itemize{
#'   \item RR: delta method on log scale. Var(log(RR)) = Var(R_exp)/R_exp^2 + Var(R_ref)/R_ref^2.
#'     CI = exp(log(RR) +/- z * SE(log(RR))).
#'   \item RD: Var(RD) = (Var(R_exp) + Var(R_ref)) * per_n^2. CI = RD +/- z * SE(RD).
#'   \item Variance from \code{\link[survival]{survfit}} (Greenwood for KM, counting-process for AJ).
#' }
#'
#' \strong{Bootstrap CI method:}
#' Percentile bootstrap with parallel computation via
#' \code{\link[parallel]{parLapply}} using L'Ecuyer-CMRG random number streams.
#' Pass an integer to \code{seed} for reproducible bootstrap CIs; by default
#' (\code{seed = NULL}) the caller's current RNG state is used and not modified.
#'
#' When \code{stratify_by} is specified (Direct Standardization):
#' \enumerate{
#'   \item Stratum weights w_k = N_k / N_total from the total population
#'   \item Within each stratum k, risk R_\{ik\} is computed for each exposure group i
#'   \item Standardized risk: Risk_std_i = Sum_k(R_\{ik\} * w_k)
#'   \item RR = Risk_std_1 / Risk_std_0; RD = Risk_std_1 - Risk_std_0
#'   \item Analytical SE: Var(Risk_std_i) = Sum_k(w_k^2 * Var(R_\{ik\}))
#'   \item Bootstrap CIs use stratified resampling with fixed original population weights
#' }
#'
#' @export
#'
#' @examples
#' csv_path <- system.file("extdata", "sample_data.csv", package = "rwetools")
#' df <- read.csv(csv_path)
#'
#' # KM with analytical CIs (fast, no bootstrap)
#' result <- estimate_rr_rd(
#'   in_df_unwt            = df,
#'   in_df_wted            = NULL,
#'   exposure_var          = "exposure",
#'   outcome_var           = "outcome",
#'   survival_time         = "follow_up_days",
#'   time_unit             = "days",
#'   rr_rd_at_timepoint    = 365,
#'   rr_rd_per_individuals = 1000,
#'   conf_int_method       = "analytical",
#'   verbose               = FALSE
#' )
#' result$estimates
#'
#' \donttest{
#' # KM with bootstrap CIs (slow) and competing risk (AJ estimator)
#' result2 <- estimate_rr_rd(
#'   in_df_unwt            = df,
#'   in_df_wted            = NULL,
#'   exposure_var          = "exposure",
#'   outcome_var           = "outcome",
#'   survival_time         = "follow_up_days",
#'   time_unit             = "days",
#'   rr_rd_at_timepoint    = 365,
#'   rr_rd_per_individuals = 1000,
#'   risk_estimator        = "AJ",
#'   competing_event_var   = "competing_event",
#'   conf_int_method       = "bootstrap",
#'   bootstrap_count       = 100,
#'   n_cores               = 1,
#'   verbose               = FALSE
#' )
#' }
estimate_rr_rd <- function(
    in_df_unwt = NULL,
    in_df_wted = NULL,
    out_xlsxpath = NULL,
    exposure_var = "exp",
    exp_value = 1,
    ref_value = 0,
    outcome_var = NULL,
    weight_var = NULL,
    survival_time = NULL,
    time_unit = c("days", "months", "years"),
    rr_rd_at_timepoint = 365,
    rr_rd_per_individuals = 1000,
    confidence_level = 0.95,
    conf_int_method = c("bootstrap", "analytical", "both"),
    risk_estimator = c("KM", "AJ"),
    competing_event_var = NULL,
    bootstrap_count = 500,
    stratify_by = NULL,
    n_cores = NULL,
    seed = NULL,
    readme_text = NULL,
    verbose = TRUE) {

  # Check required packages
  if (!requireNamespace("survival", quietly = TRUE)) {
    stop("Package 'survival' is required but not installed.")
  }
  if (!is.null(out_xlsxpath)) {
    if (!requireNamespace("openxlsx", quietly = TRUE)) {
      stop("Package 'openxlsx' is required for Excel output but not installed.")
    }
  }

  if (verbose) message("\n========================================")
  if (verbose) message("RISK RATIO & RISK DIFFERENCE ESTIMATION")
  if (!is.null(stratify_by)) {
    if (verbose) message("*** Direct Standardization Enabled ***")
  }
  if (verbose) message("========================================")

  # --- Label for per-individual denominator ---
  per_label <- format(rr_rd_per_individuals, big.mark = ",", scientific = FALSE)

  #### 1. Input Validation & Setup ########################
  # Match arguments
  time_unit <- match.arg(time_unit)
  conf_int_method <- match.arg(conf_int_method)

  # Normalize risk_estimator aliases
  risk_estimator_aliases <- c(
    "KM" = "KM", "Kaplan-Meier" = "KM", "kaplan-meier" = "KM",
    "AJ" = "AJ", "Aalen-Johansen" = "AJ", "aalen-johansen" = "AJ"
  )
  if (!risk_estimator[1] %in% names(risk_estimator_aliases)) {
    stop(paste0(
      "risk_estimator must be one of: ",
      paste(names(risk_estimator_aliases), collapse = ", "),
      ". Got: '", risk_estimator[1], "'"
    ))
  }
  risk_estimator <- risk_estimator_aliases[risk_estimator[1]]
  names(risk_estimator) <- NULL

  # Validate competing_event_var for AJ

  if (risk_estimator == "AJ") {
    if (is.null(competing_event_var)) {
      stop("competing_event_var must be specified when risk_estimator = 'AJ'.")
    }
  }

  do_boot       <- conf_int_method %in% c("bootstrap", "both")
  do_analytical  <- conf_int_method %in% c("analytical", "both")

  # Check that at least one dataset is provided
  if (is.null(in_df_unwt) && is.null(in_df_wted)) {
    stop("At least one of in_df_unwt or in_df_wted must be provided.")
  }

  # Determine which analyses to perform
  do_unweighted <- !is.null(in_df_unwt)
  do_weighted   <- !is.null(in_df_wted)

  # Validate required inputs
  if (is.null(outcome_var)) {
    stop("outcome_var must be specified.")
  }
  if (is.null(survival_time)) {
    stop("survival_time must be specified.")
  }
  if (!is.numeric(confidence_level) || length(confidence_level) != 1L ||
      confidence_level <= 0 || confidence_level >= 1) {
    stop("confidence_level must be a single numeric value between 0 and 1 (e.g., 0.95).")
  }
  if (!is.numeric(rr_rd_at_timepoint) || length(rr_rd_at_timepoint) != 1L ||
      rr_rd_at_timepoint <= 0) {
    stop("rr_rd_at_timepoint must be a single positive number.")
  }
  if (do_boot) {
    if (!is.numeric(bootstrap_count) || length(bootstrap_count) != 1L || bootstrap_count < 1) {
      stop("bootstrap_count must be a positive integer (>= 1).")
    }
    bootstrap_count <- as.integer(bootstrap_count)
  }
  if (!is.numeric(rr_rd_per_individuals) || length(rr_rd_per_individuals) != 1L ||
      rr_rd_per_individuals <= 0) {
    stop("rr_rd_per_individuals must be a single positive number (e.g., 100, 1000, 10000, 100000).")
  }

  # Resolve n_cores (only relevant for bootstrap)
  if (do_boot) {
    if (is.null(n_cores)) {
      n_cores <- max(1L, parallel::detectCores() - 1L)
    }
    n_cores <- as.integer(max(1L, n_cores))
  }

  # Check variables exist in provided datasets
  if (do_unweighted) {
    required_vars <- c(outcome_var, exposure_var, survival_time)
    if (!is.null(stratify_by)) required_vars <- c(required_vars, stratify_by)
    if (risk_estimator == "AJ") required_vars <- c(required_vars, competing_event_var)
    missing_vars <- required_vars[!required_vars %in% names(in_df_unwt)]
    if (length(missing_vars) > 0L) {
      stop(paste("Variables not found in unweighted dataset:",
                 paste(missing_vars, collapse = ", ")))
    }
    # Validate competing_event_var is binary in unweighted data
    if (risk_estimator == "AJ") {
      cev_vals <- unique(in_df_unwt[[competing_event_var]])
      cev_vals <- cev_vals[!is.na(cev_vals)]
      if (!all(cev_vals %in% c(0, 1))) {
        stop(paste0("competing_event_var ('", competing_event_var,
                    "') must be binary (0/1) in unweighted dataset. Found: ",
                    paste(sort(unique(cev_vals)), collapse = ", ")))
      }
    }
  }

  if (do_weighted) {
    required_vars <- c(outcome_var, exposure_var, survival_time)
    if (!is.null(stratify_by)) required_vars <- c(required_vars, stratify_by)
    if (risk_estimator == "AJ") required_vars <- c(required_vars, competing_event_var)
    missing_vars <- required_vars[!required_vars %in% names(in_df_wted)]
    if (length(missing_vars) > 0L) {
      stop(paste("Variables not found in weighted dataset:",
                 paste(missing_vars, collapse = ", ")))
    }
    if (is.null(weight_var)) {
      stop("weight_var must be specified when in_df_wted is provided.")
    }
    if (!weight_var %in% names(in_df_wted)) {
      stop(paste("Weight variable", weight_var, "not found in weighted dataset."))
    }
    # Validate competing_event_var is binary in weighted data
    if (risk_estimator == "AJ") {
      cev_vals <- unique(in_df_wted[[competing_event_var]])
      cev_vals <- cev_vals[!is.na(cev_vals)]
      if (!all(cev_vals %in% c(0, 1))) {
        stop(paste0("competing_event_var ('", competing_event_var,
                    "') must be binary (0/1) in weighted dataset. Found: ",
                    paste(sort(unique(cev_vals)), collapse = ", ")))
      }
    }
  }

  # Validate stratify_by if specified
  if (!is.null(stratify_by)) {
    if (do_unweighted) {
      n_strata_unwt <- length(unique(in_df_unwt[[stratify_by]]))
      if (n_strata_unwt < 2L) {
        stop(paste("stratify_by variable", stratify_by,
                   "must have at least 2 levels in unweighted data."))
      }
    }
    if (do_weighted) {
      n_strata_wted <- length(unique(in_df_wted[[stratify_by]]))
      if (n_strata_wted < 2L) {
        stop(paste("stratify_by variable", stratify_by,
                   "must have at least 2 levels in weighted data."))
      }
    }
  }

  # Print analysis parameters
  estimator_label <- ifelse(risk_estimator == "AJ", "Aalen-Johansen", "Kaplan-Meier")
  if (verbose) {
    message(sprintf("Outcome: %s", outcome_var))
    message(sprintf("Exposure: %s (Exposed=%s, Reference=%s)",
                    exposure_var, exp_value, ref_value))
    message(sprintf("Survival Time: %s (%s)", survival_time, time_unit))
    message(sprintf("Timepoint for RR/RD: %s %s", rr_rd_at_timepoint, time_unit))
    message(sprintf("Risk Estimator: %s (%s)", risk_estimator, estimator_label))
    if (risk_estimator == "AJ") {
      message(sprintf("Competing Event Variable: %s", competing_event_var))
    }
    message(sprintf("Risk expressed per: %s individuals", per_label))
    message(sprintf("Confidence Level: %.0f%%", confidence_level * 100))
    message(sprintf("CI Method: %s", conf_int_method))
    if (do_boot) {
      message(sprintf("Bootstrap Iterations: %d", bootstrap_count))
      message(sprintf("Parallel cores: %d", n_cores))
    }
    message(sprintf("Unweighted Analysis: %s", ifelse(do_unweighted, "Yes", "No")))
    message(sprintf("Weighted Analysis: %s", ifelse(do_weighted, "Yes", "No")))
    if (do_weighted) message(sprintf("Weight Variable: %s", weight_var))
    if (!is.null(stratify_by)) {
      message(sprintf("Stratification Variable: %s (Direct Standardization)", stratify_by))
    }
    message("")
  }

  # z-value for analytical CIs
  z_val <- stats::qnorm(1 - (1 - confidence_level) / 2)

  # Initialize results storage
  results <- list()
  all_estimates <- data.frame()
  all_cuminc_data <- data.frame()
  all_stratum_details <- data.frame()

  #### 3. Run Analysis (common logic for both unweighted and weighted) ########
  # Run a single analysis branch (unweighted or weighted)
  # Handles both unstratified and stratified (standardized) analysis
  # Returns: List with estimates row, cumulative incidence data, stratum details, bootstrap matrix
  run_analysis <- function(data, analysis_label, wt_var = NULL) {

    if (verbose) {
      message(sprintf("  [%s] Analysis", analysis_label))
      message(sprintf("  %s", paste(rep("-", nchar(analysis_label) + 12), collapse = "")))
    }

    # --- Recode exposure if needed ---
    unique_exp_vals <- unique(data[[exposure_var]])
    if (!(all(unique_exp_vals %in% c(0, 1)))) {
      if (verbose) message("  Recoding exposure variable...")
      data[[exposure_var]] <- ifelse(data[[exposure_var]] == exp_value, 1L,
                                     ifelse(data[[exposure_var]] == ref_value, 0L, NA_integer_))
      exp_val_int <- 1L
      ref_val_int <- 0L
    } else {
      exp_val_int <- exp_value
      ref_val_int <- ref_value
    }

    # --- Remove missing values ---
    cols_to_check <- c(survival_time, outcome_var, exposure_var)
    if (!is.null(wt_var)) cols_to_check <- c(cols_to_check, wt_var)
    if (!is.null(stratify_by)) cols_to_check <- c(cols_to_check, stratify_by)
    if (risk_estimator == "AJ") cols_to_check <- c(cols_to_check, competing_event_var)

    n_before <- nrow(data)
    complete_mask <- stats::complete.cases(data[, cols_to_check, drop = FALSE])
    data <- data[complete_mask, , drop = FALSE]
    n_after <- nrow(data)
    if (n_before > n_after && verbose) {
      message(sprintf("  Removed %d observations with missing values", n_before - n_after))
    }

    # --- Slim data: keep only needed columns for bootstrap efficiency ---
    keep_cols <- unique(c(survival_time, outcome_var, exposure_var, wt_var, stratify_by,
                          if (risk_estimator == "AJ") competing_event_var))
    data <- data[, keep_cols, drop = FALSE]

    if (verbose) {
      message(sprintf("  Sample size: %d", nrow(data)))
      message(sprintf("  Events: %d", sum(data[[outcome_var]])))
    }
    if (!is.null(wt_var) && verbose) {
      message(sprintf("  Weight summary: Min=%.3f, Median=%.3f, Max=%.3f",
                      min(data[[wt_var]]),
                      stats::median(data[[wt_var]]),
                      max(data[[wt_var]])))
    }

    ###### Stratification info ################
    if (!is.null(stratify_by)) {
      strata_tab <- table(data[[stratify_by]])
      if (verbose) {
        message(sprintf("  Stratification by: %s (%d strata)",
                        stratify_by, length(strata_tab)))
        for (s in names(strata_tab)) {
          message(sprintf("    Stratum '%s': N = %d (%.1f%%)",
                          s, strata_tab[s], 100 * strata_tab[s] / sum(strata_tab)))
        }
      }
    }

    # --- Common: count N and events per group ---
    # Events_* = total events across all follow-up.
    # EventsTP_* = events occurring within the user-specified timepoint window
    # (survival_time <= rr_rd_at_timepoint, same unit convention used by survival::Surv).
    is_ref <- data[[exposure_var]] == ref_val_int
    is_exp <- data[[exposure_var]] == exp_val_int
    in_tp  <- data[[survival_time]] <= rr_rd_at_timepoint

    if (!is.null(wt_var)) {
      wt <- data[[wt_var]]
      n_ref <- sum(wt[is_ref])
      n_exp <- sum(wt[is_exp])
      events_ref <- sum(data[[outcome_var]][is_ref] * wt[is_ref])
      events_exp <- sum(data[[outcome_var]][is_exp] * wt[is_exp])
      eventsTP_ref <- sum(data[[outcome_var]][is_ref & in_tp] * wt[is_ref & in_tp])
      eventsTP_exp <- sum(data[[outcome_var]][is_exp & in_tp] * wt[is_exp & in_tp])
    } else {
      n_ref <- sum(is_ref)
      n_exp <- sum(is_exp)
      events_ref <- sum(data[[outcome_var]][is_ref])
      events_exp <- sum(data[[outcome_var]][is_exp])
      eventsTP_ref <- sum(data[[outcome_var]][is_ref & in_tp])
      eventsTP_exp <- sum(data[[outcome_var]][is_exp & in_tp])
    }

    # ========================================================
    # BRANCH A: Stratified (Standardized) Analysis
    # ========================================================
    if (!is.null(stratify_by)) {

      if (verbose) message("  Calculating standardized point estimates...")

      data_by_stratum <- split(data, data[[stratify_by]])

      # Point estimates with analytical SEs
      std_result <- calc_standardized_rr_rd(
        data = data,
        time_var = survival_time,
        event_var = outcome_var,
        exp_var = exposure_var,
        strat_var = stratify_by,
        weight_var = wt_var,
        timepoint = rr_rd_at_timepoint,
        exp_val = exp_val_int,
        ref_val = ref_val_int,
        per_n = rr_rd_per_individuals,
        pre_split = data_by_stratum,
        risk_estimator = risk_estimator,
        competing_event_var = competing_event_var
      )

      strat_detail_df <- std_result$strat_details
      strat_detail_df$Analysis <- analysis_label
      strat_detail_df$Risk_Estimator <- risk_estimator

      if (verbose) {
        message("\n  Stratum-level risks:")
        for (i in seq_len(nrow(strat_detail_df))) {
          message(sprintf("    Stratum '%s' (w=%.3f): Risk_exp=%.4f, Risk_ref=%.4f",
                          strat_detail_df$stratum[i],
                          strat_detail_df$w_k[i],
                          strat_detail_df$risk_exp[i],
                          strat_detail_df$risk_ref[i]))
        }
        message(sprintf("  Standardized Risk_exp: %.4f, Risk_ref: %.4f",
                        std_result$risk_std_exp, std_result$risk_std_ref))
      }

      # Point estimates
      RR_point <- std_result$RR
      RD_point <- std_result$RD

      # Build cumulative incidence summary
      cuminc_data <- data.frame(
        Exposure_Value = c(as.character(ref_val_int), as.character(exp_val_int)),
        N_Individuals = c(n_ref, n_exp),
        N_Events = c(events_ref, events_exp),
        Timepoint = rr_rd_at_timepoint,
        Risk_SE = c(std_result$se_risk_std_ref, std_result$se_risk_std_exp),
        RiskperN = c(std_result$risk_std_ref * rr_rd_per_individuals,
                     std_result$risk_std_exp * rr_rd_per_individuals),
        Risk = c(std_result$risk_std_ref, std_result$risk_std_exp),
        Risk_Var = c(std_result$se_risk_std_ref^2, std_result$se_risk_std_exp^2),
        Risk_Estimator = risk_estimator,
        Analysis = paste0(analysis_label, " (Standardized)"),
        stringsAsFactors = FALSE
      )

      # --- Analytical CI (stratified) ---
      if (do_analytical) {
        # RR: log-scale CI
        if (!is.na(std_result$se_log_RR) && std_result$se_log_RR > 0) {
          RR_lci_ana <- exp(log(RR_point) - z_val * std_result$se_log_RR)
          RR_uci_ana <- exp(log(RR_point) + z_val * std_result$se_log_RR)
        } else {
          RR_lci_ana <- NA_real_
          RR_uci_ana <- NA_real_
        }
        RR_se_ana <- std_result$lnRR_se_analytical

        # RD: normal approximation
        RD_se_ana <- std_result$RD_se_analytical
        RD_lci_ana <- RD_point - z_val * RD_se_ana
        RD_uci_ana <- RD_point + z_val * RD_se_ana
      }

      # --- Bootstrap CI (stratified) ---
      boot_results <- NULL
      if (do_boot) {
        fixed_strat_weights <- std_result$strat_weights

        if (verbose) {
          message(sprintf("  Running %d stratified bootstrap iterations (%d cores)...",
                          bootstrap_count, n_cores))
        }

        boot_results <- run_parallel_bootstrap(
          n_cores = n_cores,
          bootstrap_count = bootstrap_count,
          boot_fn = boot_standardized_rr_rd_single,
          seed = seed,
          data_by_stratum = data_by_stratum,
          strat_var = stratify_by,
          strat_weights = fixed_strat_weights,
          time_var = survival_time,
          event_var = outcome_var,
          exp_var = exposure_var,
          weight_var = wt_var,
          timepoint = rr_rd_at_timepoint,
          exp_val = exp_val_int,
          ref_val = ref_val_int,
          per_n = rr_rd_per_individuals,
          risk_estimator = risk_estimator,
          competing_event_var = competing_event_var
        )
      }

      # Risk per n from standardized estimates
      risk_ref_per_n <- std_result$risk_std_ref * rr_rd_per_individuals
      risk_exp_per_n <- std_result$risk_std_exp * rr_rd_per_individuals

      analysis_label_out <- paste0(analysis_label, " (Standardized)")

    } else {

      # ========================================================
      # BRANCH B: Unstratified Analysis
      # ========================================================
      if (verbose) message("  Calculating point estimates...")

      cuminc <- calc_km_cumulative_incidence(
        data = data,
        time_var = survival_time,
        event_var = outcome_var,
        exp_var = exposure_var,
        weight_var = wt_var,
        timepoint = rr_rd_at_timepoint,
        per_n = rr_rd_per_individuals,
        risk_estimator = risk_estimator,
        competing_event_var = competing_event_var
      )

      rr_rd_vals <- calc_rr_rd(cuminc, exp_val_int, ref_val_int, rr_rd_per_individuals)
      RR_point <- rr_rd_vals["RR"]
      RD_point <- rr_rd_vals["RD"]

      # Update N from cuminc for unweighted case
      if (is.null(wt_var)) {
        n_ref <- cuminc$N_Individuals[cuminc$Exposure_Value == as.character(ref_val_int)]
        n_exp <- cuminc$N_Individuals[cuminc$Exposure_Value == as.character(exp_val_int)]
        events_ref <- cuminc$N_Events[cuminc$Exposure_Value == as.character(ref_val_int)]
        events_exp <- cuminc$N_Events[cuminc$Exposure_Value == as.character(exp_val_int)]
      }

      # Risk per n from unstratified estimates
      risk_ref_per_n <- rr_rd_vals["risk_ref"] * rr_rd_per_individuals
      risk_exp_per_n <- rr_rd_vals["risk_exp"] * rr_rd_per_individuals

      cuminc$Analysis <- analysis_label
      cuminc$Risk_Estimator <- risk_estimator
      cuminc_data <- cuminc
      strat_detail_df <- NULL

      # --- Analytical CI (unstratified) ---
      if (do_analytical) {
        ana_ci <- calc_analytical_se_unstratified(
          cuminc_data = cuminc,
          exp_val = exp_val_int,
          ref_val = ref_val_int,
          per_n = rr_rd_per_individuals,
          confidence_level = confidence_level
        )
        RR_se_ana  <- ana_ci$lnRR_se
        RR_lci_ana <- ana_ci$RR_lci
        RR_uci_ana <- ana_ci$RR_uci
        RD_se_ana  <- ana_ci$RD_se_per_n
        RD_lci_ana <- ana_ci$RD_lci_per_n
        RD_uci_ana <- ana_ci$RD_uci_per_n
      }

      # --- Bootstrap CI (unstratified) ---
      boot_results <- NULL
      if (do_boot) {
        if (verbose) {
          message(sprintf("  Running %d bootstrap iterations (%d cores)...",
                          bootstrap_count, n_cores))
        }

        boot_results <- run_parallel_bootstrap(
          n_cores = n_cores,
          bootstrap_count = bootstrap_count,
          boot_fn = boot_rr_rd_single,
          seed = seed,
          data = data,
          time_var = survival_time,
          event_var = outcome_var,
          exp_var = exposure_var,
          weight_var = wt_var,
          timepoint = rr_rd_at_timepoint,
          exp_val = exp_val_int,
          ref_val = ref_val_int,
          per_n = rr_rd_per_individuals,
          risk_estimator = risk_estimator,
          competing_event_var = competing_event_var
        )
      }

      analysis_label_out <- analysis_label
    }

    # ========================================================
    # Assemble estimates row (common for both branches)
    # ========================================================
    alpha <- 1 - confidence_level
    ci_probs <- c(alpha / 2, 1 - alpha / 2)

    # --- Bootstrap quantities ---
    if (do_boot && !is.null(boot_results)) {
      RR_lci_boot    <- stats::quantile(boot_results[, "RR"], probs = ci_probs[1], na.rm = TRUE)
      RR_uci_boot    <- stats::quantile(boot_results[, "RR"], probs = ci_probs[2], na.rm = TRUE)
      RR_se_boot     <- stats::sd(boot_results[, "RR"], na.rm = TRUE)
      RD_lci_boot    <- stats::quantile(boot_results[, "RD"], probs = ci_probs[1], na.rm = TRUE)
      RD_uci_boot    <- stats::quantile(boot_results[, "RD"], probs = ci_probs[2], na.rm = TRUE)
      RD_se_boot     <- stats::sd(boot_results[, "RD"], na.rm = TRUE)
      boot_n         <- bootstrap_count
      boot_valid_n   <- sum(!is.na(boot_results[, "RR"]))
    } else {
      RR_lci_boot  <- NA_real_
      RR_uci_boot  <- NA_real_
      RR_se_boot   <- NA_real_
      RD_lci_boot  <- NA_real_
      RD_uci_boot  <- NA_real_
      RD_se_boot   <- NA_real_
      boot_n       <- NA_integer_
      boot_valid_n <- NA_integer_
    }

    # --- Analytical quantities (set NA if not computed) ---
    if (!do_analytical) {
      RR_se_ana  <- NA_real_
      RR_lci_ana <- NA_real_
      RR_uci_ana <- NA_real_
      RD_se_ana  <- NA_real_
      RD_lci_ana <- NA_real_
      RD_uci_ana <- NA_real_
    }

    # --- Formatted CI strings ---
    fmt_ci <- function(est, lcl, ucl, digits = 3) {
      if (is.na(est) || is.na(lcl) || is.na(ucl)) return(NA_character_)
      sprintf("%.*f (%.*f, %.*f)", digits, est, digits, lcl, digits, ucl)
    }

    est_row <- data.frame(
      Analysis         = analysis_label_out,
      Risk_Estimator   = risk_estimator,
      Stratified_By    = ifelse(!is.null(stratify_by), stratify_by, NA_character_),
      Timepoint        = rr_rd_at_timepoint,
      Time_Unit        = time_unit,
      Per_Individuals  = rr_rd_per_individuals,
      N_Ref            = n_ref,
      N_Exp            = n_exp,
      Events_Ref       = events_ref,
      Events_Exp       = events_exp,
      EventsTP_Ref     = eventsTP_ref,
      EventsTP_Exp     = eventsTP_exp,
      RiskperN_Ref     = as.numeric(risk_ref_per_n),
      RiskperN_Exp     = as.numeric(risk_exp_per_n),
      # --- Bootstrap: RR ---
      RR_Boot          = as.numeric(RR_point),
      RR_LCI_Boot      = as.numeric(RR_lci_boot),
      RR_UCI_Boot      = as.numeric(RR_uci_boot),
      RR_SE_Boot       = as.numeric(RR_se_boot),
      RR_CI_Boot       = fmt_ci(RR_point, RR_lci_boot, RR_uci_boot, 2),
      # --- Analytical: RR ---
      RR_Analytical     = as.numeric(RR_point),
      RR_LCI_Analytical = as.numeric(RR_lci_ana),
      RR_UCI_Analytical = as.numeric(RR_uci_ana),
      lnRR_SE_Analytical  = as.numeric(RR_se_ana),
      RR_CI_Analytical  = fmt_ci(RR_point, RR_lci_ana, RR_uci_ana, 2),
      # --- Bootstrap: RD ---
      RDperN_Boot           = as.numeric(RD_point),
      RDperN_LCI_Boot       = as.numeric(RD_lci_boot),
      RDperN_UCI_Boot       = as.numeric(RD_uci_boot),
      RDperN_SE_Boot        = as.numeric(RD_se_boot),
      RDperN_CI_Boot        = fmt_ci(RD_point, RD_lci_boot, RD_uci_boot, 2),
      # --- Analytical: RD ---
      RDperN_Analytical     = as.numeric(RD_point),
      RDperN_LCI_Analytical = as.numeric(RD_lci_ana),
      RDperN_UCI_Analytical = as.numeric(RD_uci_ana),
      RDperN_SE_Analytical  = as.numeric(RD_se_ana),
      RDperN_CI_Analytical  = fmt_ci(RD_point, RD_lci_ana, RD_uci_ana, 2),
      # --- Bootstrap metadata ---
      Bootstrap_N       = boot_n,
      Bootstrap_Valid_N = boot_valid_n,
      stringsAsFactors  = FALSE
    )
    rownames(est_row) <- NULL

    # --- Print summary ---
    if (verbose) {
      if (do_boot && !is.na(RR_lci_boot)) {
        message(sprintf("\n  RR (bootstrap): %s", fmt_ci(RR_point, RR_lci_boot, RR_uci_boot, 2)))
        message(sprintf("  RD per %s (bootstrap): %s",
                        per_label, fmt_ci(RD_point, RD_lci_boot, RD_uci_boot, 2)))
      }
      if (do_analytical && !is.na(RR_lci_ana)) {
        message(sprintf("  RR (analytical): %s", fmt_ci(RR_point, RR_lci_ana, RR_uci_ana, 2)))
        message(sprintf("  RD per %s (analytical): %s",
                        per_label, fmt_ci(RD_point, RD_lci_ana, RD_uci_ana, 2)))
      }
    }

    # Log-transformed CI for Risk (delta method): Risk * exp(∓ z * Risk_SE / Risk)
    cuminc_data$Risk_LCI <- ifelse(
      cuminc_data$Risk > 0,
      cuminc_data$Risk * exp(-z_val * cuminc_data$Risk_SE / cuminc_data$Risk),
      NA_real_
    )
    cuminc_data$Risk_UCI <- ifelse(
      cuminc_data$Risk > 0,
      cuminc_data$Risk * exp(+z_val * cuminc_data$Risk_SE / cuminc_data$Risk),
      NA_real_
    )
    # RiskperN CI: scale Risk CI bounds by per_n
    cuminc_data$RiskperN_LCI <- cuminc_data$Risk_LCI * rr_rd_per_individuals
    cuminc_data$RiskperN_UCI <- cuminc_data$Risk_UCI * rr_rd_per_individuals

    # Ensure consistent column order for rbind compatibility
    cuminc_col_order <- c("Exposure_Value", "N_Individuals", "N_Events", "Timepoint",
                           "Risk", "Risk_LCI", "Risk_UCI", "Risk_SE", "Risk_Var",
                           "RiskperN", "RiskperN_LCI", "RiskperN_UCI",
                           "Risk_Estimator", "Analysis")
    cuminc_data <- cuminc_data[, cuminc_col_order, drop = FALSE]

    return(list(
      est_row = est_row,
      cuminc_data = cuminc_data,
      strat_details = strat_detail_df,
      boot_results = boot_results
    ))
  }

  ##### 3a. Process Unweighted Analysis ########################
  if (do_unweighted) {
    if (verbose) {
      message("Step 1a: Unweighted Analysis")
      message("----------------------------")
    }

    unwt_result <- run_analysis(
      data = in_df_unwt,
      analysis_label = "Unweighted",
      wt_var = NULL
    )

    all_estimates <- rbind(all_estimates, unwt_result$est_row)
    all_cuminc_data <- rbind(all_cuminc_data, unwt_result$cuminc_data)
    if (!is.null(unwt_result$strat_details)) {
      all_stratum_details <- rbind(all_stratum_details, unwt_result$strat_details)
    }
    results$unweighted_bootstrap <- unwt_result$boot_results
  }

  ##### 3b. Process Weighted Analysis ####################
  if (do_weighted) {
    if (verbose) {
      message("Step 1b: Weighted Analysis")
      message("--------------------------")
    }

    wted_result <- run_analysis(
      data = in_df_wted,
      analysis_label = "Weighted",
      wt_var = weight_var
    )

    all_estimates <- rbind(all_estimates, wted_result$est_row)
    all_cuminc_data <- rbind(all_cuminc_data, wted_result$cuminc_data)
    if (!is.null(wted_result$strat_details)) {
      all_stratum_details <- rbind(all_stratum_details, wted_result$strat_details)
    }
    results$weighted_bootstrap <- wted_result$boot_results
  }

  #### 5. Print Summary ################################

  if (verbose) {
    message("\n========================================")
    message("SUMMARY OF RESULTS")
    message(sprintf("(Risk Difference per %s)", per_label))
    message("========================================")

    # Display columns based on conf_int_method
    display_cols <- c("Analysis", "Stratified_By", "Timepoint", "Time_Unit", "Per_Individuals")
    if (do_boot)       display_cols <- c(display_cols, "RR_CI_Boot", "RDperN_CI_Boot")
    if (do_analytical)  display_cols <- c(display_cols, "RR_CI_Analytical", "RDperN_CI_Analytical")

    message("Risk Ratio and Risk Difference Estimates:")
    print(all_estimates[, display_cols, drop = FALSE], row.names = FALSE)
    message("")

    if (nrow(all_stratum_details) > 0L) {
      message("Stratum-Level Details:")
      print(all_stratum_details, row.names = FALSE)
      message("")
    }
  }

  #### 6. Save Results ############################

  if (verbose) {
    message("Step 2: Saving Results")
    message("----------------------")
  }

  if (!is.null(out_xlsxpath)) {
    wb <- openxlsx::createWorkbook()

    # README Sheet
    openxlsx::addWorksheet(wb, "README")

    ci_method_desc <- switch(
      conf_int_method,
      "bootstrap"  = "Percentile bootstrap",
      "analytical" = "Analytical (delta method / Greenwood)",
      "both"       = "Both bootstrap and analytical"
    )

    risk_method_desc <- ifelse(
      risk_estimator == "AJ",
      "Cumulative incidence via Aalen-Johansen estimator (competing risks).",
      "Cumulative incidence via Kaplan-Meier (1 - survival)."
    )

    method_detail <- paste0(
      risk_method_desc, " ",
      ifelse(!is.null(stratify_by),
             paste0("Standardized via direct standardization (total population as standard). ",
                    "Stratum weights w_k = N_k/N_total. ",
                    "Risk_std = Sum(w_k * Risk_k). "),
             ""),
      ifelse(do_analytical,
             paste0("Analytical SE: ",
                    ifelse(risk_estimator == "AJ",
                           "Counting-process variance (Aalen, 1978). ",
                           "Greenwood-based Var(Risk). "),
                    "RR CI via delta method on log scale. ",
                    "RD CI via normal approximation. ",
                    "Risk CI via log transformation (delta method). "),
             ""),
      ifelse(do_boot,
             paste0("Bootstrap CIs via percentile method. ",
                    ifelse(!is.null(stratify_by),
                           "Stratified resampling with fixed population weights. ",
                           ""),
                    sprintf("Parallel RNG: L'Ecuyer-CMRG (seed=%s). ",
                            ifelse(is.null(seed),
                                   "user RNG state",
                                   as.character(seed)))),
             ""),
      sprintf("Risk expressed per %s individuals.", per_label)
    )

    analysis_type_desc <- paste0(
      "Risk Ratio & Risk Difference (", estimator_label,
      ifelse(!is.null(stratify_by), ", Direct Standardization", ""),
      ifelse(risk_estimator == "AJ", ", Competing Risks", ""),
      ")"
    )

    readme_content <- data.frame(
      Item = c(
        "Analysis Type",
        "Date Generated",
        "Outcome Variable",
        "Exposure Variable",
        if (risk_estimator == "AJ") "Competing Event Variable" else character(0),
        "Survival Time Variable",
        "Time Unit",
        "Timepoint for RR/RD",
        "Risk Estimator",
        "Risk Expressed Per",
        "Confidence Level",
        "CI Method",
        if (do_boot) c("Bootstrap Iterations", "Parallel Cores") else character(0),
        "Weight Variable",
        "Stratification Variable",
        "",
        "Method Description",
        "",
        "Notes"
      ),
      Value = c(
        analysis_type_desc,
        as.character(Sys.time()),
        outcome_var,
        paste0(exposure_var, " (Exposed=", exp_value, ", Reference=", ref_value, ")"),
        if (risk_estimator == "AJ") competing_event_var else character(0),
        survival_time,
        time_unit,
        as.character(rr_rd_at_timepoint),
        paste0(risk_estimator, " (", estimator_label, ")"),
        paste0(per_label, " individuals"),
        paste0(confidence_level * 100, "%"),
        ci_method_desc,
        if (do_boot) c(as.character(bootstrap_count), as.character(n_cores)) else character(0),
        ifelse(do_weighted, weight_var, "N/A"),
        ifelse(!is.null(stratify_by), stratify_by, "N/A (Unstratified)"),
        "",
        method_detail,
        "",
        ifelse(is.null(readme_text), "", readme_text)
      ),
      stringsAsFactors = FALSE
    )
    openxlsx::writeData(wb, "README", readme_content)
    openxlsx::setColWidths(wb, "README", cols = 1:2, widths = c(30, 100))

    # Analysis Summary Sheet
    openxlsx::addWorksheet(wb, "Analysis Summary")
    openxlsx::writeDataTable(wb, "Analysis Summary", all_estimates)
    openxlsx::setColWidths(wb, "Analysis Summary", cols = seq_len(ncol(all_estimates)),
                           widths = "auto")

    # Cumulative Incidence Sheet
    openxlsx::addWorksheet(wb, "Cumulative Incidence")
    openxlsx::writeDataTable(wb, "Cumulative Incidence", all_cuminc_data)
    openxlsx::setColWidths(wb, "Cumulative Incidence", cols = seq_len(ncol(all_cuminc_data)),
                           widths = "auto")

    # Stratum Details Sheet (only if stratified)
    if (nrow(all_stratum_details) > 0L) {
      openxlsx::addWorksheet(wb, "Stratum Details")
      openxlsx::writeDataTable(wb, "Stratum Details", all_stratum_details)
      openxlsx::setColWidths(wb, "Stratum Details", cols = seq_len(ncol(all_stratum_details)),
                             widths = "auto")
    }

    # Check if output directory exists
    out_dir <- dirname(out_xlsxpath)
    if (out_dir != "." && !dir.exists(out_dir)) {
      dir.create(out_dir, recursive = TRUE)
      if (verbose) message(sprintf("  Created output directory: %s", out_dir))
    }

    openxlsx::saveWorkbook(wb, out_xlsxpath, overwrite = TRUE)
    if (verbose) message(sprintf("\u2713 Results saved to: %s", out_xlsxpath))
  }

  if (verbose) {
    message("\n========================================")
    message("ANALYSIS COMPLETE")
    message("========================================")
  }

  # Store results
  results$estimates <- all_estimates
  results$cumulative_incidence <- all_cuminc_data
  if (nrow(all_stratum_details) > 0L) {
    results$stratum_details <- all_stratum_details
  }

  return(invisible(results))
}