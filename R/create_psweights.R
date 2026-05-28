# MW OW IPTW STABIPTW ##########
#' 
#' Calculate propensity score weights and optionally generate diagnostic reports.
#' 
#' This function calculates various types of propensity score weights (matching weights,
#' overlap weights, IPTW, stabilized IPTW) and adds them to a dataset that already 
#' contains propensity scores. Optionally generates comprehensive diagnostic reports
#' including assumption checks, balance tables, and plots.
#'
#' @param in_df Data frame containing the input data with PS already calculated (optional if in_csvpath provided)
#' @param in_csvpath Character string. Path to input CSV file with PS already calculated (optional if in_df provided)
#' @param out_csvpath Character string. Path for output CSV file (optional)
#' @param out_xlsxpath_report Character string. Path for output Excel diagnostic report (optional).
#'                            If NULL, no diagnostic report is generated.
#' @param out_dir_plots Character string. Directory to save plot files (optional).
#'                      If NULL, no plots are saved.
#' @param exposure_var Character string. Name of the binary exposure/treatment variable column (default: "exp")
#' @param exp_value Value representing the exposed/treated group (default: 1)
#' @param ref_value Value representing the reference/control group (default: 0)
#' @param ps_var Character string. Name of the PS variable column in the data (default: "ps")
#' @param weight_method Character string. Type of weight to calculate: "mw" (matching weights),
#'                      "ow" (overlap weights), "iptw" (inverse probability weights), 
#'                      "stabiptw" (stabilized IPTW). Only ONE method can be specified.
#' @param weight_var Character string. Name for the weight variable column (default: "psweight")
#' @param estimand Character string. Target estimand: "ATT" (average treatment effect on treated),
#'                 "ATE" (average treatment effect), or "ATO" (average treatment effect in overlap population).
#'                 Only used for IPTW methods. For MW/OW, estimand is determined by the method itself.
#' @param make_unwt_wt_table1 Logical. Whether to generate unweighted and weighted Table 1 (default: FALSE).
#'                            Only used if out_xlsxpath_report is provided.
#' @param table1_cont_vars Character vector. Names of continuous variables for Table 1 (optional, auto-detected if NULL)
#' @param table1_binary_vars Character vector. Names of binary variables for Table 1 (optional, auto-detected if NULL)
#' @param table1_cat_vars Character vector. Names of categorical variables for Table 1 (optional, auto-detected if NULL)
#' @param std_diff_threshold Numeric. Threshold for acceptable standardized difference (default: 0.1)
#' @param readme_text Character string. Optional message to include in README sheet of Excel report
#' @param verbose Logical. Print progress messages (default TRUE).
#'
#' @return A data frame with the weight column added. Also saves to CSV if path specified.
#'
#' @section Side Effects:
#' \itemize{
#'   \item Writes a CSV file when \code{out_csvpath} is provided.
#'   \item Creates directories, writes an Excel diagnostic report, and saves
#'     PNG plot files when the corresponding path arguments are supplied.
#' }
#'
#' @export
#'
#' @examples
#' csv_path <- system.file("extdata", "sample_data.csv", package = "rwetools")
#' df_ps <- estimate_ps(
#'   in_df        = read.csv(csv_path),
#'   exposure_var = "exposure",
#'   class_vars   = c("cat1", "cat2", "cat3", "cat4"),
#'   cont_vars    = c("cont1", "cont2", "cont3"),
#'   verbose      = FALSE
#' )
#'
#' # Add matching weights
#' df_mw <- create_ps_weights(
#'   in_df         = df_ps,
#'   exposure_var  = "exposure",
#'   ps_var        = "ps",
#'   weight_method = "mw",
#'   weight_var    = "mw_wt",
#'   verbose       = FALSE
#' )
#' summary(df_mw$mw_wt)
#'
#' \donttest{
#' # Add IPTW weights with Excel diagnostic report (requires openxlsx, ggplot2)
#' if (requireNamespace("openxlsx", quietly = TRUE) &&
#'     requireNamespace("ggplot2", quietly = TRUE)) {
#'   out_xlsx <- tempfile(fileext = ".xlsx")
#'   out_dir  <- tempdir()
#'   df_iptw  <- create_ps_weights(
#'     in_df               = df_ps,
#'     exposure_var        = "exposure",
#'     ps_var              = "ps",
#'     weight_method       = "iptw",
#'     weight_var          = "iptw_wt",
#'     estimand            = "ATE",
#'     out_xlsxpath_report = out_xlsx,
#'     make_unwt_wt_table1 = TRUE,
#'     out_dir_plots       = out_dir,
#'     verbose             = FALSE
#'   )
#' }
#' }
create_ps_weights <- function(
    in_df = NULL,
    in_csvpath = NULL,
    out_csvpath = NULL,
    out_xlsxpath_report = NULL,
    out_dir_plots = NULL,
    exposure_var = "exp",
    exp_value = 1,
    ref_value = 0,
    ps_var = "ps",
    weight_method = c("mw", "ow", "iptw", "stabiptw"),
    weight_var = "psweight",
    estimand = c("ATT", "ATE", "ATO"),
    make_unwt_wt_table1 = FALSE,
    table1_cont_vars = NULL,
    table1_binary_vars = NULL,
    table1_cat_vars = NULL,
    std_diff_threshold = 0.1,
    readme_text = NULL,
    verbose = TRUE) {
  
  # INPUT VALIDATION ===============================
  # Check data input
  if (is.null(in_df) && is.null(in_csvpath)) {
    stop("Either in_df or in_csvpath must be provided")
  }
  
  if (!is.null(in_df) && !is.null(in_csvpath)) {
    warning("Both in_df and in_csvpath provided. Using in_df.")
  }
  if (!is.null(in_df) && !is.data.frame(in_df)) {
    stop("in_df must be a data frame")
  }
  
  # Validate weight_method - only ONE allowed
  weight_method <- tolower(weight_method[1])  # Take first if multiple provided
  valid_methods <- c("mw", "ow", "iptw", "stabiptw")
  if (!weight_method %in% valid_methods) {
    stop(paste("Invalid weight_method. Must be one of:", paste(valid_methods, collapse = ", ")))
  }
  if (exp_value == ref_value) {
    stop("exp_value and ref_value must be different")
  }
  
  # Validate estimand
  estimand <- match.arg(estimand)
  
  if (verbose) message("\n========================================")
  if (verbose) message("PROPENSITY SCORE WEIGHTING")
  if (verbose) message("========================================")
  
  # STEP 1: LOAD DATA          ######################  
  if (verbose) message("Step 1: Loading Data")
  if (verbose) message("----------------------------------------")
  
  if (!is.null(in_df)) {
    data_work <- in_df
    if (verbose) message("Data loaded from in_df object")
  } else {
    data_work <- utils::read.csv(in_csvpath, stringsAsFactors = FALSE)
    if (verbose) message(sprintf("Data loaded from: %s", in_csvpath))
  }
  
  if (verbose) message(sprintf("Total observations: %d", nrow(data_work)))
  
  # Check if PS variable exists
  if (!ps_var %in% names(data_work)) {
    stop(paste("PS variable", ps_var, "not found in the dataset"))
  }
  if (!is.numeric(data_work[[ps_var]])) {
    stop("PS variable must be numeric")
  }
  if (any(data_work[[ps_var]] < 0 | data_work[[ps_var]] > 1, na.rm = TRUE)) {
    warning("PS values outside [0, 1] detected - results may be unreliable")
  }
  
  # Check if exposure variable exists
  if (!exposure_var %in% names(data_work)) {
    stop(paste("Exposure variable", exposure_var, "not found in the dataset"))
  }
  
  # STEP 2: PREPARE DATA           ######################  
  if (verbose) message("\nStep 2: Data Preparation")
  if (verbose) message("----------------------------------------")
  
  # Recode exposure if needed
  unique_exposure_vals <- unique(data_work[[exposure_var]])
  
  if (!(all(unique_exposure_vals %in% c(0, 1)))) {
    if (verbose) message(sprintf("Recoding exposure: %s -> 1, %s -> 0", exp_value, ref_value))
    data_work[[exposure_var]] <- ifelse(data_work[[exposure_var]] == exp_value, 1,
                                        ifelse(data_work[[exposure_var]] == ref_value, 0, NA))
  }
  
  # Get PS and exposure vectors
  ps <- data_work[[ps_var]]
  exp <- data_work[[exposure_var]]
  
  # Show sample sizes
  n_exp <- sum(exp == 1, na.rm = TRUE)
  n_ref <- sum(exp == 0, na.rm = TRUE)
  if (verbose) message(sprintf("Exposed (n=%d), Unexposed (n=%d)", n_exp, n_ref))
  
  # PS summary by exposure
  if (verbose) message("\nPS Distribution by Exposure:")
  ps_summary <- summarize_ps_by_group(ps, exp)
  if (verbose) print(ps_summary)
  
  # Calculate C-statistic (unweighted - PS model discrimination)
  c_stat_result <- calc_c_statistic(exp, ps, label = "PS Model (unweighted)", verbose = verbose)
  c_stat <- if (!is.null(c_stat_result)) c_stat_result$c_stat else NA_real_
  c_stat_se <- if (!is.null(c_stat_result)) c_stat_result$se else NA_real_
  
  # STEP 3: CALCULATE WEIGHTS           ######################  
  if (verbose) message("\nStep 3: Calculating Weights")
  if (verbose) message("----------------------------------------")
  if (verbose) message(sprintf("Method: %s", toupper(weight_method)))
  if (verbose) message(sprintf("Estimand: %s", estimand))
  if (verbose) message(sprintf("Weight variable: %s", weight_var))
  
  if (weight_method == "mw") {
    # Matching weights (ATT)
    # MW = min(PS, 1-PS) / (exp*PS + (1-exp)*(1-PS))
    mw_numerator <- pmin(ps, 1 - ps)
    mw_denominator <- exp * ps + (1 - exp) * (1 - ps)
    data_work[[weight_var]] <- mw_numerator / mw_denominator
    if (verbose) message("  Matching weights target ATT (Average Treatment on Treated)")
    
  } else if (weight_method == "ow") {
    # Overlap weights (ATO)
    # OW = exp*(1-PS) + (1-exp)*PS
    data_work[[weight_var]] <- exp * (1 - ps) + (1 - exp) * ps
    if (verbose) message("  Overlap weights target ATO (Average Treatment in Overlap)")
    
  } else if (weight_method == "iptw") {
    # Inverse Probability of Treatment Weights
    if (estimand == "ATT") {
      # ATT: weight unexposed to look like exposed
      data_work[[weight_var]] <- ifelse(exp == 1, 1, ps / (1 - ps))
      if (verbose) message("  IPTW for ATT: Exposed get weight=1, Unexposed get PS/(1-PS)")
      
    } else if (estimand == "ATE") {
      # ATE: weight both groups to population
      data_work[[weight_var]] <- ifelse(exp == 1, 1/ps, 1/(1 - ps))
      if (verbose) message("  IPTW for ATE: Exposed get 1/PS, Unexposed get 1/(1-PS)")
      
    } else if (estimand == "ATO") {
      # ATO: weight to overlap population
      data_work[[weight_var]] <- ifelse(exp == 1, 1 - ps, ps)
      if (verbose) message("  IPTW for ATO: Exposed get (1-PS), Unexposed get PS")
    }
    
  } else if (weight_method == "stabiptw") {
    # Stabilized IPTW - need marginal probability
    marginal_prob <- mean(exp, na.rm = TRUE)
    if (verbose) message(sprintf("  Marginal probability of exposure: %.3f", marginal_prob))
    
    if (estimand == "ATT") {
      # Stabilized ATT weights
      data_work[[weight_var]] <- ifelse(exp == 1, 1, 
                                        marginal_prob * ps / ((1 - marginal_prob) * (1 - ps)))
      
    } else if (estimand == "ATE") {
      # Stabilized ATE weights
      data_work[[weight_var]] <- ifelse(exp == 1, 
                                        marginal_prob / ps,
                                        (1 - marginal_prob) / (1 - ps))
      
    } else if (estimand == "ATO") {
      # Stabilized ATO weights (less common)
      ps_overlap <- ps * (1 - ps)
      marginal_overlap <- mean(ps_overlap, na.rm = TRUE)
      data_work[[weight_var]] <- ifelse(exp == 1, 
                                        marginal_overlap * (1 - ps) / ps_overlap,
                                        marginal_overlap * ps / ps_overlap)
    }
    if (verbose) message("  Stabilized IPTW calculated")
  }
  
  # Show weight distribution
  weights <- data_work[[weight_var]]
  
  if (verbose) message("\nWeight Distribution by Exposure:")
  weight_summary <- do.call(rbind, lapply(c(0, 1), function(g) {
    w_g <- weights[exp == g & !is.na(exp)]
    data.frame(
      exposure = g,
      n = length(w_g),
      mean_weight = mean(w_g, na.rm = TRUE),
      sd_weight = stats::sd(w_g, na.rm = TRUE),
      min_weight = min(w_g, na.rm = TRUE),
      q25_weight = stats::quantile(w_g, 0.25, na.rm = TRUE),
      median_weight = stats::median(w_g, na.rm = TRUE),
      q75_weight = stats::quantile(w_g, 0.75, na.rm = TRUE),
      max_weight = max(w_g, na.rm = TRUE)
    )
  }))
  rownames(weight_summary) <- NULL
  if (verbose) print(weight_summary)
  
  # Check for extreme weights
  extreme_weights <- sum(weights > 10, na.rm = TRUE)
  if (extreme_weights > 0) {
    if (verbose) message(sprintf("\nWarning: %d observations have weight > 10", extreme_weights))
  }
  
  # Calculate effective sample size
  if (weight_method %in% c("iptw", "stabiptw")) {
    ess_exposed <- sum(weights[exp == 1], na.rm = TRUE)^2 / 
      sum(weights[exp == 1]^2, na.rm = TRUE)
    ess_unexposed <- sum(weights[exp == 0], na.rm = TRUE)^2 / 
      sum(weights[exp == 0]^2, na.rm = TRUE)
    
    if (verbose) message("\nEffective sample size:")
    if (verbose) message(sprintf("  Exposed: %.1f (%.1f%% of original)", 
                                 ess_exposed, 100 * ess_exposed / n_exp))
    if (verbose) message(sprintf("  Unexposed: %.1f (%.1f%% of original)", 
                                 ess_unexposed, 100 * ess_unexposed / n_ref))
  }
  
  # Calculate post-weighting C-statistic
  c_stat_wt_result <- calc_c_statistic(exp, ps, weights_vec = weights, label = "Post-weighting", verbose = verbose)
  c_stat_weighted <- if (!is.null(c_stat_wt_result)) c_stat_wt_result$c_stat else NA_real_
  c_stat_weighted_se <- if (!is.null(c_stat_wt_result)) c_stat_wt_result$se else NA_real_
  
  # STEP 4: GENERATE DIAGNOSTIC REPORT (if requested)         ######################
  if (!is.null(out_xlsxpath_report)) {
    
    if (verbose) message("\n========================================")
    if (verbose) message("GENERATING DIAGNOSTIC REPORT")
    if (verbose) message("========================================")
    
    # Check required packages
    if (!requireNamespace("openxlsx", quietly = TRUE)) {
      stop("Package 'openxlsx' is required for diagnostic report. Please install it.")
    }
    
    # Create output directory if needed
    out_dir <- dirname(out_xlsxpath_report)
    if (!dir.exists(out_dir) && out_dir != ".") {
      dir.create(out_dir, recursive = TRUE)
      if (verbose) message(sprintf("Created output directory: %s", out_dir))
    }
    
    # Initialize workbook
    wb <- openxlsx::createWorkbook()
    
    # README Sheet (if provided)
    if (!is.null(readme_text)) {
      add_readme_sheet(wb, readme_text, verbose = verbose)
    }
    
    # PS Assumptions Check         ######################
    if (verbose) message("  Checking PS assumptions...")
    
    ps_check <- check_ps_assumptions_internal(data_work, ps_var, exposure_var, verbose = FALSE)
    
    # Create summary dataframe
    assumptions_summary <- data.frame(
      Check = c("Perfect Separation", "Near Separation", "Common Support", 
                "Practical Positivity", "PS Discrimination"),
      Status = c(
        ifelse(ps_check$perfect_separation$detected, "WARNING", "OK"),
        ifelse(ps_check$near_separation$pct_extreme > 5, "WARNING", "OK"),
        ifelse(ps_check$positivity$has_overlap, "OK", "WARNING"),
        ifelse(ps_check$practical_positivity$pct_extreme > 10, "WARNING", "OK"),
        ifelse(abs(ps_check$ps_discrimination$std_diff) < 0.1, "WARNING", "OK")
      ),
      Details = c(
        sprintf("PS=0: %d, PS=1: %d", ps_check$perfect_separation$n_ps_0, 
                ps_check$perfect_separation$n_ps_1),
        sprintf("PS<0.001: %d (%.1f%%), PS>0.999: %d (%.1f%%)", 
                ps_check$near_separation$n_near_0,
                ps_check$near_separation$n_near_0 * 100 / nrow(data_work),
                ps_check$near_separation$n_near_1,
                ps_check$near_separation$n_near_1 * 100 / nrow(data_work)),
        sprintf("Overlap: [%.4f, %.4f], Outside: %d (%.1f%%)", 
                ps_check$positivity$overlap_region[1],
                ps_check$positivity$overlap_region[2],
                ps_check$positivity$n_outside_overlap,
                ps_check$positivity$pct_outside_overlap),
        sprintf("PS<0.1: %d (%.1f%%), PS>0.9: %d (%.1f%%)", 
                ps_check$practical_positivity$n_ps_below_0.1,
                ps_check$practical_positivity$n_ps_below_0.1 * 100 / nrow(data_work),
                ps_check$practical_positivity$n_ps_above_0.9,
                ps_check$practical_positivity$n_ps_above_0.9 * 100 / nrow(data_work)),
        sprintf("Std Diff: %.3f, Mean(Treated): %.3f, Mean(Control): %.3f", 
                ps_check$ps_discrimination$std_diff,
                ps_check$ps_discrimination$mean_ps_treated,
                ps_check$ps_discrimination$mean_ps_control)
      ),
      stringsAsFactors = FALSE
    )
    
    # Positivity details
    positivity_details <- data.frame(
      Metric = c("Treated PS Range", "Control PS Range", "Overlap Region",
                 "N Outside Overlap", "% Outside Overlap"),
      Value = c(
        sprintf("[%.4f, %.4f]", ps_check$positivity$treated_range[1], 
                ps_check$positivity$treated_range[2]),
        sprintf("[%.4f, %.4f]", ps_check$positivity$control_range[1], 
                ps_check$positivity$control_range[2]),
        sprintf("[%.4f, %.4f]", ps_check$positivity$overlap_region[1], 
                ps_check$positivity$overlap_region[2]),
        as.character(ps_check$positivity$n_outside_overlap),
        sprintf("%.1f%%", ps_check$positivity$pct_outside_overlap)
      ),
      stringsAsFactors = FALSE
    )
    
    # Add PS Assumptions sheet
    openxlsx::addWorksheet(wb, "PS_Assumptions")
    openxlsx::writeData(wb, "PS_Assumptions", "PS ASSUMPTIONS CHECK SUMMARY", startCol = 1, startRow = 1)
    openxlsx::writeDataTable(wb, "PS_Assumptions", assumptions_summary, startRow = 3)
    openxlsx::writeData(wb, "PS_Assumptions", "POSITIVITY DETAILS", startCol = 1, startRow = 10)
    openxlsx::writeDataTable(wb, "PS_Assumptions", positivity_details, startRow = 12)
    if (verbose) message("  Added PS Assumptions sheet")
    
    # C-statistics Sheet         ######################
    c_stat_df <- data.frame(
      Measure = c("PS Model C-statistic (unweighted)",
                  "Post-weighting C-statistic"),
      C_statistic = c(
        if (!is.na(c_stat)) sprintf("%.4f", c_stat) else "N/A",
        if (!is.na(c_stat_weighted)) sprintf("%.4f", c_stat_weighted) else "N/A"
      ),
      SE = c(
        if (!is.na(c_stat_se)) sprintf("%.4f", c_stat_se) else "N/A",
        if (!is.na(c_stat_weighted_se)) sprintf("%.4f", c_stat_weighted_se) else "N/A"
      ),
      Interpretation = c(
        "Discrimination of PS model (higher = better separation)",
        "Closer to 0.5 = better balance after weighting"
      ),
      stringsAsFactors = FALSE
    )
    openxlsx::addWorksheet(wb, "C_statistic")
    openxlsx::writeDataTable(wb, "C_statistic", c_stat_df)
    if (verbose) message("  Added C-statistic sheet")
    
    # Weight Distribution Sheet         ######################
    weight_dist_df <- data.frame(
      Group = c("Exposed", "Unexposed"),
      N = c(n_exp, n_ref),
      Mean_Weight = c(mean(weights[exp == 1], na.rm = TRUE),
                      mean(weights[exp == 0], na.rm = TRUE)),
      SD_Weight = c(stats::sd(weights[exp == 1], na.rm = TRUE),
                    stats::sd(weights[exp == 0], na.rm = TRUE)),
      Min_Weight = c(min(weights[exp == 1], na.rm = TRUE),
                     min(weights[exp == 0], na.rm = TRUE)),
      Median_Weight = c(stats::median(weights[exp == 1], na.rm = TRUE),
                        stats::median(weights[exp == 0], na.rm = TRUE)),
      Max_Weight = c(max(weights[exp == 1], na.rm = TRUE),
                     max(weights[exp == 0], na.rm = TRUE)),
      stringsAsFactors = FALSE
    )
    
    # Add ESS for IPTW methods
    if (weight_method %in% c("iptw", "stabiptw")) {
      weight_dist_df$Effective_N <- c(ess_exposed, ess_unexposed)
      weight_dist_df$ESS_Pct <- c(100 * ess_exposed / n_exp, 100 * ess_unexposed / n_ref)
    }
    
    openxlsx::addWorksheet(wb, "Weight_Distribution")
    openxlsx::writeData(wb, "Weight_Distribution", 
                        sprintf("Weight Method: %s, Estimand: %s", toupper(weight_method), estimand),
                        startCol = 1, startRow = 1)
    openxlsx::writeDataTable(wb, "Weight_Distribution", weight_dist_df, startRow = 3)
    if (verbose) message("  Added Weight Distribution sheet")
    
    # Table 1 (if requested)         ######################
    if (make_unwt_wt_table1) {
      if (verbose) message("  Generating Table 1...")
      
      # Auto-detect variable types if not provided
      if (is.null(table1_cont_vars) && is.null(table1_binary_vars) && is.null(table1_cat_vars)) {
        if (verbose) message("    Auto-detecting variable types...")
        if (verbose) message("    WARNING: Please verify variable classifications are correct!")
        
        # Exclude certain variables from auto-detection
        exclude_vars <- c(exposure_var, ps_var, weight_var, "id", "ID", "patient_id", "patid")
        candidate_vars <- setdiff(names(data_work), exclude_vars)
        
        # Use get_var_types function
        var_types <- get_var_types(data_work[candidate_vars], max_cat_levels = 12)
        table1_cont_vars <- var_types$cont_vars
        table1_binary_vars <- var_types$binary_vars
        table1_cat_vars <- var_types$cat_vars
        
        if (verbose) message(sprintf("    Detected: %d continuous, %d binary, %d categorical variables",
                                     length(table1_cont_vars), length(table1_binary_vars), length(table1_cat_vars)))
      }
      
      # Generate Unweighted Table 1
      if (verbose) message("    Creating unweighted Table 1...")
      unweighted_table1 <- build_table1(
        in_df = data_work,
        exposure_var = exposure_var,
        exp_value = 1,
        ref_value = 0,
        cont_vars = table1_cont_vars,
        binary_vars = table1_binary_vars,
        cat_vars = table1_cat_vars,
        use_weights = FALSE,
        verbose = verbose
      )
      
      # Generate Weighted Table 1
      if (verbose) message("    Creating weighted Table 1...")
      weighted_table1 <- build_table1(
        in_df = data_work,
        exposure_var = exposure_var,
        exp_value = 1,
        ref_value = 0,
        cont_vars = table1_cont_vars,
        binary_vars = table1_binary_vars,
        cat_vars = table1_cat_vars,
        use_weights = TRUE,
        weight_var = weight_var,
        verbose = verbose
      )
      
      # Add to workbook
      openxlsx::addWorksheet(wb, "Unweighted_Table1")
      openxlsx::writeDataTable(wb, "Unweighted_Table1", unweighted_table1)
      
      openxlsx::addWorksheet(wb, "Weighted_Table1")
      openxlsx::writeDataTable(wb, "Weighted_Table1", weighted_table1)
      
      if (verbose) message("  Added Unweighted and Weighted Table 1 sheets")
      
      # Balance Comparison Sheet=======================================================
      # Extract standardized differences from both tables
      # Find the Std_diff column (might have different names)
      std_diff_col <- grep("Std_diff|StdDiff|std_diff", names(unweighted_table1), value = TRUE)[1]
      var_col <- grep("Variable|variable", names(unweighted_table1), value = TRUE)[1]
      
      if (!is.null(std_diff_col) && !is.null(var_col)) {
        balance_comparison <- data.frame(
          Variable = unweighted_table1[[var_col]],
          Crude_Std_Diff = as.numeric(gsub("[^0-9.-]", "", unweighted_table1[[std_diff_col]])),
          Weighted_Std_Diff = as.numeric(gsub("[^0-9.-]", "", weighted_table1[[std_diff_col]])),
          stringsAsFactors = FALSE
        )
        
        balance_comparison$Improvement <- abs(balance_comparison$Weighted_Std_Diff) < 
          abs(balance_comparison$Crude_Std_Diff)
        balance_comparison$Balanced <- abs(balance_comparison$Weighted_Std_Diff) < 
          (std_diff_threshold * 100)  # Convert to percentage scale
        
        openxlsx::addWorksheet(wb, "Balance_Comparison")
        openxlsx::writeDataTable(wb, "Balance_Comparison", balance_comparison)
        if (verbose) message("  Added Balance Comparison sheet")
      }
    }
    
    # Generate Plots (if directory provided)         ######################
    if (!is.null(out_dir_plots)) {
      if (verbose) message("  Generating diagnostic plots...")
      
      if (!requireNamespace("ggplot2", quietly = TRUE)) {
        warning("Package 'ggplot2' is required for plots. Skipping plot generation.")
      } else {
        # Create plot directory
        if (!dir.exists(out_dir_plots)) {
          dir.create(out_dir_plots, recursive = TRUE)
        }
        
        # Prepare data for plotting
        plot_data <- data.frame(
          ps = ps,
          exposure = factor(exp, levels = c(0, 1), labels = c("Control", "Treated")),
          weight = weights
        )
        
        # Plot 1: Unweighted PS Distribution
        p_ps_unweighted <- ggplot2::ggplot(plot_data, ggplot2::aes(x = ps, fill = exposure)) +
          ggplot2::geom_histogram(alpha = 0.6, position = "identity", bins = 30) +
          ggplot2::facet_wrap(~exposure, ncol = 1, scales = "free_y") +
          ggplot2::labs(title = "Unweighted Propensity Score Distribution",
                        x = "Propensity Score", y = "Count", fill = "Group") +
          ggplot2::theme_minimal() +
          ggplot2::theme(legend.position = "bottom") +
          ggplot2::scale_fill_manual(values = c("Control" = "#377eb8", "Treated" = "#e41a1c"))
        
        ggplot2::ggsave(file.path(out_dir_plots, "ps_distribution_unweighted.png"),
                        plot = p_ps_unweighted, width = 10, height = 6, dpi = 150)
        
        # Plot 2: Weighted PS Distribution
        p_ps_weighted <- ggplot2::ggplot(plot_data, ggplot2::aes(x = ps, weight = weight, fill = exposure)) +
          ggplot2::geom_histogram(alpha = 0.6, position = "identity", bins = 30) +
          ggplot2::facet_wrap(~exposure, ncol = 1, scales = "free_y") +
          ggplot2::labs(title = "Weighted Propensity Score Distribution",
                        x = "Propensity Score", y = "Weighted Count", fill = "Group") +
          ggplot2::theme_minimal() +
          ggplot2::theme(legend.position = "bottom") +
          ggplot2::scale_fill_manual(values = c("Control" = "#377eb8", "Treated" = "#e41a1c"))
        
        ggplot2::ggsave(file.path(out_dir_plots, "ps_distribution_weighted.png"),
                        plot = p_ps_weighted, width = 10, height = 6, dpi = 150)
        
        # Plot 3: Weight Distribution Box Plot
        y_limit <- stats::quantile(weights, 0.99, na.rm = TRUE) * 1.1
        
        p_weight_box <- ggplot2::ggplot(plot_data, ggplot2::aes(x = exposure, y = weight, fill = exposure)) +
          ggplot2::geom_boxplot(alpha = 0.7) +
          ggplot2::geom_hline(yintercept = 1, linetype = "dashed", color = "gray50") +
          ggplot2::labs(title = "Distribution of Propensity Score Weights",
                        x = "Group", y = "Weight") +
          ggplot2::theme_minimal() +
          ggplot2::theme(legend.position = "none") +
          ggplot2::scale_fill_manual(values = c("Control" = "#377eb8", "Treated" = "#e41a1c")) +
          ggplot2::coord_cartesian(ylim = c(0, y_limit))
        
        ggplot2::ggsave(file.path(out_dir_plots, "weight_distribution_boxplot.png"),
                        plot = p_weight_box, width = 8, height = 6, dpi = 150)
        
        # Plot 4: Balance/Love Plot (only if Table 1 was generated)
        if (make_unwt_wt_table1 && exists("balance_comparison")) {
          balance_plot_data <- balance_comparison[!is.na(balance_comparison$Crude_Std_Diff) & 
                                                    !is.na(balance_comparison$Weighted_Std_Diff), ]
          if (nrow(balance_plot_data) > 0) {
            create_love_plot(
              variable_names = balance_plot_data$Variable,
              crude_std_diff = balance_plot_data$Crude_Std_Diff,
              adjusted_std_diff = balance_plot_data$Weighted_Std_Diff,
              crude_label = "Unweighted", adjusted_label = "Weighted",
              title = "Standardized Differences in Covariates",
              output_path = file.path(out_dir_plots, "balance_love_plot.png"),
              colors = c("Unweighted" = "#377eb8", "Weighted" = "#e41a1c")
            )
          }
        }
        
        if (verbose) message(sprintf("  Saved plots to: %s", out_dir_plots))
      }
    }
    
    # Save workbook
    openxlsx::saveWorkbook(wb, out_xlsxpath_report, overwrite = TRUE)
    if (verbose) message(sprintf("\nDiagnostic report saved to: %s", out_xlsxpath_report))
  }
  
  # STEP 5: SAVE OUTPUTS #############################
  if (verbose) message("\n========================================")
  if (verbose) message("SAVING OUTPUTS")
  if (verbose) message("========================================")
  
  # Save to CSV if path specified
  if (!is.null(out_csvpath)) {
    out_dir <- dirname(out_csvpath)
    if (!dir.exists(out_dir) && out_dir != ".") {
      dir.create(out_dir, recursive = TRUE)
    }
    utils::write.csv(data_work, file = out_csvpath, row.names = FALSE)
    if (verbose) message(sprintf("Data saved to CSV: %s", out_csvpath))
  }
  
  if (verbose) message("\n========================================")
  if (verbose) message("PS WEIGHTING COMPLETE")
  if (verbose) message("========================================")
  
  # Return the data with weights added
  return(invisible(data_work))
}

#### Check PS Assumptions (INTERNAL HELPER for create_ps_weights) ####
#' Check propensity score assumptions
#'
#' Internal function used by \code{\link{create_ps_weights}} and other PS
#' functions to run diagnostic checks on propensity scores, including perfect
#' separation, positivity, overlap, and covariate balance.
#'
#' @param data Data frame containing the propensity score and exposure columns.
#' @param ps_var Character string. Name of the propensity score column.
#' @param exposure Character string. Name of the binary exposure column.
#' @param verbose Logical. Print progress messages (default TRUE).
#' @return A list of assumption-check results, each element containing a
#'   \code{detected} flag and relevant summary statistics.
#' @keywords internal
check_ps_assumptions_internal <- function(data, ps_var, exposure, verbose = TRUE) {
  
  results <- list()
  
  ps_vals <- data[[ps_var]]
  exp_vals <- data[[exposure]]
  
  if (verbose) {
    if (verbose) message("\n========================================")
    if (verbose) message("PROPENSITY SCORE ASSUMPTIONS CHECK")
    if (verbose) message("========================================")
  }
  
  # 1. Check for perfect separation
  if (verbose) message("1. Perfect Separation Check:")
  
  n_ps_0 <- sum(ps_vals == 0, na.rm = TRUE)
  n_ps_1 <- sum(ps_vals == 1, na.rm = TRUE)
  
  results$perfect_separation <- list(
    detected = (n_ps_0 > 0 || n_ps_1 > 0),
    n_ps_0 = n_ps_0,
    n_ps_1 = n_ps_1
  )
  
  if (verbose) {
    if (results$perfect_separation$detected) {
      if (verbose) message(sprintf("  WARNING: Perfect separation detected!"))
      if (verbose) message(sprintf("    - %d observations with PS = 0", n_ps_0))
      if (verbose) message(sprintf("    - %d observations with PS = 1", n_ps_1))
    } else {
      if (verbose) message("  OK: No perfect separation detected")
    }
  }
  
  # 2. Check for near separation
  if (verbose) message("\n2. Near Separation Check:")
  
  n_near_0 <- sum(ps_vals < 0.001, na.rm = TRUE)
  n_near_1 <- sum(ps_vals > 0.999, na.rm = TRUE)
  
  results$near_separation <- list(
    n_near_0 = n_near_0,
    n_near_1 = n_near_1,
    pct_extreme = 100 * (n_near_0 + n_near_1) / length(ps_vals)
  )
  
  if (verbose) {
    if (n_near_0 > 0 || n_near_1 > 0) {
      if (verbose) message(sprintf("  Near separation detected:"))
      if (verbose) message(sprintf("    - %d observations with PS < 0.001 (%.1f%%)", 
                                   n_near_0, 100*n_near_0/length(ps_vals)))
      if (verbose) message(sprintf("    - %d observations with PS > 0.999 (%.1f%%)", 
                                   n_near_1, 100*n_near_1/length(ps_vals)))
    } else {
      if (verbose) message("  OK: No near separation detected")
    }
  }
  
  # 3. Check positivity assumption (common support)
  if (verbose) message("\n3. Positivity/Common Support Check:")
  
  # Split PS by exposure
  ps_treated <- ps_vals[exp_vals == 1 & !is.na(exp_vals)]
  ps_control <- ps_vals[exp_vals == 0 & !is.na(exp_vals)]
  
  # Calculate overlap region
  min_treated <- min(ps_treated, na.rm = TRUE)
  max_treated <- max(ps_treated, na.rm = TRUE)
  min_control <- min(ps_control, na.rm = TRUE)
  max_control <- max(ps_control, na.rm = TRUE)
  
  overlap_lower <- max(min_treated, min_control)
  overlap_upper <- min(max_treated, max_control)
  
  # Check if there's overlap
  has_overlap <- overlap_lower < overlap_upper
  
  # Calculate percentage outside overlap
  n_outside <- sum(ps_vals < overlap_lower | ps_vals > overlap_upper, na.rm = TRUE)
  pct_outside <- 100 * n_outside / length(ps_vals)
  
  results$positivity <- list(
    has_overlap = has_overlap,
    overlap_region = c(lower = overlap_lower, upper = overlap_upper),
    treated_range = c(min = min_treated, max = max_treated),
    control_range = c(min = min_control, max = max_control),
    n_outside_overlap = n_outside,
    pct_outside_overlap = pct_outside
  )
  
  if (verbose) {
    if (has_overlap) {
      if (verbose) message(sprintf("  OK: Common support exists"))
      if (verbose) message(sprintf("    Overlap region: [%.4f, %.4f]", overlap_lower, overlap_upper))
      if (verbose) message(sprintf("    Treated PS range: [%.4f, %.4f]", min_treated, max_treated))
      if (verbose) message(sprintf("    Control PS range: [%.4f, %.4f]", min_control, max_control))
      if (verbose) message(sprintf("    Observations outside overlap: %d (%.1f%%)", n_outside, pct_outside))
    } else {
      if (verbose) message(sprintf("  WARNING: No common support!"))
      if (verbose) message(sprintf("    Treated PS range: [%.4f, %.4f]", min_treated, max_treated))
      if (verbose) message(sprintf("    Control PS range: [%.4f, %.4f]", min_control, max_control))
    }
  }
  
  # 4. Check for practical positivity violations
  if (verbose) message("\n4. Practical Positivity Check:")
  
  # Check PS < 0.1 or > 0.9
  n_low_ps <- sum(ps_vals < 0.1, na.rm = TRUE)
  n_high_ps <- sum(ps_vals > 0.9, na.rm = TRUE)
  
  results$practical_positivity <- list(
    n_ps_below_0.1 = n_low_ps,
    n_ps_above_0.9 = n_high_ps,
    pct_extreme = 100 * (n_low_ps + n_high_ps) / length(ps_vals)
  )
  
  if (verbose) {
    if (verbose) message(sprintf("    PS < 0.1: %d observations (%.1f%%)", 
                                 n_low_ps, 100*n_low_ps/length(ps_vals)))
    if (verbose) message(sprintf("    PS > 0.9: %d observations (%.1f%%)", 
                                 n_high_ps, 100*n_high_ps/length(ps_vals)))
    if (results$practical_positivity$pct_extreme > 10) {
      if (verbose) message("  Warning: >10% of observations have extreme PS values")
    }
  }
  
  # 5. Check balance of PS between groups
  if (verbose) message("\n5. PS Discrimination Check:")
  
  mean_ps_treated <- mean(ps_treated, na.rm = TRUE)
  mean_ps_control <- mean(ps_control, na.rm = TRUE)
  sd_ps_treated <- stats::sd(ps_treated, na.rm = TRUE)
  sd_ps_control <- stats::sd(ps_control, na.rm = TRUE)
  
  # Standardized difference in PS
  pooled_sd <- sqrt((sd_ps_treated^2 + sd_ps_control^2) / 2)
  std_diff_ps <- (mean_ps_treated - mean_ps_control) / pooled_sd
  
  results$ps_discrimination <- list(
    mean_ps_treated = mean_ps_treated,
    mean_ps_control = mean_ps_control,
    std_diff = std_diff_ps
  )
  
  if (verbose) {
    if (verbose) message(sprintf("    Mean PS in Treated: %.3f (SD: %.3f)", mean_ps_treated, sd_ps_treated))
    if (verbose) message(sprintf("    Mean PS in Control: %.3f (SD: %.3f)", mean_ps_control, sd_ps_control))
    if (verbose) message(sprintf("    Standardized difference: %.3f", std_diff_ps))
    if (abs(std_diff_ps) < 0.1) {
      if (verbose) message("  Warning: PS poorly discriminates between groups (std diff < 0.1)")
    }
  }
  
  return(results)
}


# MATCHING create_ps_matched_cohort ############ 
#' 
#' Perform propensity score matching and optionally generate diagnostic reports.
#' 
#' This function performs 1:k nearest neighbor PS matching using MatchIt and 
#' generates comprehensive diagnostic reports including assumption checks, 
#' balance tables, and plots.
#'
#' @param in_df Data frame containing the input data with PS already calculated (optional if in_csvpath provided)
#' @param in_csvpath Character string. Path to input CSV file with PS already calculated (optional if in_df provided)
#' @param out_csvpath_matcheddata Character string. Path for matched cohort CSV file (optional)
#' @param out_csvpath_crudedata_w_matchindicator Character string. Path for crude data with match indicator CSV file (optional)
#' @param out_xlsxpath_report Character string. Path for output Excel diagnostic report (optional).
#' @param out_dir_plots Character string. Directory to save plot files (optional).
#' @param exposure_var Character string. Name of the binary exposure/treatment variable column (default: "exp")
#' @param exp_value Value representing the exposed/treated group (default: 1)
#' @param ref_value Value representing the reference/control group (default: 0)
#' @param ps_var Character string. Name of the PS variable column in the data (default: "ps")
#' @param ratio Integer. Matching ratio (1:k). Default is 1 for 1:1 matching.
#' @param caliper Numeric. Caliper width in standard deviations of logit(PS). Default is 0.2.
#' @param replace Logical. Whether to match with replacement. Default is FALSE.
#' @param make_crude_matched_table1 Logical. Whether to generate crude and matched Table 1 (default: FALSE).
#' @param table1_cont_vars Character vector. Names of continuous variables for Table 1 (optional, auto-detected if NULL)
#' @param table1_binary_vars Character vector. Names of binary variables for Table 1 (optional, auto-detected if NULL)
#' @param table1_cat_vars Character vector. Names of categorical variables for Table 1 (optional, auto-detected if NULL)
#' @param std_diff_threshold Numeric. Threshold for acceptable standardized difference (default: 0.1)
#' @param readme_text Character string. Optional message to include in README sheet of Excel report
#' @param verbose Logical. Print progress messages (default TRUE).
#'
#' @return A data frame containing the matched cohort only (crude data with match indicator
#'   can be saved to CSV via out_csvpath_crudedata_w_matchindicator).
#'
#' @section Side Effects:
#' \itemize{
#'   \item Writes matched-cohort and/or crude-data CSV files when the
#'     corresponding path arguments are provided.
#'   \item Creates directories, writes an Excel diagnostic report, and saves
#'     PNG plot files when the corresponding path arguments are supplied.
#' }
#'
#' @export
#'
#' @examples
#' \donttest{
#' # Requires MatchIt
#' if (requireNamespace("MatchIt", quietly = TRUE)) {
#'   csv_path <- system.file("extdata", "sample_data.csv", package = "rwetools")
#'   df_ps <- estimate_ps(
#'     in_df        = read.csv(csv_path),
#'     exposure_var = "exposure",
#'     class_vars   = c("cat1", "cat2", "cat3", "cat4"),
#'     cont_vars    = c("cont1", "cont2", "cont3"),
#'     verbose      = FALSE
#'   )
#'
#'   matched <- create_ps_matched_cohort(
#'     in_df        = df_ps,
#'     exposure_var = "exposure",
#'     ps_var       = "ps",
#'     ratio        = 1,
#'     caliper      = 0.2,
#'     verbose      = FALSE
#'   )
#'   nrow(matched)
#' }
#' }

create_ps_matched_cohort <- function(
    in_df = NULL,
    in_csvpath = NULL,
    out_csvpath_matcheddata = NULL,
    out_csvpath_crudedata_w_matchindicator = NULL,
    out_xlsxpath_report = NULL,
    out_dir_plots = NULL,
    exposure_var = "exp",
    exp_value = 1,
    ref_value = 0,
    ps_var = "ps",
    ratio = 1,
    caliper = 0.2,
    replace = FALSE,
    make_crude_matched_table1 = FALSE,
    table1_cont_vars = NULL,
    table1_binary_vars = NULL,
    table1_cat_vars = NULL,
    std_diff_threshold = 0.1,
    readme_text = NULL,
    verbose = TRUE) {
  
  # INPUT VALIDATION         ######################  
  # Check data input
  if (is.null(in_df) && is.null(in_csvpath)) {
    stop("Either in_df or in_csvpath must be provided")
  }
  
  if (!is.null(in_df) && !is.null(in_csvpath)) {
    warning("Both in_df and in_csvpath provided. Using in_df.")
  }
  if (!is.null(in_df) && !is.data.frame(in_df)) {
    stop("in_df must be a data frame")
  }
  if (!is.numeric(ratio) || length(ratio) != 1 || ratio < 1) {
    stop("ratio must be a single number >= 1")
  }
  if (!is.numeric(caliper) || length(caliper) != 1 || caliper <= 0) {
    stop("caliper must be a single positive number")
  }
  if (exp_value == ref_value) {
    stop("exp_value and ref_value must be different")
  }
  
  # Check MatchIt package
  if (!requireNamespace("MatchIt", quietly = TRUE)) {
    stop("Package 'MatchIt' is required. Please install it with: install.packages('MatchIt')")
  }
  
  if (verbose) message("\n========================================")
  if (verbose) message("PROPENSITY SCORE MATCHING")
  if (verbose) message("========================================")
  
  # STEP 1: LOAD DATA         ######################  
  if (verbose) message("Step 1: Loading Data")
  if (verbose) message("----------------------------------------")
  
  if (!is.null(in_df)) {
    data_work <- as.data.frame(in_df)
    if (verbose) message("Data loaded from in_df object")
  } else {
    data_work <- utils::read.csv(in_csvpath, stringsAsFactors = FALSE)
    if (verbose) message(sprintf("Data loaded from: %s", in_csvpath))
  }
  
  # Add row ID for tracking
  data_work$.row_id <- seq_len(nrow(data_work))
  
  if (verbose) message(sprintf("Total observations: %d", nrow(data_work)))
  
  # Check if PS variable exists
  if (!ps_var %in% names(data_work)) {
    stop(paste("PS variable", ps_var, "not found in the dataset"))
  }
  if (!is.numeric(data_work[[ps_var]])) {
    stop("PS variable must be numeric")
  }
  if (any(data_work[[ps_var]] < 0 | data_work[[ps_var]] > 1, na.rm = TRUE)) {
    warning("PS values outside [0, 1] detected \u2014 results may be unreliable")
  }
  
  # Check if exposure variable exists
  if (!exposure_var %in% names(data_work)) {
    stop(paste("Exposure variable", exposure_var, "not found in the dataset"))
  }
  
  # STEP 2: PREPARE DATA         ######################  
  if (verbose) message("\nStep 2: Data Preparation")
  if (verbose) message("----------------------------------------")
  
  # Create binary exposure variable for MatchIt
  unique_exposure_vals <- unique(data_work[[exposure_var]])
  
  # Create .treat variable (1 = treated, 0 = control)
  data_work$.treat <- ifelse(data_work[[exposure_var]] == exp_value, 1L,
                             ifelse(data_work[[exposure_var]] == ref_value, 0L, NA_integer_))
  
  # Remove observations with missing exposure
  n_missing_exp <- sum(is.na(data_work$.treat))
  if (n_missing_exp > 0) {
    if (verbose) message(sprintf("Removing %d observations with missing exposure", n_missing_exp))
    data_work <- data_work[!is.na(data_work$.treat), ]
  }
  
  # Get PS vector
  ps <- data_work[[ps_var]]
  
  # Check for missing PS
  n_missing_ps <- sum(is.na(ps))
  if (n_missing_ps > 0) {
    if (verbose) message(sprintf("Removing %d observations with missing PS", n_missing_ps))
    data_work <- data_work[!is.na(data_work[[ps_var]]), ]
    ps <- data_work[[ps_var]]
  }
  
  # Show sample sizes
  n_treated <- sum(data_work$.treat == 1)
  n_control <- sum(data_work$.treat == 0)
  if (verbose) message(sprintf("Treated (n=%d), Control (n=%d)", n_treated, n_control))
  
  # PS summary by exposure
  if (verbose) message("\nPS Distribution by Exposure (Pre-matching):")
  ps_summary_pre <- summarize_ps_by_group(ps, data_work$.treat,
                                          group_labels = c("0" = "Control", "1" = "Treated"))
  if (verbose) print(ps_summary_pre)
  
  # Calculate C-statistic (pre-matching - PS model discrimination)
  c_stat_result <- calc_c_statistic(data_work$.treat, ps, label = "PS Model (pre-matching)", verbose = verbose)
  c_stat <- if (!is.null(c_stat_result)) c_stat_result$c_stat else NA_real_
  c_stat_se <- if (!is.null(c_stat_result)) c_stat_result$se else NA_real_
  
  # STEP 3: PERFORM MATCHING         ######################
  if (verbose) message("\nStep 3: Performing PS Matching")
  if (verbose) message("----------------------------------------")
  if (verbose) message(sprintf("Method: Nearest Neighbor Matching"))
  if (verbose) message(sprintf("Ratio: 1:%d", ratio))
  if (verbose) message(sprintf("Caliper: %.2f SD of logit(PS)", caliper))
  if (verbose) message(sprintf("With replacement: %s", ifelse(replace, "Yes", "No")))
  if (verbose) message(sprintf("Estimand: ATT (Average Treatment Effect on Treated)"))
  
  # Create formula for MatchIt (using pre-computed PS via distance argument)
  match_formula <- stats::as.formula(paste(".treat ~ 1"))
  
  # Run MatchIt
  match_out <- MatchIt::matchit(
    formula = match_formula,
    data = data_work,
    method = "nearest",
    distance = data_work[[ps_var]],  # Use pre-computed PS
    ratio = ratio,
    caliper = caliper,
    std.caliper = TRUE,  # Caliper in SD units of logit(PS)
    replace = replace,
    estimand = "ATT"
  )
  
  # Get matching summary
  match_summary_obj <- summary(match_out)
  
  # Extract matched data
  matched_data <- MatchIt::match.data(match_out)
  
  # STEP 4: CREATE MATCH INDICATORS         ######################  
  if (verbose) message("\nStep 4: Creating Match Indicators")
  if (verbose) message("----------------------------------------")
  
  # Initialize match indicators in original data
  data_work$matched <- 0L
  data_work$match_id <- NA_character_
  
  # Get match matrix from MatchIt
  match_matrix <- match_out$match.matrix
  
  # Process match pairs
  if (!is.null(match_matrix)) {
    treated_ids <- rownames(match_matrix)
    
    match_counter <- 1
    for (i in seq_along(treated_ids)) {
      treated_row <- as.integer(treated_ids[i])
      control_rows <- as.integer(match_matrix[i, ])
      control_rows <- control_rows[!is.na(control_rows)]
      
      if (length(control_rows) > 0) {
        # Create unique match_id for this pair/set
        mid <- sprintf("M%05d", match_counter)
        
        # Mark treated
        data_work$matched[treated_row] <- 1L
        data_work$match_id[treated_row] <- mid
        
        # Mark controls
        for (ctrl_row in control_rows) {
          data_work$matched[ctrl_row] <- 1L
          # For replacement matching, a control might be matched multiple times
          if (is.na(data_work$match_id[ctrl_row])) {
            data_work$match_id[ctrl_row] <- mid
          } else {
            # Append additional match_id
            data_work$match_id[ctrl_row] <- paste(data_work$match_id[ctrl_row], mid, sep = ";")
          }
        }
        
        match_counter <- match_counter + 1
      }
    }
  }
  
  # Create final matched cohort (subset)
  matched_cohort <- data_work[data_work$matched == 1, ]
  
  # For matched cohort, add weights from MatchIt if available
  if ("weights" %in% names(matched_data)) {
    # Merge weights back
    weight_df <- data.frame(
      .row_id = as.integer(rownames(matched_data)),
      .match_weights = matched_data$weights
    )
    matched_cohort <- merge(matched_cohort, weight_df, by = ".row_id", all.x = TRUE)
    matched_cohort$.match_weights[is.na(matched_cohort$.match_weights)] <- 1
  } else {
    matched_cohort$.match_weights <- 1
  }
  
  # Show matching results
  n_matched_treated <- sum(matched_cohort$.treat == 1)
  n_matched_control <- sum(matched_cohort$.treat == 0)
  n_unmatched_treated <- n_treated - n_matched_treated
  n_unmatched_control <- n_control - n_matched_control
  
  if (verbose) message(sprintf("\nMatching Results:"))
  if (verbose) message(sprintf("  Treated: %d matched, %d unmatched (%.1f%% matched)", 
                               n_matched_treated, n_unmatched_treated, 100 * n_matched_treated / n_treated))
  if (verbose) message(sprintf("  Control: %d matched, %d unmatched (%.1f%% matched)", 
                               n_matched_control, n_unmatched_control, 100 * n_matched_control / n_control))
  if (verbose) message(sprintf("  Total matched cohort: %d", nrow(matched_cohort)))
  if (verbose) message(sprintf("  Unique match pairs: %d", length(unique(stats::na.omit(data_work$match_id)))))
  
  # PS summary in matched cohort
  if (verbose) message("\nPS Distribution by Exposure (Post-matching):")
  
  ps_matched <- matched_cohort[[ps_var]]
  treat_matched <- matched_cohort$.treat
  
  ps_summary_post <- summarize_ps_by_group(ps_matched, treat_matched,
                                           group_labels = c("0" = "Control", "1" = "Treated"))
  if (verbose) print(ps_summary_post)
  
  # Calculate post-matching C-statistic
  c_stat_post_result <- calc_c_statistic(treat_matched, ps_matched,
                                         weights_vec = matched_cohort$.match_weights,
                                         label = "Post-matching",
                                         verbose = verbose)
  c_stat_post <- if (!is.null(c_stat_post_result)) c_stat_post_result$c_stat else NA_real_
  c_stat_post_se <- if (!is.null(c_stat_post_result)) c_stat_post_result$se else NA_real_
  
  # STEP 5: GENERATE DIAGNOSTIC REPORT (if requested)         ######################
  
  if (!is.null(out_xlsxpath_report)) {
    
    if (verbose) message("\n========================================")
    if (verbose) message("GENERATING DIAGNOSTIC REPORT")
    if (verbose) message("========================================")
    
    # Check required packages
    if (!requireNamespace("openxlsx", quietly = TRUE)) {
      stop("Package 'openxlsx' is required for diagnostic report. Please install it.")
    }
    
    # Create output directory if needed
    out_dir <- dirname(out_xlsxpath_report)
    if (!dir.exists(out_dir) && out_dir != ".") {
      dir.create(out_dir, recursive = TRUE)
      if (verbose) message(sprintf("Created output directory: %s", out_dir))
    }
    
    # Initialize workbook
    wb <- openxlsx::createWorkbook()
    
    # README Sheet (if provided)
    if (!is.null(readme_text)) {
      add_readme_sheet(wb, readme_text, verbose = verbose)
    }
    
    # Matching Summary Sheet
    matching_summary_df <- data.frame(
      Metric = c(
        "Matching Method",
        "Ratio",
        "Caliper (SD of logit PS)",
        "With Replacement",
        "Estimand",
        "",
        "Pre-matching Treated",
        "Pre-matching Control",
        "Pre-matching Total",
        "",
        "Post-matching Treated",
        "Post-matching Control",
        "Post-matching Total",
        "",
        "Unmatched Treated",
        "Unmatched Control",
        "Percent Treated Matched",
        "Percent Control Matched",
        "",
        "Unique Match Pairs"
      ),
      Value = c(
        "Nearest Neighbor",
        sprintf("1:%d", ratio),
        sprintf("%.2f", caliper),
        ifelse(replace, "Yes", "No"),
        "ATT",
        "",
        as.character(n_treated),
        as.character(n_control),
        as.character(n_treated + n_control),
        "",
        as.character(n_matched_treated),
        as.character(n_matched_control),
        as.character(nrow(matched_cohort)),
        "",
        as.character(n_unmatched_treated),
        as.character(n_unmatched_control),
        sprintf("%.1f%%", 100 * n_matched_treated / n_treated),
        sprintf("%.1f%%", 100 * n_matched_control / n_control),
        "",
        as.character(length(unique(stats::na.omit(data_work$match_id))))
      ),
      stringsAsFactors = FALSE
    )
    
    openxlsx::addWorksheet(wb, "Matching_Summary")
    openxlsx::writeDataTable(wb, "Matching_Summary", matching_summary_df)
    if (verbose) message("  Added Matching Summary sheet")
    
    # PS Distribution Sheet
    ps_dist_df <- rbind(
      cbind(Sample = "Pre-matching", ps_summary_pre),
      cbind(Sample = "Post-matching", ps_summary_post)
    )
    
    openxlsx::addWorksheet(wb, "PS_Distribution")
    openxlsx::writeDataTable(wb, "PS_Distribution", ps_dist_df)
    if (verbose) message("  Added PS Distribution sheet")
    
    # C-statistic Sheet
    c_stat_df <- data.frame(
      Measure = c("PS Model C-statistic (pre-matching)",
                  "Post-matching C-statistic"),
      C_statistic = c(
        if (!is.na(c_stat)) sprintf("%.4f", c_stat) else "N/A",
        if (!is.na(c_stat_post)) sprintf("%.4f", c_stat_post) else "N/A"
      ),
      SE = c(
        if (!is.na(c_stat_se)) sprintf("%.4f", c_stat_se) else "N/A",
        if (!is.na(c_stat_post_se)) sprintf("%.4f", c_stat_post_se) else "N/A"
      ),
      Interpretation = c(
        "Discrimination of PS model (higher = better separation)",
        "Closer to 0.5 = better balance after matching"
      ),
      stringsAsFactors = FALSE
    )
    openxlsx::addWorksheet(wb, "C_statistic")
    openxlsx::writeDataTable(wb, "C_statistic", c_stat_df)
    if (verbose) message("  Added C-statistic sheet")
    
    # Table 1 (if requested)
    if (make_crude_matched_table1) {
      if (verbose) message("  Generating Table 1...")
      
      # Check if build_table1 and get_var_types functions exist
      if (!exists("build_table1") || !exists("get_var_types")) {
        warning("build_table1 or get_var_types function not found. Skipping Table 1 generation.")
      } else {
        
        # Auto-detect variable types if not provided
        if (is.null(table1_cont_vars) && is.null(table1_binary_vars) && is.null(table1_cat_vars)) {
          if (verbose) message("    Auto-detecting variable types...")
          if (verbose) message("    WARNING: Please verify variable classifications are correct!")
          
          # Exclude certain variables from auto-detection
          exclude_vars <- c(exposure_var, ps_var, ".treat", ".row_id", "matched", 
                            "match_id", ".match_weights", "weights", "subclass",
                            "id", "ID", "patient_id", "patid")
          candidate_vars <- setdiff(names(data_work), exclude_vars)
          
          # Use get_var_types function
          var_types <- get_var_types(data_work[candidate_vars], max_cat_levels = 12)
          table1_cont_vars <- var_types$cont_vars
          table1_binary_vars <- var_types$binary_vars
          table1_cat_vars <- var_types$cat_vars
          
          if (verbose) message(sprintf("    Detected: %d continuous, %d binary, %d categorical variables",
                                       length(table1_cont_vars), length(table1_binary_vars), length(table1_cat_vars)))
        }
        
        # Prepare data for Table 1
        # Crude: use original exposure variable
        crude_table1_data <- data_work
        crude_table1_data[[exposure_var]] <- ifelse(crude_table1_data$.treat == 1, exp_value, ref_value)
        
        # Matched: use matched cohort
        matched_table1_data <- matched_cohort
        matched_table1_data[[exposure_var]] <- ifelse(matched_table1_data$.treat == 1, exp_value, ref_value)
        
        # Generate Crude Table 1
        if (verbose) message("    Creating crude (pre-matching) Table 1...")
        crude_table1 <- tryCatch({
          build_table1(
            in_df = crude_table1_data,
            exposure_var = exposure_var,
            exp_value = exp_value,
            ref_value = ref_value,
            cont_vars = table1_cont_vars,
            binary_vars = table1_binary_vars,
            cat_vars = table1_cat_vars,
            use_weights = FALSE,
            verbose = verbose
          )
        }, error = function(e) {
          warning(paste("Error creating crude Table 1:", e$message))
          NULL
        })
        
        # Generate Matched Table 1
        if (verbose) message("    Creating matched (post-matching) Table 1...")
        matched_table1 <- tryCatch({
          build_table1(
            in_df = matched_table1_data,
            exposure_var = exposure_var,
            exp_value = exp_value,
            ref_value = ref_value,
            cont_vars = table1_cont_vars,
            binary_vars = table1_binary_vars,
            cat_vars = table1_cat_vars,
            use_weights = FALSE,
            verbose = verbose
          )
        }, error = function(e) {
          warning(paste("Error creating matched Table 1:", e$message))
          NULL
        })
        
        # Add to workbook
        if (!is.null(crude_table1)) {
          openxlsx::addWorksheet(wb, "Crude_Table1")
          openxlsx::writeDataTable(wb, "Crude_Table1", crude_table1)
          if (verbose) message("  Added Crude Table 1 sheet")
        }
        
        if (!is.null(matched_table1)) {
          openxlsx::addWorksheet(wb, "Matched_Table1")
          openxlsx::writeDataTable(wb, "Matched_Table1", matched_table1)
          if (verbose) message("  Added Matched Table 1 sheet")
        }
        
        # Balance Comparison Sheet
        if (!is.null(crude_table1) && !is.null(matched_table1)) {
          # Extract standardized differences from both tables
          std_diff_col <- grep("Std_diff|StdDiff|std_diff", names(crude_table1), value = TRUE)[1]
          var_col <- grep("Variable|variable", names(crude_table1), value = TRUE)[1]
          
          if (!is.null(std_diff_col) && !is.null(var_col)) {
            balance_comparison <- data.frame(
              Variable = crude_table1[[var_col]],
              Crude_Std_Diff = as.numeric(gsub("[^0-9.-]", "", crude_table1[[std_diff_col]])),
              Matched_Std_Diff = as.numeric(gsub("[^0-9.-]", "", matched_table1[[std_diff_col]])),
              stringsAsFactors = FALSE
            )
            
            balance_comparison$Improvement <- abs(balance_comparison$Matched_Std_Diff) < 
              abs(balance_comparison$Crude_Std_Diff)
            balance_comparison$Balanced <- abs(balance_comparison$Matched_Std_Diff) < 
              (std_diff_threshold * 100)  # Convert to percentage scale if needed
            
            openxlsx::addWorksheet(wb, "Balance_Comparison")
            openxlsx::writeDataTable(wb, "Balance_Comparison", balance_comparison)
            if (verbose) message("  Added Balance Comparison sheet")
          }
        }
      }
    }
    
    # Generate Plots (if directory provided)
    if (!is.null(out_dir_plots)) {
      if (verbose) message("  Generating diagnostic plots...")
      
      if (!requireNamespace("ggplot2", quietly = TRUE)) {
        warning("Package 'ggplot2' is required for plots. Skipping plot generation.")
      } else {
        # Create plot directory
        if (!dir.exists(out_dir_plots)) {
          dir.create(out_dir_plots, recursive = TRUE)
        }
        
        # Prepare data for plotting
        # Pre-matching (crude) data
        plot_data_crude <- data.frame(
          ps = data_work[[ps_var]],
          group = factor(data_work$.treat, levels = c(0, 1), labels = c("Control", "Treated")),
          sample = "Pre-matching"
        )
        
        # Post-matching data
        plot_data_matched <- data.frame(
          ps = matched_cohort[[ps_var]],
          group = factor(matched_cohort$.treat, levels = c(0, 1), labels = c("Control", "Treated")),
          sample = "Post-matching"
        )
        
        plot_data_combined <- rbind(plot_data_crude, plot_data_matched)
        plot_data_combined$sample <- factor(plot_data_combined$sample, 
                                            levels = c("Pre-matching", "Post-matching"))
        
        # Plot 1: PS Distribution - Pre vs Post matching (faceted)
        p_ps_comparison <- ggplot2::ggplot(plot_data_combined, 
                                           ggplot2::aes(x = ps, fill = group)) +
          ggplot2::geom_histogram(alpha = 0.6, position = "identity", bins = 30) +
          ggplot2::facet_grid(sample ~ group, scales = "free_y") +
          ggplot2::labs(title = "Propensity Score Distribution: Pre vs Post Matching",
                        x = "Propensity Score", y = "Count", fill = "Group") +
          ggplot2::theme_minimal() +
          ggplot2::theme(legend.position = "bottom") +
          ggplot2::scale_fill_manual(values = c("Control" = "#377eb8", "Treated" = "#e41a1c"))
        
        ggplot2::ggsave(file.path(out_dir_plots, "ps_distribution_comparison.png"),
                        plot = p_ps_comparison, width = 10, height = 8, dpi = 150)
        
        # Plot 2: PS Distribution - Mirrored histogram (Pre-matching)
        p_ps_mirror_pre <- ggplot2::ggplot(plot_data_crude, ggplot2::aes(x = ps, fill = group)) +
          ggplot2::geom_histogram(data = subset(plot_data_crude, group == "Treated"),
                                  ggplot2::aes(y = ggplot2::after_stat(count)),
                                  alpha = 0.7, bins = 30) +
          ggplot2::geom_histogram(data = subset(plot_data_crude, group == "Control"),
                                  ggplot2::aes(y = -ggplot2::after_stat(count)),
                                  alpha = 0.7, bins = 30) +
          ggplot2::labs(title = "PS Distribution (Pre-matching)",
                        x = "Propensity Score", y = "Count", fill = "Group") +
          ggplot2::theme_minimal() +
          ggplot2::theme(legend.position = "bottom") +
          ggplot2::scale_fill_manual(values = c("Control" = "#377eb8", "Treated" = "#e41a1c")) +
          ggplot2::geom_hline(yintercept = 0, color = "gray30")
        
        ggplot2::ggsave(file.path(out_dir_plots, "ps_distribution_pre_mirror.png"),
                        plot = p_ps_mirror_pre, width = 10, height = 6, dpi = 150)
        
        # Plot 3: PS Distribution - Mirrored histogram (Post-matching)
        p_ps_mirror_post <- ggplot2::ggplot(plot_data_matched, ggplot2::aes(x = ps, fill = group)) +
          ggplot2::geom_histogram(data = subset(plot_data_matched, group == "Treated"),
                                  ggplot2::aes(y = ggplot2::after_stat(count)),
                                  alpha = 0.7, bins = 30) +
          ggplot2::geom_histogram(data = subset(plot_data_matched, group == "Control"),
                                  ggplot2::aes(y = -ggplot2::after_stat(count)),
                                  alpha = 0.7, bins = 30) +
          ggplot2::labs(title = "PS Distribution (Post-matching)",
                        x = "Propensity Score", y = "Count", fill = "Group") +
          ggplot2::theme_minimal() +
          ggplot2::theme(legend.position = "bottom") +
          ggplot2::scale_fill_manual(values = c("Control" = "#377eb8", "Treated" = "#e41a1c")) +
          ggplot2::geom_hline(yintercept = 0, color = "gray30")
        
        ggplot2::ggsave(file.path(out_dir_plots, "ps_distribution_post_mirror.png"),
                        plot = p_ps_mirror_post, width = 10, height = 6, dpi = 150)
        
        # Plot 4: Balance/Love Plot (only if Table 1 was generated)
        if (make_crude_matched_table1 && exists("balance_comparison")) {
          balance_plot_data <- balance_comparison[!is.na(balance_comparison$Crude_Std_Diff) & 
                                                    !is.na(balance_comparison$Matched_Std_Diff), ]
          if (nrow(balance_plot_data) > 0) {
            create_love_plot(
              variable_names = balance_plot_data$Variable,
              crude_std_diff = balance_plot_data$Crude_Std_Diff,
              adjusted_std_diff = balance_plot_data$Matched_Std_Diff,
              crude_label = "Pre-matching", adjusted_label = "Post-matching",
              title = "Standardized Differences: Pre vs Post Matching",
              output_path = file.path(out_dir_plots, "balance_love_plot.png"),
              colors = c("Pre-matching" = "#377eb8", "Post-matching" = "#e41a1c")
            )
          }
        }
        
        if (verbose) message(sprintf("  Saved plots to: %s", out_dir_plots))
      }
    }
    
    # Save workbook
    openxlsx::saveWorkbook(wb, out_xlsxpath_report, overwrite = TRUE)
    if (verbose) message(sprintf("\nDiagnostic report saved to: %s", out_xlsxpath_report))
  }
  
  # STEP 6: CLEAN UP AND SAVE OUTPUTS         ######################  
  if (verbose) message("\n========================================")
  if (verbose) message("SAVING OUTPUTS")
  if (verbose) message("========================================")
  
  # Clean up internal columns from output data
  internal_cols <- c(".treat", ".row_id")
  
  # Prepare matched cohort output
  matched_output <- matched_cohort[, !names(matched_cohort) %in% internal_cols, drop = FALSE]
  
  # Prepare crude data with indicator output
  crude_output <- data_work[, !names(data_work) %in% c(".treat"), drop = FALSE]
  # Keep .row_id for reference or remove it
  crude_output <- crude_output[, !names(crude_output) %in% c(".row_id"), drop = FALSE]
  
  # Save matched data to CSV if path specified
  if (!is.null(out_csvpath_matcheddata)) {
    out_dir <- dirname(out_csvpath_matcheddata)
    if (!dir.exists(out_dir) && out_dir != ".") {
      dir.create(out_dir, recursive = TRUE)
    }
    utils::write.csv(matched_output, file = out_csvpath_matcheddata, row.names = FALSE)
    if (verbose) message(sprintf("Matched data saved to CSV: %s", out_csvpath_matcheddata))
  }
  
  # Save crude data with indicator to CSV if path specified
  if (!is.null(out_csvpath_crudedata_w_matchindicator)) {
    out_dir <- dirname(out_csvpath_crudedata_w_matchindicator)
    if (!dir.exists(out_dir) && out_dir != ".") {
      dir.create(out_dir, recursive = TRUE)
    }
    utils::write.csv(crude_output, file = out_csvpath_crudedata_w_matchindicator, row.names = FALSE)
    if (verbose) message(sprintf("Crude data with match indicator saved to CSV: %s", 
                                 out_csvpath_crudedata_w_matchindicator))
  }
  
  if (verbose) message("\n========================================")
  if (verbose) message("PS MATCHING COMPLETE")
  if (verbose) message("========================================")
  
  # Return matched cohort data frame
  return(invisible(matched_output))
}



# FS create_ps_fs_weights  ################################
#' Calculate PS Fine Stratification Weights, Trim Non-overlapping Regions,
#' and Generate Diagnostics
#'
#' Combined function that (1) creates PS fine strata, (2) calculates stratification
#' weights, (3) trims non-overlapping PS regions, and optionally (4) builds a
#' crude vs weighted balance table (Table 1) and (5) generates diagnostic plots.
#'
#' @param in_df Data frame with PS already calculated (the "crude" / untrimmed data)
#' @param out_csvpath Character. Path for output CSV of trimmed+weighted data (optional)
#' @param out_xlsxpath_report Character. Path for Excel diagnostic report (optional)
#' @param out_dir_plots Character. Directory to save diagnostic plots (NULL = skip plots)
#' @param trim_nonoverlap_region Logical. Trim non-overlapping PS regions (default TRUE)
#' @param exposure_var Character. Name of binary exposure variable (default "exp")
#' @param exp_value Value representing exposed group (default 1)
#' @param ref_value Value representing reference group (default 0)
#' @param ps_var Character. Name of PS variable in the data (required)
#' @param weight_var Character. Name for the stratification weight variable (default "ps_fs_wt")
#' @param number_of_strata Integer. Number of strata to create (default 50)
#' @param stratification_method Character. "exposure" or "cohort" (default "exposure")
#' @param estimand Character. "ATT" or "ATE" (default "ATT")
#' @param make_unwt_wt_table1 Logical. Build crude vs weighted balance table (default FALSE)
#' @param table1_cont_vars Character vector. Continuous vars for Table 1 (auto-detected if NULL)
#' @param table1_binary_vars Character vector. Binary vars for Table 1 (auto-detected if NULL)
#' @param table1_cat_vars Character vector. Categorical vars for Table 1 (auto-detected if NULL)
#' @param std_diff_threshold Numeric. Threshold for acceptable standardised difference (default 0.1)
#' @param readme_text Character. Optional message for README sheet in Excel report
#' @param verbose Logical. Print progress messages (default TRUE).
#'
#' @return Invisibly returns the trimmed+weighted data frame.
#'
#' @section Side Effects:
#' \itemize{
#'   \item Writes a CSV file when \code{out_csvpath} is provided.
#'   \item Creates directories, writes an Excel diagnostic report, and saves
#'     PNG plot files when the corresponding path arguments are supplied.
#' }
#'
#' @export
#'
#' @examples
#' csv_path <- system.file("extdata", "sample_data.csv", package = "rwetools")
#' df_ps <- estimate_ps(
#'   in_df        = read.csv(csv_path),
#'   exposure_var = "exposure",
#'   class_vars   = c("cat1", "cat2", "cat3", "cat4"),
#'   cont_vars    = c("cont1", "cont2", "cont3"),
#'   verbose      = FALSE
#' )
#'
#' result <- create_ps_fs_weights(
#'   in_df                 = df_ps,
#'   exposure_var          = "exposure",
#'   ps_var                = "ps",
#'   weight_var            = "fs_wt",
#'   number_of_strata      = 10,
#'   stratification_method = "exposure",
#'   estimand              = "ATT",
#'   verbose               = FALSE
#' )
#' summary(result$fs_wt)
create_ps_fs_weights <- function(
    in_df,
    out_csvpath = NULL,
    out_xlsxpath_report = NULL,
    out_dir_plots = NULL,
    trim_nonoverlap_region = TRUE,
    exposure_var = "exp",
    exp_value = 1,
    ref_value = 0,
    ps_var = NULL,
    weight_var = "ps_fs_wt",
    number_of_strata = 50,
    stratification_method = c("exposure", "cohort"),
    estimand = c("ATT", "ATE"),
    make_unwt_wt_table1 = FALSE,
    table1_cont_vars = NULL,
    table1_binary_vars = NULL,
    table1_cat_vars = NULL,
    std_diff_threshold = 0.1,
    readme_text = NULL,
    verbose = TRUE) {
  
  # Package loading 
    if (!is.null(out_xlsxpath_report)) {
    if (!requireNamespace("openxlsx", quietly = TRUE)) stop("Package 'openxlsx' is required for Excel diagnostic reports.")
  }
    if (!is.null(out_dir_plots)) {
    if (!requireNamespace("ggplot2", quietly = TRUE)) stop("Package 'ggplot2' is required for plots.")
  }
  
  #  Validate inputs 
  if (is.null(ps_var))                    stop("ps_var must be specified")
  if (!ps_var %in% names(in_df))          stop(paste("PS variable", ps_var, "not found in dataset"))
  if (!exposure_var %in% names(in_df))    stop(paste("Exposure variable", exposure_var, "not found in dataset"))
  if (!is.numeric(in_df[[ps_var]]))        stop("PS variable must be numeric")
  if (any(in_df[[ps_var]] < 0 | in_df[[ps_var]] > 1, na.rm = TRUE)) {
    warning("PS values outside [0, 1] detected - results may be unreliable")
  }
  
  stratification_method <- match.arg(stratification_method)
  estimand              <- match.arg(estimand)
  
  # Handle empty input 
  if (nrow(in_df) == 0) {
    warning("Input data has 0 rows. Returning empty output.")
    return(invisible(in_df))
  }
  
  if (verbose) message("\n========================================")
  if (verbose) message("PS FINE STRATIFICATION WITH TRIMMING")
  if (verbose) message("========================================")
  
  original_n <- nrow(in_df)
  data_work  <- in_df                        # working copy (will be trimmed)
  in_crude   <- in_df                        # keep untouched copy for diagnostics
  
  # Recode exposure to 0/1 if needed 
  uv <- unique(data_work[[exposure_var]])
  if (!(all(uv %in% c(0, 1)))) {
    if (verbose) message("Recoding exposure variable...")
    data_work[[exposure_var]] <- ifelse(data_work[[exposure_var]] == exp_value, 1L,
                                        ifelse(data_work[[exposure_var]] == ref_value, 0L, NA_integer_))
    in_crude[[exposure_var]]  <- ifelse(in_crude[[exposure_var]] == exp_value, 1L,
                                        ifelse(in_crude[[exposure_var]] == ref_value, 0L, NA_integer_))
  }
  
  # Create internal row-order ID (used to restore order after rbind in exposure stratification)
  # Generate a collision-safe column name
  .row_id_col <- ".ps_fs_row_order"
  while (.row_id_col %in% names(data_work)) .row_id_col <- paste0(.row_id_col, "_")
  data_work[[.row_id_col]] <- seq_len(nrow(data_work))
  
  results <- list()
  
  
  #### STEP 1: Trim Non-overlapping PS Regions ########################################################
  if (trim_nonoverlap_region) {
    if (verbose) message("Step 1: Trimming Non-overlapping PS Regions")
    if (verbose) message("--------------------------------------------")
    
    exp_idx  <- which(data_work[[exposure_var]] == 1)
    ref_idx  <- which(data_work[[exposure_var]] == 0)
    ps_vals  <- data_work[[ps_var]]
    
    min_exp <- min(ps_vals[exp_idx], na.rm = TRUE)
    max_exp <- max(ps_vals[exp_idx], na.rm = TRUE)
    min_ref <- min(ps_vals[ref_idx], na.rm = TRUE)
    max_ref <- max(ps_vals[ref_idx], na.rm = TRUE)
    
    overlap_lower <- max(min_exp, min_ref)
    overlap_upper <- min(max_exp, max_ref)
    
    if (overlap_lower >= overlap_upper) {
      warning("No overlap in propensity score distributions!")
      results$trim_summary <- list(overlap = FALSE,
                                   bounds  = c(lower = overlap_lower, upper = overlap_upper),
                                   n_trimmed = 0L)
    } else {
      in_overlap <- ps_vals >= overlap_lower & ps_vals <= overlap_upper
      
      # Build trimming summary table (base R)
      trimmed_flag <- ifelse(in_overlap, "No", "Yes")
      grp          <- data_work[[exposure_var]]
      trim_tab     <- as.data.frame(table(trimmed = trimmed_flag, exposure = grp))
      # Reshape wide
      trim_wide <- stats::reshape(trim_tab, direction = "wide", idvar = "trimmed",
                                  timevar = "exposure", v.names = "Freq")
      names(trim_wide) <- gsub("Freq\\.", "exp_", names(trim_wide))
      # ensure both columns present
      if (!"exp_0" %in% names(trim_wide)) trim_wide$exp_0 <- 0L
      if (!"exp_1" %in% names(trim_wide)) trim_wide$exp_1 <- 0L
      trim_wide$total     <- trim_wide$exp_0 + trim_wide$exp_1
      trim_wide$pct_exp_0 <- round(100 * trim_wide$exp_0 / sum(trim_wide$exp_0), 2)
      trim_wide$pct_exp_1 <- round(100 * trim_wide$exp_1 / sum(trim_wide$exp_1), 2)
      trim_wide$pct_total <- round(100 * trim_wide$total / sum(trim_wide$total), 2)
      
      if (verbose) message(sprintf("Original N: %d", nrow(data_work)))
      if (verbose) message(sprintf("Overlap bounds: [%.6f, %.6f]", overlap_lower, overlap_upper))
      if (verbose) print(trim_wide)
      
      data_trimmed <- data_work[in_overlap, , drop = FALSE]
      n_trimmed    <- nrow(data_work) - nrow(data_trimmed)
      
      if (verbose) message(sprintf("\nAfter trimming N: %d", nrow(data_trimmed)))
      if (verbose) message(sprintf("Trimmed N: %d (%.1f%%)\n", n_trimmed, 100 * n_trimmed / nrow(data_work)))
      
      results$trim_summary <- list(
        overlap    = TRUE,
        bounds     = c(lower = overlap_lower, upper = overlap_upper),
        n_trimmed  = n_trimmed,
        trim_table = trim_wide
      )
      data_work <- data_trimmed
    }
  } else {
    if (verbose) message("Step 1: Skipping trimming (trim_nonoverlap_region = FALSE)")
    results$trim_summary <- list(overlap = NA, bounds = NA, n_trimmed = 0L)
  }
  
  #### STEP 2: Create PS Strata ##########################################################
  if (verbose) message("Step 2: Creating PS Strata")
  if (verbose) message("--------------------------")
  if (verbose) message(sprintf("Method: %s", stratification_method))
  if (verbose) message(sprintf("Number of strata: %d", number_of_strata))
  
  if (stratification_method == "exposure") {
    # Stratify based on exposed group PS distribution 
    exp_mask  <- data_work[[exposure_var]] == 1
    ref_mask  <- data_work[[exposure_var]] == 0
    
    exposed_data   <- data_work[exp_mask, , drop = FALSE]
    unexposed_data <- data_work[ref_mask, , drop = FALSE]
    
    # Sort exposed by PS
    exposed_data <- exposed_data[order(exposed_data[[ps_var]]), , drop = FALSE]
    
    n_exp <- nrow(exposed_data)
    if (n_exp < number_of_strata) {
      stop(sprintf("Not enough exposed subjects (%d) for %d strata", n_exp, number_of_strata))
    }
    
    # Rank stratum (mimicking SAS PROC RANK)
    ranks <- floor(rank(exposed_data[[ps_var]], ties.method = "first") /
                     (n_exp / number_of_strata))
    ranks[ranks == number_of_strata] <- number_of_strata - 1L
    exposed_data$strata <- ranks + 1L
    
    # Strata boundaries from exposed
    strata_mins <- tapply(exposed_data[[ps_var]], exposed_data$strata, min, na.rm = TRUE)
    cut_points  <- c(0, as.numeric(strata_mins[-1]), Inf)
    
    # Assign strata to unexposed
    unexposed_data$strata <- as.integer(cut(
      unexposed_data[[ps_var]],
      breaks = cut_points,
      labels = seq_len(number_of_strata),
      right = FALSE, include.lowest = TRUE
    ))
    
    stratified_data <- rbind(exposed_data, unexposed_data)
    stratified_data <- stratified_data[order(stratified_data[[.row_id_col]]), , drop = FALSE]
    
  } else {
    # Stratify based on entire cohort distribution 
    stratified_data <- data_work[order(data_work[[ps_var]]), , drop = FALSE]
    n_total <- nrow(stratified_data)
    ranks   <- floor(rank(stratified_data[[ps_var]], ties.method = "first") /
                       (n_total / number_of_strata))
    ranks[ranks == number_of_strata] <- number_of_strata - 1L
    stratified_data$strata <- ranks + 1L
  }
  
  # Strata boundaries (final)
  strata_ps_min <- tapply(stratified_data[[ps_var]], stratified_data$strata, min, na.rm = TRUE)
  strata_ps_max <- tapply(stratified_data[[ps_var]], stratified_data$strata, max, na.rm = TRUE)
  strata_ps_mean <- tapply(stratified_data[[ps_var]], stratified_data$strata, mean, na.rm = TRUE)
  strata_n      <- tapply(stratified_data[[ps_var]], stratified_data$strata, length)
  
  if (verbose) message(sprintf("\nStrata created: %d", length(unique(stratified_data$strata))))
  
  #### STEP 3: Calculate Stratification Weights ##########################################################
  if (verbose) message("\nStep 3: Calculating Stratification Weights")
  if (verbose) message("-------------------------------------------")
  if (verbose) message(sprintf("Estimand: %s", estimand))
  
  # Stratum-specific counts (base R)
  exp_vec <- stratified_data[[exposure_var]]
  str_vec <- stratified_data$strata
  
  total_exp_by_str   <- tapply(exp_vec == 1, str_vec, sum, na.rm = TRUE)
  total_unexp_by_str <- tapply(exp_vec == 0, str_vec, sum, na.rm = TRUE)
  
  strata_summary <- data.frame(
    strata      = as.integer(names(total_exp_by_str)),
    n_exp       = as.integer(total_exp_by_str),
    n_ref       = as.integer(total_unexp_by_str),
    stringsAsFactors = FALSE
  )
  
  # Keep only strata with both groups
  complete_mask    <- strata_summary$n_exp > 0 & strata_summary$n_ref > 0
  n_dropped_strata <- sum(!complete_mask)
  if (n_dropped_strata > 0) {
    if (verbose) message(sprintf("Warning: %d strata with only one exposure group dropped", n_dropped_strata))
  }
  strata_complete <- strata_summary[complete_mask, , drop = FALSE]
  
  sum_exp   <- sum(strata_complete$n_exp)
  sum_unexp <- sum(strata_complete$n_ref)
  sum_total <- sum_exp + sum_unexp
  
  # Calculate weights per stratum
  if (estimand == "ATT") {
    strata_complete$weight_exp <- 1
    strata_complete$weight_ref <- (strata_complete$n_exp / sum_exp) /
      (strata_complete$n_ref / sum_unexp)
    strata_complete$att_weighted_unexp <- strata_complete$n_ref * strata_complete$weight_ref
  } else {
    strata_complete$weight_exp <- ((strata_complete$n_exp + strata_complete$n_ref) / sum_total) /
      (strata_complete$n_exp / sum_exp)
    strata_complete$weight_ref <- ((strata_complete$n_exp + strata_complete$n_ref) / sum_total) /
      (strata_complete$n_ref / sum_unexp)
    strata_complete$ate_weighted_exp   <- strata_complete$n_exp * strata_complete$weight_exp
    strata_complete$ate_weighted_unexp <- strata_complete$n_ref * strata_complete$weight_ref
  }
  
  # Build lookup: strata -> weight for exp & unexp
  wt_exp_lookup   <- stats::setNames(strata_complete$weight_exp, strata_complete$strata)
  wt_unexp_lookup <- stats::setNames(strata_complete$weight_ref, strata_complete$strata)
  
  # Assign weights
  str_char <- as.character(stratified_data$strata)
  stratified_data[[weight_var]] <- ifelse(
    stratified_data[[exposure_var]] == 1,
    wt_exp_lookup[str_char],
    wt_unexp_lookup[str_char]
  )
  
  # Remove observations with NA weights (from incomplete strata)
  n_before   <- nrow(stratified_data)
  final_data <- stratified_data[!is.na(stratified_data[[weight_var]]), , drop = FALSE]
  n_removed  <- n_before - nrow(final_data)
  if (n_removed > 0) {
    if (verbose) message(sprintf("\n%d observations removed (strata with only one exposure group)", n_removed))
  }
  
  # Add PS statistics to strata summary
  strata_complete$n_total  <- strata_complete$n_exp + strata_complete$n_ref
  strata_complete$mean_ps  <- as.numeric(strata_ps_mean[as.character(strata_complete$strata)])
  strata_complete$min_ps   <- as.numeric(strata_ps_min[as.character(strata_complete$strata)])
  strata_complete$max_ps   <- as.numeric(strata_ps_max[as.character(strata_complete$strata)])
  
  # Reorder columns: strata, n_total, n_exp, n_ref, mean_ps, min_ps, max_ps, weight_exp, weight_ref
  strata_complete <- strata_complete[, c("strata", "n_total", "n_exp", "n_ref", 
                                         "mean_ps", "min_ps", "max_ps", 
                                         "weight_exp", "weight_ref"), drop = FALSE]
  
  # Overall summary
  if (estimand == "ATT") {
    overall_summary <- data.frame(
      strata                = "Whole sample",
      n_total_exp           = sum_exp,
      n_total_unexp         = sum_unexp,
      total_att_weighted_unexp = round(sum(strata_complete$n_ref * strata_complete$weight_ref)),
      stringsAsFactors = FALSE
    )
  } else {
    overall_summary <- data.frame(
      strata                   = "Whole sample",
      n_total_exp              = sum_exp,
      n_total_unexp            = sum_unexp,
      total_ate_weighted_exp   = round(sum(strata_complete$n_exp * strata_complete$weight_exp)),
      total_ate_weighted_unexp = round(sum(strata_complete$n_ref * strata_complete$weight_ref)),
      stringsAsFactors = FALSE
    )
  }
  
  results$strata_summary  <- strata_complete
  results$overall_summary <- overall_summary
  
  # Print summary
  if (verbose) message("\nFinal Sample Size:")
  if (verbose) message(sprintf("  Original N: %d", original_n))
  if (trim_nonoverlap_region) {
    if (verbose) message(sprintf("  After trimming: %d", nrow(data_work)))
  }
  if (verbose) message(sprintf("  After stratification: %d", nrow(final_data)))
  if (verbose) message(sprintf("  Number of complete strata: %d", nrow(strata_complete)))
  
  if (verbose) message("\nFirst few strata:")
  if (verbose) print(utils::head(strata_complete))
  
  # Weight distribution summary (base R)
  wt_vals <- final_data[[weight_var]]
  exp_grp <- final_data[[exposure_var]]
  weight_summary <- do.call(rbind, lapply(sort(unique(exp_grp)), function(g) {
    w <- wt_vals[exp_grp == g]
    data.frame(
      exposure    = g,
      n           = length(w),
      mean_weight = mean(w, na.rm = TRUE),
      sd_weight   = stats::sd(w, na.rm = TRUE),
      min_weight  = min(w, na.rm = TRUE),
      max_weight  = max(w, na.rm = TRUE),
      stringsAsFactors = FALSE
    )
  }))
  if (verbose) message("\nWeight Distribution:")
  if (verbose) print(weight_summary)
  
  # Calculate C-statistics
  c_stat_result <- calc_c_statistic(final_data[[exposure_var]], final_data[[ps_var]],
                                    label = "PS Model (unweighted)",
                                    verbose = verbose)
  c_stat <- if (!is.null(c_stat_result)) c_stat_result$c_stat else NA_real_
  c_stat_se <- if (!is.null(c_stat_result)) c_stat_result$se else NA_real_
  
  c_stat_wt_result <- calc_c_statistic(final_data[[exposure_var]], final_data[[ps_var]],
                                       weights_vec = final_data[[weight_var]],
                                       label = "Post-weighting (FS)",
                                       verbose = verbose)
  c_stat_weighted <- if (!is.null(c_stat_wt_result)) c_stat_wt_result$c_stat else NA_real_
  c_stat_weighted_se <- if (!is.null(c_stat_wt_result)) c_stat_wt_result$se else NA_real_
  
  #### STEP 4: Save trimmed+weighted data outputs ####################
  if (verbose) message("\nStep 4: Saving Outputs")
  if (verbose) message("----------------------")
  
  if (!is.null(out_csvpath)) {
    out_dir <- dirname(out_csvpath)
    if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
    utils::write.csv(final_data, file = out_csvpath, row.names = FALSE)
    if (verbose) message(sprintf("Data saved to CSV: %s", out_csvpath))
  }
  
  #### STEP 5: Diagnostics - Balance Table (optional) ####################
  combined_table <- NULL
  n_balanced <- n_total_vars <- n_improved <- NA
  
  if (make_unwt_wt_table1) {
    if (verbose) message("\nStep 5: Creating Combined Balance Table (Crude vs Weighted)")
    if (verbose) message("------------------------------------------------------------")
    
    # Auto-detect variable types if none provided
    if (is.null(table1_cont_vars) && is.null(table1_binary_vars) && is.null(table1_cat_vars)) {
      if (verbose) message("  Auto-detecting variable types...")
      exclude_vars <- c(exposure_var, ps_var, weight_var, "strata", "id", "ID", .row_id_col)
      candidate_vars <- setdiff(names(in_crude), exclude_vars)
      candidate_vars <- intersect(candidate_vars, names(final_data))
      
      var_types          <- get_var_types(in_crude[candidate_vars], max_cat_levels = 12)
      table1_cont_vars   <- var_types$cont_vars
      table1_binary_vars <- var_types$binary_vars
      table1_cat_vars    <- var_types$cat_vars
      
      if (verbose) message(sprintf("  Detected: %d continuous, %d binary, %d categorical",
                                   length(table1_cont_vars), length(table1_binary_vars), length(table1_cat_vars)))
    }
    
    if (verbose) message("  Calculating unweighted statistics...")
    unweighted_table1 <- build_table1(
      in_df             = in_crude,
      exposure_var      = exposure_var,
      exp_value         = exp_value,
      ref_value         = ref_value,
      cont_vars         = table1_cont_vars,
      binary_vars       = table1_binary_vars,
      cat_vars          = table1_cat_vars,
      use_weights       = FALSE,
      add_n_of_patients_row = TRUE,
      verbose           = verbose
    )
    
    if (verbose) message("  Calculating weighted statistics...")
    weighted_table1 <- build_table1(
      in_df             = final_data,
      exposure_var      = exposure_var,
      exp_value         = exp_value,
      ref_value         = ref_value,
      cont_vars         = table1_cont_vars,
      binary_vars       = table1_binary_vars,
      cat_vars          = table1_cat_vars,
      weight_var        = weight_var,
      use_weights       = TRUE,
      add_n_of_patients_row = TRUE,
      verbose           = verbose
    )
    
    # Combine: crude columns + weighted columns side by side
    # Standard column names from build_table1: Variable, Type, Total, Exp, Ref, Std_diff, ...
    combined_table <- data.frame(
      Variable         = unweighted_table1$Variable,
      Type             = unweighted_table1$Type,
      Total_Crude      = unweighted_table1$Total,
      Exp_Crude        = unweighted_table1$Exp,
      Ref_Crude        = unweighted_table1$Ref,
      Std_diff_Crude   = unweighted_table1$Std_diff,
      stringsAsFactors = FALSE
    )
    
    # Match weighted rows to crude rows by Variable
    wt_idx <- match(combined_table$Variable, weighted_table1$Variable)
    combined_table$Total_Weighted    <- weighted_table1$Total[wt_idx]
    combined_table$Exp_Weighted      <- weighted_table1$Exp[wt_idx]
    combined_table$Ref_Weighted      <- weighted_table1$Ref[wt_idx]
    combined_table$Std_diff_Weighted <- weighted_table1$Std_diff[wt_idx]
    
    # Balance statistics
    valid_rows <- combined_table$Variable != "n_of_patients"
    sd_crude_num   <- suppressWarnings(as.numeric(combined_table$Std_diff_Crude[valid_rows]))
    sd_weighted_num <- suppressWarnings(as.numeric(combined_table$Std_diff_Weighted[valid_rows]))
    
    n_balanced   <- sum(abs(sd_weighted_num) < (std_diff_threshold * 100), na.rm = TRUE)
    n_total_vars <- sum(!is.na(sd_weighted_num))
    n_improved   <- sum(abs(sd_weighted_num) < abs(sd_crude_num), na.rm = TRUE)
    
    if (verbose) message(sprintf("\n  Variables balanced (|Std Diff| < %.0f%%): %d/%d",
                                 std_diff_threshold * 100, n_balanced, n_total_vars))
    if (verbose) message(sprintf("  Variables with improved balance: %d/%d", n_improved, n_total_vars))
    
    results$combined_table <- combined_table
  }
  
  #### STEP 6: Diagnostics- PS Distribution & Plots (optional) ################
  if (!is.null(out_dir_plots)) {
    if (verbose) message("\nStep 6: Creating Diagnostic Plots")
    if (verbose) message("----------------------------------")
    
    if (!requireNamespace("ggplot2", quietly = TRUE)) {
      stop(
        "PS Diagnostic Plots requires the 'ggplot2' package.\n",
        "Please install it using install.packages('ggplot2')."
      )
    }
    if (!dir.exists(out_dir_plots)) dir.create(out_dir_plots, recursive = TRUE)
    
    # Derive plot prefix from xlsx path or use generic
    plot_prefix <- if (!is.null(out_xlsxpath_report)) {
      px <- sub("\\.[^.]*$", "", basename(out_xlsxpath_report))
      gsub("_fs_diagnostic$|_diagnostic$", "", px)
    } else {
      "ps_fs"
    }
    
    # Colour palette
    fill_vals  <- c("0" = "#e41a1c", "1" = "#377eb8")
    fill_labs  <- c("0" = "Reference", "1" = "Exposure")
    
    # Unweighted data 
    if (verbose) message("  Creating unweighted PS distribution plots...")
    unwt_data <- in_crude
    unwt_data$exposure_factor <- factor(unwt_data[[exposure_var]], levels = c(0, 1), labels = c("0", "1"))
    
    # 6a. Unweighted density
    p <- ggplot2::ggplot(unwt_data, ggplot2::aes(x = .data[[ps_var]], fill = exposure_factor, color = exposure_factor)) +
      ggplot2::geom_density(alpha = 0.3, adjust = 1.5) +
      ggplot2::scale_fill_manual(values = fill_vals, labels = fill_labs) +
      ggplot2::scale_color_manual(values = fill_vals, labels = fill_labs) +
      ggplot2::labs(title = "Unweighted PS Distribution - Density", x = "PS", y = "Density") +
      ggplot2::theme_minimal() + ggplot2::theme(plot.title = ggplot2::element_text(size = 14, face = "bold", hjust = 0.5),
                                                legend.position = "top") + ggplot2::xlim(0, 1)
    ggplot2::ggsave(file.path(out_dir_plots, paste0(plot_prefix, "_ps_distr_density_unwt.png")),
                    plot = p, width = 8, height = 6, dpi = 150)
    
    # 6b. Unweighted histogram
    p <- ggplot2::ggplot(unwt_data, ggplot2::aes(x = .data[[ps_var]],
                                                 y = ggplot2::after_stat(count / tapply(count, group, sum)[group] * 100),
                                                 fill = exposure_factor)) +
      ggplot2::geom_histogram(alpha = 0.5, position = "identity", bins = 100) +
      ggplot2::scale_fill_manual(values = fill_vals, labels = fill_labs) +
      ggplot2::labs(title = "Unweighted PS Distribution - Histogram", x = "PS", y = "% patients (within group)") +
      ggplot2::theme_minimal() + ggplot2::theme(plot.title = ggplot2::element_text(size = 14, face = "bold", hjust = 0.5),
                                                legend.position = "top") + ggplot2::xlim(0, 1)
    ggplot2::ggsave(file.path(out_dir_plots, paste0(plot_prefix, "_ps_distr_histog_unwt.png")),
                    plot = p, width = 8, height = 6, dpi = 150)
    
    # 6c. Unweighted both
    p <- ggplot2::ggplot(unwt_data, ggplot2::aes(x = .data[[ps_var]])) +
      ggplot2::geom_histogram(ggplot2::aes(y = ggplot2::after_stat(density), fill = exposure_factor),
                              alpha = 0.5, position = "identity", bins = 100) +
      ggplot2::geom_density(ggplot2::aes(y = ggplot2::after_stat(density), color = exposure_factor), linewidth = 1.2, adjust = 1.5) +
      ggplot2::scale_fill_manual(values = fill_vals, labels = fill_labs) +
      ggplot2::scale_color_manual(values = fill_vals, labels = fill_labs, guide = "none") +
      ggplot2::labs(title = "Unweighted PS Distribution", x = "PS", y = "Density") +
      ggplot2::theme_minimal() + ggplot2::theme(plot.title = ggplot2::element_text(size = 14, face = "bold", hjust = 0.5),
                                                legend.position = "top") + ggplot2::xlim(0, 1)
    ggplot2::ggsave(file.path(out_dir_plots, paste0(plot_prefix, "_ps_distr_unwt.png")),
                    plot = p, width = 8, height = 6, dpi = 150)
    
    # Weighted data 
    if (verbose) message("  Creating weighted PS distribution plots...")
    wt_data <- final_data
    wt_data$exposure_factor <- factor(wt_data[[exposure_var]], levels = c(0, 1), labels = c("0", "1"))
    
    # 6d. Weighted density
    p <- ggplot2::ggplot(wt_data, ggplot2::aes(x = .data[[ps_var]], weight = .data[[weight_var]],
                                               fill = exposure_factor, color = exposure_factor)) +
      ggplot2::geom_density(alpha = 0.3, adjust = 1.5) +
      ggplot2::scale_fill_manual(values = fill_vals, labels = fill_labs) +
      ggplot2::scale_color_manual(values = fill_vals, labels = fill_labs) +
      ggplot2::labs(title = "Weighted PS Distribution - Density", x = "PS", y = "Density") +
      ggplot2::theme_minimal() + ggplot2::theme(plot.title = ggplot2::element_text(size = 14, face = "bold", hjust = 0.5),
                                                legend.position = "top") + ggplot2::xlim(0, 1)
    ggplot2::ggsave(file.path(out_dir_plots, paste0(plot_prefix, "_ps_distr_density_wt.png")),
                    plot = p, width = 8, height = 6, dpi = 150)
    
    # 6e. Weighted histogram
    p <- ggplot2::ggplot(wt_data, ggplot2::aes(x = .data[[ps_var]], weight = .data[[weight_var]],
                                               fill = exposure_factor)) +
      ggplot2::geom_histogram(ggplot2::aes(y = ggplot2::after_stat(count / tapply(count, group, sum)[group] * 100)),
                              alpha = 0.5, position = "identity", bins = 100) +
      ggplot2::scale_fill_manual(values = fill_vals, labels = fill_labs) +
      ggplot2::labs(title = "Weighted PS Distribution - Histogram", x = "PS", y = "% patients (within group)") +
      ggplot2::theme_minimal() + ggplot2::theme(plot.title = ggplot2::element_text(size = 14, face = "bold", hjust = 0.5),
                                                legend.position = "top") + ggplot2::xlim(0, 1)
    ggplot2::ggsave(file.path(out_dir_plots, paste0(plot_prefix, "_ps_distr_histog_wt.png")),
                    plot = p, width = 8, height = 6, dpi = 150)
    
    # 6f. Weighted both
    p <- ggplot2::ggplot(wt_data, ggplot2::aes(x = .data[[ps_var]], weight = .data[[weight_var]])) +
      ggplot2::geom_histogram(ggplot2::aes(y = ggplot2::after_stat(density), fill = exposure_factor),
                              alpha = 0.5, position = "identity", bins = 100) +
      ggplot2::geom_density(ggplot2::aes(y = ggplot2::after_stat(density), color = exposure_factor), linewidth = 1.2, adjust = 1.5) +
      ggplot2::scale_fill_manual(values = fill_vals, labels = fill_labs) +
      ggplot2::scale_color_manual(values = fill_vals, labels = fill_labs, guide = "none") +
      ggplot2::labs(title = "Weighted PS Distribution", x = "PS", y = "Density") +
      ggplot2::theme_minimal() + ggplot2::theme(plot.title = ggplot2::element_text(size = 14, face = "bold", hjust = 0.5),
                                                legend.position = "top") + ggplot2::xlim(0, 1)
    ggplot2::ggsave(file.path(out_dir_plots, paste0(plot_prefix, "_ps_distr_wt.png")),
                    plot = p, width = 8, height = 6, dpi = 150)
    
    # 6g. Weight distribution box plot
    if (verbose) message("  Creating weight distribution box plot...")
    box_data <- final_data
    box_data$exposure_group <- factor(box_data[[exposure_var]], levels = c(0, 1),
                                      labels = c("Reference", "Exposure"))
    y_limit <- stats::quantile(final_data[[weight_var]], 0.99, na.rm = TRUE) * 1.1
    
    p <- ggplot2::ggplot(box_data, ggplot2::aes(x = exposure_group, y = .data[[weight_var]], fill = exposure_group)) +
      ggplot2::geom_boxplot(alpha = 0.7, outlier.alpha = 0.3) +
      ggplot2::geom_hline(yintercept = 1, linetype = "dashed", color = "gray50", linewidth = 1) +
      ggplot2::labs(title = "Distribution of Fine Stratification Weights", x = "Group", y = "Weight") +
      ggplot2::theme_minimal() + ggplot2::theme(plot.title = ggplot2::element_text(size = 14, face = "bold", hjust = 0.5),
                                                legend.position = "none") +
      ggplot2::scale_fill_manual(values = c("Exposure" = "#377eb8", "Reference" = "#e41a1c")) +
      ggplot2::coord_cartesian(ylim = c(0, y_limit))
    ggplot2::ggsave(file.path(out_dir_plots, paste0(plot_prefix, "_boxplot.png")),
                    plot = p, width = 8, height = 6, dpi = 150)
    
    # 6h. Balance plot (Love plot) - only if table1 was created
    if (make_unwt_wt_table1 && !is.null(combined_table)) {
      if (verbose) message("  Creating balance plot...")
      bal_data <- combined_table[combined_table$Variable != "n_of_patients", , drop = FALSE]
      bal_data$sd_crude <- suppressWarnings(as.numeric(bal_data$Std_diff_Crude))
      bal_data$sd_wt    <- suppressWarnings(as.numeric(bal_data$Std_diff_Weighted))
      bal_data <- bal_data[!is.na(bal_data$sd_crude) | !is.na(bal_data$sd_wt), , drop = FALSE]
      
      create_love_plot(
        variable_names = bal_data$Variable,
        crude_std_diff = bal_data$sd_crude,
        adjusted_std_diff = bal_data$sd_wt,
        crude_label = "Unweighted", adjusted_label = "Weighted",
        title = "Standardized Differences: Unweighted vs Weighted",
        output_path = file.path(out_dir_plots, paste0(plot_prefix, "_balance.png")),
        colors = c("Unweighted" = "darkorange1", "Weighted" = "dodgerblue"),
        shapes = c("Unweighted" = 17, "Weighted" = 16),
        use_absolute = TRUE
      )
    }
    
    if (verbose) message(sprintf("  Saved plots to %s", out_dir_plots))
  }
  
  #### STEP 7: Diagnostics - Excel Report (optional) ################
  if (!is.null(out_xlsxpath_report)) {
    if (verbose) message("\nStep 7: Creating Excel Diagnostic Report")
    if (verbose) message("-----------------------------------------")
    
    out_dir <- dirname(out_xlsxpath_report)
    if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
    
    wb <- openxlsx::createWorkbook()
    
    # README
    if (!is.null(readme_text)) {
      add_readme_sheet(wb, readme_text, verbose = verbose)
    }
    
    # Summary sheet 
    openxlsx::addWorksheet(wb, "Summary")
    
    # PS Overlap Region
    has_overlap_bounds <- trim_nonoverlap_region &&
      !is.null(results$trim_summary$bounds) &&
      !anyNA(results$trim_summary$bounds)
    ps_overlap_min <- if (has_overlap_bounds) as.character(results$trim_summary$bounds["lower"]) else "N/A"
    ps_overlap_max <- if (has_overlap_bounds) as.character(results$trim_summary$bounds["upper"]) else "N/A"
    
    # Transpose overall_summary (1-row wide) into Parameter-Value rows (excluding "strata" column)
    os_cols <- setdiff(names(overall_summary), "strata")
    # Clean display names for overall_summary parameters
    os_display_names <- c(
      n_total_exp              = "N Exposure",
      n_total_unexp            = "N Reference",
      total_att_weighted_unexp = "ATT Weighted N Reference",
      total_ate_weighted_exp   = "ATE Weighted N Exposure",
      total_ate_weighted_unexp = "ATE Weighted N Reference"
    )
    overall_params <- data.frame(
      Parameter = ifelse(os_cols %in% names(os_display_names),
                         os_display_names[os_cols], os_cols),
      Value     = as.character(unlist(overall_summary[1, os_cols])),
      stringsAsFactors = FALSE
    )
    
    summary_info <- data.frame(
      Parameter = c("Exposure", "Exposure Value", "Reference Value",
                    "Estimand", "Stratification Method", "Number of Strata",
                    "Original N", "N after trimming", "Final N",
                    "N strata with both groups",
                    "PS Model C-statistic (unweighted)",
                    "Post-weighting C-statistic (FS)",
                    "PS Overlap Region_Min",
                    "PS Overlap Region_Max",
                    overall_params$Parameter),
      Value = c(exposure_var, as.character(exp_value), as.character(ref_value),
                estimand, stratification_method, as.character(number_of_strata),
                as.character(original_n),
                as.character(original_n - results$trim_summary$n_trimmed),
                as.character(nrow(final_data)),
                as.character(nrow(strata_complete)),
                if (!is.na(c_stat)) sprintf("%.4f", c_stat) else "N/A",
                if (!is.na(c_stat_weighted)) sprintf("%.4f", c_stat_weighted) else "N/A",
                ps_overlap_min,
                ps_overlap_max,
                overall_params$Value),
      stringsAsFactors = FALSE
    )
    openxlsx::writeDataTable(wb, "Summary", summary_info, startRow = 1)
    
    # Sample Size sheet
    openxlsx::addWorksheet(wb, "Sample Size")
    n_crude_exp <- sum(in_crude[[exposure_var]] == 1, na.rm = TRUE)
    n_crude_ref <- sum(in_crude[[exposure_var]] == 0, na.rm = TRUE)
    n_trim_exp  <- sum(final_data[[exposure_var]] == 1, na.rm = TRUE)
    n_trim_ref  <- sum(final_data[[exposure_var]] == 0, na.rm = TRUE)
    sample_size_comparison <- data.frame(
      exposure    = c(0, 1),
      Crude       = c(n_crude_ref, n_crude_exp),
      Trimmed     = c(n_trim_ref, n_trim_exp),
      n_trimmed   = c(n_crude_ref - n_trim_ref, n_crude_exp - n_trim_exp),
      pct_trimmed = round(100 * c(n_crude_ref - n_trim_ref, n_crude_exp - n_trim_exp) /
                            c(n_crude_ref, n_crude_exp), 1),
      stringsAsFactors = FALSE
    )
    openxlsx::writeDataTable(wb, "Sample Size", sample_size_comparison, startRow = 1)
    
    # PS Distribution sheet
    openxlsx::addWorksheet(wb, "PS Distribution")
    ps_comp <- do.call(rbind, lapply(c("Crude", "Trimmed"), function(ds) {
      d <- if (ds == "Crude") in_crude else final_data
      do.call(rbind, lapply(c(0, 1), function(g) {
        ps <- d[[ps_var]][d[[exposure_var]] == g]
        data.frame(dataset = ds, exposure = g, n = length(ps),
                   mean_ps = mean(ps, na.rm = TRUE), sd_ps = stats::sd(ps, na.rm = TRUE),
                   min_ps = min(ps, na.rm = TRUE), max_ps = max(ps, na.rm = TRUE),
                   stringsAsFactors = FALSE)
      }))
    }))
    openxlsx::writeDataTable(wb, "PS Distribution", ps_comp, startRow = 1)
    
    # Strata Detail
    openxlsx::addWorksheet(wb, "Strata Detail")
    openxlsx::writeDataTable(wb, "Strata Detail", strata_complete)
    
    # Weight Distribution 
    openxlsx::addWorksheet(wb, "Weight Distribution")
    openxlsx::writeDataTable(wb, "Weight Distribution", weight_summary)
    
    # --- CHANGED: Removed "Interpretation" column from c_stat_df ---
    c_stat_df <- data.frame(
      Measure = c("PS Model C-statistic (unweighted)",
                  "Post-weighting C-statistic (FS)"),
      C_statistic = c(
        if (!is.na(c_stat)) sprintf("%.4f", c_stat) else "N/A",
        if (!is.na(c_stat_weighted)) sprintf("%.4f", c_stat_weighted) else "N/A"
      ),
      SE = c(
        if (!is.na(c_stat_se)) sprintf("%.4f", c_stat_se) else "N/A",
        if (!is.na(c_stat_weighted_se)) sprintf("%.4f", c_stat_weighted_se) else "N/A"
      ),
      stringsAsFactors = FALSE
    )
    openxlsx::addWorksheet(wb, "C_statistic")
    openxlsx::writeDataTable(wb, "C_statistic", c_stat_df)
    if (verbose) message("  Added C-statistic sheet")
    
    # Balance Table (if created) 
    if (make_unwt_wt_table1 && !is.null(combined_table)) {
      openxlsx::addWorksheet(wb, "Balance Table")
      openxlsx::writeData(wb, "Balance Table",
                          "COMBINED BALANCE TABLE: Crude vs Fine Stratification Weighted",
                          startRow = 1)
      openxlsx::writeData(wb, "Balance Table",
                          sprintf("Variables balanced (|Std Diff| < %.0f%%): %d/%d | Improved: %d/%d",
                                  std_diff_threshold * 100, n_balanced, n_total_vars, n_improved, n_total_vars),
                          startRow = 2)
      openxlsx::writeDataTable(wb, "Balance Table", combined_table, startRow = 4)
    }
    
    openxlsx::saveWorkbook(wb, out_xlsxpath_report, overwrite = TRUE)
    if (verbose) message(sprintf("Diagnostics saved to: %s", out_xlsxpath_report))
  }
  
  # Clean up internal columns
  final_data[[.row_id_col]] <- NULL
  
  # DONE
  if (verbose) message("\n========================================")
  if (verbose) message("PS FINE STRATIFICATION COMPLETE")
  if (verbose) message("========================================")
  
  return(invisible(final_data))
}
