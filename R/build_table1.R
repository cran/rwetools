# BUILD TABLE 1 ################################################
#' Build a Table 1 comparing baseline characteristics between groups
#'
#' Generates a descriptive Table 1 with counts, percentages, means, SDs,
#' crude differences, standardized mean differences (SMD), and missingness.
#' Supports inverse probability weighting via a user-supplied weight column
#' and optionally saves the output to an Excel file.
#'
#' @param in_df Data frame containing the analytic cohort.
#' @param out_xlsxpath Character string. File path for Excel output, or
#'   \code{NULL} to skip (default \code{NULL}).
#' @param exposure_var Character string. Name of the binary exposure column.
#' @param exp_value Value in \code{exposure_var} representing the exposed group
#'   (default 1).
#' @param ref_value Value in \code{exposure_var} representing the reference group
#'   (default 0).
#' @param use_weights Logical. Apply weights? (default FALSE).
#' @param weight_var Character string. Name of the weight column (default
#'   \code{"psweight"}).
#' @param cont_vars Character vector of continuous variable names.
#' @param cat_vars Character vector or named list of categorical variable names.
#'   If a named list, each element should be a character vector of factor levels.
#' @param binary_vars Character vector of binary variable names.
#' @param drop_vars Character vector of variable names to exclude.
#' @param drop_varpattern Character vector of regex patterns; matching variables
#'   are excluded.
#' @param Var_colname,Vartype_colname,Total_colname,Exp_colname,Ref_colname
#'   Column-name overrides for the output table.
#' @param CrudeDiff_colname,StdDiff_colname Column-name overrides for
#'   difference columns.
#' @param MissingTotal_colname,MissingExp_colname,MissingRef_colname
#'   Column-name overrides for missingness columns.
#' @param ExistingTotal_colname,ExistingExp_colname,ExistingRef_colname
#'   Column-name overrides for non-missing-count columns. Set to \code{NULL}
#'   to omit (default for Exp/Ref).
#' @param mean_decimal,sd_decimal Integer. Decimal places for mean and SD.
#' @param n_decimal Integer. Decimal places for counts (default 0).
#' @param pct_decimal Integer. Decimal places for percentages (default 1).
#' @param use_absolute_values_for_diff Logical. Show absolute differences?
#'   (default FALSE).
#' @param add_n_of_patients_row Logical. Prepend an \code{"n_of_patients"} row?
#'   (default TRUE).
#' @param verbose Logical. Print progress messages (default TRUE).
#'
#' @return Invisibly returns the Table 1 data frame.
#'
#' @section Side Effects:
#' When \code{out_xlsxpath} is not \code{NULL}, creates the output directory
#' (if needed) and writes an Excel workbook.
#'
#' @examples
#' csv_path <- system.file("extdata", "sample_data.csv", package = "rwetools")
#' df <- read.csv(csv_path)
#'
#' # Unweighted Table 1
#' tbl <- build_table1(
#'   in_df        = df,
#'   exposure_var = "exposure",
#'   cont_vars    = c("cont1", "cont2", "cont3"),
#'   binary_vars  = c("binary1", "binary2"),
#'   cat_vars     = c("cat1", "cat2"),
#'   verbose      = FALSE
#' )
#' head(tbl)
#'
#' \donttest{
#' # Weighted Table 1 with Excel output (requires openxlsx)
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
#'   tbl_wt <- build_table1(
#'     in_df        = df_wt,
#'     out_xlsxpath = out_xlsx,
#'     exposure_var = "exposure",
#'     use_weights  = TRUE,
#'     weight_var   = "mw_wt",
#'     cont_vars    = c("cont1", "cont2", "cont3"),
#'     binary_vars  = c("binary1", "binary2"),
#'     cat_vars     = c("cat1", "cat2"),
#'     verbose      = FALSE
#'   )
#' }
#' }
#'
#' @export
build_table1 <- function(in_df,
                         out_xlsxpath = NULL,
                         exposure_var,
                         exp_value = 1,
                         ref_value = 0,
                         use_weights = FALSE,
                         weight_var = "psweight",
                         cont_vars = NULL,
                         cat_vars = NULL,
                         binary_vars = NULL,
                         drop_vars = NULL,
                         drop_varpattern = NULL,
                         Var_colname = "Variable",
                         Vartype_colname = "Type",
                         Total_colname = "Total",
                         Exp_colname = "Exp",
                         Ref_colname = "Ref",
                         CrudeDiff_colname = "Crude_diff",
                         StdDiff_colname = "Std_diff",
                         MissingTotal_colname = "Missing_Total_N_Pct",
                         MissingExp_colname = "Missing_Exp_N_Pct",
                         MissingRef_colname = "Missing_Ref_N_Pct",
                         ExistingTotal_colname = "Existing_Total_N_denom",
                         ExistingExp_colname = NULL,
                         ExistingRef_colname = NULL,
                         mean_decimal = 2,
                         sd_decimal = 2,
                         n_decimal = 0,
                         pct_decimal = 1,
                         use_absolute_values_for_diff = FALSE,
                         add_n_of_patients_row = TRUE,
                         verbose = TRUE){

  if (verbose) message("\U0001f680 Starting Table 1 generation...")

  # ================================================================
  # Input validation
  # ================================================================
  if (is.null(in_df) || !is.data.frame(in_df)) {
    stop("in_df must be a non-NULL data frame")
  }
  if (is.null(Var_colname) || Var_colname == "") {
    stop("Var_colname cannot be NULL or empty - this column is required!")
  }
  if (is.null(Exp_colname) || Exp_colname == "") {
    stop("Exp_colname cannot be NULL or empty - this column is required!")
  }
  if (is.null(Ref_colname) || Ref_colname == "") {
    stop("Ref_colname cannot be NULL or empty - this column is required!")
  }

  # Check if output directory exists, create if it doesn't
  if (!is.null(out_xlsxpath)) {
    out_directory <- dirname(out_xlsxpath)
    if (!dir.exists(out_directory)) {
      dir.create(out_directory, recursive = TRUE)
      if (verbose) message(sprintf("  Created output directory: %s", out_directory))
    }
  }

  # Determine if we need to calculate existing (non-missing) counts
  calculate_existing <- !is.null(ExistingTotal_colname) ||
    !is.null(ExistingExp_colname) ||
    !is.null(ExistingRef_colname)

  if (calculate_existing) {
    if (verbose) message("  Will calculate existing (non-missing) patient counts")
  }

  # Make a copy to avoid modifying original data
  dat <- in_df

  # ================================================================
  # Exposure validation: warn and drop rows with NA in exposure_var
  # ================================================================
  if (!exposure_var %in% names(dat)) stop("exposure column not found.")

  n_exp_na <- sum(is.na(dat[[exposure_var]]))
  if (n_exp_na > 0) {
    warning(sprintf(
      "exposure_var '%s' contains %d missing value(s) (%.1f%% of %d rows). These rows will have NA group assignment and are excluded from group-level calculations.",
      exposure_var, n_exp_na, 100 * n_exp_na / nrow(dat), nrow(dat)
    ))
  }

  # Recode exposure to 0/1
  g <- dat[[exposure_var]]
  if (verbose) message(sprintf("  Recoding exposure: '%s' (ref) -> 0, '%s' (exp) -> 1", ref_value, exp_value))
  g_recoded <- ifelse(g == ref_value, 0L,
                      ifelse(g == exp_value, 1L, NA_integer_))

  # Check if all values are either ref_value or exp_value
  unique_vals <- unique(g[!is.na(g)])
  if (!all(unique_vals %in% c(ref_value, exp_value))) {
    stop(sprintf("Exposure variable contains values other than ref_value (%s) and exp_value (%s): %s",
                 ref_value, exp_value, paste(setdiff(unique_vals, c(ref_value, exp_value)), collapse = ", ")))
  }

  dat$.grp <- g_recoded

  n_exp <- sum(dat$.grp == 1, na.rm = TRUE)
  n_ref <- sum(dat$.grp == 0, na.rm = TRUE)
  n_missing <- sum(is.na(dat$.grp))
  if (verbose) message(sprintf("  Exposure groups: Ref (n=%d), Exp (n=%d), Missing (n=%d)", n_ref, n_exp, n_missing))

  # ================================================================
  # Drop variables by exact name
  # ================================================================
  if (!is.null(drop_vars) && length(drop_vars) > 0) {
    if (verbose) message(sprintf("  Dropping %d variables by exact name", length(drop_vars)))
    cont_vars   <- setdiff(cont_vars, drop_vars)
    binary_vars <- setdiff(binary_vars, drop_vars)

    if (is.list(cat_vars) && !is.null(names(cat_vars))) {
      cat_vars <- cat_vars[!names(cat_vars) %in% drop_vars]
    } else {
      cat_vars <- setdiff(cat_vars, drop_vars)
    }
  }

  # Drop variables by pattern matching
  if (!is.null(drop_varpattern) && length(drop_varpattern) > 0) {
    if (verbose) message(sprintf("  Dropping variables matching %d pattern(s)", length(drop_varpattern)))

    pattern_regex <- paste(drop_varpattern, collapse = "|")

    if (!is.null(cont_vars) && length(cont_vars) > 0) {
      matched_cont <- cont_vars[grepl(pattern_regex, cont_vars)]
      if (length(matched_cont) > 0 && verbose) {
        message(sprintf("    Dropping %d continuous variables: %s",
                        length(matched_cont), paste(matched_cont, collapse = ", ")))
      }
      cont_vars <- cont_vars[!grepl(pattern_regex, cont_vars)]
    }

    if (!is.null(binary_vars) && length(binary_vars) > 0) {
      matched_binary <- binary_vars[grepl(pattern_regex, binary_vars)]
      if (length(matched_binary) > 0 && verbose) {
        message(sprintf("    Dropping %d binary variables: %s",
                        length(matched_binary), paste(matched_binary, collapse = ", ")))
      }
      binary_vars <- binary_vars[!grepl(pattern_regex, binary_vars)]
    }

    if (!is.null(cat_vars) && length(cat_vars) > 0) {
      if (is.list(cat_vars) && !is.null(names(cat_vars))) {
        matched_cat <- names(cat_vars)[grepl(pattern_regex, names(cat_vars))]
        if (length(matched_cat) > 0 && verbose) {
          message(sprintf("    Dropping %d categorical variables: %s",
                          length(matched_cat), paste(matched_cat, collapse = ", ")))
        }
        cat_vars <- cat_vars[!grepl(pattern_regex, names(cat_vars))]
      } else {
        matched_cat <- cat_vars[grepl(pattern_regex, cat_vars)]
        if (length(matched_cat) > 0 && verbose) {
          message(sprintf("    Dropping %d categorical variables: %s",
                          length(matched_cat), paste(matched_cat, collapse = ", ")))
        }
        cat_vars <- cat_vars[!grepl(pattern_regex, cat_vars)]
      }
    }
  }

  # ================================================================
  # Prepare weights (with NA validation when use_weights = TRUE)
  # ================================================================
  if (verbose) message("Step 1: Preparing weights...")
  w <- rep(1, nrow(dat))
  using_w <- FALSE

  if (use_weights && !is.null(weight_var) && weight_var %in% names(dat)) {
    # Warn about NA in weight_var (but do NOT drop rows — substitute with 1, same as original)
    n_wt_na <- sum(is.na(dat[[weight_var]]))
    if (n_wt_na > 0) {
      warning(sprintf(
        "weight_var '%s' contains %d missing value(s) (%.1f%% of %d rows). These will be assigned weight = 1.",
        weight_var, n_wt_na, 100 * n_wt_na / nrow(dat), nrow(dat)
      ))
    }

    w_try <- suppressWarnings(as.numeric(dat[[weight_var]]))
    w <- ifelse(is.na(w_try) | w_try <= 0, 1, w_try)
    using_w <- any(w > 1)
    if (verbose) message(sprintf("  Weight summary: min=%.2f, max=%.2f, mean=%.2f",
                                 min(w, na.rm=TRUE), max(w, na.rm=TRUE), mean(w, na.rm=TRUE)))
  } else {
    if (verbose) message("  Using unweighted analysis (all weights = 1)")
  }
  dat$..w <- w

  # ================================================================
  # Pre-compute shared vectors and group indices (ONCE)
  # ================================================================
  w_vec   <- dat$..w
  grp_vec <- dat$.grp
  idx_ref <- which(grp_vec == 0)
  idx_tx  <- which(grp_vec == 1)

  # Table parameters
  if (verbose) message("Step 2: Creating survey design objects...")
  N_all <- nrow(dat)
  N_ref <- length(idx_ref)
  N_tx  <- length(idx_tx)

  # Create survey design objects ONCE for all continuous variables
  des_all <- survey::svydesign(ids = ~1, weights = ~..w, data = dat)
  des_ref <- subset(des_all, .grp == 0)
  des_tx  <- subset(des_all, .grp == 1)

  # ================================================================
  # Create n_of_patients row if requested
  # ================================================================
  n_patients_row <- NULL
  if (add_n_of_patients_row) {
    if (verbose) message("Step 3: Creating n_of_patients row...")

    if (use_weights && using_w) {
      sum_w_exp <- sum(w_vec[idx_tx], na.rm = TRUE)
      sum_w_ref <- sum(w_vec[idx_ref], na.rm = TRUE)
      sum_w_total <- sum_w_exp + sum_w_ref

      n_patients_row <- list(
        Variable = "n_of_patients",
        Type = "Other",
        Total = formatC(round(sum_w_total, n_decimal), format = "f", digits = n_decimal),
        Exp = formatC(round(sum_w_exp, n_decimal), format = "f", digits = n_decimal),
        Ref = formatC(round(sum_w_ref, n_decimal), format = "f", digits = n_decimal),
        Crude_diff = formatC(0, format = "f", digits = mean_decimal),
        Std_diff = formatC(0, format = "f", digits = 3),
        Missing_Total_N_Pct = fmt_n_pct(0, 1, n_decimal, pct_decimal),
        Missing_Exp_N_Pct = fmt_n_pct(0, 1, n_decimal, pct_decimal),
        Missing_Ref_N_Pct = fmt_n_pct(0, 1, n_decimal, pct_decimal)
      )
    } else {
      n_patients_row <- list(
        Variable = "n_of_patients",
        Type = "Other",
        Total = formatC(N_all, format = "f", digits = n_decimal),
        Exp = formatC(N_tx, format = "f", digits = n_decimal),
        Ref = formatC(N_ref, format = "f", digits = n_decimal),
        Crude_diff = formatC(0, format = "f", digits = mean_decimal),
        Std_diff = formatC(0, format = "f", digits = 3),
        Missing_Total_N_Pct = fmt_n_pct(0, 1, n_decimal, pct_decimal),
        Missing_Exp_N_Pct = fmt_n_pct(0, 1, n_decimal, pct_decimal),
        Missing_Ref_N_Pct = fmt_n_pct(0, 1, n_decimal, pct_decimal)
      )
    }

    if (calculate_existing) {
      n_patients_row$Existing_Total_N <- n_patients_row$Total
      n_patients_row$Existing_Exp_N <- n_patients_row$Exp
      n_patients_row$Existing_Ref_N <- n_patients_row$Ref
    }
  }

  # ================================================================
  # Process Continuous Variables (reusing pre-built survey designs)
  # ================================================================
  if (verbose) message("Step 4: Processing Continuous Variables...")

  cont_tbl <- if (!is.null(cont_vars) && length(cont_vars) > 0) {
    do.call(rbind, lapply(cont_vars, function(v) {
      if (!v %in% names(dat)) return(NULL)
      if (verbose) message(sprintf("  Processing variable: %s", v))

      x_vec <- dat[[v]]
      miss_all <- sum(is.na(x_vec))
      miss_ref <- sum(is.na(x_vec[idx_ref]))
      miss_tx  <- sum(is.na(x_vec[idx_tx]))

      existing_all <- existing_ref <- existing_tx <- NULL
      if (calculate_existing) {
        existing_all <- sum(!is.na(x_vec))
        existing_ref <- sum(!is.na(x_vec[idx_ref]))
        existing_tx  <- sum(!is.na(x_vec[idx_tx]))
      }

      # Reuse pre-built design objects (no svydesign/subset per variable)
      stats_all <- get_weighted_stats(des_all, v)
      stats_ref <- get_weighted_stats(des_ref, v)
      stats_tx  <- get_weighted_stats(des_tx, v)

      m_all <- unname(stats_all["mean"]); s_all <- unname(stats_all["sd"])
      m_ref <- unname(stats_ref["mean"]); s_ref <- unname(stats_ref["sd"])
      m_tx  <- unname(stats_tx["mean"]);  s_tx  <- unname(stats_tx["sd"])

      pooled_var <- (s_ref^2 + s_tx^2) / 2
      smd <- if (is.na(pooled_var) || pooled_var <= 0) NA_real_ else (m_tx - m_ref) / sqrt(pooled_var)

      result <- list(
        Variable = v,
        Type = "Continuous",
        Total = fmt_mean_sd(m_all, s_all, mean_decimal, sd_decimal),
        Exp = fmt_mean_sd(m_tx, s_tx, mean_decimal, sd_decimal),
        Ref = fmt_mean_sd(m_ref, s_ref, mean_decimal, sd_decimal),
        Crude_diff = fmt_num(m_tx - m_ref, digits = mean_decimal, use_abs = use_absolute_values_for_diff),
        Std_diff = fmt_num(smd, digits = 3, use_abs = use_absolute_values_for_diff),
        Missing_Total_N_Pct = fmt_n_pct(miss_all, N_all, n_decimal, pct_decimal),
        Missing_Exp_N_Pct = fmt_n_pct(miss_tx, N_tx, n_decimal, pct_decimal),
        Missing_Ref_N_Pct = fmt_n_pct(miss_ref, N_ref, n_decimal, pct_decimal)
      )

      if (calculate_existing) {
        result$Existing_Total_N <- formatC(existing_all, format = "f", digits = n_decimal)
        result$Existing_Exp_N <- formatC(existing_tx, format = "f", digits = n_decimal)
        result$Existing_Ref_N <- formatC(existing_ref, format = "f", digits = n_decimal)
      }

      result
    }))
  } else {
    NULL
  }

  # ================================================================
  # Process Binary Variables (using vectorized xtabs)
  # ================================================================
  if (verbose) message("Step 5: Processing Binary Variables...")
  binary_tbl <- if (!is.null(binary_vars) && length(binary_vars) > 0) {
    do.call(rbind, lapply(binary_vars, function(v) {
      process_catbin_direct(v, dat,
                            levels_to_show = "non_zero",
                            N_all = N_all, N_ref = N_ref, N_tx = N_tx,
                            idx_ref = idx_ref, idx_tx = idx_tx,
                            w_vec = w_vec, grp_vec = grp_vec,
                            using_w = using_w,
                            n_decimal = n_decimal,
                            pct_decimal = pct_decimal,
                            use_abs = use_absolute_values_for_diff,
                            calculate_existing = calculate_existing,
                            verbose = verbose)
    }))
  } else {
    NULL
  }

  # ================================================================
  # Process Categorical Variables (using vectorized xtabs)
  # ================================================================
  if (verbose) message("Step 6: Processing Categorical Variables...")
  cat_tbl <- if (!is.null(cat_vars) && length(cat_vars) > 0) {
    is_named_list <- is.list(cat_vars) &&
      !is.null(names(cat_vars)) &&
      all(names(cat_vars) != "") &&
      length(names(cat_vars)) == length(cat_vars)

    if (is_named_list) {
      if (verbose) message("  Using predefined levels from cat_vars list")
      do.call(rbind, lapply(names(cat_vars), function(var_name) {
        predefined_levels <- as.character(cat_vars[[var_name]])

        process_catbin_direct(
          var_name, dat,
          levels_to_show = "all",
          predefined_levels = predefined_levels,
          N_all = N_all, N_ref = N_ref, N_tx = N_tx,
          idx_ref = idx_ref, idx_tx = idx_tx,
          w_vec = w_vec, grp_vec = grp_vec,
          using_w = using_w,
          n_decimal = n_decimal,
          pct_decimal = pct_decimal,
          use_abs = use_absolute_values_for_diff,
          calculate_existing = calculate_existing,
          verbose = verbose
        )
      }))
    } else {
      if (verbose) message("  Using data-driven levels (cat_vars is character vector)")
      do.call(rbind, lapply(cat_vars, function(v) {
        process_catbin_direct(
          v, dat,
          levels_to_show = "all",
          predefined_levels = NULL,
          N_all = N_all, N_ref = N_ref, N_tx = N_tx,
          idx_ref = idx_ref, idx_tx = idx_tx,
          w_vec = w_vec, grp_vec = grp_vec,
          using_w = using_w,
          n_decimal = n_decimal,
          pct_decimal = pct_decimal,
          use_abs = use_absolute_values_for_diff,
          calculate_existing = calculate_existing,
          verbose = verbose
        )
      }))
    }
  } else {
    NULL
  }

  # ================================================================
  # Combine all results
  # ================================================================
  if (verbose) message("Step 7: Finalizing and combining results...")

  all_rows <- list()
  if (add_n_of_patients_row) {
    all_rows <- c(list(n_patients_row), all_rows)
  }
  if (!is.null(cont_tbl) && nrow(cont_tbl) > 0) {
    all_rows <- c(all_rows, lapply(seq_len(nrow(cont_tbl)), function(i) as.list(cont_tbl[i, ])))
  }
  if (!is.null(binary_tbl) && nrow(binary_tbl) > 0) {
    all_rows <- c(all_rows, lapply(seq_len(nrow(binary_tbl)), function(i) as.list(binary_tbl[i, ])))
  }
  if (!is.null(cat_tbl) && nrow(cat_tbl) > 0) {
    all_rows <- c(all_rows, lapply(seq_len(nrow(cat_tbl)), function(i) as.list(cat_tbl[i, ])))
  }

  # ================================================================
  # Build final table with only requested columns
  # ================================================================
  if (verbose) message("Step 8: Building final table with requested columns...")

  col_mapping <- list(
    Variable = Var_colname,
    Type = Vartype_colname,
    Total = Total_colname,
    Exp = Exp_colname,
    Ref = Ref_colname,
    Crude_diff = CrudeDiff_colname,
    Std_diff = StdDiff_colname,
    Missing_Total_N_Pct = MissingTotal_colname,
    Missing_Exp_N_Pct = MissingExp_colname,
    Missing_Ref_N_Pct = MissingRef_colname,
    Existing_Total_N = ExistingTotal_colname,
    Existing_Exp_N = ExistingExp_colname,
    Existing_Ref_N = ExistingRef_colname
  )

  active_cols <- col_mapping[!sapply(col_mapping, is.null)]

  col_list <- lapply(names(active_cols), function(internal_name) {
    col_data <- sapply(all_rows, function(row) row[[internal_name]])
    col_data
  })
  final_table <- as.data.frame(col_list, stringsAsFactors = FALSE)
  names(final_table) <- vapply(active_cols, identity, character(1))

  # ================================================================
  # Save to Excel if requested
  # ================================================================
  if (!is.null(out_xlsxpath)) {
    if (verbose) message("Step 10: Saving output to Excel file")
    wb <- openxlsx::createWorkbook()
    openxlsx::addWorksheet(wb, "Table1")
    header_style <- openxlsx::createStyle(textDecoration = "bold", halign = "center",
                                          valign = "center", border = "TopBottomLeftRight")
    openxlsx::writeData(wb, "Table1", final_table, headerStyle = header_style)
    openxlsx::setColWidths(wb, "Table1", cols = 1:ncol(final_table), widths = "auto")
    openxlsx::setColWidths(wb, "Table1", cols = 1, widths = 35)
    openxlsx::freezePane(wb, "Table1", firstActiveRow = 2)
    openxlsx::saveWorkbook(wb, out_xlsxpath, overwrite = TRUE)
    if (verbose) message("\u2705 Successfully saved Table 1 to: ", out_xlsxpath)
  }

  if (verbose) message("\u2705 Table 1 generation complete!")
  invisible(final_table)
}
