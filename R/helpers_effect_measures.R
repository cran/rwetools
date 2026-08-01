######### Shared input layer ##########################################
# HELPER: shared data-block validation for the effect functions ########
#' Validates the in_df_crude / in_df_weight / in_df_match data blocks and
#' their if_* companions (API v2 directive): weight XOR match; crude-less
#' weight/match calls allowed; orphan if_* arguments are ignored with a
#' message; required variables must exist in every provided block.
#'
#' @param in_df_crude,in_df_weight,in_df_match Data blocks (or NULL).
#' @param if_weight_weight_var Weight column name for in_df_weight, or NULL.
#' @param if_match_match_id Match-set id column name for in_df_match, or NULL.
#' @param required_vars Character vector of columns every provided block must
#'   contain (exposure/outcome/time and, when used, stratification etc.).
#' @param verbose Logical; emit the ignore-messages.
#' @return List: blocks (named list Crude/Weighted/Matched, only the provided
#'   ones), weight_var (validated or NULL), match_id (validated or NULL).
#' @noRd
validate_effect_blocks <- function(in_df_crude = NULL,
                                   in_df_weight = NULL,
                                   in_df_match = NULL,
                                   if_weight_weight_var = NULL,
                                   if_match_match_id = NULL,
                                   required_vars = character(0),
                                   verbose = TRUE) {

  if (is.null(in_df_crude) && is.null(in_df_weight) && is.null(in_df_match)) {
    stop("At least one of in_df_crude, in_df_weight, in_df_match must be provided.")
  }
  if (!is.null(in_df_weight) && !is.null(in_df_match)) {
    stop("Provide in_df_weight OR in_df_match, not both: a weighted and a ",
         "matched analysis are separate calls in this API.")
  }

  # --- weight block ---
  weight_var <- NULL
  if (!is.null(in_df_weight)) {
    if (is.null(if_weight_weight_var)) {
      stop("if_weight_weight_var must be specified when in_df_weight is provided.")
    }
    if (!if_weight_weight_var %in% names(in_df_weight)) {
      stop(paste0("Weight variable '", if_weight_weight_var,
                  "' not found in in_df_weight."))
    }
    weight_var <- if_weight_weight_var
  } else if (!is.null(if_weight_weight_var)) {
    if (verbose) message("if_weight_weight_var is ignored because in_df_weight was not provided.")
  }

  # --- match block ---
  match_id <- NULL
  if (!is.null(in_df_match)) {
    if (is.null(if_match_match_id)) {
      if (verbose) message("No if_match_match_id supplied for in_df_match: no pair ",
                           "clustering (robust SE without cluster; row-resample bootstrap).")
    } else {
      if (!if_match_match_id %in% names(in_df_match)) {
        stop(paste0("Match id variable '", if_match_match_id,
                    "' not found in in_df_match."))
      }
      match_id <- if_match_match_id
    }
  } else if (!is.null(if_match_match_id)) {
    if (verbose) message("if_match_match_id is ignored because in_df_match was not provided.")
  }

  # --- required variables per provided block ---
  blocks <- list(Crude = in_df_crude, Weighted = in_df_weight, Matched = in_df_match)
  blocks <- blocks[!vapply(blocks, is.null, logical(1))]
  for (bl in names(blocks)) {
    missing_vars <- required_vars[!required_vars %in% names(blocks[[bl]])]
    if (length(missing_vars) > 0L) {
      stop(paste0("Variables not found in in_df_",
                  c(Crude = "crude", Weighted = "weight", Matched = "match")[bl],
                  ": ", paste(missing_vars, collapse = ", ")))
    }
  }

  list(blocks = blocks, weight_var = weight_var, match_id = match_id)
}


######### HRIR Helper #################
# HELPER: Rename IR/IRD columns with dynamic IR base name ##########
#' Renames the generic intermediate column names to final names with
#' the dynamic ir_col_base and the analysis suffix.
#' @noRd
rename_ir_ird_cols <- function(df, ir_base, suffix) {
  rename_map <- c(
    "n_subjects"    = paste0("N_Subjects_", suffix),
    "n_events"      = paste0("N_Events_", suffix),
    "person_years"  = paste0("PersonYears_", suffix),
    "ir_value"      = paste0(ir_base, "_", suffix),
    "ir_lci"        = paste0("IR_LCI_", suffix),
    "ir_uci"        = paste0("IR_UCI_", suffix),
    "IR_SE"         = paste0("IR_SE_", suffix),
    "ird"           = paste0("IRD_per_Npy_", suffix),
    "ird_lci"       = paste0("IRD_LCI_", suffix),
    "ird_uci"       = paste0("IRD_UCI_", suffix),
    "IRD_SE"        = paste0("IRD_SE_", suffix)
  )
  for (old_name in names(rename_map)) {
    names(df)[names(df) == old_name] <- rename_map[old_name]
  }
  return(df)
}


# HELPER: IR/IRD calculation (Unweighted data, No stratification) ##########
#' @noRd
calc_crude_ir_ird <- function(df, exp_var, out_var, py_var, conf_level, z, suffix, multiplier, ir_base) {

  summary_by_group <- do.call(rbind, lapply(split(df, df[[exp_var]]), function(g) {
    data.frame(
      exp_val      = g[[exp_var]][1],
      n_subjects   = nrow(g),
      n_events     = sum(g[[out_var]], na.rm = TRUE),
      person_years = sum(g[[py_var]], na.rm = TRUE),
      stringsAsFactors = FALSE
    )
  }))
  names(summary_by_group)[1] <- exp_var
  rownames(summary_by_group) <- NULL
  summary_by_group$ir_value    <- multiplier * summary_by_group$n_events / summary_by_group$person_years
  summary_by_group$ir_lci      <- multiplier * stats::qgamma((1 - conf_level) / 2, summary_by_group$n_events) / summary_by_group$person_years
  summary_by_group$ir_uci      <- multiplier * stats::qgamma(1 - (1 - conf_level) / 2, summary_by_group$n_events + 1) / summary_by_group$person_years
  summary_by_group$IR_SE <- multiplier * sqrt(summary_by_group$n_events) / summary_by_group$person_years

  ref_row <- summary_by_group[summary_by_group[[exp_var]] == 0, ]
  exp_row <- summary_by_group[summary_by_group[[exp_var]] == 1, ]

  ird_value <- exp_row$ir_value - ref_row$ir_value
  se_ird <- multiplier * sqrt((ref_row$n_events / ref_row$person_years^2) +
                                (exp_row$n_events / exp_row$person_years^2))
  ird_lci_value <- ird_value - z * se_ird
  ird_uci_value <- ird_value + z * se_ird

  summary_by_group$ird          <- ird_value
  summary_by_group$ird_lci      <- ird_lci_value
  summary_by_group$ird_uci      <- ird_uci_value
  summary_by_group$IRD_SE <- se_ird

  # Total row
  total_row <- data.frame(
    exp_val      = 99,
    n_subjects   = nrow(df),
    n_events     = sum(df[[out_var]], na.rm = TRUE),
    person_years = sum(df[[py_var]], na.rm = TRUE),
    stringsAsFactors = FALSE
  )
  names(total_row)[1] <- exp_var
  total_row$ir_value    <- multiplier * total_row$n_events / total_row$person_years
  total_row$ir_lci      <- multiplier * stats::qgamma((1 - conf_level) / 2, total_row$n_events) / total_row$person_years
  total_row$ir_uci      <- multiplier * stats::qgamma(1 - (1 - conf_level) / 2, total_row$n_events + 1) / total_row$person_years
  total_row$IR_SE <- multiplier * sqrt(total_row$n_events) / total_row$person_years
  total_row$ird          <- ird_value
  total_row$ird_lci      <- ird_lci_value
  total_row$ird_uci      <- ird_uci_value
  total_row$IRD_SE <- se_ird

  result <- rbind(summary_by_group, total_row)
  result <- rename_ir_ird_cols(result, ir_base, suffix)

  return(result)
}

# HELPER: IR/IRD calculation (Weighted data, No Stratification) #####
#' @noRd
calc_crude_ir_ird_weighted <- function(df, exp_var, out_var, py_var, wt_var, conf_level, z, suffix, multiplier, ir_base) {

  summary_by_group <- do.call(rbind, lapply(split(df, df[[exp_var]]), function(g) {
    data.frame(
      exp_val      = g[[exp_var]][1],
      n_subjects   = sum(g[[wt_var]], na.rm = TRUE),
      n_events     = sum(g[[out_var]] * g[[wt_var]], na.rm = TRUE),
      person_years = sum(g[[py_var]] * g[[wt_var]], na.rm = TRUE),
      stringsAsFactors = FALSE
    )
  }))
  names(summary_by_group)[1] <- exp_var
  rownames(summary_by_group) <- NULL
  summary_by_group$ir_value <- multiplier * summary_by_group$n_events / summary_by_group$person_years

  # Robust (design-based sandwich) arm SE/CI via the joint cell engine
  # (matrix r4; P1): the previous rounded-qgamma path treated the weighted
  # event sum as a true Poisson count and was anti-conservative under IPTW.
  cell_fit <- fit_weighted_rate_cells(df, exp_var, out_var, py_var, wt_var)
  if (!is.null(cell_fit)) {
    ctr <- wtd_rate_contrasts(cell_fit, z, multiplier)
    ord <- match(as.character(summary_by_group[[exp_var]]), ctr$ir$exp_val)
    summary_by_group$ir_lci <- ctr$ir$lci[ord]
    summary_by_group$ir_uci <- ctr$ir$uci[ord]
    summary_by_group$IR_SE  <- ctr$ir$se[ord]
  } else {
    warning("Weighted IR robust variance fit failed; returning NA CIs. ",
            "Consider the bootstrap.")
    summary_by_group$ir_lci <- NA_real_
    summary_by_group$ir_uci <- NA_real_
    summary_by_group$IR_SE  <- NA_real_
  }

  ref_row <- summary_by_group[summary_by_group[[exp_var]] == 0, ]
  exp_row <- summary_by_group[summary_by_group[[exp_var]] == 1, ]

  # IRD point keeps the v0.3.0 arithmetic; SE is the joint sandwich/delta
  # (matrix r8).
  ird_value <- exp_row$ir_value - ref_row$ir_value
  se_ird <- if (!is.null(cell_fit)) ctr$ird$se else NA_real_
  ird_lci_value <- ird_value - z * se_ird
  ird_uci_value <- ird_value + z * se_ird

  summary_by_group$ird          <- ird_value
  summary_by_group$ird_lci      <- ird_lci_value
  summary_by_group$ird_uci      <- ird_uci_value
  summary_by_group$IRD_SE <- se_ird

  # Total row: crude weighted rate with a robust (intercept-only cell) CI
  total_ne <- sum(df[[out_var]] * df[[wt_var]], na.rm = TRUE)
  total_py <- sum(df[[py_var]] * df[[wt_var]], na.rm = TRUE)
  total_row <- data.frame(
    exp_val      = 99,
    n_subjects   = sum(df[[wt_var]], na.rm = TRUE),
    n_events     = total_ne,
    person_years = total_py,
    stringsAsFactors = FALSE
  )
  names(total_row)[1] <- exp_var
  total_row$ir_value <- multiplier * total_ne / total_py

  df_tot <- df
  df_tot$.const <- 0L
  tot_fit <- fit_weighted_rate_cells(df_tot, ".const", out_var, py_var, wt_var)
  if (!is.null(tot_fit)) {
    tot_ctr <- wtd_rate_contrasts_total(tot_fit, z, multiplier)
    total_row$ir_lci <- tot_ctr$lci
    total_row$ir_uci <- tot_ctr$uci
    total_row$IR_SE  <- tot_ctr$se
  } else {
    total_row$ir_lci <- NA_real_
    total_row$ir_uci <- NA_real_
    total_row$IR_SE  <- NA_real_
  }
  total_row$ird          <- ird_value
  total_row$ird_lci      <- ird_lci_value
  total_row$ird_uci      <- ird_uci_value
  total_row$IRD_SE <- se_ird

  result <- rbind(summary_by_group, total_row)
  result <- rename_ir_ird_cols(result, ir_base, suffix)

  return(result)
}



# HELPER: IR/IRD calculation (Unweighted data, Direct Standardization) ###############
#' @noRd
calc_standardized_ir_ird <- function(df, exp_var, out_var, py_var, strat_var, conf_level, z, suffix, multiplier, ir_base) {

  # Stratum-level statistics
  strat_stats <- do.call(rbind, lapply(
    split(df, list(df[[strat_var]], df[[exp_var]]), drop = TRUE),
    function(g) {
      data.frame(
        strat_val    = g[[strat_var]][1],
        exp_val      = g[[exp_var]][1],
        n_subjects   = nrow(g),
        n_events     = sum(g[[out_var]], na.rm = TRUE),
        person_years = sum(g[[py_var]], na.rm = TRUE),
        stringsAsFactors = FALSE
      )
    }
  ))
  names(strat_stats)[1:2] <- c(strat_var, exp_var)
  rownames(strat_stats) <- NULL
  strat_stats$ir_value <- multiplier * strat_stats$n_events / strat_stats$person_years

  # Total person-time per stratum (across both exposure groups)
  strat_total_pt <- stats::aggregate(
    df[[py_var]], by = list(df[[strat_var]]),
    FUN = function(x) sum(x, na.rm = TRUE)
  )
  names(strat_total_pt) <- c(strat_var, "pt_total_stratum")
  total_pt_all <- sum(strat_total_pt$pt_total_stratum)
  strat_total_pt$std_weight <- strat_total_pt$pt_total_stratum / total_pt_all

  # Merge weights into stratum stats
  strat_stats <- merge(strat_stats, strat_total_pt, by = strat_var, all.x = TRUE)

  # Standardized IR per exposure group: point/SE arithmetic unchanged; the
  # arm CI is Fay-Feuer gamma (matrix r3; collapses to Garwood qgamma for a
  # single stratum) instead of the previous log-normal interval [P2].
  std_ir_by_group <- do.call(rbind, lapply(split(strat_stats, strat_stats[[exp_var]]), function(g) {
    ff <- fay_feuer_ci(events = g$n_events, py = g$person_years,
                       std_w = g$std_weight, conf_level = conf_level)
    data.frame(
      exp_val            = g[[exp_var]][1],
      ir_value           = sum(g$ir_value * g$std_weight),
      var_ir_std         = sum(g$std_weight^2 * g$n_events / g$person_years^2) * multiplier^2,
      ir_lci             = multiplier * ff$lci,
      ir_uci             = multiplier * ff$uci,
      n_subjects         = sum(g$n_subjects),
      n_events           = sum(g$n_events),
      person_years       = sum(g$person_years),
      stringsAsFactors   = FALSE
    )
  }))
  names(std_ir_by_group)[1] <- exp_var
  rownames(std_ir_by_group) <- NULL
  std_ir_by_group$IR_SE <- sqrt(std_ir_by_group$var_ir_std)
  std_ir_by_group$var_ir_std  <- NULL

  # Standardized IRD
  ref_ir <- std_ir_by_group$ir_value[std_ir_by_group[[exp_var]] == 0]
  exp_ir <- std_ir_by_group$ir_value[std_ir_by_group[[exp_var]] == 1]
  ref_se <- std_ir_by_group$IR_SE[std_ir_by_group[[exp_var]] == 0]
  exp_se <- std_ir_by_group$IR_SE[std_ir_by_group[[exp_var]] == 1]

  ird_value <- exp_ir - ref_ir
  se_ird <- sqrt(ref_se^2 + exp_se^2)
  ird_lci_value <- ird_value - z * se_ird
  ird_uci_value <- ird_value + z * se_ird

  # Build summary table
  summary_by_group <- std_ir_by_group
  summary_by_group$ird          <- ird_value
  summary_by_group$ird_lci      <- ird_lci_value
  summary_by_group$ird_uci      <- ird_uci_value
  summary_by_group$IRD_SE <- se_ird

  # Total row (crude, not standardized)
  total_row <- data.frame(
    exp_val      = 99,
    n_subjects   = nrow(df),
    n_events     = sum(df[[out_var]], na.rm = TRUE),
    person_years = sum(df[[py_var]], na.rm = TRUE),
    stringsAsFactors = FALSE
  )
  names(total_row)[1] <- exp_var
  total_row$ir_value    <- multiplier * total_row$n_events / total_row$person_years
  total_row$ir_lci      <- multiplier * stats::qgamma((1 - conf_level) / 2, total_row$n_events) / total_row$person_years
  total_row$ir_uci      <- multiplier * stats::qgamma(1 - (1 - conf_level) / 2, total_row$n_events + 1) / total_row$person_years
  total_row$IR_SE <- multiplier * sqrt(total_row$n_events) / total_row$person_years
  total_row$ird          <- ird_value
  total_row$ird_lci      <- ird_lci_value
  total_row$ird_uci      <- ird_uci_value
  total_row$IRD_SE <- se_ird

  result <- rbind(summary_by_group, total_row)
  result <- rename_ir_ird_cols(result, ir_base, suffix)

  # Build stratum detail table
  strat_stats$exposure_group <- ifelse(
    strat_stats[[exp_var]] == 1, "Exposed",
    ifelse(strat_stats[[exp_var]] == 0, "Reference", as.character(strat_stats[[exp_var]]))
  )
  names(strat_stats)[names(strat_stats) == "ir_value"] <- ir_base
  stratum_detail <- strat_stats[, c(strat_var, "exposure_group", "n_subjects", "n_events", "person_years",
                                    ir_base, "pt_total_stratum", "std_weight"), drop = FALSE]

  return(list(summary = result, stratum_detail = stratum_detail))
}

# HELPER: IR/IRD calculation (Weighted data, Direct Standardization) ####################
#' @noRd
calc_standardized_ir_ird_weighted <- function(df, exp_var, out_var, py_var, wt_var, strat_var, conf_level, z, suffix, multiplier, ir_base) {

  # Stratum-level weighted statistics
  strat_stats <- do.call(rbind, lapply(
    split(df, list(df[[strat_var]], df[[exp_var]]), drop = TRUE),
    function(g) {
      data.frame(
        strat_val    = g[[strat_var]][1],
        exp_val      = g[[exp_var]][1],
        n_subjects   = sum(g[[wt_var]], na.rm = TRUE),
        n_events     = sum(g[[out_var]] * g[[wt_var]], na.rm = TRUE),
        person_years = sum(g[[py_var]] * g[[wt_var]], na.rm = TRUE),
        stringsAsFactors = FALSE
      )
    }
  ))
  names(strat_stats)[1:2] <- c(strat_var, exp_var)
  rownames(strat_stats) <- NULL
  strat_stats$ir_value <- multiplier * strat_stats$n_events / strat_stats$person_years

  # Total weighted person-time per stratum
  strat_total_pt <- stats::aggregate(
    df[[py_var]] * df[[wt_var]], by = list(df[[strat_var]]),
    FUN = function(x) sum(x, na.rm = TRUE)
  )
  names(strat_total_pt) <- c(strat_var, "pt_total_stratum")
  total_pt_all <- sum(strat_total_pt$pt_total_stratum)
  strat_total_pt$std_weight <- strat_total_pt$pt_total_stratum / total_pt_all

  strat_stats <- merge(strat_stats, strat_total_pt, by = strat_var, all.x = TRUE)

  # Joint sandwich/delta over arm x stratum cells (matrix r5/r9): the whole
  # weighted standardized analysis shares one robust vcov [P1/P2].
  std_w <- strat_total_pt$std_weight
  names(std_w) <- as.character(strat_total_pt[[strat_var]])
  cell_fit <- fit_weighted_rate_cells(df, exp_var, out_var, py_var, wt_var,
                                      strat_var = strat_var)
  ctr <- if (!is.null(cell_fit)) {
    wtd_std_rate_contrasts(cell_fit, std_w = std_w, z = z, multiplier = multiplier)
  } else {
    warning("Weighted standardized IR robust variance fit failed; returning ",
            "NA CIs. Consider the bootstrap.")
    NULL
  }

  # Standardized IR per exposure group: point arithmetic unchanged; SE/CI
  # from the joint engine (robust log-scale CI).
  std_ir_by_group <- do.call(rbind, lapply(split(strat_stats, strat_stats[[exp_var]]), function(g) {
    a <- as.character(g[[exp_var]][1])
    if (!is.null(ctr)) {
      sel <- ctr$ir$exp_val == a
      ir_se  <- ctr$ir$se[sel]
      ir_lci <- ctr$ir$lci[sel]
      ir_uci <- ctr$ir$uci[sel]
    } else {
      ir_se <- NA_real_; ir_lci <- NA_real_; ir_uci <- NA_real_
    }
    data.frame(
      exp_val      = g[[exp_var]][1],
      ir_value     = sum(g$ir_value * g$std_weight),
      IR_SE        = ir_se,
      ir_lci       = ir_lci,
      ir_uci       = ir_uci,
      n_subjects   = sum(g$n_subjects),
      n_events     = sum(g$n_events),
      person_years = sum(g$person_years),
      stringsAsFactors = FALSE
    )
  }))
  names(std_ir_by_group)[1] <- exp_var
  rownames(std_ir_by_group) <- NULL

  #  Standardized IRD: point unchanged; SE from the joint gradient.
  ref_ir <- std_ir_by_group$ir_value[std_ir_by_group[[exp_var]] == 0]
  exp_ir <- std_ir_by_group$ir_value[std_ir_by_group[[exp_var]] == 1]

  ird_value <- exp_ir - ref_ir
  se_ird <- if (!is.null(ctr)) ctr$ird$se else NA_real_
  ird_lci_value <- ird_value - z * se_ird
  ird_uci_value <- ird_value + z * se_ird

  summary_by_group <- std_ir_by_group
  summary_by_group$ird          <- ird_value
  summary_by_group$ird_lci      <- ird_lci_value
  summary_by_group$ird_uci      <- ird_uci_value
  summary_by_group$IRD_SE <- se_ird

  # Total row (crude weighted): robust intercept-only cell (rounding removed)
  total_ne <- sum(df[[out_var]] * df[[wt_var]], na.rm = TRUE)
  total_py <- sum(df[[py_var]] * df[[wt_var]], na.rm = TRUE)
  total_row <- data.frame(
    exp_val      = 99,
    n_subjects   = sum(df[[wt_var]], na.rm = TRUE),
    n_events     = total_ne,
    person_years = total_py,
    stringsAsFactors = FALSE
  )
  names(total_row)[1] <- exp_var
  total_row$ir_value <- multiplier * total_ne / total_py

  df_tot <- df
  df_tot$.const <- 0L
  tot_fit <- fit_weighted_rate_cells(df_tot, ".const", out_var, py_var, wt_var)
  if (!is.null(tot_fit)) {
    tot_ctr <- wtd_rate_contrasts_total(tot_fit, z, multiplier)
    total_row$ir_lci <- tot_ctr$lci
    total_row$ir_uci <- tot_ctr$uci
    total_row$IR_SE  <- tot_ctr$se
  } else {
    total_row$ir_lci <- NA_real_
    total_row$ir_uci <- NA_real_
    total_row$IR_SE  <- NA_real_
  }
  total_row$ird          <- ird_value
  total_row$ird_lci      <- ird_lci_value
  total_row$ird_uci      <- ird_uci_value
  total_row$IRD_SE <- se_ird

  result <- rbind(summary_by_group, total_row)
  result <- rename_ir_ird_cols(result, ir_base, suffix)

  strat_stats$exposure_group <- ifelse(
    strat_stats[[exp_var]] == 1, "Exposed",
    ifelse(strat_stats[[exp_var]] == 0, "Reference", as.character(strat_stats[[exp_var]]))
  )
  names(strat_stats)[names(strat_stats) == "ir_value"] <- ir_base
  stratum_detail <- strat_stats[, c(strat_var, "exposure_group", "n_subjects", "n_events", "person_years",
                                    ir_base, "pt_total_stratum", "std_weight"), drop = FALSE]

  return(list(summary = result, stratum_detail = stratum_detail))
}


# RRRD Helper Functions ####################

#' Calculate cumulative incidence at specified timepoint (unstratified)
#' Supports Kaplan-Meier (KM) and Aalen-Johansen (AJ) estimators.
#' KM: risk = 1 - S(t) via standard KM survival.
#' AJ: cumulative incidence function accounting for competing risks via
#' multi-state \code{survival::survfit} (requires \code{survival} >= 3.1).
#' Returns data frame with: Exposure_Value, N_Individuals, N_Events, Timepoint,
#' Risk_SE, RiskperN, Risk, Risk_Var. No Survival column is exposed.
#' @noRd
calc_km_cumulative_incidence <- function(data, time_var, event_var, exp_var,
                                         weight_var = NULL, timepoint, per_n,
                                         risk_estimator = "KM",
                                         competing_event_var = NULL) {

  # Set up weights
  if (!is.null(weight_var)) {
    wts <- data[[weight_var]]
  } else {
    wts <- rep(1, nrow(data))
  }

  # Count outcome events per exposure group (same for KM and AJ)
  n_events_vec <- vapply(
    split(data[[event_var]], data[[exp_var]]),
    sum, numeric(1)
  )

  if (risk_estimator == "AJ") {
    # ----- Aalen-Johansen estimator (competing risks) -----
    # Construct multi-state status:
    #   outcome_var==1                             -> 1 (event of interest)
    #   outcome_var==0 & competing_event_var==1    -> 2 (competing event)
    #   outcome_var==0 & competing_event_var==0    -> 0 (censored)
    ms_status <- ifelse(
      data[[event_var]] == 1, 1L,
      ifelse(data[[competing_event_var]] == 1, 2L, 0L)
    )
    ms_factor <- factor(ms_status, levels = c(0L, 1L, 2L))

    surv_obj <- survival::Surv(data[[time_var]], ms_factor)
    fit <- survival::survfit(surv_obj ~ data[[exp_var]], weights = wts)
    fit_summary <- summary(fit, times = timepoint, extend = TRUE)

    strata_names <- names(fit$strata)

    # Identify column index for event-of-interest state ("1")
    event_col <- which(fit$states == "1")

    # Extract CIF and SE for event state across strata
    # pstate/std.err are matrices: rows = strata, cols = states
    risk_values  <- fit_summary$pstate[, event_col]
    risk_se_values <- fit_summary$std.err[, event_col]

    result <- data.frame(
      Exposure_Value = gsub(paste0("data\\[\\[exp_var\\]\\]="), "", strata_names),
      N_Individuals = fit$n,
      N_Events = n_events_vec,
      Timepoint = timepoint,
      Risk_SE = risk_se_values,
      stringsAsFactors = FALSE
    )

    result$RiskperN <- risk_values * per_n
    result$Risk       <- risk_values
    result$Risk_Var   <- risk_se_values^2

  } else {
    # ----- Kaplan-Meier estimator -----
    surv_obj <- survival::Surv(data[[time_var]], data[[event_var]])
    km_fit <- survival::survfit(surv_obj ~ data[[exp_var]], weights = wts)
    km_summary <- summary(km_fit, times = timepoint, extend = TRUE)

    strata_names <- names(km_fit$strata)

    result <- data.frame(
      Exposure_Value = gsub(paste0("data\\[\\[exp_var\\]\\]="), "", strata_names),
      N_Individuals = km_fit$n,
      N_Events = n_events_vec,
      Timepoint = timepoint,
      Risk_SE = km_summary$std.err,
      stringsAsFactors = FALSE
    )

    result$Risk       <- 1 - km_summary$surv
    result$RiskperN <- result$Risk * per_n
    result$Risk_Var   <- result$Risk_SE^2
  }

  return(result)
}

#' Calculate cumulative incidence within a single stratum for one exposure group
#' Supports Kaplan-Meier (KM) and Aalen-Johansen (AJ) estimators.
#' @noRd
calc_km_risk_single_group <- function(data, time_var, event_var,
                                      weight_var = NULL, timepoint,
                                      risk_estimator = "KM",
                                      competing_event_var = NULL) {

  n <- nrow(data)
  n_events <- sum(data[[event_var]])

  # Handle edge case: no observations or no events
  if (n == 0) {
    return(list(risk = NA_real_, risk_var = NA_real_, n = 0L, n_events = 0L,
                survival = NA_real_, survival_se = NA_real_))
  }

  # Set up weights
  if (!is.null(weight_var)) {
    wts <- data[[weight_var]]
  } else {
    wts <- rep(1, n)
  }

  if (risk_estimator == "AJ") {
    # ----- Aalen-Johansen estimator (competing risks) -----
    ms_status <- ifelse(
      data[[event_var]] == 1, 1L,
      ifelse(data[[competing_event_var]] == 1, 2L, 0L)
    )
    ms_factor <- factor(ms_status, levels = c(0L, 1L, 2L))

    surv_obj <- survival::Surv(data[[time_var]], ms_factor)
    fit <- survival::survfit(surv_obj ~ 1, weights = wts)
    fit_summary <- summary(fit, times = timepoint, extend = TRUE)

    event_col <- which(fit$states == "1")

    # pstate may be a matrix (1-row) or a named vector depending on context
    if (is.matrix(fit_summary$pstate)) {
      risk_val <- fit_summary$pstate[1, event_col]
      risk_se  <- fit_summary$std.err[1, event_col]
    } else {
      risk_val <- fit_summary$pstate[event_col]
      risk_se  <- fit_summary$std.err[event_col]
    }

    risk_var_val <- risk_se^2
    surv    <- 1 - risk_val
    surv_se <- risk_se

  } else {
    # ----- Kaplan-Meier estimator -----
    surv_obj <- survival::Surv(data[[time_var]], data[[event_var]])
    km_fit <- survival::survfit(surv_obj ~ 1, weights = wts)
    km_summary <- summary(km_fit, times = timepoint, extend = TRUE)

    surv    <- km_summary$surv
    surv_se <- km_summary$std.err

    risk_val     <- 1 - surv
    risk_var_val <- surv_se^2
  }

  return(list(
    risk = risk_val,
    risk_var = risk_var_val,
    n = n,
    n_events = n_events,
    survival = surv,
    survival_se = surv_se
  ))
}

#' Calculate analytical SE and CI for unstratified RR and RD
#'
#' Uses Greenwood-based variance from KM survival estimates.
#' RR CI via delta method on log scale: exp(log(RR) +/- z * SE(log(RR))).
#' RD CI via normal approximation: RD +/- z * SE(RD).
#'
#' @param cuminc_data Data frame from \code{calc_km_cumulative_incidence},
#'   must contain columns: Exposure_Value, Risk, Risk_Var.
#' @param exp_val Exposed group value (as used in exposure_value column).
#' @param ref_val Reference group value (as used in exposure_value column).
#' @param per_n Denominator for risk expression (e.g., 1000).
#' @param confidence_level Numeric in (0, 1) for CI width.
#'
#' @return Named list with elements:
#'   \describe{
#'     \item{lnRR_se}{SE of log(RR) (via delta method on log scale)}
#'     \item{RR_lci, RR_uci}{CI bounds for RR (log-scale back-transformed)}
#'     \item{RD_se_per_n}{SE of RD scaled by per_n}
#'     \item{RD_lci_per_n, RD_uci_per_n}{CI bounds for RD scaled by per_n}
#'   }
#' @noRd
calc_analytical_se_unstratified <- function(cuminc_data, exp_val, ref_val,
                                            per_n, confidence_level) {

  exp_char <- as.character(exp_val)
  ref_char <- as.character(ref_val)

  risk_exp <- cuminc_data$Risk[cuminc_data$Exposure_Value == exp_char]
  risk_ref <- cuminc_data$Risk[cuminc_data$Exposure_Value == ref_char]
  var_exp  <- cuminc_data$Risk_Var[cuminc_data$Exposure_Value == exp_char]
  var_ref  <- cuminc_data$Risk_Var[cuminc_data$Exposure_Value == ref_char]

  na_result <- list(
    lnRR_se = NA_real_, RR_lci = NA_real_, RR_uci = NA_real_,
    RD_se_per_n = NA_real_, RD_lci_per_n = NA_real_, RD_uci_per_n = NA_real_
  )

  if (length(risk_exp) == 0L || length(risk_ref) == 0L) {
    return(na_result)
  }

  z <- stats::qnorm(1 - (1 - confidence_level) / 2)

  # --- RR: delta method on log scale ---
  # Var(log(RR)) = Var(R_exp)/R_exp^2 + Var(R_ref)/R_ref^2
  if (risk_exp > 0 && risk_ref > 0) {
    var_log_rr <- var_exp / risk_exp^2 + var_ref / risk_ref^2
    se_log_rr  <- sqrt(var_log_rr)
    RR         <- risk_exp / risk_ref
    lnRR_se    <- se_log_rr
    log_rr     <- log(RR)
    RR_lci     <- exp(log_rr - z * se_log_rr)
    RR_uci     <- exp(log_rr + z * se_log_rr)
  } else {
    lnRR_se <- NA_real_
    RR_lci <- NA_real_
    RR_uci <- NA_real_
  }

  # --- RD: normal approximation ---
  # Var(RD_per_n) = (Var(R_exp) + Var(R_ref)) * per_n^2
  RD_per_n     <- (risk_exp - risk_ref) * per_n
  var_rd       <- (var_exp + var_ref) * per_n^2
  RD_se_per_n  <- sqrt(var_rd)
  RD_lci_per_n <- RD_per_n - z * RD_se_per_n
  RD_uci_per_n <- RD_per_n + z * RD_se_per_n

  return(list(
    lnRR_se      = lnRR_se,
    RR_lci       = RR_lci,
    RR_uci       = RR_uci,
    RD_se_per_n  = RD_se_per_n,
    RD_lci_per_n = RD_lci_per_n,
    RD_uci_per_n = RD_uci_per_n
  ))
}

#' Calculate standardized risk, RR, and RD via direct standardization
#' @noRd
calc_standardized_rr_rd <- function(data, time_var, event_var, exp_var, strat_var,
                                    weight_var = NULL, timepoint, exp_val, ref_val,
                                    per_n, strat_weights = NULL,
                                    pre_split = NULL,
                                    risk_estimator = "KM",
                                    competing_event_var = NULL) {

  # Get strata
  strata_levels <- sort(unique(data[[strat_var]]))

  # Calculate stratum weights from data if not provided (original population proportions)
  if (is.null(strat_weights)) {
    N_total <- nrow(data)
    strat_weights <- vapply(strata_levels, function(k) {
      sum(data[[strat_var]] == k) / N_total
    }, numeric(1))
    names(strat_weights) <- as.character(strata_levels)
  }

  # Use pre-split data if provided, otherwise subset on the fly
  if (!is.null(pre_split)) {
    data_by_stratum <- pre_split
  } else {
    data_by_stratum <- split(data, data[[strat_var]])
  }

  # Accumulate stratum details in a list (avoid iterative rbind)
  strat_rows <- vector("list", length(strata_levels))

  for (j in seq_along(strata_levels)) {
    k <- strata_levels[j]
    k_char <- as.character(k)
    w_k <- strat_weights[k_char]

    # Subset to this stratum
    data_k <- data_by_stratum[[k_char]]
    data_k_exp <- data_k[data_k[[exp_var]] == exp_val, , drop = FALSE]
    data_k_ref <- data_k[data_k[[exp_var]] == ref_val, , drop = FALSE]

    # KM risk for exposed in stratum k
    km_exp <- calc_km_risk_single_group(
      data = data_k_exp,
      time_var = time_var,
      event_var = event_var,
      weight_var = weight_var,
      timepoint = timepoint,
      risk_estimator = risk_estimator,
      competing_event_var = competing_event_var
    )

    # KM risk for reference in stratum k
    km_ref <- calc_km_risk_single_group(
      data = data_k_ref,
      time_var = time_var,
      event_var = event_var,
      weight_var = weight_var,
      timepoint = timepoint,
      risk_estimator = risk_estimator,
      competing_event_var = competing_event_var
    )

    strat_rows[[j]] <- data.frame(
      stratum = k_char,
      w_k = w_k,
      n_total = nrow(data_k),
      n_exp = nrow(data_k_exp),
      n_ref = nrow(data_k_ref),
      events_exp = km_exp$n_events,
      events_ref = km_ref$n_events,
      risk_exp = km_exp$risk,
      risk_ref = km_ref$risk,
      risk_var_exp = km_exp$risk_var,
      risk_var_ref = km_ref$risk_var,
      stringsAsFactors = FALSE
    )
  }

  strat_details <- do.call(rbind, strat_rows)

  # Standardized risks: Risk_std = Sum(w_k * R_k)
  risk_std_exp <- sum(strat_details$w_k * strat_details$risk_exp, na.rm = TRUE)
  risk_std_ref <- sum(strat_details$w_k * strat_details$risk_ref, na.rm = TRUE)

  # RR and RD from standardized risks
  RR <- risk_std_exp / risk_std_ref
  RD <- (risk_std_exp - risk_std_ref) * per_n

  # Analytical SE via Greenwood: Var(Risk_std) = Sum(w_k^2 * Var(R_k))
  var_risk_std_exp <- sum(strat_details$w_k^2 * strat_details$risk_var_exp, na.rm = TRUE)
  var_risk_std_ref <- sum(strat_details$w_k^2 * strat_details$risk_var_ref, na.rm = TRUE)
  se_risk_std_exp <- sqrt(var_risk_std_exp)
  se_risk_std_ref <- sqrt(var_risk_std_ref)

  # Delta method SE for RR (log scale):
  # Var(log(RR)) ~ Var(R_exp)/R_exp^2 + Var(R_ref)/R_ref^2
  if (risk_std_exp > 0 && risk_std_ref > 0) {
    var_log_RR <- var_risk_std_exp / risk_std_exp^2 + var_risk_std_ref / risk_std_ref^2
    se_log_RR <- sqrt(var_log_RR)
    lnRR_se_analytical <- se_log_RR
  } else {
    se_log_RR <- NA_real_
    lnRR_se_analytical <- NA_real_
  }

  # SE for RD: Var(RD) = Var(R_exp) + Var(R_ref) (scaled by per_n)
  var_RD <- (var_risk_std_exp + var_risk_std_ref) * per_n^2
  RD_se_analytical <- sqrt(var_RD)

  return(list(
    risk_std_exp = risk_std_exp,
    risk_std_ref = risk_std_ref,
    RR = RR,
    RD = RD,
    lnRR_se_analytical = lnRR_se_analytical,
    RD_se_analytical = RD_se_analytical,
    se_log_RR = se_log_RR,
    se_risk_std_exp = se_risk_std_exp,
    se_risk_std_ref = se_risk_std_ref,
    var_risk_std_exp = var_risk_std_exp,
    var_risk_std_ref = var_risk_std_ref,
    strat_details = strat_details,
    strat_weights = strat_weights
  ))
}

#' Calculate RR and RD from cumulative incidence data (unstratified)
#' @noRd
calc_rr_rd <- function(cuminc_data, exp_val, ref_val, per_n) {

  risk_exp <- cuminc_data$Risk[cuminc_data$Exposure_Value == as.character(exp_val)]
  risk_ref <- cuminc_data$Risk[cuminc_data$Exposure_Value == as.character(ref_val)]

  if (length(risk_exp) == 0L || length(risk_ref) == 0L) {
    return(c(RR = NA_real_, RD = NA_real_, risk_exp = NA_real_, risk_ref = NA_real_))
  }

  RR <- risk_exp / risk_ref
  RD <- (risk_exp - risk_ref) * per_n

  return(c(RR = RR, RD = RD, risk_exp = risk_exp, risk_ref = risk_ref))
}

#' Single bootstrap iteration for unstratified RR/RD
#' match_ids (optional): match-set ids aligned to data rows; when supplied,
#' whole matched sets are resampled together (pair-level bootstrap).
#' @noRd
boot_rr_rd_single <- function(data, time_var, event_var, exp_var,
                              weight_var, timepoint, exp_val, ref_val, per_n,
                              risk_estimator = "KM",
                              competing_event_var = NULL,
                              match_ids = NULL) {

  if (is.null(match_ids)) {
    n <- nrow(data)
    boot_idx <- sample.int(n, size = n, replace = TRUE)
  } else {
    boot_idx <- pair_resample_index(match_ids)
  }
  boot_data <- data[boot_idx, , drop = FALSE]

  tryCatch({
    cuminc <- calc_km_cumulative_incidence(
      data = boot_data,
      time_var = time_var,
      event_var = event_var,
      exp_var = exp_var,
      weight_var = weight_var,
      timepoint = timepoint,
      per_n = per_n,
      risk_estimator = risk_estimator,
      competing_event_var = competing_event_var
    )
    rr_rd <- calc_rr_rd(cuminc, exp_val, ref_val, per_n)
    return(c(rr_rd["RR"], rr_rd["RD"]))
  }, error = function(e) {
    return(c(RR = NA_real_, RD = NA_real_))
  })
}

#' Single bootstrap iteration for stratified (standardized) RR/RD
#' Resamples within each stratum independently, applies fixed population
#' weights. match_ids_by_stratum (optional): list parallel to
#' data_by_stratum with the match-set ids of each stratum's rows; when
#' supplied, matched sets are resampled together within each stratum
#' (matched sets are assumed to be nested within strata).
#' @noRd
boot_standardized_rr_rd_single <- function(data_by_stratum, strat_var, strat_weights,
                                           time_var, event_var, exp_var,
                                           weight_var, timepoint, exp_val, ref_val, per_n,
                                           risk_estimator = "KM",
                                           competing_event_var = NULL,
                                           match_ids_by_stratum = NULL) {

  tryCatch({
    # Stratified resampling: resample within each stratum independently
    boot_data_list <- lapply(seq_along(data_by_stratum), function(k) {
      data_k <- data_by_stratum[[k]]
      if (is.null(match_ids_by_stratum)) {
        n_k <- nrow(data_k)
        boot_idx <- sample.int(n_k, size = n_k, replace = TRUE)
      } else {
        boot_idx <- pair_resample_index(match_ids_by_stratum[[k]])
      }
      data_k[boot_idx, , drop = FALSE]
    })
    boot_data <- do.call(rbind, boot_data_list)

    # Calculate standardized RR/RD using fixed original weights
    boot_split <- split(boot_data, boot_data[[strat_var]])

    std_result <- calc_standardized_rr_rd(
      data = boot_data,
      time_var = time_var,
      event_var = event_var,
      exp_var = exp_var,
      strat_var = strat_var,
      weight_var = weight_var,
      timepoint = timepoint,
      exp_val = exp_val,
      ref_val = ref_val,
      per_n = per_n,
      strat_weights = strat_weights,
      pre_split = boot_split,
      risk_estimator = risk_estimator,
      competing_event_var = competing_event_var
    )

    return(c(RR = std_result$RR, RD = std_result$RD))

  }, error = function(e) {
    return(c(RR = NA_real_, RD = NA_real_))
  })
}

#' Run parallel bootstrap and return results matrix
#'
#' If \code{seed} is \code{NULL} (default), no seed is set: the sequential path
#' (\code{n_cores = 1}) draws from and advances the caller's current RNG state
#' (standard R behaviour), while the parallel path seeds independent
#' L'Ecuyer-CMRG worker streams. If \code{seed} is supplied, it seeds the
#' parallel RNG stream (L'Ecuyer-CMRG) or, in the sequential fallback,
#' \code{set.seed()}; in the sequential case the caller's \code{.Random.seed} is
#' saved and restored on exit so the function has no side effect on the global
#' RNG state. (When two bootstraps run in one call -- e.g. the unweighted and
#' weighted fits -- a \code{NULL} seed keeps their sequential draws independent.)
#'
#' \code{export_varlist} names the helper objects exported to parallel workers
#' (defaults to the RR/RD helpers); pass \code{character(0)} when the bootstrap
#' function only uses namespaced calls (e.g. the Cox HR bootstrap).
#' \code{result_colnames} sets the column names of the returned matrix (defaults
#' to \code{c("RR", "RD")}). Both defaults preserve the original RR/RD behavior.
#' @noRd
run_parallel_bootstrap <- function(n_cores, bootstrap_count, boot_fn,
                                   seed = NULL,
                                   export_varlist = c(
                                     "calc_km_cumulative_incidence",
                                     "calc_km_risk_single_group",
                                     "calc_standardized_rr_rd",
                                     "calc_rr_rd"
                                   ),
                                   result_colnames = c("RR", "RD"),
                                   ...) {

  args <- list(...)

  if (n_cores > 1L) {
    cl <- parallel::makeCluster(n_cores)
    on.exit(parallel::stopCluster(cl), add = TRUE)

    # Export helper functions needed by boot_fn to the workers (if any).
    if (length(export_varlist) > 0L) {
      parallel::clusterExport(
        cl,
        varlist = export_varlist,
        envir = environment(boot_fn)
      )
    }

    # Parallel RNG (L'Ecuyer-CMRG). Only seed when the user supplies one.
    if (!is.null(seed)) {
      parallel::clusterSetRNGStream(cl, iseed = seed)
    } else {
      parallel::clusterSetRNGStream(cl)
    }

    boot_list <- parallel::parLapply(
      cl,
      X = seq_len(bootstrap_count),
      fun = function(b) {
        do.call(boot_fn, args)
      }
    )
  } else {
    # Sequential fallback. Save and restore the caller's RNG state when we
    # touch the global seed, so the function does not have side effects.
    if (!is.null(seed)) {
      if (exists(".Random.seed", envir = globalenv(), inherits = FALSE)) {
        old_seed <- get(".Random.seed", envir = globalenv(), inherits = FALSE)
        on.exit(assign(".Random.seed", old_seed, envir = globalenv()),
                add = TRUE)
      } else {
        on.exit(rm(".Random.seed", envir = globalenv()), add = TRUE)
      }
      set.seed(seed)
    }
    boot_list <- lapply(seq_len(bootstrap_count), function(b) {
      do.call(boot_fn, args)
    })
  }

  boot_results <- do.call(rbind, boot_list)
  colnames(boot_results) <- result_colnames
  return(boot_results)
}


# IRR Helper #########################################################
#' Fit the always-on marginal incidence-rate-ratio model for one block.
#'
#' Fixed per-block method (matrix main-text rule, no user choice):
#' unweighted (crude/matched) blocks use a marginal Poisson rate model with a
#' log(person-time) offset — the model SE equals the matrix r10-A log-ratio
#' Poisson approximation sqrt(1/D1 + 1/D0) exactly for this saturated
#' marginal fit. Weighted blocks use survey::svyglm with a quasi-Poisson
#' family, i.e. the design-based robust sandwich SE (matrix r12-A).
#'
#' @param data Data frame with recoded (0/1) exposure, outcome, person-time.
#' @param outcome_var,exposure_var,py_var Column names.
#' @param weight_var Weight column name, or NULL for an unweighted fit.
#' @param confidence_level Numeric in (0, 1).
#' @param base_label Character label for the Analysis column.
#' @return One-row data frame: Analysis, IRR_CI, IRR, IRR_LCI, IRR_UCI, lnIRR,
#'   lnIRR_SE, Pvalue, Model, SE_type.
#' @noRd
fit_irr_model <- function(data, outcome_var, exposure_var, py_var = "person_years",
                          weight_var = NULL,
                          confidence_level = 0.95, base_label = "Crude") {
  is_weighted <- !is.null(weight_var)
  z <- stats::qnorm(1 - (1 - confidence_level) / 2)

  fml <- stats::as.formula(
    paste0(outcome_var, " ~ ", exposure_var, " + offset(log(", py_var, "))")
  )

  if (is_weighted) {
    des <- survey::svydesign(
      ids = ~1,
      weights = stats::as.formula(paste0("~", weight_var)),
      data = data
    )
    fit <- survey::svyglm(fml, design = des, family = stats::quasipoisson())
    se_type     <- "robust (svyglm)"
    model_label <- "quasipoisson (weighted)"
  } else {
    fit <- stats::glm(fml, data = data, family = stats::poisson())
    se_type     <- "model (Poisson, r10-A)"
    model_label <- "poisson"
  }

  s   <- summary(fit)
  cf  <- s$coefficients
  ridx <- which(rownames(cf) == exposure_var)
  ln_irr <- unname(cf[ridx, 1L])
  ln_se  <- unname(cf[ridx, 2L])
  p_val  <- unname(cf[ridx, 4L])
  irr    <- exp(ln_irr)
  lci    <- exp(ln_irr - z * ln_se)
  uci    <- exp(ln_irr + z * ln_se)

  data.frame(
    Analysis   = base_label,
    IRR_CI     = sprintf("%.2f (%.2f, %.2f)", irr, lci, uci),
    IRR        = irr,
    IRR_LCI    = lci,
    IRR_UCI    = uci,
    lnIRR      = ln_irr,
    lnIRR_SE   = ln_se,
    Pvalue     = p_val,
    Model      = model_label,
    SE_type    = se_type,
    stringsAsFactors = FALSE,
    row.names  = NULL
  )
}


# IRR Helper (stratified paths) ######################################
#' Stratified marginal IRR for one block (always-on; matrix r11-A / r13-A).
#' Unweighted (crude/matched): direct-standardized IRR from per-stratum
#' Poisson cells via direct_std_irr(). Weighted: joint-sandwich
#' log-ratio delta via fit_weighted_rate_cells() + wtd_std_rate_contrasts().
#' Standardization weights = total (weighted) person-time share per stratum,
#' the same W_s used by the standardized IR/IRD calculators.
#'
#' @return One-row data frame with the same columns as fit_irr_model().
#' @noRd
fit_irr_model_stratified <- function(data, outcome_var, exposure_var,
                                     py_var = "person_years", strat_var,
                                     weight_var = NULL,
                                     confidence_level = 0.95,
                                     base_label = "Crude") {
  z <- stats::qnorm(1 - (1 - confidence_level) / 2)
  is_weighted <- !is.null(weight_var)

  if (is_weighted) {
    pt <- data[[py_var]] * data[[weight_var]]
    pt_s <- vapply(split(pt, data[[strat_var]]), function(x) sum(x, na.rm = TRUE),
                   numeric(1))
    std_w <- pt_s / sum(pt_s)
    cell_fit <- fit_weighted_rate_cells(data, exposure_var, outcome_var,
                                        py_var, weight_var, strat_var = strat_var)
    if (is.null(cell_fit)) {
      irr <- list(est = NA_real_, lci = NA_real_, uci = NA_real_,
                  ln = NA_real_, ln_se = NA_real_)
    } else {
      irr <- wtd_std_rate_contrasts(cell_fit, std_w = std_w, z = z,
                                    multiplier = 1)$irr
    }
    model_label <- "direct-standardized (weighted joint sandwich)"
    se_type     <- "robust joint sandwich / log-ratio delta"
  } else {
    pt_s <- vapply(split(data[[py_var]], data[[strat_var]]),
                   function(x) sum(x, na.rm = TRUE), numeric(1))
    std_w <- pt_s / sum(pt_s)
    lev <- names(std_w)
    ev <- function(a, s) sum(data[[outcome_var]][data[[exposure_var]] == a &
                                                   data[[strat_var]] == s], na.rm = TRUE)
    py <- function(a, s) sum(data[[py_var]][data[[exposure_var]] == a &
                                              data[[strat_var]] == s], na.rm = TRUE)
    irr <- direct_std_irr(
      events1 = vapply(lev, function(s) ev(1, s), numeric(1)),
      py1     = vapply(lev, function(s) py(1, s), numeric(1)),
      events0 = vapply(lev, function(s) ev(0, s), numeric(1)),
      py0     = vapply(lev, function(s) py(0, s), numeric(1)),
      std_w   = std_w, z = z
    )
    model_label <- "direct-standardized (Poisson strata)"
    se_type     <- "log-ratio delta (per-stratum Poisson)"
  }

  p_val <- if (is.finite(irr$ln) && is.finite(irr$ln_se) && irr$ln_se > 0) {
    2 * stats::pnorm(-abs(irr$ln / irr$ln_se))
  } else {
    NA_real_
  }

  data.frame(
    Analysis   = base_label,
    IRR_CI     = if (is.na(irr$est)) NA_character_ else
      sprintf("%.2f (%.2f, %.2f)", irr$est, irr$lci, irr$uci),
    IRR        = irr$est,
    IRR_LCI    = irr$lci,
    IRR_UCI    = irr$uci,
    lnIRR      = irr$ln,
    lnIRR_SE   = irr$ln_se,
    Pvalue     = p_val,
    Model      = model_label,
    SE_type    = se_type,
    stringsAsFactors = FALSE,
    row.names  = NULL
  )
}


# Fine-Gray Helper ###################################################
#' Fit a Fine-Gray subdistribution-hazard regression for one analysis block
#'
#' Builds a competing-risks multi-state status (1 = event of interest,
#' 2 = competing event, 0 = censored), expands the data with
#' \code{survival::finegray} (passing the propensity-score weights so that the
#' returned \code{fgwt} already equals IPCW x PS weight), and fits the
#' subdistribution Cox model with \code{weights = fgwt} and a robust SE clustered
#' on the original subject. The original subject id and the stratification
#' variable are carried as ordinary right-hand-side terms so they survive the
#' \code{finegray} expansion (an \code{id} argument alone would not). When
#' \code{strat_var} is supplied, \code{strata()} is used in BOTH \code{finegray}
#' (stratum-specific censoring distributions) and \code{coxph} (stratified
#' subdistribution baseline).
#'
#' @param data Data frame with recoded (0/1) exposure, outcome, competing-event
#'   indicator, time, optional weight and stratification columns.
#' @param time_var,event_var,competing_var,exp_var Column names (character).
#' @param weight_var Weight column name, or NULL (unweighted).
#' @param strat_var Stratification column name, or NULL.
#' @param cluster_id_var Optional cluster-id column (e.g. a match id) for the
#'   robust SE; defaults to the subject-level expansion id (.fg_id).
#' @param confidence_level Numeric in (0, 1).
#' @param base_label Character label for the Analysis column.
#' @return One-row data frame: Analysis, SHR_CI, SHR, SHR_LCI, SHR_UCI, lnSHR,
#'   lnSHR_SE, Pvalue. Returns NA estimates (with a warning) if the fit fails.
#' @noRd
fit_finegray_shr <- function(data, time_var, event_var, competing_var, exp_var,
                             weight_var = NULL, strat_var = NULL,
                             cluster_id_var = NULL,
                             confidence_level = 0.95, base_label = "Unweighted") {
  is_weighted <- !is.null(weight_var)
  z <- stats::qnorm(1 - (1 - confidence_level) / 2)

  na_row <- data.frame(
    Analysis = base_label, SHR_CI = NA_character_, SHR = NA_real_,
    SHR_LCI = NA_real_, SHR_UCI = NA_real_, lnSHR = NA_real_,
    lnSHR_SE = NA_real_, Pvalue = NA_real_,
    stringsAsFactors = FALSE, row.names = NULL
  )

  tryCatch({
    work <- data
    work[[".ms"]] <- factor(
      ifelse(work[[event_var]] == 1, 1L,
             ifelse(work[[competing_var]] == 1, 2L, 0L)),
      levels = c(0L, 1L, 2L), labels = c("censor", "event", "compete")
    )
    work[[".fg_id"]] <- seq_len(nrow(work))

    # Carry .fg_id (+ strat_var) as plain RHS terms so they survive the
    # finegray expansion; strata() sets stratum-specific censoring weights.
    rhs_terms <- c(exp_var, ".fg_id")
    if (!is.null(cluster_id_var)) rhs_terms <- c(rhs_terms, cluster_id_var)
    if (!is.null(strat_var)) {
      rhs_terms <- c(rhs_terms, strat_var, paste0("survival::strata(", strat_var, ")"))
    }
    fg_fml <- stats::as.formula(
      paste0("survival::Surv(", time_var, ", .ms) ~ ", paste(rhs_terms, collapse = " + "))
    )

    fg <- if (is_weighted) {
      survival::finegray(fg_fml, data = work, weights = work[[weight_var]], etype = "event")
    } else {
      survival::finegray(fg_fml, data = work, etype = "event")
    }

    cox_rhs <- exp_var
    if (!is.null(strat_var)) {
      cox_rhs <- paste0(cox_rhs, " + survival::strata(", strat_var, ")")
    }
    cox_fml <- stats::as.formula(
      paste0("survival::Surv(fgstart, fgstop, fgstatus) ~ ", cox_rhs)
    )
    cl <- if (!is.null(cluster_id_var)) fg[[cluster_id_var]] else fg$.fg_id
    fit <- survival::coxph(cox_fml, data = fg, weights = fg$fgwt,
                           cluster = cl, robust = TRUE)

    s    <- summary(fit, conf.int = confidence_level)
    cf   <- s$coefficients
    ci   <- s$conf.int
    ridx <- which(rownames(cf) == exp_var)
    se_col <- if ("robust se" %in% colnames(cf)) "robust se" else "se(coef)"

    shr    <- unname(ci[ridx, "exp(coef)"])
    lci    <- unname(ci[ridx, 3])
    uci    <- unname(ci[ridx, 4])
    ln_se  <- unname(cf[ridx, se_col])
    p_val  <- unname(cf[ridx, "Pr(>|z|)"])

    data.frame(
      Analysis = base_label,
      SHR_CI   = sprintf("%.2f (%.2f, %.2f)", shr, lci, uci),
      SHR      = shr,
      SHR_LCI  = lci,
      SHR_UCI  = uci,
      lnSHR    = log(shr),
      lnSHR_SE = ln_se,
      Pvalue   = p_val,
      stringsAsFactors = FALSE, row.names = NULL
    )
  }, error = function(e) {
    warning(sprintf("Fine-Gray regression failed for '%s': %s", base_label, conditionMessage(e)))
    na_row
  })
}