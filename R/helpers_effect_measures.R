######### Shared input layer ##########################################
# HELPER: shared data-block validation for the effect functions ########
#' Validates the in_df_crude / in_df_weight / in_df_match data blocks and
#' their if_* companions (API v2 directive): weight XOR match; crude-less
#' weight/match calls allowed; orphan if_* arguments are ignored with a
#' message; required variables must exist in every provided block. Also
#' rejects matched-cohort designs the matched block cannot analyse:
#' with-replacement match ids (backlog Item 6) and non-unit matching weights
#' (Item 5).
#'
#' @param in_df_crude,in_df_weight,in_df_match Data blocks (or NULL).
#' @param if_weight_weight_var Weight column name for in_df_weight, or NULL.
#' @param if_match_match_id Match-set id column name for in_df_match, or NULL.
#' @param required_vars Character vector of columns every provided block must
#'   contain (exposure/outcome/time and, when used, strata etc.).
#' @param verbose Logical; emit the ignore-messages.
#' @note Matched-set composition (strict 1:1) is NOT checked here: it needs the
#'   canonicalized exposure and the complete-case rows, so it lives in
#'   \code{prep_effect_block()} via \code{check_matched_1to1_support()}.
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

      # Item 6: with-replacement matching. create_ps_matched_cohort() joins the
      # ids of every set a reused control belongs to with ";" (e.g.
      # "M00001;M00007"). Such a value is not a cluster: cluster(match_id) and
      # the pair-resample bootstrap would both treat it as a separate third
      # set, which is neither the Abadie-Imbens variance nor a valid
      # sandwich cluster.
      if (any(grepl(";", in_df_match[[if_match_match_id]], fixed = TRUE),
              na.rm = TRUE)) {
        stop("in_df_match contains match ids that name more than one matched ",
             "set (e.g. 'M00001;M00007'), which happens with ",
             "create_ps_matched_cohort(replace = TRUE). With-replacement ",
             "matched cohorts are not supported in the matched block: the ",
             "pair-clustered SE and the pair-resample bootstrap both require ",
             "each row to belong to exactly one set. Either match without ",
             "replacement, or pass the cohort as in_df_weight with ",
             "if_weight_weight_var = '.match_weights' -- which gives a ",
             "control-reuse-weighted analysis with a weight-based robust SE, ",
             "not the Abadie-Imbens with-replacement variance.")
      }
    }
  } else if (!is.null(if_match_match_id)) {
    if (verbose) message("if_match_match_id is ignored because in_df_match was not provided.")
  }

  # Item 5: the matched block fits unweighted models.  It therefore accepts
  # matching weights only when every retained row has unit weight.  The older
  # within-arm-constant rule was sufficient for arm-level rates/risks, but not
  # for the Cox partial likelihood when the two arms carried different
  # constants.  Reject that ambiguity instead of silently choosing an
  # estimand/model the caller did not request.
  if (!is.null(in_df_match) && ".match_weights" %in% names(in_df_match)) {
    mw <- in_df_match$.match_weights
    invalid <- !is.numeric(mw) || anyNA(mw) || any(!is.finite(mw)) ||
      any(abs(mw - 1) > 1e-8)
    if (invalid) {
      stop("in_df_match is analysed as an unweighted strict 1:1 matched ",
           "cohort, so every .match_weights value must equal 1. Non-unit ",
           "weights can change at least one supported estimand (including ",
           "the Cox partial likelihood). For a weight-defined analysis, pass ",
           "the cohort as in_df_weight with if_weight_weight_var = ",
           "'.match_weights'. That route honours the weights but does not use ",
           "matched-set clustering, so it is not a substitute for ",
           "design-aware inference for variable-ratio, subclass, or ",
           "with-replacement matching.")
    }
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


# HELPER: per-block preparation, shared by all effect estimators #############
#' Prepares one analysis block: canonicalizes exposure so that 1 is always
#' \code{exp_value} and 0 always \code{ref_value}, optionally derives
#' \code{person_years} from the follow-up time, drops incomplete rows on the
#' columns the analysis actually uses, requires both arms to survive, and
#' (optionally) enforces an unambiguous competing-event coding.
#'
#' The exposure recode is UNCONDITIONAL. Every downstream calculator selects
#' arms with a hard-coded \code{== 1} / \code{== 0}, so skipping the recode
#' when the column already looks like 0/1 -- as 0.4.0 did -- silently ignored
#' \code{exp_value = 0, ref_value = 1} and analysed the wrong arm (backlog
#' Item 8).
#'
#' Promoted out of the 0.4.0 \code{estimate_hr_ir()} local \code{prep_block()}
#' when that function was split, so both public functions share one
#' preparation path. The complete-case column set is deliberately
#' function-specific: \code{estimate_ir()} passes no \code{strata_var} and no
#' \code{competing_var}, so it keeps rows that \code{estimate_hr()} would drop
#' when those columns contain NA.
#'
#' @param df Data frame for one block.
#' @param block Block label ("Crude", "Weighted", "Matched") -- used for
#'   messages and to decide which optional columns are required.
#' @param exposure_var,exp_value,ref_value Exposure column and the two values
#'   identifying the exposed and reference arms. Must differ. Rows matching
#'   neither are dropped with a warning.
#' @param outcome_var,followuptime_var Column names (character).
#' @param time_unit One of "days", "months", "years" to derive
#'   \code{person_years}, or NULL to skip that derivation.
#' @param weight_var Weight column (required for the Weighted block).
#' @param match_id Match-set id column (Matched block), or NULL.
#' @param strata_var Stratification column to require complete, or NULL.
#' @param competing_var Competing-event column to require complete, or NULL.
#' @param check_competing Logical. When TRUE, \code{competing_var} must be
#'   binary and no row may be both the event of interest and a competing
#'   event.
#' @param competing_context Label used in competing-event errors. One of
#'   "Fine-Gray" or "Aalen-Johansen".
#' @param verbose Logical.
#' @return The prepared data frame.
#' @noRd
prep_effect_block <- function(df, block, exposure_var, exp_value, ref_value,
                              outcome_var, followuptime_var,
                              time_unit = NULL,
                              weight_var = NULL, match_id = NULL,
                              strata_var = NULL, competing_var = NULL,
                              check_competing = FALSE,
                              competing_context = "Fine-Gray",
                              verbose = TRUE) {
  # --- canonicalize exposure to 1 = exp_value, 0 = ref_value ---------------
  # ALWAYS recoded (backlog Item 8). Up to 0.4.0 the recode was skipped when
  # the column already looked like 0/1, so exp_value = 0 / ref_value = 1 on an
  # already-0/1 column was silently ignored and every downstream calculator --
  # which hard-codes == 1 / == 0 -- analysed the wrong arm, inverting the
  # hazard/rate ratio and flipping the sign of the rate difference without a
  # word. Recoding unconditionally is idempotent in the default case.
  ev <- df[[exposure_var]]
  already_canonical <- is.numeric(ev) && all(ev[!is.na(ev)] %in% c(0, 1)) &&
    isTRUE(exp_value == 1) && isTRUE(ref_value == 0)
  exposure_01 <- .canonicalize_exposure(
    ev, exp_value, ref_value, exposure_var = exposure_var,
    require_both = FALSE, warn_unmatched = FALSE
  )
  n_unmatched <- exposure_01$n_unmatched
  df[[exposure_var]] <- exposure_01$value
  if (!already_canonical && verbose) {
    message(sprintf("Recoding exposure (%s block): '%s' -> 1 (exposed), '%s' -> 0 (reference).",
                    block, exp_value, ref_value))
  }
  if (n_unmatched > 0L) {
    warning(sprintf(paste0("%d row(s) in the %s block match neither exp_value ",
                           "('%s') nor ref_value ('%s') and were dropped."),
                    n_unmatched, block, exp_value, ref_value))
  }

  # person-years
  used <- c(outcome_var, exposure_var, followuptime_var)
  if (!is.null(time_unit)) {
    df$person_years <- switch(time_unit,
      days   = as.numeric(df[[followuptime_var]]) / 365.25,
      months = as.numeric(df[[followuptime_var]]) / 12,
      years  = as.numeric(df[[followuptime_var]])
    )
    used <- c(used, "person_years")
  }

  # complete cases on the columns this analysis uses
  if (block == "Weighted") used <- c(used, weight_var)
  if (block == "Matched" && !is.null(match_id)) used <- c(used, match_id)
  if (!is.null(strata_var)) used <- c(used, strata_var)
  if (!is.null(competing_var)) used <- c(used, competing_var)
  n0 <- nrow(df)
  df <- df[stats::complete.cases(df[, used, drop = FALSE]), , drop = FALSE]
  if (nrow(df) < n0) {
    warning(sprintf("%d observation(s) with missing values removed from the %s block",
                    n0 - nrow(df), block))
  }

  # both arms must survive, or every contrast below is undefined
  arm_n <- table(factor(df[[exposure_var]], levels = c(0L, 1L)))
  if (any(arm_n == 0L)) {
    empty <- if (arm_n[["1"]] == 0L) {
      sprintf("exposed (exp_value = '%s')", exp_value)
    } else {
      sprintf("reference (ref_value = '%s')", ref_value)
    }
    stop(sprintf(paste0("The %s block has no %s rows after recoding and ",
                        "complete-case removal; a contrast cannot be estimated."),
                 block, empty))
  }

  # competing risks need an unambiguous event type
  if (check_competing) {
    if (!competing_context %in% c("Fine-Gray", "Aalen-Johansen")) {
      stop("competing_context must be 'Fine-Gray' or 'Aalen-Johansen'.")
    }
    cev <- df[[competing_var]]
    cev_vals <- unique(cev[!is.na(cev)])
    competing_arg <- if (competing_context == "Aalen-Johansen") {
      "if_aj_competing_event_var"
    } else {
      "if_fg_competing_event_var"
    }
    if (!all(cev_vals %in% c(0, 1))) {
      stop(paste0(competing_arg, " ('", competing_var,
                   "') must be binary (0/1) in the ", block, " block."))
    }
    n_both <- sum(df[[outcome_var]] == 1 & cev == 1, na.rm = TRUE)
    if (n_both > 0L) {
      stop(sprintf(paste0("%s: %d row(s) in the %s block have both ",
                          "outcome_var == 1 and %s == 1, ",
                          "which is ambiguous. Resolve these first."),
                   competing_context, n_both, block, competing_arg))
    }
  }

  if (block == "Matched" && !is.null(match_id)) {
    check_matched_1to1_support(df, exposure_var, match_id, block)
  }

  df
}


# HELPER: strict matched-set composition #####################################
#' Require each retained match id to identify exactly one exposed/reference
#' pair. Runs after recoding and complete-case removal, so a pair broken by a
#' missing analysis value is rejected rather than silently analysed as an
#' independent or malformed cluster.
#' @noRd
check_matched_1to1_support <- function(df, exp_var, match_id, block = "Matched") {
  tab <- table(as.character(df[[match_id]]),
               factor(df[[exp_var]], levels = c(0L, 1L)))
  bad <- rownames(tab)[rowSums(tab) != 2L | tab[, "0"] != 1L | tab[, "1"] != 1L]
  if (length(bad) > 0L) {
    shown <- utils::head(bad, 5L)
    detail <- vapply(shown, function(id) {
      sprintf("'%s' (%d reference, %d exposed)", id, tab[id, "0"], tab[id, "1"])
    }, character(1))
    suffix <- if (length(bad) > length(shown)) "; ..." else ""
    stop(sprintf(paste0(
      "The %s block supports strict 1:1 matching only: every match id must ",
      "contain exactly one reference and one exposed row after complete-case ",
      "removal. Invalid set(s): %s%s"
    ), block, paste(detail, collapse = "; "), suffix))
  }
  invisible(TRUE)
}


# HELPER: arm x stratum support for a stratified hazard model ###########
#' Requires every level of \code{strata_var} to contain rows from BOTH
#' exposure arms.
#'
#' A \code{survival::strata()} level in which only one arm appears contributes
#' nothing to the stratified partial likelihood -- there is no within-stratum
#' comparison to make -- and \code{coxph()} absorbs it silently. The reported
#' conditional common (s)HR is then estimated on a subset of the data the
#' caller never asked to restrict to, with no indication in the output. Since
#' a stratum missing an entire arm almost always signals a data problem
#' (e.g. cohort-entry-year groups where the reference drug was not yet
#' marketed), this is an error rather than a warning.
#'
#' Runs AFTER exposure recoding and complete-case removal, so it sees the rows
#' the model will actually use.
#'
#' @param df One prepared analysis block (exposure already 0/1).
#' @param exp_var,strat_var Column names (character).
#' @param block Block label, for the message.
#' @return \code{invisible(TRUE)}; called for the side effect.
#' @noRd
check_arm_stratum_support <- function(df, exp_var, strat_var, block) {
  tab <- table(as.character(df[[strat_var]]),
               factor(df[[exp_var]], levels = c(0L, 1L)))
  empty <- rownames(tab)[tab[, "0"] == 0L | tab[, "1"] == 0L]
  if (length(empty) > 0L) {
    detail <- vapply(empty, function(lv) sprintf("'%s' (%d reference, %d exposed)",
                                                 lv, tab[lv, "0"], tab[lv, "1"]),
                     character(1))
    stop(sprintf(paste0("strata_var '%s' has %d level(s) in the %s block with ",
                        "no rows in one exposure arm: %s. survival::strata() ",
                        "would drop these levels from the partial likelihood ",
                        "without saying so, so the reported hazard ratio would ",
                        "silently describe only the remaining strata. Collapse ",
                        "or exclude these levels first."),
                 strat_var, length(empty), block,
                 paste(detail, collapse = "; ")))
  }
  invisible(TRUE)
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

# HELPER: design-based IR/IRD (analysis weights and/or clustering) #####
#' IR / IRD for a block whose variance needs a design-based sandwich: an
#' analysis-weighted block (wt_var, e.g. IPTW -- matrix r4/r8) and/or a matched
#' block whose rows are clustered in matched sets (cluster_var).
#'
#' Point estimates are the plain (weighted) sums either way; only the variance
#' uses the joint cell fit. With wt_var = NULL the sums reduce exactly to
#' calc_crude_ir_ird()'s, so routing a matched block here changes the intervals
#' and nothing else. Matched-set clustering replaces the independent-Poisson
#' Garwood / Wald intervals with cluster-robust log-Wald and joint-delta ones,
#' which may be narrower OR wider than the independence intervals depending on
#' the sign of the within-set covariance.
#' @noRd
calc_crude_ir_ird_design <- function(df, exp_var, out_var, py_var,
                                     wt_var = NULL, cluster_var = NULL,
                                     conf_level, z, suffix, multiplier, ir_base) {
  w <- if (is.null(wt_var)) rep(1, nrow(df)) else df[[wt_var]]
  df$.w_ <- w

  summary_by_group <- do.call(rbind, lapply(split(df, df[[exp_var]]), function(g) {
    data.frame(
      exp_val      = g[[exp_var]][1],
      n_subjects   = sum(g$.w_, na.rm = TRUE),
      n_events     = sum(g[[out_var]] * g$.w_, na.rm = TRUE),
      person_years = sum(g[[py_var]] * g$.w_, na.rm = TRUE),
      stringsAsFactors = FALSE
    )
  }))
  names(summary_by_group)[1] <- exp_var
  rownames(summary_by_group) <- NULL
  summary_by_group$ir_value <- multiplier * summary_by_group$n_events / summary_by_group$person_years

  # Robust (design-based sandwich) arm SE/CI via the joint cell engine
  # (matrix r4; P1): the previous rounded-qgamma path treated the weighted
  # event sum as a true Poisson count and was anti-conservative under IPTW.
  cell_fit <- fit_weighted_rate_cells(df, exp_var, out_var, py_var, wt_var,
                                      cluster_var = cluster_var)
  if (!is.null(cell_fit)) {
    ctr <- wtd_rate_contrasts(cell_fit, z, multiplier)
    ord <- match(as.character(summary_by_group[[exp_var]]), ctr$ir$exp_val)
    summary_by_group$ir_lci <- ctr$ir$lci[ord]
    summary_by_group$ir_uci <- ctr$ir$uci[ord]
    summary_by_group$IR_SE  <- ctr$ir$se[ord]
  } else {
    warning("Design-based IR robust variance fit failed; returning NA CIs. ",
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
  total_ne <- sum(df[[out_var]] * w, na.rm = TRUE)
  total_py <- sum(df[[py_var]] * w, na.rm = TRUE)
  total_row <- data.frame(
    exp_val      = 99,
    n_subjects   = sum(w, na.rm = TRUE),
    n_events     = total_ne,
    person_years = total_py,
    stringsAsFactors = FALSE
  )
  names(total_row)[1] <- exp_var
  total_row$ir_value <- multiplier * total_ne / total_py

  df_tot <- df
  df_tot$.const <- 0L
  tot_fit <- fit_weighted_rate_cells(df_tot, ".const", out_var, py_var, wt_var,
                                     cluster_var = cluster_var)
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
                          weight_var = NULL, cluster_var = NULL,
                          confidence_level = 0.95, base_label = "Crude") {
  is_weighted  <- !is.null(weight_var)
  is_clustered <- !is.null(cluster_var)
  z <- stats::qnorm(1 - (1 - confidence_level) / 2)

  fml <- stats::as.formula(
    paste0(outcome_var, " ~ ", exposure_var, " + offset(log(", py_var, "))")
  )

  if (is_weighted || is_clustered) {
    ids_f <- if (is_clustered) {
      stats::as.formula(paste0("~`", cluster_var, "`"))
    } else {
      ~1
    }
    wt_f <- if (is_weighted) {
      stats::as.formula(paste0("~`", weight_var, "`"))
    } else {
      data$.unit_w_ <- 1     # explicit unit weights: silences svydesign's note
      ~.unit_w_
    }
    des <- survey::svydesign(ids = ids_f, weights = wt_f, data = data)
    fit <- survey::svyglm(fml, design = des, family = stats::quasipoisson())
    se_type <- if (is_weighted && is_clustered) {
      "robust (svyglm, weighted + match-id cluster)"
    } else if (is_weighted) {
      "robust (svyglm)"
    } else {
      "robust (svyglm, match-id cluster)"
    }
    model_label <- if (is_weighted) {
      "quasipoisson (weighted)"
    } else {
      "quasipoisson (clustered)"
    }
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
