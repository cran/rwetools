# RR Risk Ratios and Risk Differences ##############
#' Estimate Marginal Risks, Risk Ratios, and Risk Differences
#'
#' Estimates arm-specific cumulative incidence, an exposed-versus-reference
#' risk ratio (RR), and an exposed-minus-reference risk difference (RD) at one
#' time point. A call can contain a crude cohort plus either a weighted cohort
#' or a matched cohort. The weighted and matched blocks are mutually exclusive,
#' and the crude block is optional.
#'
#' Direct standardization was removed in rwetools 0.5.0. This function has no
#' stratification argument: every reported risk measure is marginal over the
#' rows in its analysis block. To address a baseline stratifier, include it in
#' the propensity-score model; alternatively, run separate analyses within
#' levels and combine them explicitly in caller code under a prespecified
#' target-population rule.
#'
#' `risk_estimator = "KM"` reports `1 - S(t)` from Kaplan-Meier.
#' `risk_estimator = "AJ"` reports the Aalen-Johansen cumulative incidence
#' function and requires a separate binary competing-event indicator.
#'
#' @section Exposure coding and analysis rows:
#' `exp_value` and `ref_value` are always recoded internally to 1 and 0,
#' including when the source column is already coded 0/1. They must be distinct
#' scalar, non-missing values. Rows matching neither value are reported and
#' removed; rows incomplete on variables used by this risk analysis are also
#' removed. Both arms must remain. For a logical exposure column, use
#' `exp_value = TRUE` and `ref_value = FALSE`, not numeric 1 and 0.
#'
#' For Aalen-Johansen, the competing-event indicator must be binary
#' after complete-case removal, and a row cannot have both indicators
#' equal to 1.
#'
#' @section Analytical confidence intervals:
#' Methods are fixed; there are no CI-method arguments. Risk/CIF intervals use
#' a complementary log-log transformation of the estimated risk and its
#' Greenwood (KM) or counting-process (AJ) SE. At a risk of exactly 0 or 1 the
#' transform is undefined, so the bounds are `NA`. RR uses a log-scale delta
#' interval and RD uses a normal-Wald interval.
#'
#' The analytical variance in a matched block treats rows as independent. It
#' does not model within-pair covariance, even when `if_match_match_id` is
#' supplied, and is typically mildly conservative for 1:1 matching. Request a
#' bootstrap and supply `if_match_match_id` for pair-aware resampling.
#'
#' @section Matched designs and matching weights:
#' The matched block is an unweighted matched analysis. If a `.match_weights`
#' column is present, every retained value must equal 1. With
#' `if_match_match_id`, every complete-case set must contain exactly one
#' exposed and one reference row; matching with replacement is rejected.
#'
#' For subclass, variable-ratio, fixed ratios greater than 1, or
#' with-replacement matching, pass the cohort as `in_df_weight` and set
#' `if_weight_weight_var = ".match_weights"` to honour the point-estimate
#' weights. This route still uses fixed-case-weight Greenwood or
#' Aalen-Johansen variance for risk measures and does not use matched-set
#' clustering. It therefore cannot express weights and set clustering
#' together and is not design-aware inference for variable-ratio, subclass, or
#' replacement matching. The package does not produce a
#' one-row-per-unit-per-pair replacement dataset.
#'
#' @section Bootstrap:
#' `if_bootstrap_count` adds percentile intervals. Supplied analysis weights
#' are frozen, so propensity-score estimation uncertainty is not propagated.
#' A matched block with `if_match_match_id` resamples whole matched sets;
#' without the id it resamples rows. The seed controls this function only and
#' does not establish shared draws with [estimate_hr()] or [estimate_ir()].
#'
#' @section Migration from rwetools 0.4.0:
#' `estimate_rr_rd()` is now `estimate_risk()`;
#' `rr_rd_at_timepoint` is now `risk_at_timepoint`. The old
#' `stratification_var` and direct-standardization output, including the
#' `Stratified_By` column and stratum-details sheet, were removed.
#'
#' @param in_df_crude Data frame for a crude block, or `NULL`.
#' @param in_df_weight Data frame for a weighted block, or `NULL`. Requires
#'   `if_weight_weight_var` and cannot be combined with `in_df_match`.
#' @param in_df_match Data frame for a matched block, or `NULL`. Cannot be
#'   combined with `in_df_weight`.
#' @param out_xlsxpath Character path for an Excel workbook, or `NULL`.
#' @param exposure_var Character name of the exposure column. The default is
#'   `"exp"`.
#' @param exp_value,ref_value Values identifying the exposed and reference
#'   arms in `exposure_var`.
#' @param outcome_var Character name of the event indicator, coded 1 for the
#'   event of interest.
#' @param followuptime_var Character name of the follow-up time column.
#' @param time_unit Unit shared by `followuptime_var` and
#'   `risk_at_timepoint`: `"days"`, `"months"`, or `"years"`.
#' @param risk_at_timepoint Single positive time at which risk, RR, and RD are
#'   evaluated. The default is 365.
#' @param risk_per_individuals Positive reporting multiplier for arm risks and
#'   RD. The default is 1000.
#' @param confidence_level Single number strictly between 0 and 1. The default
#'   is 0.95.
#' @param risk_estimator `"KM"` or `"Kaplan-Meier"`; or `"AJ"` or
#'   `"Aalen-Johansen"` for competing risks.
#' @param if_aj_competing_event_var Character name of a binary competing-event
#'   indicator, coded 1 for a competing event. Required for AJ and ignored with
#'   a message for KM.
#' @param if_weight_weight_var Character name of the analysis-weight column in
#'   `in_df_weight`. Required for a weighted block and ignored with a message
#'   otherwise.
#' @param if_match_match_id Character name of the set id in `in_df_match`.
#'   It defines the pair-bootstrap unit but does not alter analytical risk,
#'   RR, or RD variance. Without it, bootstrap resampling is by row. Ignored
#'   with a message without a matched block.
#' @param if_bootstrap_count Positive number of percentile-bootstrap
#'   replicates, or `NULL` to omit bootstrap CIs. Weights are frozen at their
#'   supplied values.
#' @param if_bootstrap_n_cores Number of bootstrap workers, or `NULL` to use
#'   all detected cores minus one. Ignored without `if_bootstrap_count`.
#' @param if_bootstrap_seed Integer RNG seed, or `NULL`. Ignored without
#'   `if_bootstrap_count`.
#' @param readme_text Optional text for the Excel README worksheet.
#' @param verbose Logical; print progress, exclusions, and method messages.
#'
#' @return Invisibly returns a list with `estimates` (one row per block) and
#'   `cumulative_incidence` (one row per arm and block). When requested,
#'   per-block bootstrap matrices are returned as `crude_bootstrap`,
#'   `weighted_bootstrap`, or `matched_bootstrap` for the blocks present.
#'
#' @section Side effects:
#' With `out_xlsxpath`, creates its parent directory if needed and writes an
#' Excel workbook containing README, analysis-summary, and cumulative-incidence
#' sheets. With `verbose = TRUE`, prints progress and method messages.
#'
#' @seealso [estimate_hr()], [estimate_ir()], [create_ps_matched_cohort()]
#'
#' @export
#' @examples
#' df <- read.csv(system.file("extdata", "sample_data.csv",
#'                            package = "rwetools"))
#' res <- estimate_risk(
#'   in_df_crude = df, exposure_var = "exposure",
#'   outcome_var = "outcome", followuptime_var = "follow_up_days",
#'   time_unit = "days", risk_at_timepoint = 365,
#'   risk_per_individuals = 1000, verbose = FALSE
#' )
#' res$estimates
#'
#' \donttest{
#' # Aalen-Johansen competing-risk estimates with frozen-weight bootstrap CIs
#' res_aj <- estimate_risk(
#'   in_df_crude = df, exposure_var = "exposure",
#'   outcome_var = "outcome", followuptime_var = "follow_up_days",
#'   risk_at_timepoint = 365, risk_estimator = "AJ",
#'   if_aj_competing_event_var = "competing_event",
#'   if_bootstrap_count = 100, if_bootstrap_n_cores = 1,
#'   if_bootstrap_seed = 2026, verbose = FALSE
#' )
#' res_aj$estimates
#' }
estimate_risk <- function(
    in_df_crude = NULL,
    in_df_weight = NULL,
    in_df_match = NULL,
    out_xlsxpath = NULL,
    exposure_var = "exp",
    exp_value = 1,
    ref_value = 0,
    outcome_var = NULL,
    followuptime_var = NULL,
    time_unit = c("days", "months", "years"),
    risk_at_timepoint = 365,
    risk_per_individuals = 1000,
    confidence_level = 0.95,
    risk_estimator = c("KM", "AJ"),
    if_aj_competing_event_var = NULL,
    if_weight_weight_var = NULL,
    if_match_match_id = NULL,
    if_bootstrap_count = NULL,
    if_bootstrap_n_cores = NULL,
    if_bootstrap_seed = NULL,
    readme_text = NULL,
    verbose = TRUE) {

  # --- package checks ---
  if (!requireNamespace("survival", quietly = TRUE)) {
    stop("Package 'survival' is required but not installed.")
  }
  if (!is.null(out_xlsxpath) && !requireNamespace("openxlsx", quietly = TRUE)) {
    stop("Package 'openxlsx' is required for Excel output but not installed.")
  }

  if (verbose) message("\n========================================")
  if (verbose) message("RISK RATIO & RISK DIFFERENCE ESTIMATION")
  if (verbose) message("========================================")

  per_label <- format(risk_per_individuals, big.mark = ",", scientific = FALSE)

  #### 1. Input Validation & Setup ########################
  time_unit <- match.arg(time_unit)

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
  risk_estimator <- unname(risk_estimator_aliases[risk_estimator[1]])

  # if_aj_competing_event_var is AJ-only
  if (risk_estimator == "AJ") {
    if (is.null(if_aj_competing_event_var)) {
      stop("if_aj_competing_event_var must be specified when risk_estimator = 'AJ'.")
    }
  } else if (!is.null(if_aj_competing_event_var)) {
    if (verbose) message("if_aj_competing_event_var is ignored because risk_estimator = 'KM'.")
    if_aj_competing_event_var <- NULL
  }
  need_competing <- risk_estimator == "AJ"

  # --- bootstrap settings (opt-in; ALWAYS fixed-weight/frozen) ---
  do_boot <- !is.null(if_bootstrap_count)
  do_analytical <- TRUE
  if (!do_boot) {
    if (!is.null(if_bootstrap_n_cores) && verbose) {
      message("if_bootstrap_n_cores is ignored because if_bootstrap_count is NULL.")
    }
    if (!is.null(if_bootstrap_seed) && verbose) {
      message("if_bootstrap_seed is ignored because if_bootstrap_count is NULL.")
    }
  } else {
    if (!is.numeric(if_bootstrap_count) || length(if_bootstrap_count) != 1L ||
        if_bootstrap_count < 1) {
      stop("if_bootstrap_count must be a positive integer (>= 1) or NULL.")
    }
    if_bootstrap_count <- as.integer(if_bootstrap_count)
    if (is.null(if_bootstrap_n_cores)) {
      if_bootstrap_n_cores <- max(1L, parallel::detectCores() - 1L)
    }
    if_bootstrap_n_cores <- as.integer(max(1L, if_bootstrap_n_cores))
    message("Bootstrap method: fixed-weight (frozen) - weights are held at ",
            "their original values; PS-estimation uncertainty is NOT ",
            "propagated. Percentile CIs.")
  }

  # --- scalar validations ---
  if (is.null(outcome_var)) stop("outcome_var must be specified.")
  if (is.null(followuptime_var)) stop("followuptime_var must be specified.")
  if (!is.numeric(confidence_level) || length(confidence_level) != 1L ||
      confidence_level <= 0 || confidence_level >= 1) {
    stop("confidence_level must be a single numeric value between 0 and 1 (e.g., 0.95).")
  }
  if (!is.numeric(risk_at_timepoint) || length(risk_at_timepoint) != 1L ||
      risk_at_timepoint <= 0) {
    stop("risk_at_timepoint must be a single positive number.")
  }
  if (!is.numeric(risk_per_individuals) || length(risk_per_individuals) != 1L ||
      risk_per_individuals <= 0) {
    stop("risk_per_individuals must be a single positive number (e.g., 100, 1000, 10000, 100000).")
  }

  # --- block validation (shared with estimate_hr / estimate_ir) ---
  required_vars <- c(outcome_var, exposure_var, followuptime_var)
  if (need_competing) required_vars <- c(required_vars, if_aj_competing_event_var)

  vb <- validate_effect_blocks(
    in_df_crude = in_df_crude, in_df_weight = in_df_weight,
    in_df_match = in_df_match,
    if_weight_weight_var = if_weight_weight_var,
    if_match_match_id = if_match_match_id,
    required_vars = required_vars,
    verbose = verbose
  )
  blocks   <- vb$blocks
  weight_var <- vb$weight_var
  match_id   <- vb$match_id

  estimator_label <- ifelse(risk_estimator == "AJ", "Aalen-Johansen", "Kaplan-Meier")
  if (verbose) {
    message(sprintf("Outcome: %s", outcome_var))
    message(sprintf("Exposure: %s (Exposed=%s, Reference=%s)",
                    exposure_var, exp_value, ref_value))
    message(sprintf("Follow-up time: %s (%s)", followuptime_var, time_unit))
    message(sprintf("Timepoint for RR/RD: %s %s", risk_at_timepoint, time_unit))
    message(sprintf("Risk Estimator: %s (%s)", risk_estimator, estimator_label))
    if (need_competing) {
      message(sprintf("Competing Event Variable: %s", if_aj_competing_event_var))
    }
    message(sprintf("Risk expressed per: %s individuals", per_label))
    message(sprintf("Confidence Level: %.0f%%", confidence_level * 100))
    message("Analytical CIs (always): Risk/CIF = cloglog; RD = normal-Wald; RR = log-delta.")
    if (do_boot) {
      message(sprintf("Bootstrap Iterations: %d", if_bootstrap_count))
      message(sprintf("Parallel cores: %d", if_bootstrap_n_cores))
    }
    message(sprintf("Blocks: %s", paste(names(blocks), collapse = " + ")))
    if (!is.null(weight_var)) message(sprintf("Weight Variable: %s", weight_var))
    message("")
  }

  z_val <- stats::qnorm(1 - (1 - confidence_level) / 2)

  # Use the same exposure recode, complete-case removal, arm support, weight,
  # match-id and competing-event preprocessing as estimate_hr()/estimate_ir().
  blocks <- Map(function(df, bl) {
    prep_effect_block(
      df, bl, exposure_var = exposure_var, exp_value = exp_value,
      ref_value = ref_value, outcome_var = outcome_var,
      followuptime_var = followuptime_var, time_unit = NULL,
      weight_var = weight_var, match_id = match_id,
      competing_var = if (need_competing) if_aj_competing_event_var else NULL,
      check_competing = need_competing,
      competing_context = "Aalen-Johansen", verbose = verbose
    )
  }, blocks, names(blocks))

  results <- list()
  all_estimates <- data.frame()
  all_cuminc_data <- data.frame()

  #### 3. Run Analysis (common logic for all blocks) ########
  run_analysis <- function(data, analysis_label, wt_var = NULL, pair_col = NULL) {

    if (verbose) {
      message(sprintf("  [%s] Analysis", analysis_label))
      message(sprintf("  %s", paste(rep("-", nchar(analysis_label) + 12), collapse = "")))
    }

    # prep_effect_block() has already canonicalized the arms to 1/0 and
    # removed incomplete rows on exactly the columns used by this block.
    exp_val_int <- 1L
    ref_val_int <- 0L

    # --- Slim data: keep only needed columns for bootstrap efficiency ---
    keep_cols <- unique(c(followuptime_var, outcome_var, exposure_var, wt_var,
                          pair_col,
                          if (need_competing) if_aj_competing_event_var))
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

    # --- Common: count N and events per group ---
    is_ref <- data[[exposure_var]] == ref_val_int
    is_exp <- data[[exposure_var]] == exp_val_int
    in_tp  <- data[[followuptime_var]] <= risk_at_timepoint

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

    # --- point estimates, analytical CI, optional bootstrap ---
    if (verbose) message("  Calculating point estimates...")

    cuminc <- calc_km_cumulative_incidence(
      data = data,
      time_var = followuptime_var,
      event_var = outcome_var,
      exp_var = exposure_var,
      weight_var = wt_var,
      timepoint = risk_at_timepoint,
      per_n = risk_per_individuals,
      risk_estimator = risk_estimator,
      competing_event_var = if_aj_competing_event_var
    )

    rr_rd_vals <- calc_rr_rd(cuminc, exp_val_int, ref_val_int, risk_per_individuals)
    RR_point <- rr_rd_vals["RR"]
    RD_point <- rr_rd_vals["RD"]

    if (is.null(wt_var)) {
      n_ref <- cuminc$N_Individuals[cuminc$Exposure_Value == as.character(ref_val_int)]
      n_exp <- cuminc$N_Individuals[cuminc$Exposure_Value == as.character(exp_val_int)]
      events_ref <- cuminc$N_Events[cuminc$Exposure_Value == as.character(ref_val_int)]
      events_exp <- cuminc$N_Events[cuminc$Exposure_Value == as.character(exp_val_int)]
    }

    risk_ref_per_n <- rr_rd_vals["risk_ref"] * risk_per_individuals
    risk_exp_per_n <- rr_rd_vals["risk_exp"] * risk_per_individuals

    cuminc$Analysis <- analysis_label
    cuminc$Risk_Estimator <- risk_estimator
    cuminc_data <- cuminc

    # --- Analytical CI: RR log-delta, RD normal-Wald ---
    ana_ci <- calc_analytical_se_unstratified(
      cuminc_data = cuminc,
      exp_val = exp_val_int,
      ref_val = ref_val_int,
      per_n = risk_per_individuals,
      confidence_level = confidence_level
    )
    RR_se_ana  <- ana_ci$lnRR_se
    RR_lci_ana <- ana_ci$RR_lci
    RR_uci_ana <- ana_ci$RR_uci
    RD_se_ana  <- ana_ci$RD_se_per_n
    RD_lci_ana <- ana_ci$RD_lci_per_n
    RD_uci_ana <- ana_ci$RD_uci_per_n

    # --- Bootstrap CI (unstratified; matched = pair resample) ---
    boot_results <- NULL
    if (do_boot) {
      if (verbose) {
        message(sprintf("  Running %d bootstrap iterations (%d cores)...",
                        if_bootstrap_count, if_bootstrap_n_cores))
      }
      boot_results <- run_parallel_bootstrap(
        n_cores = if_bootstrap_n_cores,
        bootstrap_count = if_bootstrap_count,
        boot_fn = boot_rr_rd_single,
        seed = if_bootstrap_seed,
        export_varlist = c("calc_km_cumulative_incidence",
                           "calc_km_risk_single_group", "calc_rr_rd",
                           "pair_resample_index"),
        data = data,
        time_var = followuptime_var,
        event_var = outcome_var,
        exp_var = exposure_var,
        weight_var = wt_var,
        timepoint = risk_at_timepoint,
        exp_val = exp_val_int,
        ref_val = ref_val_int,
        per_n = risk_per_individuals,
        risk_estimator = risk_estimator,
        competing_event_var = if_aj_competing_event_var,
        match_ids = if (!is.null(pair_col)) data[[pair_col]] else NULL
      )
    }

    analysis_label_out <- analysis_label

    # ========================================================
    # Assemble estimates row (common for both branches)
    # ========================================================
    alpha <- 1 - confidence_level
    ci_probs <- c(alpha / 2, 1 - alpha / 2)

    if (do_boot && !is.null(boot_results)) {
      RR_lci_boot    <- stats::quantile(boot_results[, "RR"], probs = ci_probs[1], na.rm = TRUE)
      RR_uci_boot    <- stats::quantile(boot_results[, "RR"], probs = ci_probs[2], na.rm = TRUE)
      RR_se_boot     <- stats::sd(boot_results[, "RR"], na.rm = TRUE)
      RD_lci_boot    <- stats::quantile(boot_results[, "RD"], probs = ci_probs[1], na.rm = TRUE)
      RD_uci_boot    <- stats::quantile(boot_results[, "RD"], probs = ci_probs[2], na.rm = TRUE)
      RD_se_boot     <- stats::sd(boot_results[, "RD"], na.rm = TRUE)
      boot_n         <- if_bootstrap_count
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

    fmt_ci <- function(est, lcl, ucl, digits = 3) {
      if (is.na(est) || is.na(lcl) || is.na(ucl)) return(NA_character_)
      sprintf("%.*f (%.*f, %.*f)", digits, est, digits, lcl, digits, ucl)
    }

    est_row <- data.frame(
      Analysis         = analysis_label_out,
      Risk_Estimator   = risk_estimator,
      Timepoint        = risk_at_timepoint,
      Time_Unit        = time_unit,
      Per_Individuals  = risk_per_individuals,
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

    if (verbose) {
      if (do_boot && !is.na(RR_lci_boot)) {
        message(sprintf("\n  RR (bootstrap): %s", fmt_ci(RR_point, RR_lci_boot, RR_uci_boot, 2)))
        message(sprintf("  RD per %s (bootstrap): %s",
                        per_label, fmt_ci(RD_point, RD_lci_boot, RD_uci_boot, 2)))
      }
      if (!is.na(RR_lci_ana)) {
        message(sprintf("  RR (analytical): %s", fmt_ci(RR_point, RR_lci_ana, RR_uci_ana, 2)))
        message(sprintf("  RD per %s (analytical): %s",
                        per_label, fmt_ci(RD_point, RD_lci_ana, RD_uci_ana, 2)))
      }
    }

    # cloglog CI for the 0-1 Risk/CIF (matrix r22-r29; single site - the
    # crude/weighted/matched cases all funnel through here). Risk on {0,1}:
    # cloglog is undefined -> NA bounds + a bootstrap-recommending message.
    cl_ci <- cloglog_risk_ci(cuminc_data$Risk, cuminc_data$Risk_SE, z_val)
    cuminc_data$Risk_LCI <- cl_ci$lci
    cuminc_data$Risk_UCI <- cl_ci$uci
    if (any(!is.na(cuminc_data$Risk) &
              (cuminc_data$Risk <= 0 | cuminc_data$Risk >= 1))) {
      message("Risk estimate of exactly 0 or 1: the cloglog CI is undefined ",
              "(NA returned). Consider the bootstrap (if_bootstrap_count) ",
              "for an interval.")
    }
    cuminc_data$RiskperN_LCI <- cuminc_data$Risk_LCI * risk_per_individuals
    cuminc_data$RiskperN_UCI <- cuminc_data$Risk_UCI * risk_per_individuals

    cuminc_col_order <- c("Exposure_Value", "N_Individuals", "N_Events", "Timepoint",
                           "Risk", "Risk_LCI", "Risk_UCI", "Risk_SE", "Risk_Var",
                           "RiskperN", "RiskperN_LCI", "RiskperN_UCI",
                           "Risk_Estimator", "Analysis")
    cuminc_data <- cuminc_data[, cuminc_col_order, drop = FALSE]

    return(list(
      est_row = est_row,
      cuminc_data = cuminc_data,
      boot_results = boot_results
    ))
  }

  ##### 3a-c. Process the blocks ########################
  for (bl in names(blocks)) {
    if (verbose) {
      message(sprintf("Step 1: %s Analysis", bl))
      message("----------------------------")
    }
    pair_col <- if (bl == "Matched") match_id else NULL
    if (bl == "Matched" && is.null(match_id) && do_boot && verbose) {
      message("Matched block bootstrap: row resampling (no if_match_match_id).")
    }
    bl_result <- run_analysis(
      data = blocks[[bl]],
      analysis_label = bl,
      wt_var = if (bl == "Weighted") weight_var else NULL,
      pair_col = pair_col
    )
    all_estimates <- rbind(all_estimates, bl_result$est_row)
    all_cuminc_data <- rbind(all_cuminc_data, bl_result$cuminc_data)
    results[[paste0(tolower(bl), "_bootstrap")]] <- bl_result$boot_results
  }

  #### 5. Print Summary ################################
  if (verbose) {
    message("\n========================================")
    message("SUMMARY OF RESULTS")
    message(sprintf("(Risk Difference per %s)", per_label))
    message("========================================")

    display_cols <- c("Analysis", "Timepoint", "Time_Unit", "Per_Individuals")
    if (do_boot) display_cols <- c(display_cols, "RR_CI_Boot", "RDperN_CI_Boot")
    display_cols <- c(display_cols, "RR_CI_Analytical", "RDperN_CI_Analytical")

    message("Risk Ratio and Risk Difference Estimates:")
    print(all_estimates[, display_cols, drop = FALSE], row.names = FALSE)
    message("")
  }

  #### 6. Save Results ############################
  if (!is.null(out_xlsxpath)) {
    wb <- openxlsx::createWorkbook()

    openxlsx::addWorksheet(wb, "README")

    risk_method_desc <- ifelse(
      risk_estimator == "AJ",
      "Cumulative incidence via Aalen-Johansen estimator (competing risks).",
      "Cumulative incidence via Kaplan-Meier (1 - survival)."
    )

    method_detail <- paste0(
      risk_method_desc, " ",
      "Analytical SE: ",
      ifelse(risk_estimator == "AJ",
             "Counting-process variance (Aalen, 1978). ",
             "Greenwood-based Var(Risk). "),
      "Risk/CIF CI via complementary log-log (cloglog) transformation ",
      "(Kalbfleisch & Prentice 2002); Risk of exactly 0/1 -> NA. ",
      "RR CI via delta method on log scale. RD CI via normal approximation. ",
      ifelse("Matched" %in% names(blocks),
             paste0("Matched-block analytical variance treats rows as independent ",
                    "and does not model within-set correlation; when a match id ",
                    "and bootstrap are supplied, whole matched sets are resampled. "),
             ""),
      ifelse(do_boot,
             paste0("Bootstrap CIs via percentile method, fixed-weight (frozen); ",
                    ifelse(!is.null(match_id),
                           "matched sets resampled by match id; ",
                           ""),
                    sprintf("parallel RNG: L'Ecuyer-CMRG (seed=%s). ",
                            ifelse(is.null(if_bootstrap_seed),
                                   "user RNG state",
                                   as.character(if_bootstrap_seed)))),
             ""),
      sprintf("Risk expressed per %s individuals.", per_label)
    )

    analysis_type_desc <- paste0(
      "Risk Ratio & Risk Difference (", estimator_label,
      ifelse(risk_estimator == "AJ", ", Competing Risks", ""),
      ")"
    )

    readme_content <- data.frame(
      Item = c(
        "Analysis Type",
        "Date Generated",
        "Outcome Variable",
        "Exposure Variable",
        if (need_competing) "Competing Event Variable" else character(0),
        "Follow-up Time Variable",
        "Time Unit",
        "Timepoint for RR/RD",
        "Risk Estimator",
        "Risk Expressed Per",
        "Confidence Level",
        "CI Method",
        if (do_boot) c("Bootstrap Iterations", "Parallel Cores") else character(0),
        "Blocks",
        "Weight Variable",
        "Match Id",
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
        if (need_competing) if_aj_competing_event_var else character(0),
        followuptime_var,
        time_unit,
        as.character(risk_at_timepoint),
        paste0(risk_estimator, " (", estimator_label, ")"),
        paste0(per_label, " individuals"),
        paste0(confidence_level * 100, "%"),
        paste0("Analytical (always)",
               ifelse(do_boot, " + percentile bootstrap (fixed-weight)", "")),
        if (do_boot) c(as.character(if_bootstrap_count),
                       as.character(if_bootstrap_n_cores)) else character(0),
        paste(names(blocks), collapse = " + "),
        ifelse(is.null(weight_var), "N/A", weight_var),
        ifelse(is.null(match_id), "N/A", match_id),
        "",
        method_detail,
        "",
        ifelse(is.null(readme_text), "", readme_text)
      ),
      stringsAsFactors = FALSE
    )
    openxlsx::writeData(wb, "README", readme_content)
    openxlsx::setColWidths(wb, "README", cols = 1:2, widths = c(30, 100))

    openxlsx::addWorksheet(wb, "Analysis Summary")
    openxlsx::writeDataTable(wb, "Analysis Summary", all_estimates)
    openxlsx::setColWidths(wb, "Analysis Summary", cols = seq_len(ncol(all_estimates)),
                           widths = "auto")

    openxlsx::addWorksheet(wb, "Cumulative Incidence")
    openxlsx::writeDataTable(wb, "Cumulative Incidence", all_cuminc_data)
    openxlsx::setColWidths(wb, "Cumulative Incidence", cols = seq_len(ncol(all_cuminc_data)),
                           widths = "auto")

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

  results$estimates <- all_estimates
  results$cumulative_incidence <- all_cuminc_data

  return(invisible(results))
}
