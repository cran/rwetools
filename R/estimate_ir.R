# IR Incidence Rates, Rate Differences and Rate Ratios ###############
#' Estimate Marginal Incidence Rates, Rate Differences, and Rate Ratios
#'
#' Calculates incidence rates (IR), exposed-minus-reference rate differences
#' (IRD), and exposed-versus-reference rate ratios (IRR). A call can contain a
#' crude cohort plus either a weighted cohort or a matched cohort. The weighted
#' and matched blocks are mutually exclusive, and the crude block is optional.
#'
#' Direct standardization was removed in rwetools 0.5.0. This function has no
#' stratification argument: every reported rate measure is marginal over the
#' rows in its analysis block. To address a baseline stratifier, include it in
#' the propensity-score model; alternatively, run separate analyses within
#' levels and combine them explicitly in caller code under a prespecified
#' target-population rule.
#'
#' @section Exposure coding and analysis rows:
#' `exp_value` and `ref_value` are always recoded internally to 1 and 0,
#' including when the source column is already coded 0/1. They must be distinct
#' scalar, non-missing values. Rows matching neither value are reported and
#' removed; rows incomplete on variables used by this rate analysis are also
#' removed. Both arms must remain. For a logical exposure column, use
#' `exp_value = TRUE` and `ref_value = FALSE`, not numeric 1 and 0.
#'
#' Follow-up is converted to person-years using 365.25 days per year, 12 months
#' per year, or the supplied years directly. The IR and IRD are multiplied by
#' `ir_per_pyears`; the IRR is unitless.
#'
#' @section Analytical confidence intervals:
#' Methods are fixed by block; there are no CI-method arguments.
#' \itemize{
#'   \item Crude: the arm-specific IR uses a Garwood exact-Poisson interval;
#'     IRD uses the independent-Poisson variance with a normal-Wald interval;
#'     IRR uses the model SE from a marginal Poisson rate model.
#'   \item Weighted: IR, IRD, and IRR are contrasts from one saturated weighted
#'     quasi-Poisson cell model with a design-based sandwich covariance matrix.
#'   \item Matched with `if_match_match_id`: the same cell model declares the
#'     matched set as the sampling unit, so all three intervals are
#'     cluster-robust. The point estimates do not change. The intervals may be
#'     narrower or wider than independence intervals, depending on the sign of
#'     within-set covariance.
#'   \item Matched without `if_match_match_id`: rows are treated as independent
#'     and the crude analytical methods are used.
#' }
#'
#' @section Matched designs and matching weights:
#' The matched block is an unweighted matched analysis. If a `.match_weights`
#' column is present, every retained value must equal 1. With
#' `if_match_match_id`, every complete-case set must contain exactly one
#' exposed and one reference row; matching with replacement is rejected.
#'
#' For subclass, variable-ratio, fixed ratios greater than 1, or
#' with-replacement matching, pass the cohort as `in_df_weight` and set
#' `if_weight_weight_var = ".match_weights"`. This honours the weights and uses
#' weight-based robust inference, but it does not retain matched-set clustering.
#' In particular, variable-ratio matching would require weights and set
#' clustering together, which this block API cannot express, and the weighted
#' route is not the Abadie-Imbens variance for replacement matching. The
#' package does not produce a one-row-per-unit-per-pair replacement dataset.
#'
#' @section Migration from rwetools 0.4.0:
#' `estimate_hr_ir()` was split into [estimate_hr()] and `estimate_ir()`.
#' Hazard-ratio arguments moved to [estimate_hr()]. `stratification_var` was
#' removed from the rate API together with direct standardization. Fine-Gray
#' competing-event missingness no longer narrows the incidence-rate rows.
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
#' @param time_unit Unit of `followuptime_var`: `"days"`, `"months"`, or
#'   `"years"`.
#' @param ir_per_pyears Reporting multiplier for IR and IRD. Must be one of 1,
#'   100, 1000, 10000, or 100000.
#' @param confidence_level Single number strictly between 0 and 1. The default
#'   is 0.95.
#' @param if_weight_weight_var Character name of the analysis-weight column in
#'   `in_df_weight`. Required for a weighted block and ignored with a message
#'   otherwise.
#' @param if_match_match_id Character name of the set id in `in_df_match`.
#'   When supplied, it defines the analytical cluster and pair-level bootstrap
#'   unit. Without it, matched rows are treated as independent analytically and
#'   resampled by row. Ignored with a message without a matched block.
#' @param if_bootstrap_count Positive number of percentile-bootstrap
#'   replicates, or `NULL` to omit bootstrap CIs. Weights are frozen at their
#'   supplied values, so PS-estimation uncertainty is not propagated.
#' @param if_bootstrap_n_cores Number of bootstrap workers, or `NULL` to use
#'   all detected cores minus one. Ignored without `if_bootstrap_count`.
#' @param if_bootstrap_seed Integer RNG seed, or `NULL`. It makes this
#'   function reproducible, but does not coordinate draws with [estimate_hr()].
#'   Ignored without `if_bootstrap_count`.
#' @param readme_text Optional text for a README worksheet when Excel output is
#'   requested.
#' @param verbose Logical; print progress, exclusions, and method messages.
#'
#' @return Invisibly returns a list with `incidence_rates`,
#'   `incidence_rate_ratios`, and `ir_per_pyears`. Per-block bootstrap draw
#'   matrices are returned in `bootstrap` when requested.
#'
#' @section Side effects:
#' With `out_xlsxpath`, creates its parent directory if needed and writes an
#' Excel workbook. With `verbose = TRUE`, prints progress and method messages.
#'
#' @seealso [estimate_hr()], [estimate_risk()], [create_ps_matched_cohort()]
#'
#' @export
#' @examples
#' df <- read.csv(system.file("extdata", "sample_data.csv",
#'                            package = "rwetools"))
#' res <- estimate_ir(
#'   in_df_crude      = df,
#'   exposure_var     = "exposure",
#'   outcome_var      = "outcome",
#'   followuptime_var = "follow_up_days",
#'   time_unit        = "days",
#'   verbose          = FALSE
#' )
#' res$incidence_rates
#' res$incidence_rate_ratios
#'
#' \donttest{
#' # Crude and weighted blocks with frozen-weight bootstrap CIs
#' df_ps <- estimate_ps(
#'   in_df = df, exposure_var = "exposure",
#'   class_vars = c("cat1", "cat2", "cat3", "cat4"),
#'   cont_vars = c("cont1", "cont2", "cont3"), verbose = FALSE
#' )
#' df_wt <- create_iptw(
#'   in_df = df_ps, exposure_var = "exposure", ps_var = "ps",
#'   weight_var = "iptw", verbose = FALSE
#' )
#' res_boot <- estimate_ir(
#'   in_df_crude = df, in_df_weight = df_wt,
#'   if_weight_weight_var = "iptw", exposure_var = "exposure",
#'   outcome_var = "outcome", followuptime_var = "follow_up_days",
#'   if_bootstrap_count = 200, if_bootstrap_n_cores = 1,
#'   if_bootstrap_seed = 2026, verbose = FALSE
#' )
#' res_boot$incidence_rates
#' }
estimate_ir <- function(
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
    ir_per_pyears = 1000,
    confidence_level = 0.95,
    if_weight_weight_var = NULL,
    if_match_match_id = NULL,
    if_bootstrap_count = NULL,
    if_bootstrap_n_cores = NULL,
    if_bootstrap_seed = NULL,
    readme_text = NULL,
    verbose = TRUE) {

  # --- package checks ---
  if (!is.null(out_xlsxpath) && !requireNamespace("openxlsx", quietly = TRUE)) {
    stop("Package 'openxlsx' is required for Excel output but not installed.")
  }

  if (verbose) message("\n========================================")
  if (verbose) message("INCIDENCE RATE & RATE-RATIO ESTIMATION")
  if (verbose) message("========================================")

  # --- scalar validations ---
  if (is.null(outcome_var)) stop("outcome_var must be specified")
  if (is.null(followuptime_var)) stop("followuptime_var must be specified")
  allowed_multipliers <- c(1, 100, 1000, 10000, 100000)
  if (!ir_per_pyears %in% allowed_multipliers) {
    stop(paste0("ir_per_pyears must be one of: ",
                paste(allowed_multipliers, collapse = ", "),
                ". Got: ", ir_per_pyears))
  }
  if (!is.numeric(confidence_level) || length(confidence_level) != 1 ||
      confidence_level <= 0 || confidence_level >= 1) {
    stop("confidence_level must be a single numeric value between 0 and 1 (e.g., 0.95)")
  }
  time_unit <- match.arg(time_unit)

  # --- bootstrap settings (opt-in; ALWAYS fixed-weight/frozen) ---
  do_boot <- !is.null(if_bootstrap_count)
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
  }

  # --- block validation (shared with estimate_hr / estimate_risk) ---
  required_vars <- c(outcome_var, exposure_var, followuptime_var)

  vb <- validate_effect_blocks(
    in_df_crude = in_df_crude, in_df_weight = in_df_weight,
    in_df_match = in_df_match,
    if_weight_weight_var = if_weight_weight_var,
    if_match_match_id = if_match_match_id,
    required_vars = required_vars,
    verbose = verbose
  )
  blocks     <- vb$blocks
  weight_var <- vb$weight_var
  match_id   <- vb$match_id

  ir_col_base <- "IR_per_Npy"
  z_score <- stats::qnorm(1 - (1 - confidence_level) / 2)

  if (verbose) {
    message(sprintf("Outcome: %s", outcome_var))
    message(sprintf("Follow-up time: %s (%s)", followuptime_var, time_unit))
    message(sprintf("Blocks: %s", paste(names(blocks), collapse = " + ")))
    message(sprintf("Confidence level: %.0f%%", confidence_level * 100))
    message("Analytical CIs (fixed per block):")
    message("  Crude: IR = exact Poisson (Garwood), IRD = Wald, IRR = Poisson model SE.")
    message("  Weighted: IR/IRD/IRR = design-based robust sandwich (joint delta).")
    if (!is.null(match_id)) {
      message("  Matched: IR/IRD/IRR = cluster-robust sandwich on the match id; ",
              "intervals may be narrower or wider than independence intervals.")
    } else if ("Matched" %in% names(blocks)) {
      message("  Matched (no if_match_match_id): rows treated as independent -- ",
              "IR = Garwood, IRD = Wald, IRR = Poisson model SE.")
    }
    message("")
  }

  # --- per-block preparation ---
  blocks <- Map(function(df, bl) {
    prep_effect_block(
      df, bl, exposure_var = exposure_var, exp_value = exp_value,
      ref_value = ref_value, outcome_var = outcome_var,
      followuptime_var = followuptime_var, time_unit = time_unit,
      weight_var = weight_var, match_id = match_id, verbose = verbose
    )
  }, blocks, names(blocks))

  results <- list()

  # STEP 1: incidence rates ###########################################
  if (verbose) message("Step 1: Incidence rates & rate differences")
  incidence_rates <- NULL

  for (bl in names(blocks)) {
    df <- blocks[[bl]]
    wv <- if (bl == "Weighted") weight_var else NULL
    cv <- if (bl == "Matched") match_id else NULL
    block_summary <- if (is.null(wv) && is.null(cv)) {
      # crude, or matched without a match id: independent-Poisson intervals
      calc_crude_ir_ird(df, exposure_var, outcome_var, "person_years",
                        confidence_level, z_score, bl, ir_per_pyears, ir_col_base)
    } else {
      # weighted (analysis weights) or matched with a match id (clustered):
      # design-based sandwich intervals (backlog Item 2)
      calc_crude_ir_ird_design(df, exposure_var, outcome_var, "person_years",
                               wt_var = wv, cluster_var = cv,
                               confidence_level, z_score, bl,
                               ir_per_pyears, ir_col_base)
    }
    if (is.null(incidence_rates)) {
      incidence_rates <- block_summary
    } else {
      add_cols <- setdiff(names(block_summary), exposure_var)
      incidence_rates <- cbind(incidence_rates, block_summary[, add_cols, drop = FALSE])
    }
  }

  incidence_rates$Exposure_Group <- ifelse(
    incidence_rates[[exposure_var]] == 99, "Total",
    ifelse(incidence_rates[[exposure_var]] == 1, "Exposed",
           ifelse(incidence_rates[[exposure_var]] == 0, "Reference",
                  as.character(incidence_rates[[exposure_var]])))
  )
  col_order <- c("Exposure_Group",
                 setdiff(names(incidence_rates), c("Exposure_Group", exposure_var)))
  incidence_rates <- incidence_rates[, col_order]
  if (verbose) print(incidence_rates, digits = 3)

  # STEP 2: incidence rate ratios (always on) #########################
  if (verbose) message("\nStep 2: Incidence rate ratios (marginal, always reported)")
  incidence_rate_ratios <- NULL
  for (bl in names(blocks)) {
    df <- blocks[[bl]]
    wv <- if (bl == "Weighted") weight_var else NULL
    if (any(df$person_years <= 0, na.rm = TRUE)) {
      stop(sprintf("IRR rate model requires positive person-time; non-positive person_years in the %s block.", bl))
    }
    irr_row <- fit_irr_model(df, outcome_var, exposure_var, "person_years",
                             weight_var = wv,
                             cluster_var = if (bl == "Matched") match_id else NULL,
                             confidence_level = confidence_level,
                             base_label = bl)
    incidence_rate_ratios <- rbind(incidence_rate_ratios, irr_row)
  }
  if (verbose) print(incidence_rate_ratios, digits = 3)

  # STEP 3: bootstrap layer (opt-in, fixed-weight) ####################
  if (do_boot) {
    if (verbose) {
      message(sprintf("\nStep 3: Bootstrap (%d replicates, %d core(s))",
                      if_bootstrap_count, if_bootstrap_n_cores))
    }
    message("Bootstrap method: fixed-weight (frozen) - weights are held at their ",
            "original values; PS-estimation uncertainty is NOT propagated. ",
            "Percentile CIs.")

    boot_stats <- list()
    for (bl in names(blocks)) {
      df <- blocks[[bl]]
      wv <- if (bl == "Weighted") weight_var else NULL
      use_pairs <- (bl == "Matched") && !is.null(match_id)
      if (bl == "Matched" && is.null(match_id) && verbose) {
        message("Matched block bootstrap: row resampling (no if_match_match_id).")
      }

      boot_fn <- local({
        df_b <- df; wv_b <- wv; use_pairs_b <- use_pairs
        match_col <- if (use_pairs) df[[match_id]] else NULL
        exposure_b <- exposure_var; outcome_b <- outcome_var
        mult_b <- ir_per_pyears
        function() {
          idx <- if (use_pairs_b) {
            pair_resample_index(match_col)
          } else {
            sample.int(nrow(df_b), nrow(df_b), replace = TRUE)
          }
          d <- df_b[idx, , drop = FALSE]
          out <- c(IR1 = NA_real_, IR0 = NA_real_,
                   IRD = NA_real_, IRR = NA_real_)
          tryCatch({
            w <- if (is.null(wv_b)) rep(1, nrow(d)) else d[[wv_b]]
            arm <- function(a) {
              sel <- d[[exposure_b]] == a
              c(ev = sum(d[[outcome_b]][sel] * w[sel]),
                py = sum(d$person_years[sel] * w[sel]))
            }
            a1 <- arm(1); a0 <- arm(0)
            ir1 <- a1["ev"] / a1["py"]; ir0 <- a0["ev"] / a0["py"]
            out["IR1"] <- mult_b * ir1
            out["IR0"] <- mult_b * ir0
            out["IRD"] <- mult_b * (ir1 - ir0)
            out["IRR"] <- ir1 / ir0
            out
          }, error = function(e) out)
        }
      })

      boot_stats[[bl]] <- run_parallel_bootstrap(
        n_cores = if_bootstrap_n_cores,
        bootstrap_count = if_bootstrap_count,
        boot_fn = boot_fn,
        seed = if_bootstrap_seed,
        export_varlist = c("pair_resample_index"),
        result_colnames = c("IR1", "IR0", "IRD", "IRR")
      )
    }
    results$bootstrap <- boot_stats

    # summarize percentile CIs / SEs into the output tables
    alpha <- 1 - confidence_level
    qs <- function(x) c(lci = unname(stats::quantile(x, alpha / 2, na.rm = TRUE)),
                        uci = unname(stats::quantile(x, 1 - alpha / 2, na.rm = TRUE)),
                        se  = stats::sd(x, na.rm = TRUE))
    # pre-initialize boot columns (assigning into a subset of a nonexistent
    # data.frame column is not safe)
    incidence_rate_ratios$IRR_BootLCI <- NA_real_
    incidence_rate_ratios$IRR_BootUCI <- NA_real_
    incidence_rate_ratios$IRR_BootSE  <- NA_real_
    for (bl in names(boot_stats)) {
      bm <- boot_stats[[bl]]
      s_ir1 <- qs(bm[, "IR1"]); s_ir0 <- qs(bm[, "IR0"])
      s_ird <- qs(bm[, "IRD"]); s_irr <- qs(bm[, "IRR"])

      grp <- incidence_rates$Exposure_Group
      incidence_rates[[paste0("IR_BootLCI_", bl)]] <-
        ifelse(grp == "Exposed", s_ir1["lci"],
               ifelse(grp == "Reference", s_ir0["lci"], NA_real_))
      incidence_rates[[paste0("IR_BootUCI_", bl)]] <-
        ifelse(grp == "Exposed", s_ir1["uci"],
               ifelse(grp == "Reference", s_ir0["uci"], NA_real_))
      incidence_rates[[paste0("IR_BootSE_", bl)]] <-
        ifelse(grp == "Exposed", s_ir1["se"],
               ifelse(grp == "Reference", s_ir0["se"], NA_real_))
      incidence_rates[[paste0("IRD_BootLCI_", bl)]] <- s_ird["lci"]
      incidence_rates[[paste0("IRD_BootUCI_", bl)]] <- s_ird["uci"]
      incidence_rates[[paste0("IRD_BootSE_",  bl)]] <- s_ird["se"]

      irr_sel <- incidence_rate_ratios$Analysis == bl
      incidence_rate_ratios$IRR_BootLCI[irr_sel] <- s_irr["lci"]
      incidence_rate_ratios$IRR_BootUCI[irr_sel] <- s_irr["uci"]
      incidence_rate_ratios$IRR_BootSE[irr_sel]  <- s_irr["se"]
    }
  }

  # STEP 4: Excel output ##############################################
  if (!is.null(out_xlsxpath)) {
    wb <- openxlsx::createWorkbook()
    if (!is.null(readme_text)) add_readme_sheet(wb, readme_text, verbose = verbose)

    openxlsx::addWorksheet(wb, "Analysis Summary")
    summary_info <- data.frame(
      Parameter = c("Outcome", "Exposure", "Follow-up Time Variable",
                    "Time Unit", "IR per Person-Years Multiplier",
                    "Confidence Level", "Blocks", "Weight Variable",
                    "Match Id", "Analytical CI methods", "Bootstrap"),
      Value = c(outcome_var, exposure_var, followuptime_var,
                time_unit, format(ir_per_pyears, big.mark = ","),
                paste0(confidence_level * 100, "%"),
                paste(names(blocks), collapse = " + "),
                ifelse(is.null(weight_var), "N/A", weight_var),
                ifelse(is.null(match_id), "N/A", match_id),
                paste0("Crude: IR=Garwood, IRD=Wald, IRR=Poisson model SE. ",
                       "Weighted: IR/IRD/IRR=design-based robust sandwich. ",
                       "Matched: cluster-robust sandwich on the match id when ",
                       "if_match_match_id is supplied (intervals may be narrower ",
                       "or wider), otherwise the ",
                       "independent-row methods."),
                ifelse(do_boot,
                       sprintf("fixed-weight (frozen), %d replicates, seed=%s",
                               if_bootstrap_count,
                               ifelse(is.null(if_bootstrap_seed), "user RNG",
                                      as.character(if_bootstrap_seed))),
                       "off"))
    )
    openxlsx::writeDataTable(wb, "Analysis Summary", summary_info, startRow = 1)

    openxlsx::addWorksheet(wb, "Incidence Rates")
    openxlsx::writeDataTable(wb, "Incidence Rates", incidence_rates, startRow = 1)

    openxlsx::addWorksheet(wb, "Incidence Rate Ratios")
    openxlsx::writeDataTable(wb, "Incidence Rate Ratios", incidence_rate_ratios, startRow = 1)

    out_dir <- dirname(out_xlsxpath)
    if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
    openxlsx::saveWorkbook(wb, out_xlsxpath, overwrite = TRUE)
    if (verbose) message(sprintf("\u2713 Results saved to: %s", out_xlsxpath))
  }

  # --- assemble results ---
  results$incidence_rates <- incidence_rates
  results$incidence_rate_ratios <- incidence_rate_ratios
  results$ir_per_pyears <- ir_per_pyears

  if (verbose) message("\nANALYSIS COMPLETE")
  return(invisible(results))
}
