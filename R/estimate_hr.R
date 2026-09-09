# HR Hazard Ratios ###################################################
#' Estimate a Cox Hazard Ratio or Fine-Gray Subdistribution Hazard Ratio
#'
#' Fits either a Cox proportional-hazards model or a Fine-Gray model for a
#' binary exposure. A call can contain a crude cohort plus either a weighted
#' cohort or a matched cohort. The weighted and matched blocks are mutually
#' exclusive, and the crude block is optional.
#'
#' @section Stratification is not standardization:
#' `strata_var` adds `survival::strata()` to the hazard model. It allows a
#' separate baseline hazard in each level while estimating one conditional
#' common exposure (s)HR. This is not a marginal standardized rate or risk.
#' Every retained level must contain both exposure arms; the function stops
#' otherwise because a single-arm level contributes no within-level comparison
#' to the partial likelihood.
#'
#' Do not automatically pass the `strata` column produced by
#' [create_ps_fs_weights()] as `strata_var`. Fine-stratification weights and a
#' stratified hazard model are different adjustments; using both is a separate
#' analysis choice that changes the fitted hazard model.
#'
#' @section Exposure coding and analysis rows:
#' `exp_value` and `ref_value` are always recoded internally to 1 and 0,
#' including when the source column is already coded 0/1. They must be distinct
#' scalar, non-missing values. Rows matching neither value are reported and
#' removed; rows incomplete on variables used by the requested model are also
#' removed. Both arms must remain. For a logical exposure column, use
#' `exp_value = TRUE` and `ref_value = FALSE`, not numeric 1 and 0.
#'
#' @section Analytical confidence intervals:
#' Methods are fixed by block; there are no CI-method arguments.
#' \itemize{
#'   \item Cox: crude uses the model SE, weighted uses a robust sandwich SE,
#'     and matched uses a robust SE clustered on `if_match_match_id` when the
#'     id is supplied.
#'   \item Fine-Gray: the expanded-data Cox fit uses robust inference clustered
#'     on subject for crude and weighted analyses, or on match id for a matched
#'     analysis. Weighted expansion rows carry the product of IPCW and the
#'     supplied analysis weight.
#' }
#'
#' @section Matched designs and matching weights:
#' The matched block is an unweighted matched analysis. If a `.match_weights`
#' column is present, every retained value must equal 1. With
#' `if_match_match_id`, every complete-case set must contain exactly one
#' exposed and one reference row; matching with replacement is rejected because
#' a reused unit belongs to multiple sets. Without a match id, inference is not
#' clustered by matched set and the bootstrap resamples rows.
#'
#' For subclass, variable-ratio, fixed ratios greater than 1, or
#' with-replacement matching, pass the cohort as `in_df_weight` and set
#' `if_weight_weight_var = ".match_weights"` to honour the matching weights.
#' That route provides a weight-based robust SE but does not use matched-set
#' clustering. It therefore cannot express weights and set clustering
#' together, and it is not the Abadie-Imbens variance for matching with
#' replacement. The package also does not produce the one-row-per-unit-per-pair
#' representation required for design-aware replacement inference.
#'
#' @section Migration from rwetools 0.4.0:
#' `estimate_hr_ir()` was split into `estimate_hr()` and [estimate_ir()]. The
#' old `stratification_var` argument is now `strata_var` here only. `time_unit`
#' and `ir_per_pyears` belong to [estimate_ir()] and are not arguments to this
#' function. With `hr_model = "Fine-Gray"`, missing values in the competing
#' event variable affect this hazard analysis but no longer remove rows from a
#' separate incidence-rate analysis.
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
#' @param followuptime_var Character name of the follow-up time column. The
#'   Cox/Fine-Gray coefficient is invariant to a constant change of time unit.
#' @param confidence_level Single number strictly between 0 and 1. The default
#'   is 0.95.
#' @param strata_var Character name of a model-stratification column, or
#'   `NULL`. See **Stratification is not standardization**.
#' @param hr_model `"Cox"` for a cause-specific HR or `"Fine-Gray"` for a
#'   subdistribution HR.
#' @param if_fg_competing_event_var Character name of a binary competing-event
#'   indicator, coded 1 for a competing event. Required for Fine-Gray and
#'   ignored with a message for Cox.
#' @param if_weight_weight_var Character name of the analysis-weight column in
#'   `in_df_weight`. Required for a weighted block and ignored with a message
#'   otherwise.
#' @param if_match_match_id Character name of the set id in `in_df_match`.
#'   When supplied, it defines the robust-SE cluster and pair-level bootstrap
#'   unit. Ignored with a message without a matched block.
#' @param if_bootstrap_count Positive number of percentile-bootstrap
#'   replicates, or `NULL` to omit bootstrap CIs. Weights are frozen at their
#'   supplied values, so PS-estimation uncertainty is not propagated.
#' @param if_bootstrap_n_cores Number of bootstrap workers, or `NULL` to use
#'   all detected cores minus one. Ignored without `if_bootstrap_count`.
#' @param if_bootstrap_seed Integer RNG seed, or `NULL`. It makes this
#'   function reproducible, but does not coordinate draws with [estimate_ir()].
#'   Ignored without `if_bootstrap_count`.
#' @param readme_text Optional text for a README worksheet when Excel output is
#'   requested.
#' @param verbose Logical; print progress, exclusions, and method messages.
#'
#' @return Invisibly returns a list. Cox analysis supplies `hazard_ratios` and
#'   per-block `models`; Fine-Gray supplies `subdist_hazard`. Bootstrap draws
#'   are in `bootstrap` when requested, and `strata_var` is returned when used.
#'
#' @section Side effects:
#' With `out_xlsxpath`, creates its parent directory if needed and writes an
#' Excel workbook. With `verbose = TRUE`, prints progress and method messages.
#'
#' @seealso [estimate_ir()], [estimate_risk()], [create_ps_matched_cohort()]
#'
#' @export
#' @examples
#' df <- read.csv(system.file("extdata", "sample_data.csv",
#'                            package = "rwetools"))
#' res <- estimate_hr(
#'   in_df_crude      = df,
#'   exposure_var     = "exposure",
#'   outcome_var      = "outcome",
#'   followuptime_var = "follow_up_days",
#'   verbose          = FALSE
#' )
#' res$hazard_ratios
#'
#' \donttest{
#' # Crude and weighted Fine-Gray blocks
#' df_ps <- estimate_ps(
#'   in_df = df, exposure_var = "exposure",
#'   class_vars = c("cat1", "cat2", "cat3", "cat4"),
#'   cont_vars = c("cont1", "cont2", "cont3"), verbose = FALSE
#' )
#' df_wt <- create_iptw(
#'   in_df = df_ps, exposure_var = "exposure", ps_var = "ps",
#'   weight_var = "iptw", verbose = FALSE
#' )
#' res_fg <- estimate_hr(
#'   in_df_crude = df, in_df_weight = df_wt,
#'   if_weight_weight_var = "iptw", exposure_var = "exposure",
#'   outcome_var = "outcome", followuptime_var = "follow_up_days",
#'   hr_model = "Fine-Gray",
#'   if_fg_competing_event_var = "competing_event", verbose = FALSE
#' )
#' res_fg$subdist_hazard
#' }
estimate_hr <- function(
    in_df_crude = NULL,
    in_df_weight = NULL,
    in_df_match = NULL,
    out_xlsxpath = NULL,
    exposure_var = "exp",
    exp_value = 1,
    ref_value = 0,
    outcome_var = NULL,
    followuptime_var = NULL,
    confidence_level = 0.95,
    strata_var = NULL,
    hr_model = c("Cox", "Fine-Gray"),
    if_fg_competing_event_var = NULL,
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
  if (verbose) message("HAZARD RATIO ESTIMATION")
  if (verbose) message("========================================")

  # --- scalar validations ---
  if (is.null(outcome_var)) stop("outcome_var must be specified")
  if (is.null(followuptime_var)) stop("followuptime_var must be specified")
  if (!is.numeric(confidence_level) || length(confidence_level) != 1 ||
      confidence_level <= 0 || confidence_level >= 1) {
    stop("confidence_level must be a single numeric value between 0 and 1 (e.g., 0.95)")
  }
  hr_model <- match.arg(hr_model)

  # --- Fine-Gray settings ---
  if (hr_model == "Fine-Gray") {
    if (is.null(if_fg_competing_event_var)) {
      stop("if_fg_competing_event_var must be specified when hr_model = 'Fine-Gray'.")
    }
  } else if (!is.null(if_fg_competing_event_var)) {
    if (verbose) message("if_fg_competing_event_var is ignored because hr_model = 'Cox'.")
    if_fg_competing_event_var <- NULL
  }

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

  # --- block validation (shared with estimate_ir / estimate_risk) ---
  required_vars <- c(outcome_var, exposure_var, followuptime_var)
  if (!is.null(strata_var)) required_vars <- c(required_vars, strata_var)
  if (!is.null(if_fg_competing_event_var)) required_vars <- c(required_vars, if_fg_competing_event_var)

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
  do_stratify <- !is.null(strata_var)

  if (verbose) {
    message(sprintf("Outcome: %s", outcome_var))
    message(sprintf("Follow-up time: %s", followuptime_var))
    message(sprintf("Hazard model: %s", hr_model))
    message(sprintf("Blocks: %s", paste(names(blocks), collapse = " + ")))
    if (do_stratify) {
      message(sprintf("Strata: %s (survival::strata(); conditional common %s)",
                      strata_var,
                      ifelse(hr_model == "Cox", "HR", "sHR")))
    }
    message(sprintf("Confidence level: %.0f%%", confidence_level * 100))
    message(sprintf("Analytical CIs (fixed per block): %s Crude = model SE, Weighted = robust, Matched = robust%s.",
                    ifelse(hr_model == "Cox", "HR", "Fine-Gray sHR"),
                    ifelse(!is.null(match_id), " + cluster(match id)", "")))
    message("")
  }

  # --- per-block preparation ---
  match_id_missing_msged <- FALSE
  blocks <- Map(function(df, bl) {
    prep_effect_block(
      df, bl, exposure_var = exposure_var, exp_value = exp_value,
      ref_value = ref_value, outcome_var = outcome_var,
      followuptime_var = followuptime_var, time_unit = NULL,
      weight_var = weight_var, match_id = match_id,
      strata_var = strata_var,
      competing_var = if_fg_competing_event_var,
      check_competing = (hr_model == "Fine-Gray"), verbose = verbose
    )
  }, blocks, names(blocks))

  # every strata() level needs both arms, or it silently leaves the model
  if (do_stratify) {
    for (bl in names(blocks)) {
      check_arm_stratum_support(blocks[[bl]], exposure_var, strata_var, bl)
    }
  }

  results <- list()
  results$models <- list()

  # STEP 1: hazard model ##############################################
  if (verbose) message(sprintf("Step 1: %s model", hr_model))

  fit_cox_block <- function(df, block) {
    lhs <- paste0("survival::Surv(", followuptime_var, ", ", outcome_var, ")")
    rhs <- exposure_var
    if (do_stratify) rhs <- paste0(rhs, " + survival::strata(", strata_var, ")")

    if (block == "Crude") {
      fml <- stats::as.formula(paste(lhs, "~", rhs))
      model <- survival::coxph(fml, data = df, robust = FALSE)
      se_col <- "se(coef)"
      label <- "Crude"
    } else if (block == "Weighted") {
      fml <- stats::as.formula(paste(lhs, "~", rhs))
      model <- survival::coxph(fml, data = df, weights = df[[weight_var]],
                               robust = TRUE)
      se_col <- "robust se"
      label <- "Weighted (robust SE)"
    } else { # Matched
      rhs_fit <- rhs
      if (!is.null(match_id)) {
        rhs_fit <- paste0(rhs_fit, " + cluster(", match_id, ")")
        label <- "Matched (robust SE, match-id cluster)"
      } else {
        if (verbose && !match_id_missing_msged) {
          message("Matched block without if_match_match_id: robust SE without pair clustering.")
          match_id_missing_msged <<- TRUE
        }
        label <- "Matched (robust SE)"
      }
      fml <- stats::as.formula(paste(lhs, "~", rhs_fit))
      model <- survival::coxph(fml, data = df, robust = TRUE)
      se_col <- "robust se"
    }
    if (do_stratify) label <- paste0(label, " (stratified by ", strata_var, ")")

    s    <- summary(model, conf.int = confidence_level)
    cf   <- s$coefficients
    ci   <- s$conf.int
    ridx <- which(rownames(cf) == exposure_var)
    row_out <- data.frame(
      Analysis = label,
      HR_CI    = sprintf("%.2f (%.2f, %.2f)", unname(ci[ridx, "exp(coef)"]),
                         unname(ci[ridx, 3]), unname(ci[ridx, 4])),
      HR       = unname(ci[ridx, "exp(coef)"]),
      HR_LCI   = unname(ci[ridx, 3]),
      HR_UCI   = unname(ci[ridx, 4]),
      lnHR     = log(unname(ci[ridx, "exp(coef)"])),
      lnHR_SE  = unname(cf[ridx, se_col]),
      Pvalue   = unname(cf[ridx, "Pr(>|z|)"]),
      stringsAsFactors = FALSE, row.names = NULL
    )
    list(row = row_out, model = model)
  }

  hazard_ratios <- NULL
  subdist_hazard <- NULL
  for (bl in names(blocks)) {
    df <- blocks[[bl]]
    if (hr_model == "Cox") {
      fit <- fit_cox_block(df, bl)
      hazard_ratios <- rbind(hazard_ratios, fit$row)
      results$models[[bl]] <- fit$model
    } else {
      lab <- bl
      if (do_stratify) lab <- paste0(lab, " (stratified by ", strata_var, ")")
      fg_row <- fit_finegray_shr(
        data = df, time_var = followuptime_var, event_var = outcome_var,
        competing_var = if_fg_competing_event_var, exp_var = exposure_var,
        weight_var = if (bl == "Weighted") weight_var else NULL,
        strat_var = strata_var,
        cluster_id_var = if (bl == "Matched") match_id else NULL,
        confidence_level = confidence_level, base_label = lab
      )
      subdist_hazard <- rbind(subdist_hazard, fg_row)
    }
  }
  if (verbose && !is.null(hazard_ratios))  print(hazard_ratios,  digits = 3)
  if (verbose && !is.null(subdist_hazard)) print(subdist_hazard, digits = 3)

  # STEP 2: bootstrap layer (opt-in, fixed-weight) ####################
  if (do_boot) {
    if (verbose) {
      message(sprintf("\nStep 2: Bootstrap (%d replicates, %d core(s))",
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
        ftime_b <- followuptime_var; strat_b <- strata_var
        do_strat_b <- do_stratify
        hr_model_b <- hr_model; fg_var_b <- if_fg_competing_event_var
        function() {
          idx <- if (use_pairs_b) {
            pair_resample_index(match_col)
          } else {
            sample.int(nrow(df_b), nrow(df_b), replace = TRUE)
          }
          d <- df_b[idx, , drop = FALSE]
          out <- c(HR = NA_real_)
          tryCatch({
            lhs <- paste0("survival::Surv(", ftime_b, ", ", outcome_b, ")")
            rhs <- exposure_b
            if (do_strat_b) rhs <- paste0(rhs, " + survival::strata(", strat_b, ")")
            fml <- stats::as.formula(paste(lhs, "~", rhs))
            environment(fml) <- globalenv()
            if (hr_model_b == "Cox") {
              m <- if (is.null(wv_b)) {
                survival::coxph(fml, data = d)
              } else {
                survival::coxph(fml, data = d, weights = d[[wv_b]])
              }
              out["HR"] <- unname(exp(stats::coef(m)[exposure_b]))
            } else {
              fg_row <- fit_finegray_shr(
                data = d, time_var = ftime_b, event_var = outcome_b,
                competing_var = fg_var_b, exp_var = exposure_b,
                weight_var = wv_b, strat_var = strat_b,
                confidence_level = 0.95, base_label = "boot"
              )
              out["HR"] <- fg_row$SHR
            }
            out
          }, error = function(e) out)
        }
      })

      boot_stats[[bl]] <- run_parallel_bootstrap(
        n_cores = if_bootstrap_n_cores,
        bootstrap_count = if_bootstrap_count,
        boot_fn = boot_fn,
        seed = if_bootstrap_seed,
        export_varlist = c("pair_resample_index", "fit_finegray_shr"),
        result_colnames = c("HR")
      )
    }
    results$bootstrap <- boot_stats

    # summarize percentile CIs / SEs into the output table
    alpha <- 1 - confidence_level
    qs <- function(x) c(lci = unname(stats::quantile(x, alpha / 2, na.rm = TRUE)),
                        uci = unname(stats::quantile(x, 1 - alpha / 2, na.rm = TRUE)),
                        se  = stats::sd(x, na.rm = TRUE))
    if (hr_model == "Cox") {
      hazard_ratios$HR_BootLCI <- NA_real_
      hazard_ratios$HR_BootUCI <- NA_real_
      hazard_ratios$HR_BootSE  <- NA_real_
    } else {
      subdist_hazard$SHR_BootLCI <- NA_real_
      subdist_hazard$SHR_BootUCI <- NA_real_
      subdist_hazard$SHR_BootSE  <- NA_real_
    }
    for (bl in names(boot_stats)) {
      s_hr <- qs(boot_stats[[bl]][, "HR"])
      if (hr_model == "Cox") {
        hr_sel <- grepl(paste0("^", bl), hazard_ratios$Analysis)
        hazard_ratios$HR_BootLCI[hr_sel] <- s_hr["lci"]
        hazard_ratios$HR_BootUCI[hr_sel] <- s_hr["uci"]
        hazard_ratios$HR_BootSE[hr_sel]  <- s_hr["se"]
      } else {
        fg_sel <- grepl(paste0("^", bl), subdist_hazard$Analysis)
        subdist_hazard$SHR_BootLCI[fg_sel] <- s_hr["lci"]
        subdist_hazard$SHR_BootUCI[fg_sel] <- s_hr["uci"]
        subdist_hazard$SHR_BootSE[fg_sel]  <- s_hr["se"]
      }
    }
  }

  # STEP 3: Excel output ##############################################
  if (!is.null(out_xlsxpath)) {
    wb <- openxlsx::createWorkbook()
    if (!is.null(readme_text)) add_readme_sheet(wb, readme_text, verbose = verbose)

    openxlsx::addWorksheet(wb, "Analysis Summary")
    summary_info <- data.frame(
      Parameter = c("Hazard Model", "Outcome", "Exposure",
                    "Follow-up Time Variable", "Confidence Level",
                    "Blocks", "Weight Variable", "Match Id",
                    "Strata Variable", "Analytical CI methods", "Bootstrap"),
      Value = c(hr_model, outcome_var, exposure_var, followuptime_var,
                paste0(confidence_level * 100, "%"),
                paste(names(blocks), collapse = " + "),
                ifelse(is.null(weight_var), "N/A", weight_var),
                ifelse(is.null(match_id), "N/A", match_id),
                ifelse(do_stratify, strata_var, "N/A"),
                paste0("HR/sHR crude=model SE, weighted=robust, ",
                       "matched=robust+cluster(match id)."),
                ifelse(do_boot,
                       sprintf("fixed-weight (frozen), %d replicates, seed=%s",
                               if_bootstrap_count,
                               ifelse(is.null(if_bootstrap_seed), "user RNG",
                                      as.character(if_bootstrap_seed))),
                       "off"))
    )
    openxlsx::writeDataTable(wb, "Analysis Summary", summary_info, startRow = 1)

    if (!is.null(hazard_ratios)) {
      openxlsx::addWorksheet(wb, "Hazard Ratios")
      openxlsx::writeDataTable(wb, "Hazard Ratios", hazard_ratios, startRow = 1)
    }
    if (!is.null(subdist_hazard)) {
      openxlsx::addWorksheet(wb, "Fine-Gray SHR")
      openxlsx::writeDataTable(wb, "Fine-Gray SHR", subdist_hazard, startRow = 1)
    }

    out_dir <- dirname(out_xlsxpath)
    if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
    openxlsx::saveWorkbook(wb, out_xlsxpath, overwrite = TRUE)
    if (verbose) message(sprintf("\u2713 Results saved to: %s", out_xlsxpath))
  }

  # --- assemble results ---
  if (!is.null(hazard_ratios))  results$hazard_ratios  <- hazard_ratios
  if (!is.null(subdist_hazard)) results$subdist_hazard <- subdist_hazard
  if (do_stratify) results$strata_var <- strata_var

  if (verbose) message("\nANALYSIS COMPLETE")
  return(invisible(results))
}
