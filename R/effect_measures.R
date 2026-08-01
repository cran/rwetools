# HR IR Incidence Rates and Hazard Ratios ############################################
#' Estimate Incidence Rates, Hazard Ratios and Incidence Rate Ratios
#'
#' Estimates incidence rates (IR), incidence rate differences (IRD), a hazard
#' model (Cox HR or Fine-Gray subdistribution HR), and marginal incidence
#' rate ratios (IRR) for a time-to-event outcome, for up to two analysis
#' blocks: a crude block (\code{in_df_crude}) plus either a weighted block
#' (\code{in_df_weight}, e.g. IPTW) or a matched block (\code{in_df_match}).
#'
#' Analytical confidence intervals are always computed with a fixed,
#' per-block method (there are no CI-method arguments):
#' \itemize{
#'   \item Crude/matched IR: Garwood exact-Poisson (gamma); standardized IR
#'     arm CIs: Fay--Feuer gamma (single stratum collapses to Garwood).
#'   \item Weighted IR/IRD: design-based robust (sandwich) variance via a
#'     joint weighted quasi-Poisson cell model (Lumley 2004); standardized
#'     variants use the joint sandwich/delta method.
#'   \item HR: crude block uses the Cox model SE; weighted block the robust
#'     sandwich SE; matched block the robust SE clustered on
#'     \code{if_match_match_id} (Lin--Wei 1989; Austin 2016).
#'   \item Fine-Gray sHR (\code{hr_model = "Fine-Gray"}): robust SE clustered
#'     on the subject (crude/weighted) or the match id (matched); weighted
#'     expansion weights are IPCW x PS weight (Fine & Gray 1999).
#'   \item IRR (always reported): crude/matched blocks use a marginal Poisson
#'     rate model (model SE, equal to sqrt(1/D1 + 1/D0)); the weighted block
#'     uses survey::svyglm quasi-Poisson (robust SE). With
#'     \code{stratification_var}, the IRR is direct-standardized.
#' }
#'
#' @param in_df_crude Data frame for the crude (unadjusted) block, or NULL.
#' @param in_df_weight Data frame for the weighted block (requires
#'   \code{if_weight_weight_var}), or NULL. Cannot be combined with
#'   \code{in_df_match}.
#' @param in_df_match Data frame for the matched block, or NULL. Cannot be
#'   combined with \code{in_df_weight}.
#' @param out_xlsxpath Character path for an Excel output file, or NULL.
#' @param exposure_var Character. Binary exposure variable name.
#' @param exp_value,ref_value Values of \code{exposure_var} for the exposed
#'   and reference groups.
#' @param outcome_var Character. Event indicator (0 = censored, 1 = event).
#' @param followuptime_var Character. Follow-up time variable (required).
#' @param time_unit "days", "months" or "years" -- the unit of
#'   \code{followuptime_var}; person-years are derived as days/365.25,
#'   months/12, or years as-is.
#' @param ir_per_pyears Multiplier for IR/IRD (1, 100, 1000, 10000, 100000).
#' @param confidence_level Confidence level in (0, 1). Default 0.95.
#' @param stratification_var Character or NULL. When supplied, the hazard
#'   model is stratified via \code{strata()} and IR/IRD/IRR are
#'   direct-standardized with the total (weighted) person-time distribution
#'   as the standard.
#' @param hr_model "Cox" (cause-specific HR) or "Fine-Gray" (subdistribution
#'   HR; requires \code{if_fg_competing_event_var}).
#' @param if_fg_competing_event_var Character or NULL. Competing-event
#'   indicator (1 = competing event); required for \code{hr_model =
#'   "Fine-Gray"}, ignored (with a message) for "Cox".
#' @param if_weight_weight_var Character or NULL. Weight column in
#'   \code{in_df_weight}; required with that block, ignored (message)
#'   otherwise.
#' @param if_match_match_id Character or NULL. Match-set id column in
#'   \code{in_df_match}; enables pair clustering (HR/Fine-Gray SE) and
#'   pair-level bootstrap resampling. Without it the matched block uses a
#'   non-clustered robust SE and row resampling (a message is emitted).
#'   Ignored (message) without \code{in_df_match}.
#' @param if_bootstrap_count Integer or NULL. When supplied, adds percentile
#'   bootstrap CIs (columns \code{*_Boot}) for IR, IRD, IRR and the hazard
#'   estimate. The bootstrap is ALWAYS fixed-weight ("frozen"): weights and
#'   standardization shares are held at their original values, so
#'   PS-estimation uncertainty is not propagated (a message states this).
#' @param if_bootstrap_n_cores Integer or NULL. Cores for the bootstrap
#'   (NULL = all minus one); ignored (message) without
#'   \code{if_bootstrap_count}.
#' @param if_bootstrap_seed Integer or NULL. Bootstrap RNG seed; ignored
#'   (message) without \code{if_bootstrap_count}.
#' @param readme_text Optional README text for the Excel output.
#' @param verbose Logical. Print progress and method messages (default TRUE).
#'
#' @return Invisibly, a list with \code{incidence_rates} (one row set with
#'   per-block column suffixes \code{_Crude} / \code{_Weighted} /
#'   \code{_Matched}), \code{hazard_ratios} (Cox) or \code{subdist_hazard}
#'   (Fine-Gray), \code{incidence_rate_ratios}, per-block model objects
#'   (\code{models}), stratum details when standardized, bootstrap matrices
#'   when run, and \code{ir_per_pyears}.
#'
#' @section Side Effects:
#' Creates the output directory and writes an Excel workbook when
#' \code{out_xlsxpath} is provided.
#'
#' @export
#'
#' @examples
#' csv_path <- system.file("extdata", "sample_data.csv", package = "rwetools")
#' df <- read.csv(csv_path)
#'
#' # Crude block only
#' res <- estimate_hr_ir(
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
#' # Crude + weighted blocks, Fine-Gray hazard model, bootstrap CIs
#' df_ps <- estimate_ps(
#'   in_df = df, exposure_var = "exposure",
#'   class_vars = c("cat1", "cat2", "cat3", "cat4"),
#'   cont_vars  = c("cont1", "cont2", "cont3"),
#'   verbose = FALSE
#' )
#' df_wt <- create_iptw(
#'   in_df = df_ps, exposure_var = "exposure",
#'   ps_var = "ps", weight_var = "iptw", verbose = FALSE
#' )
#' res_fg <- estimate_hr_ir(
#'   in_df_crude               = df,
#'   in_df_weight              = df_wt,
#'   if_weight_weight_var      = "iptw",
#'   exposure_var              = "exposure",
#'   outcome_var               = "outcome",
#'   followuptime_var          = "follow_up_days",
#'   hr_model                  = "Fine-Gray",
#'   if_fg_competing_event_var = "competing_event",
#'   verbose                   = FALSE
#' )
#' res_fg$subdist_hazard
#'
#' res_boot <- estimate_hr_ir(
#'   in_df_crude          = df,
#'   exposure_var         = "exposure",
#'   outcome_var          = "outcome",
#'   followuptime_var     = "follow_up_days",
#'   if_bootstrap_count   = 200,
#'   if_bootstrap_n_cores = 1,
#'   if_bootstrap_seed    = 2026,
#'   verbose              = FALSE
#' )
#' }
estimate_hr_ir <- function(
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
    stratification_var = NULL,
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
  if (verbose) message("INCIDENCE RATE, HAZARD & RATE-RATIO ESTIMATION")
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
  hr_model  <- match.arg(hr_model)

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

  # --- block validation (shared with estimate_rr_rd) ---
  required_vars <- c(outcome_var, exposure_var, followuptime_var)
  if (!is.null(stratification_var)) required_vars <- c(required_vars, stratification_var)
  if (!is.null(if_fg_competing_event_var)) required_vars <- c(required_vars, if_fg_competing_event_var)

  vb <- validate_effect_blocks(
    in_df_crude = in_df_crude, in_df_weight = in_df_weight,
    in_df_match = in_df_match,
    if_weight_weight_var = if_weight_weight_var,
    if_match_match_id = if_match_match_id,
    required_vars = required_vars, verbose = verbose
  )
  blocks    <- vb$blocks
  weight_var <- vb$weight_var
  match_id   <- vb$match_id
  do_stratify <- !is.null(stratification_var)

  ir_col_base <- "IR_per_Npy"
  z_score <- stats::qnorm(1 - (1 - confidence_level) / 2)

  if (verbose) {
    message(sprintf("Outcome: %s", outcome_var))
    message(sprintf("Follow-up time: %s (%s)", followuptime_var, time_unit))
    message(sprintf("Hazard model: %s", hr_model))
    message(sprintf("Blocks: %s", paste(names(blocks), collapse = " + ")))
    if (do_stratify) {
      message(sprintf("Stratification: %s (Cox strata(); IR/IRD/IRR direct standardization)",
                      stratification_var))
    }
    message(sprintf("Confidence level: %.0f%%", confidence_level * 100))
    message("Analytical CIs (fixed per block): IR crude/matched = exact Poisson (Garwood),")
    message("  standardized = Fay-Feuer; weighted IR/IRD = robust sandwich (joint delta);")
    message(sprintf("  %s: Crude = model SE, Weighted = robust, Matched = robust%s;",
                    ifelse(hr_model == "Cox", "HR", "Fine-Gray sHR"),
                    ifelse(!is.null(match_id), " + cluster(match id)", "")))
    message("  IRR: crude/matched = marginal Poisson (model SE), weighted = svyglm robust.")
    message("")
  }

  # --- per-block preparation ---
  match_id_missing_msged <- FALSE
  prep_block <- function(df, block) {
    # recode exposure
    uv <- unique(df[[exposure_var]])
    if (!(all(uv %in% c(0, 1)))) {
      if (verbose) message(sprintf("Recoding exposure variable (%s block)...", block))
      df[[exposure_var]] <- ifelse(df[[exposure_var]] == exp_value, 1,
                                   ifelse(df[[exposure_var]] == ref_value, 0, NA))
    }
    # person-years
    df$person_years <- switch(time_unit,
      days   = as.numeric(df[[followuptime_var]]) / 365.25,
      months = as.numeric(df[[followuptime_var]]) / 12,
      years  = as.numeric(df[[followuptime_var]])
    )
    # complete cases on used columns
    used <- c(outcome_var, exposure_var, followuptime_var, "person_years")
    if (block == "Weighted") used <- c(used, weight_var)
    if (block == "Matched" && !is.null(match_id)) used <- c(used, match_id)
    if (do_stratify) used <- c(used, stratification_var)
    if (!is.null(if_fg_competing_event_var)) used <- c(used, if_fg_competing_event_var)
    n0 <- nrow(df)
    df <- df[stats::complete.cases(df[, used, drop = FALSE]), , drop = FALSE]
    if (nrow(df) < n0) {
      warning(sprintf("%d observation(s) with missing values removed from the %s block",
                      n0 - nrow(df), block))
    }
    # Fine-Gray needs an unambiguous event type
    if (hr_model == "Fine-Gray") {
      cev <- df[[if_fg_competing_event_var]]
      cev_vals <- unique(cev[!is.na(cev)])
      if (!all(cev_vals %in% c(0, 1))) {
        stop(paste0("if_fg_competing_event_var ('", if_fg_competing_event_var,
                    "') must be binary (0/1) in the ", block, " block."))
      }
      n_both <- sum(df[[outcome_var]] == 1 & cev == 1, na.rm = TRUE)
      if (n_both > 0L) {
        stop(sprintf(paste0("Fine-Gray: %d row(s) in the %s block have both ",
                            "outcome_var == 1 and if_fg_competing_event_var == 1, ",
                            "which is ambiguous. Resolve these first."), n_both, block))
      }
    }
    df
  }
  blocks <- Map(prep_block, blocks, names(blocks))

  results <- list()
  results$models <- list()

  # STEP 1: incidence rates ###########################################
  if (verbose) message("Step 1: Incidence rates & rate differences")
  incidence_rates <- NULL
  stratum_details <- list()

  for (bl in names(blocks)) {
    df <- blocks[[bl]]
    wv <- if (bl == "Weighted") weight_var else NULL
    if (do_stratify) {
      ir_res <- if (is.null(wv)) {
        calc_standardized_ir_ird(df, exposure_var, outcome_var, "person_years",
                                 stratification_var, confidence_level, z_score,
                                 bl, ir_per_pyears, ir_col_base)
      } else {
        calc_standardized_ir_ird_weighted(df, exposure_var, outcome_var, "person_years",
                                          wv, stratification_var, confidence_level,
                                          z_score, bl, ir_per_pyears, ir_col_base)
      }
      block_summary <- ir_res$summary
      stratum_details[[bl]] <- ir_res$stratum_detail
    } else {
      block_summary <- if (is.null(wv)) {
        calc_crude_ir_ird(df, exposure_var, outcome_var, "person_years",
                          confidence_level, z_score, bl, ir_per_pyears, ir_col_base)
      } else {
        calc_crude_ir_ird_weighted(df, exposure_var, outcome_var, "person_years",
                                   wv, confidence_level, z_score, bl,
                                   ir_per_pyears, ir_col_base)
      }
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

  # STEP 2: hazard model ##############################################
  if (verbose) message(sprintf("\nStep 2: %s model", hr_model))

  fit_cox_block <- function(df, block) {
    lhs <- paste0("survival::Surv(", followuptime_var, ", ", outcome_var, ")")
    rhs <- exposure_var
    if (do_stratify) rhs <- paste0(rhs, " + survival::strata(", stratification_var, ")")

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
    if (do_stratify) label <- paste0(label, " (stratified by ", stratification_var, ")")

    s    <- summary(model, conf.int = confidence_level)
    cf   <- s$coefficients
    ci   <- s$conf.int
    ridx <- which(rownames(cf) == exposure_var)
    data.frame(
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
    ) -> row_out
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
      if (do_stratify) lab <- paste0(lab, " (stratified by ", stratification_var, ")")
      fg_row <- fit_finegray_shr(
        data = df, time_var = followuptime_var, event_var = outcome_var,
        competing_var = if_fg_competing_event_var, exp_var = exposure_var,
        weight_var = if (bl == "Weighted") weight_var else NULL,
        strat_var = stratification_var,
        cluster_id_var = if (bl == "Matched") match_id else NULL,
        confidence_level = confidence_level, base_label = lab
      )
      subdist_hazard <- rbind(subdist_hazard, fg_row)
    }
  }
  if (verbose && !is.null(hazard_ratios))  print(hazard_ratios,  digits = 3)
  if (verbose && !is.null(subdist_hazard)) print(subdist_hazard, digits = 3)

  # STEP 3: incidence rate ratios (always on) #########################
  if (verbose) message("\nStep 3: Incidence rate ratios (marginal, always reported)")
  incidence_rate_ratios <- NULL
  for (bl in names(blocks)) {
    df <- blocks[[bl]]
    wv <- if (bl == "Weighted") weight_var else NULL
    if (any(df$person_years <= 0, na.rm = TRUE)) {
      stop(sprintf("IRR rate model requires positive person-time; non-positive person_years in the %s block.", bl))
    }
    irr_row <- if (do_stratify) {
      fit_irr_model_stratified(df, outcome_var, exposure_var, "person_years",
                               strat_var = stratification_var, weight_var = wv,
                               confidence_level = confidence_level, base_label = bl)
    } else {
      fit_irr_model(df, outcome_var, exposure_var, "person_years",
                    weight_var = wv, confidence_level = confidence_level,
                    base_label = bl)
    }
    incidence_rate_ratios <- rbind(incidence_rate_ratios, irr_row)
  }
  if (verbose) print(incidence_rate_ratios, digits = 3)

  # STEP 4: bootstrap layer (opt-in, fixed-weight) ####################
  if (do_boot) {
    if (verbose) {
      message(sprintf("\nStep 4: Bootstrap (%d replicates, %d core(s))",
                      if_bootstrap_count, if_bootstrap_n_cores))
    }
    message("Bootstrap method: fixed-weight (frozen) - weights and standardization ",
            "shares are held at their original values; PS-estimation uncertainty ",
            "is NOT propagated. Percentile CIs.")

    boot_stats <- list()
    for (bl in names(blocks)) {
      df <- blocks[[bl]]
      wv <- if (bl == "Weighted") weight_var else NULL
      use_pairs <- (bl == "Matched") && !is.null(match_id)
      if (bl == "Matched" && is.null(match_id) && verbose) {
        message("Matched block bootstrap: row resampling (no if_match_match_id).")
      }

      # fixed standardization shares (frozen, like the weights)
      std_w_fixed <- NULL
      if (do_stratify) {
        pt <- df$person_years * (if (is.null(wv)) 1 else df[[wv]])
        pt_s <- vapply(split(pt, df[[stratification_var]]),
                       function(x) sum(x, na.rm = TRUE), numeric(1))
        std_w_fixed <- pt_s / sum(pt_s)
      }

      boot_fn <- local({
        df_b <- df; wv_b <- wv; use_pairs_b <- use_pairs
        match_col <- if (use_pairs) df[[match_id]] else NULL
        exposure_b <- exposure_var; outcome_b <- outcome_var
        ftime_b <- followuptime_var; strat_b <- stratification_var
        std_w_b <- std_w_fixed; do_strat_b <- do_stratify
        hr_model_b <- hr_model; fg_var_b <- if_fg_competing_event_var
        mult_b <- ir_per_pyears
        function() {
          idx <- if (use_pairs_b) {
            pair_resample_index(match_col)
          } else {
            sample.int(nrow(df_b), nrow(df_b), replace = TRUE)
          }
          d <- df_b[idx, , drop = FALSE]
          out <- c(HR = NA_real_, IR1 = NA_real_, IR0 = NA_real_,
                   IRD = NA_real_, IRR = NA_real_)
          tryCatch({
            w <- if (is.null(wv_b)) rep(1, nrow(d)) else d[[wv_b]]
            if (do_strat_b) {
              lev <- names(std_w_b)
              ir_arm <- function(a) {
                sum(vapply(lev, function(s) {
                  sel <- d[[exposure_b]] == a & as.character(d[[strat_b]]) == s
                  ev <- sum(d[[outcome_b]][sel] * w[sel])
                  py <- sum(d$person_years[sel] * w[sel])
                  if (py > 0) std_w_b[[s]] * ev / py else 0
                }, numeric(1)))
              }
              ir1 <- ir_arm(1); ir0 <- ir_arm(0)
            } else {
              arm <- function(a) {
                sel <- d[[exposure_b]] == a
                c(ev = sum(d[[outcome_b]][sel] * w[sel]),
                  py = sum(d$person_years[sel] * w[sel]))
              }
              a1 <- arm(1); a0 <- arm(0)
              ir1 <- a1["ev"] / a1["py"]; ir0 <- a0["ev"] / a0["py"]
            }
            out["IR1"] <- mult_b * ir1
            out["IR0"] <- mult_b * ir0
            out["IRD"] <- mult_b * (ir1 - ir0)
            out["IRR"] <- ir1 / ir0
            # hazard point estimate
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

      boot_mat <- run_parallel_bootstrap(
        n_cores = if_bootstrap_n_cores,
        bootstrap_count = if_bootstrap_count,
        boot_fn = boot_fn,
        seed = if_bootstrap_seed,
        export_varlist = c("pair_resample_index", "fit_finegray_shr"),
        result_colnames = c("HR", "IR1", "IR0", "IRD", "IRR")
      )
      boot_stats[[bl]] <- boot_mat
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
      bm <- boot_stats[[bl]]
      s_ir1 <- qs(bm[, "IR1"]); s_ir0 <- qs(bm[, "IR0"])
      s_ird <- qs(bm[, "IRD"]); s_irr <- qs(bm[, "IRR"]); s_hr <- qs(bm[, "HR"])

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

  # STEP 5: Excel output ##############################################
  if (!is.null(out_xlsxpath)) {
    wb <- openxlsx::createWorkbook()
    if (!is.null(readme_text)) add_readme_sheet(wb, readme_text, verbose = verbose)

    openxlsx::addWorksheet(wb, "Analysis Summary")
    summary_info <- data.frame(
      Parameter = c("Hazard Model", "Outcome", "Exposure", "Follow-up Time Variable",
                    "Time Unit", "IR per Person-Years Multiplier", "Confidence Level",
                    "Blocks", "Weight Variable", "Match Id",
                    "Stratification Variable",
                    "Analytical CI methods",
                    "Bootstrap"),
      Value = c(hr_model, outcome_var, exposure_var, followuptime_var,
                time_unit, format(ir_per_pyears, big.mark = ","),
                paste0(confidence_level * 100, "%"),
                paste(names(blocks), collapse = " + "),
                ifelse(is.null(weight_var), "N/A", weight_var),
                ifelse(is.null(match_id), "N/A", match_id),
                ifelse(do_stratify, stratification_var, "N/A"),
                paste0("IR crude/matched=Garwood; standardized=Fay-Feuer; ",
                       "weighted=robust sandwich. HR/sHR crude=model SE, ",
                       "weighted=robust, matched=robust+cluster. ",
                       "IRR crude/matched=Poisson model SE, weighted=svyglm robust."),
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

    if (!is.null(hazard_ratios)) {
      openxlsx::addWorksheet(wb, "Hazard Ratios")
      openxlsx::writeDataTable(wb, "Hazard Ratios", hazard_ratios, startRow = 1)
    }
    if (!is.null(subdist_hazard)) {
      openxlsx::addWorksheet(wb, "Fine-Gray SHR")
      openxlsx::writeDataTable(wb, "Fine-Gray SHR", subdist_hazard, startRow = 1)
    }
    openxlsx::addWorksheet(wb, "Incidence Rate Ratios")
    openxlsx::writeDataTable(wb, "Incidence Rate Ratios", incidence_rate_ratios, startRow = 1)

    for (bl in names(stratum_details)) {
      sheet_name <- sprintf("Stratum Detail (%s Std)", bl)
      openxlsx::addWorksheet(wb, sheet_name)
      openxlsx::writeDataTable(wb, sheet_name, stratum_details[[bl]], startRow = 1)
    }

    out_dir <- dirname(out_xlsxpath)
    if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
    openxlsx::saveWorkbook(wb, out_xlsxpath, overwrite = TRUE)
    if (verbose) message(sprintf("\u2713 Results saved to: %s", out_xlsxpath))
  }

  # --- assemble results ---
  results$incidence_rates <- incidence_rates
  if (!is.null(hazard_ratios))  results$hazard_ratios  <- hazard_ratios
  if (!is.null(subdist_hazard)) results$subdist_hazard <- subdist_hazard
  results$incidence_rate_ratios <- incidence_rate_ratios
  results$ir_per_pyears <- ir_per_pyears
  if (do_stratify) {
    results$stratification_var <- stratification_var
    for (bl in names(stratum_details)) {
      results[[paste0("stratum_details_", tolower(bl))]] <- stratum_details[[bl]]
    }
  }

  if (verbose) message("\nANALYSIS COMPLETE")
  return(invisible(results))
}


# RR Risk Ratios and Risk Differences ##############
#' Estimate Risk Ratios and Risk Differences using Cumulative Incidence
#'
#' Estimates risks (cumulative incidence), risk ratios (RR) and risk
#' differences (RD) at a specified timepoint for up to two analysis blocks:
#' a crude block (\code{in_df_crude}) plus either a weighted block
#' (\code{in_df_weight}, e.g. IPTW) or a matched block (\code{in_df_match}).
#' Two estimators are supported: Kaplan-Meier (\code{"KM"}, 1 - S(t)) and
#' Aalen-Johansen (\code{"AJ"}, the cumulative incidence function under
#' competing risks; requires \code{if_aj_competing_event_var}).
#'
#' Analytical confidence intervals are always computed with fixed methods:
#' \itemize{
#'   \item Risk / CIF: complementary log-log (cloglog) interval on the final
#'     (standardized/weighted) estimate with its combined SE
#'     (Kalbfleisch & Prentice 2002) -- for a crude KM risk this reproduces
#'     \code{survfit(conf.type = "log-log")} exactly. A risk of exactly 0 or
#'     1 has no defined cloglog CI: NA is returned and a message recommends
#'     the bootstrap.
#'   \item RD: normal-Wald with Greenwood (KM) / counting-process (AJ)
#'     variance. RR: delta method on the log scale.
#'   \item Standardization (with \code{stratification_var}): stratum weights
#'     w_k = N_k / N_total; \code{Var(Risk_std) = Sum(w_k^2 Var(R_k))}.
#' }
#' The percentile bootstrap is opt-in via \code{if_bootstrap_count} and is
#' ALWAYS fixed-weight ("frozen"): analysis weights and stratum shares are
#' held at their original values (PS-estimation uncertainty is not
#' propagated); matched blocks resample matched sets by
#' \code{if_match_match_id}.
#'
#' @param in_df_crude Data frame for the crude block, or NULL.
#' @param in_df_weight Data frame for the weighted block (requires
#'   \code{if_weight_weight_var}); cannot be combined with \code{in_df_match}.
#' @param in_df_match Data frame for the matched block; cannot be combined
#'   with \code{in_df_weight}.
#' @param out_xlsxpath Character path for an Excel output file, or NULL.
#' @param exposure_var Character. Binary exposure variable name.
#' @param exp_value,ref_value Values of \code{exposure_var} for the exposed
#'   and reference groups.
#' @param outcome_var Character. Event indicator (0 = censored, 1 = event).
#' @param followuptime_var Character. Follow-up time variable (required).
#' @param time_unit "days", "months" or "years" -- the unit of
#'   \code{followuptime_var} and \code{rr_rd_at_timepoint}.
#' @param rr_rd_at_timepoint Numeric timepoint for the cumulative incidence
#'   (same unit as \code{followuptime_var}). Default 365.
#' @param risk_per_individuals Denominator for expressing risks and risk
#'   differences (default 1000).
#' @param confidence_level Confidence level in (0, 1). Default 0.95.
#' @param risk_estimator "KM" (Kaplan-Meier) or "AJ" (Aalen-Johansen;
#'   requires \code{if_aj_competing_event_var}). Aliases "Kaplan-Meier" /
#'   "Aalen-Johansen" are accepted.
#' @param if_aj_competing_event_var Character or NULL. Competing-event
#'   indicator (1 = competing event) for the AJ estimator; required with
#'   \code{risk_estimator = "AJ"}, ignored (message) with "KM".
#' @param stratification_var Character or NULL. Direct standardization
#'   variable (stratum weights from the total population).
#' @param if_weight_weight_var Character or NULL. Weight column in
#'   \code{in_df_weight}; required with that block, ignored (message)
#'   otherwise.
#' @param if_match_match_id Character or NULL. Match-set id column in
#'   \code{in_df_match}; enables pair-level bootstrap resampling. Without it
#'   the matched block uses row resampling (a message is emitted). Ignored
#'   (message) without \code{in_df_match}.
#' @param if_bootstrap_count Integer or NULL. When supplied, adds percentile
#'   bootstrap CIs (the \code{*_Boot} columns). Fixed-weight (frozen) only.
#' @param if_bootstrap_n_cores Integer or NULL. Cores for the bootstrap
#'   (NULL = all minus one); ignored (message) without
#'   \code{if_bootstrap_count}.
#' @param if_bootstrap_seed Integer or NULL. Bootstrap RNG seed; ignored
#'   (message) without \code{if_bootstrap_count}.
#' @param readme_text Optional README text for the Excel output.
#' @param verbose Logical. Print progress messages (default TRUE).
#'
#' @return Invisibly, a list with \code{estimates} (one row per block; both
#'   analytical and \code{*_Boot} columns), \code{cumulative_incidence}
#'   (per-arm risks with cloglog CIs), \code{stratum_details} (when
#'   standardized), and per-block bootstrap matrices when run.
#'
#' @section Side Effects:
#' Creates the output directory and writes an Excel workbook when
#' \code{out_xlsxpath} is provided.
#'
#' @export
#'
#' @examples
#' csv_path <- system.file("extdata", "sample_data.csv", package = "rwetools")
#' df <- read.csv(csv_path)
#'
#' # Crude KM risks with analytical (cloglog / Wald / log-delta) CIs
#' res <- estimate_rr_rd(
#'   in_df_crude          = df,
#'   exposure_var         = "exposure",
#'   outcome_var          = "outcome",
#'   followuptime_var     = "follow_up_days",
#'   time_unit            = "days",
#'   rr_rd_at_timepoint   = 365,
#'   risk_per_individuals = 1000,
#'   verbose              = FALSE
#' )
#' res$estimates
#'
#' \donttest{
#' # AJ (competing risks) + bootstrap CIs
#' res2 <- estimate_rr_rd(
#'   in_df_crude               = df,
#'   exposure_var              = "exposure",
#'   outcome_var               = "outcome",
#'   followuptime_var          = "follow_up_days",
#'   rr_rd_at_timepoint        = 365,
#'   risk_estimator            = "AJ",
#'   if_aj_competing_event_var = "competing_event",
#'   if_bootstrap_count        = 100,
#'   if_bootstrap_n_cores      = 1,
#'   if_bootstrap_seed         = 2026,
#'   verbose                   = FALSE
#' )
#' res2$estimates
#' }
estimate_rr_rd <- function(
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
    rr_rd_at_timepoint = 365,
    risk_per_individuals = 1000,
    confidence_level = 0.95,
    risk_estimator = c("KM", "AJ"),
    if_aj_competing_event_var = NULL,
    stratification_var = NULL,
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
  if (!is.null(stratification_var)) {
    if (verbose) message("*** Direct Standardization Enabled ***")
  }
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
    message("Bootstrap method: fixed-weight (frozen) - weights and stratum ",
            "shares are held at their original values; PS-estimation ",
            "uncertainty is NOT propagated. Percentile CIs.")
  }

  # --- scalar validations ---
  if (is.null(outcome_var)) stop("outcome_var must be specified.")
  if (is.null(followuptime_var)) stop("followuptime_var must be specified.")
  if (!is.numeric(confidence_level) || length(confidence_level) != 1L ||
      confidence_level <= 0 || confidence_level >= 1) {
    stop("confidence_level must be a single numeric value between 0 and 1 (e.g., 0.95).")
  }
  if (!is.numeric(rr_rd_at_timepoint) || length(rr_rd_at_timepoint) != 1L ||
      rr_rd_at_timepoint <= 0) {
    stop("rr_rd_at_timepoint must be a single positive number.")
  }
  if (!is.numeric(risk_per_individuals) || length(risk_per_individuals) != 1L ||
      risk_per_individuals <= 0) {
    stop("risk_per_individuals must be a single positive number (e.g., 100, 1000, 10000, 100000).")
  }

  # --- block validation (shared with estimate_hr_ir) ---
  required_vars <- c(outcome_var, exposure_var, followuptime_var)
  if (!is.null(stratification_var)) required_vars <- c(required_vars, stratification_var)
  if (need_competing) required_vars <- c(required_vars, if_aj_competing_event_var)

  vb <- validate_effect_blocks(
    in_df_crude = in_df_crude, in_df_weight = in_df_weight,
    in_df_match = in_df_match,
    if_weight_weight_var = if_weight_weight_var,
    if_match_match_id = if_match_match_id,
    required_vars = required_vars, verbose = verbose
  )
  blocks   <- vb$blocks
  weight_var <- vb$weight_var
  match_id   <- vb$match_id

  # competing-event indicator must be binary in every provided block
  if (need_competing) {
    for (bl in names(blocks)) {
      cev_vals <- unique(blocks[[bl]][[if_aj_competing_event_var]])
      cev_vals <- cev_vals[!is.na(cev_vals)]
      if (!all(cev_vals %in% c(0, 1))) {
        stop(paste0("if_aj_competing_event_var ('", if_aj_competing_event_var,
                    "') must be binary (0/1) in the ", bl, " block. Found: ",
                    paste(sort(unique(cev_vals)), collapse = ", ")))
      }
    }
  }

  # stratification variable needs >= 2 levels per block
  if (!is.null(stratification_var)) {
    for (bl in names(blocks)) {
      if (length(unique(blocks[[bl]][[stratification_var]])) < 2L) {
        stop(paste("stratification_var", stratification_var,
                   "must have at least 2 levels in the", bl, "block."))
      }
    }
  }

  estimator_label <- ifelse(risk_estimator == "AJ", "Aalen-Johansen", "Kaplan-Meier")
  if (verbose) {
    message(sprintf("Outcome: %s", outcome_var))
    message(sprintf("Exposure: %s (Exposed=%s, Reference=%s)",
                    exposure_var, exp_value, ref_value))
    message(sprintf("Follow-up time: %s (%s)", followuptime_var, time_unit))
    message(sprintf("Timepoint for RR/RD: %s %s", rr_rd_at_timepoint, time_unit))
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
    if (!is.null(stratification_var)) {
      message(sprintf("Stratification Variable: %s (Direct Standardization)",
                      stratification_var))
    }
    message("")
  }

  z_val <- stats::qnorm(1 - (1 - confidence_level) / 2)

  results <- list()
  all_estimates <- data.frame()
  all_cuminc_data <- data.frame()
  all_stratum_details <- data.frame()

  #### 3. Run Analysis (common logic for all blocks) ########
  run_analysis <- function(data, analysis_label, wt_var = NULL, pair_col = NULL) {

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
    cols_to_check <- c(followuptime_var, outcome_var, exposure_var)
    if (!is.null(wt_var)) cols_to_check <- c(cols_to_check, wt_var)
    if (!is.null(pair_col)) cols_to_check <- c(cols_to_check, pair_col)
    if (!is.null(stratification_var)) cols_to_check <- c(cols_to_check, stratification_var)
    if (need_competing) cols_to_check <- c(cols_to_check, if_aj_competing_event_var)

    n_before <- nrow(data)
    complete_mask <- stats::complete.cases(data[, cols_to_check, drop = FALSE])
    data <- data[complete_mask, , drop = FALSE]
    n_after <- nrow(data)
    if (n_before > n_after && verbose) {
      message(sprintf("  Removed %d observations with missing values", n_before - n_after))
    }

    # --- Slim data: keep only needed columns for bootstrap efficiency ---
    keep_cols <- unique(c(followuptime_var, outcome_var, exposure_var, wt_var,
                          pair_col, stratification_var,
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

    ###### Stratification info ################
    if (!is.null(stratification_var)) {
      strata_tab <- table(data[[stratification_var]])
      if (verbose) {
        message(sprintf("  Stratification by: %s (%d strata)",
                        stratification_var, length(strata_tab)))
      }
    }

    # --- Common: count N and events per group ---
    is_ref <- data[[exposure_var]] == ref_val_int
    is_exp <- data[[exposure_var]] == exp_val_int
    in_tp  <- data[[followuptime_var]] <= rr_rd_at_timepoint

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
    if (!is.null(stratification_var)) {

      if (verbose) message("  Calculating standardized point estimates...")

      data_by_stratum <- split(data, data[[stratification_var]])

      std_result <- calc_standardized_rr_rd(
        data = data,
        time_var = followuptime_var,
        event_var = outcome_var,
        exp_var = exposure_var,
        strat_var = stratification_var,
        weight_var = wt_var,
        timepoint = rr_rd_at_timepoint,
        exp_val = exp_val_int,
        ref_val = ref_val_int,
        per_n = risk_per_individuals,
        pre_split = data_by_stratum,
        risk_estimator = risk_estimator,
        competing_event_var = if_aj_competing_event_var
      )

      strat_detail_df <- std_result$strat_details
      strat_detail_df$Analysis <- analysis_label
      strat_detail_df$Risk_Estimator <- risk_estimator

      RR_point <- std_result$RR
      RD_point <- std_result$RD

      cuminc_data <- data.frame(
        Exposure_Value = c(as.character(ref_val_int), as.character(exp_val_int)),
        N_Individuals = c(n_ref, n_exp),
        N_Events = c(events_ref, events_exp),
        Timepoint = rr_rd_at_timepoint,
        Risk_SE = c(std_result$se_risk_std_ref, std_result$se_risk_std_exp),
        RiskperN = c(std_result$risk_std_ref * risk_per_individuals,
                     std_result$risk_std_exp * risk_per_individuals),
        Risk = c(std_result$risk_std_ref, std_result$risk_std_exp),
        Risk_Var = c(std_result$se_risk_std_ref^2, std_result$se_risk_std_exp^2),
        Risk_Estimator = risk_estimator,
        Analysis = paste0(analysis_label, " (Standardized)"),
        stringsAsFactors = FALSE
      )

      # --- Analytical CI (stratified): RR log-delta, RD normal-Wald ---
      if (!is.na(std_result$se_log_RR) && std_result$se_log_RR > 0) {
        RR_lci_ana <- exp(log(RR_point) - z_val * std_result$se_log_RR)
        RR_uci_ana <- exp(log(RR_point) + z_val * std_result$se_log_RR)
      } else {
        RR_lci_ana <- NA_real_
        RR_uci_ana <- NA_real_
      }
      RR_se_ana <- std_result$lnRR_se_analytical

      RD_se_ana <- std_result$RD_se_analytical
      RD_lci_ana <- RD_point - z_val * RD_se_ana
      RD_uci_ana <- RD_point + z_val * RD_se_ana

      # --- Bootstrap CI (stratified; fixed strat weights) ---
      boot_results <- NULL
      if (do_boot) {
        fixed_strat_weights <- std_result$strat_weights
        if (verbose) {
          message(sprintf("  Running %d stratified bootstrap iterations (%d cores)...",
                          if_bootstrap_count, if_bootstrap_n_cores))
        }
        boot_results <- run_parallel_bootstrap(
          n_cores = if_bootstrap_n_cores,
          bootstrap_count = if_bootstrap_count,
          boot_fn = boot_standardized_rr_rd_single,
          seed = if_bootstrap_seed,
          export_varlist = c("calc_km_cumulative_incidence",
                             "calc_km_risk_single_group",
                             "calc_standardized_rr_rd", "calc_rr_rd",
                             "pair_resample_index"),
          data_by_stratum = data_by_stratum,
          strat_var = stratification_var,
          strat_weights = fixed_strat_weights,
          time_var = followuptime_var,
          event_var = outcome_var,
          exp_var = exposure_var,
          weight_var = wt_var,
          timepoint = rr_rd_at_timepoint,
          exp_val = exp_val_int,
          ref_val = ref_val_int,
          per_n = risk_per_individuals,
          risk_estimator = risk_estimator,
          competing_event_var = if_aj_competing_event_var,
          match_ids_by_stratum = if (!is.null(pair_col)) {
            lapply(data_by_stratum, function(d) d[[pair_col]])
          } else {
            NULL
          }
        )
      }

      risk_ref_per_n <- std_result$risk_std_ref * risk_per_individuals
      risk_exp_per_n <- std_result$risk_std_exp * risk_per_individuals
      analysis_label_out <- paste0(analysis_label, " (Standardized)")

    } else {

      # ========================================================
      # BRANCH B: Unstratified Analysis
      # ========================================================
      if (verbose) message("  Calculating point estimates...")

      cuminc <- calc_km_cumulative_incidence(
        data = data,
        time_var = followuptime_var,
        event_var = outcome_var,
        exp_var = exposure_var,
        weight_var = wt_var,
        timepoint = rr_rd_at_timepoint,
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
      strat_detail_df <- NULL

      # --- Analytical CI (unstratified): RR log-delta, RD normal-Wald ---
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
                             "calc_km_risk_single_group",
                             "calc_standardized_rr_rd", "calc_rr_rd",
                             "pair_resample_index"),
          data = data,
          time_var = followuptime_var,
          event_var = outcome_var,
          exp_var = exposure_var,
          weight_var = wt_var,
          timepoint = rr_rd_at_timepoint,
          exp_val = exp_val_int,
          ref_val = ref_val_int,
          per_n = risk_per_individuals,
          risk_estimator = risk_estimator,
          competing_event_var = if_aj_competing_event_var,
          match_ids = if (!is.null(pair_col)) data[[pair_col]] else NULL
        )
      }

      analysis_label_out <- analysis_label
    }

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
      Stratified_By    = ifelse(!is.null(stratification_var), stratification_var, NA_character_),
      Timepoint        = rr_rd_at_timepoint,
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

    # ★ cloglog CI for the 0-1 Risk/CIF (matrix r22-r29; single site — all
    # four cases funnel through here). Applied to the FINAL (standardized/
    # weighted) estimate with its combined SE. Risk on {0,1}: cloglog is
    # undefined -> NA bounds + a bootstrap-recommending message.
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
      strat_details = strat_detail_df,
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
    if (!is.null(bl_result$strat_details)) {
      all_stratum_details <- rbind(all_stratum_details, bl_result$strat_details)
    }
    results[[paste0(tolower(bl), "_bootstrap")]] <- bl_result$boot_results
  }

  #### 5. Print Summary ################################
  if (verbose) {
    message("\n========================================")
    message("SUMMARY OF RESULTS")
    message(sprintf("(Risk Difference per %s)", per_label))
    message("========================================")

    display_cols <- c("Analysis", "Stratified_By", "Timepoint", "Time_Unit", "Per_Individuals")
    if (do_boot) display_cols <- c(display_cols, "RR_CI_Boot", "RDperN_CI_Boot")
    display_cols <- c(display_cols, "RR_CI_Analytical", "RDperN_CI_Analytical")

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
      ifelse(!is.null(stratification_var),
             paste0("Standardized via direct standardization (total population as standard). ",
                    "Stratum weights w_k = N_k/N_total. ",
                    "Risk_std = Sum(w_k * Risk_k). "),
             ""),
      "Analytical SE: ",
      ifelse(risk_estimator == "AJ",
             "Counting-process variance (Aalen, 1978). ",
             "Greenwood-based Var(Risk). "),
      "Risk/CIF CI via complementary log-log (cloglog) transformation ",
      "(Kalbfleisch & Prentice 2002); Risk of exactly 0/1 -> NA. ",
      "RR CI via delta method on log scale. RD CI via normal approximation. ",
      ifelse(do_boot,
             paste0("Bootstrap CIs via percentile method, fixed-weight (frozen); ",
                    ifelse(!is.null(stratification_var),
                           "stratified resampling with fixed population weights; ",
                           ""),
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
      ifelse(!is.null(stratification_var), ", Direct Standardization", ""),
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
        if (need_competing) if_aj_competing_event_var else character(0),
        followuptime_var,
        time_unit,
        as.character(rr_rd_at_timepoint),
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
        ifelse(!is.null(stratification_var), stratification_var, "N/A (Unstratified)"),
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

    if (nrow(all_stratum_details) > 0L) {
      openxlsx::addWorksheet(wb, "Stratum Details")
      openxlsx::writeDataTable(wb, "Stratum Details", all_stratum_details)
      openxlsx::setColWidths(wb, "Stratum Details", cols = seq_len(ncol(all_stratum_details)),
                             widths = "auto")
    }

    out_dir <- dirname(out_xlsxpath)
    if (out_dir != "." && !dir.exists(out_dir)) {
      dir.create(out_dir, recursive = TRUE)
      if (verbose) message(sprintf("  Created output directory: %s", out_dir))
    }

    openxlsx::saveWorkbook(wb, out_xlsxpath, overwrite = TRUE)
    if (verbose) message(sprintf("\\u2713 Results saved to: %s", out_xlsxpath))
  }

  if (verbose) {
    message("\n========================================")
    message("ANALYSIS COMPLETE")
    message("========================================")
  }

  results$estimates <- all_estimates
  results$cumulative_incidence <- all_cuminc_data
  if (nrow(all_stratum_details) > 0L) {
    results$stratum_details <- all_stratum_details
  }

  return(invisible(results))
}
