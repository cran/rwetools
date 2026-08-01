# Effect-estimation statistical engines ###############################
# Lowest layer of the effect-estimation stack: pure statistical engines
# (CI formulas, weighted rate-cell fits, delta-method contrasts, resampling)
# that return raw estimates and know nothing about the reporting format.
# Called by the calculators in helpers_effect_measures.R and the public
# functions in effect_measures.R. (Introduced by the 0.4.0 API v2 redesign,
# EFFECT_CI_V2_PLAN.)

# HELPER: Fay-Feuer gamma CI for a directly standardized rate ##########
#' Fay-Feuer (1997) gamma confidence interval for a directly standardized
#' incidence rate IR_std = sum_s(std_w_s * d_s / py_s), with
#' Var(IR_std) = sum_s(std_w_s^2 * d_s / py_s^2). For a single stratum with
#' std_w = 1 the interval collapses exactly to the Garwood exact-Poisson
#' interval qgamma(alpha/2, d)/py, qgamma(1 - alpha/2, d + 1)/py (matrix r3;
#' regression requirement).
#'
#' @param events Numeric vector of per-stratum event counts d_s.
#' @param py Numeric vector of per-stratum person-times py_s (> 0).
#' @param std_w Numeric vector of standardization weights W_s (sum to 1).
#' @param conf_level Numeric in (0, 1).
#' @return List: rate, var, se, lci, uci (all on the per-1-person-time scale;
#'   the caller applies the reporting multiplier).
#' @noRd
fay_feuer_ci <- function(events, py, std_w, conf_level = 0.95) {
  if (any(py <= 0, na.rm = TRUE) || anyNA(events) || anyNA(py) || anyNA(std_w)) {
    return(list(rate = NA_real_, var = NA_real_, se = NA_real_,
                lci = NA_real_, uci = NA_real_))
  }
  alpha <- 1 - conf_level

  x  <- sum(std_w * events / py)        # standardized rate
  v  <- sum(std_w^2 * events / py^2)    # its variance estimate
  wm <- max(std_w / py)                 # largest single-event contribution

  # Lower: gamma with shape x^2/v, scale v/x (0 when no events anywhere)
  lci <- if (x > 0) {
    stats::qgamma(alpha / 2, shape = x^2 / v, scale = v / x)
  } else {
    0
  }

  # Upper: shift by the maximal weight-to-py ratio (Fay-Feuer 1997)
  x_up <- x + wm
  v_up <- v + wm^2
  uci  <- stats::qgamma(1 - alpha / 2, shape = x_up^2 / v_up, scale = v_up / x_up)

  list(rate = x, var = v, se = sqrt(v), lci = lci, uci = uci)
}


# HELPER: cloglog CI for a 0-1 risk / CIF ##############################
#' Complementary log-log confidence interval for a cumulative risk (KM
#' 1 - S(t)) or CIF (Aalen-Johansen), applied to the FINAL estimate with its
#' combined SE (matrix r22-r29): g = log(-log(1 - R)),
#' SE_g = SE_R / ((1 - R) * (-log(1 - R))), CI = 1 - exp(-exp(g +/- z*SE_g)).
#' For a crude KM risk this reproduces survfit(conf.type = "log-log") exactly.
#' Boundary policy: Risk on {0, 1} (cloglog undefined) returns NA bounds; the
#' caller emits the bootstrap-recommending message.
#'
#' @param risk Numeric vector of risks/CIFs in [0, 1].
#' @param se Numeric vector of SEs of risk (natural 0-1 scale).
#' @param z Numeric z-value for the confidence level.
#' @return List of numeric vectors lci, uci (NA where risk is on the boundary
#'   or inputs are missing).
#' @noRd
cloglog_risk_ci <- function(risk, se, z) {
  ok <- is.finite(risk) & is.finite(se) & risk > 0 & risk < 1

  g    <- rep(NA_real_, length(risk))
  se_g <- rep(NA_real_, length(risk))
  g[ok]    <- log(-log(1 - risk[ok]))
  se_g[ok] <- se[ok] / ((1 - risk[ok]) * (-log(1 - risk[ok])))

  lci <- 1 - exp(-exp(g - z * se_g))
  uci <- 1 - exp(-exp(g + z * se_g))
  list(lci = lci, uci = uci)
}


# HELPER: joint-sandwich engine for weighted rate cells ################
#' Fit the saturated weighted quasi-Poisson rate model over analysis cells
#' (arm, or arm x stratum) and return the cell log-rates with their joint
#' robust (design-based sandwich) variance-covariance matrix
#' (survey::svyglm; Lumley 2004, Binder 1983). This is the shared engine for
#' the weighted IR / IRD / IRR analytical methods (matrix r4/r8/r12-A and,
#' with strata, r5/r9/r13-A): every downstream quantity is a delta-method
#' contrast on these cell log-rates.
#'
#' Cells are exposure x stratum combinations; exp(log-rate) of a saturated
#' cell equals the weighted event sum / weighted person-time of that cell.
#' A cell with zero weighted events keeps the fit alive (large negative
#' coefficient); its own CI is set to NA downstream.
#'
#' @param df Data frame (already recoded exposure 0/1, person-time column).
#' @param exp_var,out_var,py_var,wt_var Column names (character).
#' @param strat_var Optional stratification column name, or NULL.
#' @return List: cells (data.frame with exp_val, strat_val, log_rate,
#'   n_events_w, py_w, zero_events flag), vcov (joint robust vcov, cells in
#'   row order of `cells`), or NULL invisibly if the fit fails.
#' @noRd
fit_weighted_rate_cells <- function(df, exp_var, out_var, py_var, wt_var,
                                    strat_var = NULL) {
  work <- df
  if (is.null(strat_var)) {
    work$.cell <- factor(work[[exp_var]])
  } else {
    work$.cell <- interaction(work[[exp_var]], work[[strat_var]],
                              drop = TRUE, sep = "|", lex.order = TRUE)
  }

  fml <- stats::as.formula(
    paste0(out_var, " ~ 0 + .cell + offset(log(", py_var, "))")
  )

  fit <- tryCatch({
    des <- survey::svydesign(
      ids = ~1,
      weights = stats::as.formula(paste0("~", wt_var)),
      data = work
    )
    survey::svyglm(fml, design = des, family = stats::quasipoisson())
  }, error = function(e) NULL)
  if (is.null(fit)) return(NULL)

  cf <- stats::coef(fit)
  V  <- stats::vcov(fit)

  cell_levels <- levels(work$.cell)
  coef_names  <- paste0(".cell", cell_levels)
  # Guard against dropped coefficients (should not happen for saturated cells)
  if (!all(coef_names %in% names(cf))) return(NULL)
  cf <- cf[coef_names]
  V  <- V[coef_names, coef_names, drop = FALSE]

  # Weighted event sums / person-time per cell (for reporting & zero flags)
  ne_w <- vapply(split(work[[out_var]] * work[[wt_var]], work$.cell),
                 function(x) sum(x, na.rm = TRUE), numeric(1))[cell_levels]
  py_w <- vapply(split(work[[py_var]] * work[[wt_var]], work$.cell),
                 function(x) sum(x, na.rm = TRUE), numeric(1))[cell_levels]

  if (is.null(strat_var)) {
    exp_val   <- cell_levels
    strat_val <- NA_character_
  } else {
    parts     <- strsplit(cell_levels, "|", fixed = TRUE)
    exp_val   <- vapply(parts, `[`, character(1), 1L)
    strat_val <- vapply(parts, `[`, character(1), 2L)
  }

  cells <- data.frame(
    exp_val     = exp_val,
    strat_val   = strat_val,
    log_rate    = unname(cf),
    n_events_w  = unname(ne_w),
    py_w        = unname(py_w),
    zero_events = unname(ne_w) <= 0,
    stringsAsFactors = FALSE
  )

  list(cells = cells, vcov = V)
}


# HELPER: weighted IR/IRD/IRR contrasts from rate cells (unstratified) #
#' Delta-method contrasts on the joint cell fit for the UNSTRATIFIED weighted
#' analysis: per-arm IR with robust log-scale Wald CI (r4), IRD with identity
#' -scale joint-delta Wald CI (r8), IRR with log-ratio delta CI (r12-A;
#' numerically identical to the saturated svyglm exposure coefficient).
#'
#' @param cell_fit Result of fit_weighted_rate_cells (strat_var = NULL).
#' @param z z-value; multiplier = reporting multiplier (per N person-years).
#' @return List: ir (data.frame per arm: exp_val, ir, lci, uci, se),
#'   ird (est, lci, uci, se), irr (est, lci, uci, ln, ln_se).
#'   Arms with zero weighted events get NA CIs.
#' @noRd
wtd_rate_contrasts <- function(cell_fit, z, multiplier) {
  cells <- cell_fit$cells
  V     <- cell_fit$vcov
  i0 <- which(cells$exp_val == "0")
  i1 <- which(cells$exp_val == "1")

  l  <- cells$log_rate
  ir <- exp(l)

  # Per-arm IR: robust log-scale Wald
  se_log <- sqrt(diag(V))
  ir_lci <- exp(l - z * se_log)
  ir_uci <- exp(l + z * se_log)
  ir_se  <- ir * se_log                       # natural-scale delta SE
  bad <- cells$zero_events
  ir_lci[bad] <- NA_real_; ir_uci[bad] <- NA_real_; ir_se[bad] <- NA_real_

  ir_df <- data.frame(
    exp_val = cells$exp_val,
    ir  = multiplier * ir,
    lci = multiplier * ir_lci,
    uci = multiplier * ir_uci,
    se  = multiplier * ir_se,
    stringsAsFactors = FALSE
  )

  # IRD: identity-scale joint delta
  ird     <- ir[i1] - ir[i0]
  var_ird <- ir[i1]^2 * V[i1, i1] + ir[i0]^2 * V[i0, i0] -
    2 * ir[i1] * ir[i0] * V[i1, i0]
  se_ird  <- sqrt(var_ird)
  ird_out <- list(
    est = multiplier * ird,
    lci = multiplier * (ird - z * se_ird),
    uci = multiplier * (ird + z * se_ird),
    se  = multiplier * se_ird
  )

  # IRR: log-ratio delta (== exposure coefficient of the 2-parameter fit)
  ln_irr    <- l[i1] - l[i0]
  var_lnirr <- V[i1, i1] + V[i0, i0] - 2 * V[i1, i0]
  se_lnirr  <- sqrt(var_lnirr)
  irr_out <- list(
    est   = exp(ln_irr),
    lci   = exp(ln_irr - z * se_lnirr),
    uci   = exp(ln_irr + z * se_lnirr),
    ln    = ln_irr,
    ln_se = se_lnirr
  )
  if (any(cells$zero_events)) {
    irr_out[c("lci", "uci")] <- list(NA_real_, NA_real_)
  }

  list(ir = ir_df, ird = ird_out, irr = irr_out)
}


# HELPER: single-cell robust IR CI (Total rows of weighted analyses) ###
#' Log-scale Wald CI for the single (intercept-only) cell of a weighted
#' rate fit; used for the crude weighted Total row (replaces the rounded
#' qgamma remnant).
#' @noRd
wtd_rate_contrasts_total <- function(cell_fit, z, multiplier) {
  l  <- cell_fit$cells$log_rate[1]
  se <- sqrt(cell_fit$vcov[1, 1])
  if (cell_fit$cells$zero_events[1]) {
    return(list(lci = NA_real_, uci = NA_real_, se = NA_real_))
  }
  list(
    lci = multiplier * exp(l - z * se),
    uci = multiplier * exp(l + z * se),
    se  = multiplier * exp(l) * se
  )
}


# HELPER: weighted standardized IR/IRD/IRR from rate cells (stratified) #
#' Joint sandwich/delta contrasts for the WEIGHTED STANDARDIZED analysis
#' (matrix r5/r9/r13-A): per-arm IR_wstd = sum_s W_s * exp(l_as) with robust
#' log-scale CI, IRD_wstd via the joint gradient, IRR_wstd via the log-ratio
#' joint gradient. Stratum shares std_w (named by stratum level) are treated
#' as fixed.
#'
#' @param cell_fit Result of fit_weighted_rate_cells (with strat_var).
#' @param std_w Named numeric vector of standardization weights W_s
#'   (names = stratum levels as character, sum to 1).
#' @param z z-value; multiplier = reporting multiplier.
#' @return List: ir (per arm: exp_val, ir, lci, uci, se), ird, irr
#'   (same shapes as wtd_rate_contrasts).
#' @noRd
wtd_std_rate_contrasts <- function(cell_fit, std_w, z, multiplier) {
  cells <- cell_fit$cells
  V     <- cell_fit$vcov
  n     <- nrow(cells)

  w_s <- std_w[cells$strat_val]              # W_s aligned to cell order
  ir_cell <- exp(cells$log_rate)

  arm_stats <- lapply(c("0", "1"), function(a) {
    sel <- cells$exp_val == a
    est <- sum(w_s[sel] * ir_cell[sel])
    grad <- numeric(n)
    grad[sel] <- w_s[sel] * ir_cell[sel]     # d(est)/d(l_c)
    var <- as.numeric(t(grad) %*% V %*% grad)
    list(est = est, grad = grad, var = var)
  })
  names(arm_stats) <- c("0", "1")

  # A zero-event CELL contributes (numerically) nothing to the standardized
  # estimate or its gradient, mirroring the unweighted standardized variance
  # (and Fay-Feuer), so the arm CI stays defined as long as the standardized
  # estimate itself is positive.
  ir_df <- do.call(rbind, lapply(c("0", "1"), function(a) {
    s <- arm_stats[[a]]
    se_nat <- sqrt(s$var)
    if (s$est > 0) {
      se_log <- se_nat / s$est                # robust log-scale CI (r5-I)
      lci <- exp(log(s$est) - z * se_log)
      uci <- exp(log(s$est) + z * se_log)
    } else {
      lci <- NA_real_; uci <- NA_real_
    }
    data.frame(exp_val = a,
               ir  = multiplier * s$est,
               lci = multiplier * lci,
               uci = multiplier * uci,
               se  = multiplier * se_nat,
               stringsAsFactors = FALSE)
  }))

  # IRD_wstd: joint gradient (+ arm1 cells, - arm0 cells)
  grad_ird <- arm_stats[["1"]]$grad - arm_stats[["0"]]$grad
  ird      <- arm_stats[["1"]]$est - arm_stats[["0"]]$est
  se_ird   <- sqrt(as.numeric(t(grad_ird) %*% V %*% grad_ird))
  ird_out <- list(
    est = multiplier * ird,
    lci = multiplier * (ird - z * se_ird),
    uci = multiplier * (ird + z * se_ird),
    se  = multiplier * se_ird
  )

  # IRR_wstd: log-ratio joint gradient
  e1 <- arm_stats[["1"]]$est
  e0 <- arm_stats[["0"]]$est
  if (e1 > 0 && e0 > 0) {
    grad_ln  <- arm_stats[["1"]]$grad / e1 - arm_stats[["0"]]$grad / e0
    se_lnirr <- sqrt(as.numeric(t(grad_ln) %*% V %*% grad_ln))
    ln_irr   <- log(e1) - log(e0)
    irr_out <- list(
      est   = exp(ln_irr),
      lci   = exp(ln_irr - z * se_lnirr),
      uci   = exp(ln_irr + z * se_lnirr),
      ln    = ln_irr,
      ln_se = se_lnirr
    )
  } else {
    irr_out <- list(est = NA_real_, lci = NA_real_, uci = NA_real_,
                    ln = NA_real_, ln_se = NA_real_)
  }

  list(ir = ir_df, ird = ird_out, irr = irr_out)
}


# HELPER: direct-standardized IRR (unweighted, matrix r11-A) ###########
#' Direct-standardized IRR for the unweighted stratified analysis:
#' IR_a_std = sum_s(W_s * d_as / py_as), Var(IR_a_std) =
#' sum_s(W_s^2 * d_as / py_as^2) (independent Poisson strata), then
#' log-ratio delta: Var(log IRR_std) = Var(IR1)/IR1^2 + Var(IR0)/IR0^2.
#'
#' @param events1,py1 Exposed-arm per-stratum events and person-times.
#' @param events0,py0 Reference-arm per-stratum events and person-times
#'   (same stratum order as arm 1).
#' @param std_w Standardization weights W_s (same order, sum to 1).
#' @param z z-value for the confidence level.
#' @return List: est, lci, uci, ln, ln_se (NA CIs when an arm rate is 0).
#' @noRd
direct_std_irr <- function(events1, py1, events0, py0, std_w, z) {
  ir1 <- sum(std_w * events1 / py1)
  ir0 <- sum(std_w * events0 / py0)
  v1  <- sum(std_w^2 * events1 / py1^2)
  v0  <- sum(std_w^2 * events0 / py0^2)

  if (ir1 > 0 && ir0 > 0) {
    ln_irr <- log(ir1) - log(ir0)
    ln_se  <- sqrt(v1 / ir1^2 + v0 / ir0^2)
    list(est = exp(ln_irr),
         lci = exp(ln_irr - z * ln_se),
         uci = exp(ln_irr + z * ln_se),
         ln = ln_irr, ln_se = ln_se)
  } else {
    list(est = if (ir0 > 0) ir1 / ir0 else NA_real_,
         lci = NA_real_, uci = NA_real_,
         ln = NA_real_, ln_se = NA_real_)
  }
}


# HELPER: pair-level bootstrap resample index ##########################
#' Row indices for one matched-cohort bootstrap replicate: matched sets
#' (match ids) are resampled with replacement and every row of each sampled
#' set is kept together (pair-level resampling; used when if_match_match_id
#' is available).
#'
#' @param match_ids Vector of match-set ids, one per data row.
#' @return Integer vector of row indices for the replicate.
#' @noRd
pair_resample_index <- function(match_ids) {
  ids       <- unique(match_ids)
  sampled   <- sample(ids, size = length(ids), replace = TRUE)
  idx_by_id <- split(seq_along(match_ids), match_ids)
  unlist(idx_by_id[as.character(sampled)], use.names = FALSE)
}
