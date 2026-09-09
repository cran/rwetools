# Effect-estimation statistical engines ###############################
# Lowest layer of the effect-estimation stack: pure statistical engines
# (CI formulas, weighted rate-cell fits, delta-method contrasts, resampling)
# that return raw estimates and know nothing about the reporting format.
# Called by the calculators in helpers_effect_measures.R and the public
# functions in effect_measures.R. (Introduced by the 0.4.0 API v2 redesign,
# EFFECT_CI_V2_PLAN.)


# HELPER: cloglog CI for a 0-1 risk / CIF ##############################
#' Complementary log-log confidence interval for a cumulative risk (KM
#' 1 - S(t)) or CIF (Aalen-Johansen), applied to the FINAL estimate with its
#' combined SE (matrix r22-r29): g = log(-log(1 - R)),
#' SE_g = SE_R / ((1 - R) * (-log(1 - R))), CI = 1 - exp(-exp(g +/- z*SE_g)).
#' For a crude KM risk this reproduces survfit(conf.type = "log-log") exactly.
#' Boundary policy: Risk on {0, 1} (cloglog undefined) returns NA bounds; the
#' caller emits the bootstrap-recommending message.
#'
#' @param risk Numeric vector of risks/CIFs in `[0, 1]`.
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
#' Fit the saturated weighted quasi-Poisson rate model over the analysis arms
#' and return the cell log-rates with their joint robust (design-based
#' sandwich) variance-covariance matrix (survey::svyglm; Lumley 2004,
#' Binder 1983). This is the shared engine for the weighted IR / IRD / IRR
#' analytical methods (matrix r4/r8/r12-A): every downstream quantity is a
#' delta-method contrast on these cell log-rates.
#'
#' One cell per exposure arm; exp(log-rate) of a saturated cell equals the
#' weighted event sum / weighted person-time of that cell. A cell with zero
#' weighted events keeps the fit alive (large negative coefficient); its own
#' CI is set to NA downstream.
#'
#' Two design features may be declared, independently:
#' \code{wt_var} attaches analysis weights (the IPTW / weighted block), and
#' \code{cluster_var} declares the sampling unit (a matched-set id), which
#' turns the sandwich into a cluster-robust one. With neither, the fit is an
#' ordinary unweighted saturated quasi-Poisson and the sandwich is the plain
#' heteroskedasticity-robust one.
#'
#' @param df Data frame (already recoded exposure 0/1, person-time column).
#' @param exp_var,out_var,py_var Column names (character).
#' @param wt_var Analysis-weight column, or NULL for unweighted.
#' @param cluster_var Cluster-id column (e.g. a matched-set id), or NULL for
#'   independent rows.
#' @return List: cells (data.frame with exp_val, log_rate, n_events_w, py_w,
#'   zero_events flag), vcov (joint robust vcov, cells in row order of
#'   `cells`), or NULL invisibly if the fit fails.
#' @noRd
fit_weighted_rate_cells <- function(df, exp_var, out_var, py_var,
                                    wt_var = NULL, cluster_var = NULL) {
  work <- df
  work$.cell <- factor(work[[exp_var]])
  w <- if (is.null(wt_var)) rep(1, nrow(work)) else work[[wt_var]]
  cell_levels <- levels(work$.cell)

  # A single-level .cell (the intercept-only "total row" call, where exp_var is
  # a constant column) cannot go through `~ 0 + .cell`: a one-level factor has
  # no contrasts, and model.matrix() errors with "contrasts can be applied only
  # to factors with 2 or more levels". Up to 0.4.0 that error was swallowed by
  # the caller's tryCatch, so the weighted block's Total-row IR CI was silently
  # NA in every release. Use a plain intercept in that case; the coefficient is
  # the same log-rate and the sandwich the same variance.
  single_cell <- length(cell_levels) == 1L
  fml <- if (single_cell) {
    stats::as.formula(paste0(out_var, " ~ 1 + offset(log(", py_var, "))"))
  } else {
    stats::as.formula(paste0(out_var, " ~ 0 + .cell + offset(log(", py_var, "))"))
  }

  fit <- tryCatch({
    ids_f <- if (is.null(cluster_var)) {
      ~1
    } else {
      stats::as.formula(paste0("~`", cluster_var, "`"))
    }
    # an explicit unit-weight column keeps svydesign from warning about
    # "no weights or probabilities supplied"; it is numerically the same design
    wt_f <- if (is.null(wt_var)) {
      work$.unit_w_ <- 1
      ~.unit_w_
    } else {
      stats::as.formula(paste0("~`", wt_var, "`"))
    }
    des <- survey::svydesign(ids = ids_f, weights = wt_f, data = work)
    survey::svyglm(fml, design = des, family = stats::quasipoisson())
  }, error = function(e) NULL)
  if (is.null(fit)) return(NULL)

  cf <- stats::coef(fit)
  V  <- stats::vcov(fit)

  coef_names <- if (single_cell) "(Intercept)" else paste0(".cell", cell_levels)
  # Guard against dropped coefficients (should not happen for saturated cells)
  if (!all(coef_names %in% names(cf))) return(NULL)
  cf <- cf[coef_names]
  V  <- V[coef_names, coef_names, drop = FALSE]

  # Weighted event sums / person-time per cell (for reporting & zero flags)
  ne_w <- vapply(split(work[[out_var]] * w, work$.cell),
                 function(x) sum(x, na.rm = TRUE), numeric(1))[cell_levels]
  py_w <- vapply(split(work[[py_var]] * w, work$.cell),
                 function(x) sum(x, na.rm = TRUE), numeric(1))[cell_levels]

  cells <- data.frame(
    exp_val     = cell_levels,
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
#' @param cell_fit Result of fit_weighted_rate_cells.
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
