# Unit tests for the statistical engines in R/helpers_effect_engines.R.
# Internal helpers are exercised via rwetools::: .

# --- shared toy data ---------------------------------------------------
make_toy_df <- function(n = 400, seed = 42) {
  set.seed(seed)
  data.frame(
    exposure     = rep(c(0L, 1L), each = n / 2),
    outcome      = rbinom(n, 1, 0.25),
    person_years = stats::rgamma(n, shape = 2, scale = 1.5),
    iptw         = stats::runif(n, 0.5, 3),
    stratum      = sample(c("a", "b", "c"), n, replace = TRUE),
    stringsAsFactors = FALSE
  )
}

# Fay-Feuer ==============================================================

test_that("fay_feuer_ci collapses to Garwood qgamma for a single stratum", {
  for (d in c(1, 5, 40)) {
    for (py in c(10, 123.4)) {
      for (cl in c(0.90, 0.95, 0.99)) {
        ff <- rwetools:::fay_feuer_ci(events = d, py = py, std_w = 1,
                                      conf_level = cl)
        a <- 1 - cl
        expect_equal(ff$lci, stats::qgamma(a / 2, d) / py, tolerance = 1e-12)
        expect_equal(ff$uci, stats::qgamma(1 - a / 2, d + 1) / py,
                     tolerance = 1e-12)
      }
    }
  }
})

test_that("fay_feuer_ci handles zero events (single and multi-stratum)", {
  ff0 <- rwetools:::fay_feuer_ci(events = 0, py = 50, std_w = 1)
  expect_identical(ff0$lci, 0)
  expect_equal(ff0$uci, stats::qgamma(0.975, 1) / 50, tolerance = 1e-12)

  # one empty stratum among several: finite interval containing the rate
  ff <- rwetools:::fay_feuer_ci(events = c(4, 0, 7), py = c(100, 80, 150),
                                std_w = c(0.3, 0.2, 0.5))
  expect_true(ff$lci > 0 && ff$uci > ff$lci)
  expect_true(ff$lci < ff$rate && ff$rate < ff$uci)
})

test_that("multi-stratum Fay-Feuer differs from the log-normal CI", {
  ev <- c(3, 12); py <- c(60, 400); w <- c(0.35, 0.65)
  ff <- rwetools:::fay_feuer_ci(ev, py, w)
  se_log <- ff$se / ff$rate
  ln_lci <- exp(log(ff$rate) - stats::qnorm(0.975) * se_log)
  expect_false(isTRUE(all.equal(ff$lci, ln_lci, tolerance = 1e-6)))
})

# cloglog ================================================================

test_that("cloglog_risk_ci reproduces survfit conf.type = 'log-log' exactly", {
  df <- make_toy_df(seed = 7)
  df$time <- stats::rexp(nrow(df), rate = 0.4)
  fit <- survival::survfit(survival::Surv(time, outcome) ~ 1, data = df,
                           conf.type = "log-log")
  tp <- stats::median(df$time)
  s  <- summary(fit, times = tp, extend = TRUE)

  risk <- 1 - s$surv
  ci <- rwetools:::cloglog_risk_ci(risk = risk, se = s$std.err,
                                   z = stats::qnorm(0.975))
  # survival CI for S maps to risk CI with bounds swapped
  expect_equal(ci$lci, 1 - s$upper, tolerance = 1e-10)
  expect_equal(ci$uci, 1 - s$lower, tolerance = 1e-10)
})

test_that("cloglog_risk_ci stays inside (0, 1) and brackets the risk", {
  set.seed(11)
  risk <- stats::runif(50, 0.02, 0.90)
  se   <- stats::runif(50, 0.001, 0.08)
  ci <- rwetools:::cloglog_risk_ci(risk, se, z = stats::qnorm(0.975))
  expect_true(all(ci$lci > 0 & ci$uci < 1))
  expect_true(all(ci$lci < risk & risk < ci$uci))

  # extreme risk/SE combos may saturate to the boundary by underflow, but
  # never escape [0,1] and never invert
  rx <- c(0.001, 0.995); sx <- c(0.2, 0.2)
  cx <- rwetools:::cloglog_risk_ci(rx, sx, z = 1.96)
  expect_true(all(cx$lci >= 0 & cx$uci <= 1))
  expect_true(all(cx$lci <= rx & rx <= cx$uci))
})

test_that("cloglog_risk_ci returns NA on the 0/1 boundary and for NA input", {
  ci <- rwetools:::cloglog_risk_ci(c(0, 1, NA, 0.3), c(0.01, 0.01, 0.01, NA),
                                   z = 1.96)
  expect_true(all(is.na(ci$lci[1:4])))
  expect_true(all(is.na(ci$uci[1:4])))
})

# joint-sandwich rate-cell engine =======================================

test_that("saturated cell log-rates equal weighted sums, IRR matches svyglm", {
  df <- make_toy_df(seed = 42)
  fitc <- rwetools:::fit_weighted_rate_cells(
    df, exp_var = "exposure", out_var = "outcome",
    py_var = "person_years", wt_var = "iptw"
  )
  expect_false(is.null(fitc))

  # exp(log-rate) == weighted events / weighted person-time per cell
  for (a in c("0", "1")) {
    sel <- df$exposure == as.integer(a)
    rate_hand <- sum(df$outcome[sel] * df$iptw[sel]) /
      sum(df$person_years[sel] * df$iptw[sel])
    rate_fit <- exp(fitc$cells$log_rate[fitc$cells$exp_val == a])
    expect_equal(rate_fit, rate_hand, tolerance = 1e-8)
  }

  # log-ratio delta == exposure coefficient of the 2-parameter svyglm
  ctr <- rwetools:::wtd_rate_contrasts(fitc, z = stats::qnorm(0.975),
                                       multiplier = 1)
  des <- survey::svydesign(ids = ~1, weights = ~iptw, data = df)
  ref <- survey::svyglm(outcome ~ exposure + offset(log(person_years)),
                        design = des, family = stats::quasipoisson())
  expect_equal(ctr$irr$ln, unname(stats::coef(ref)["exposure"]),
               tolerance = 1e-8)
  expect_equal(ctr$irr$ln_se,
               unname(sqrt(diag(stats::vcov(ref))["exposure"])),
               tolerance = 1e-6)
})

test_that("weighted IRD joint delta matches hand-computed gradient algebra", {
  df <- make_toy_df(seed = 3)
  fitc <- rwetools:::fit_weighted_rate_cells(
    df, exp_var = "exposure", out_var = "outcome",
    py_var = "person_years", wt_var = "iptw"
  )
  ctr <- rwetools:::wtd_rate_contrasts(fitc, z = stats::qnorm(0.975),
                                       multiplier = 1000)
  i0 <- which(fitc$cells$exp_val == "0")
  i1 <- which(fitc$cells$exp_val == "1")
  r  <- exp(fitc$cells$log_rate)
  V  <- fitc$vcov
  var_hand <- r[i1]^2 * V[i1, i1] + r[i0]^2 * V[i0, i0] -
    2 * r[i1] * r[i0] * V[i1, i0]
  expect_equal(ctr$ird$se, 1000 * sqrt(var_hand), tolerance = 1e-10)
  expect_equal(ctr$ird$est, 1000 * (r[i1] - r[i0]), tolerance = 1e-10)
})

test_that("stratified engine matches hand-standardized rates (clean data)", {
  df <- make_toy_df(seed = 42)
  fitc <- rwetools:::fit_weighted_rate_cells(
    df, exp_var = "exposure", out_var = "outcome",
    py_var = "person_years", wt_var = "iptw", strat_var = "stratum"
  )
  expect_false(is.null(fitc))
  expect_identical(nrow(fitc$cells), 6L)

  # standardization weights: weighted person-time share per stratum
  pt_s <- vapply(split(df$person_years * df$iptw, df$stratum), sum, numeric(1))
  w_s  <- pt_s / sum(pt_s)

  ctr <- rwetools:::wtd_std_rate_contrasts(fitc, std_w = w_s,
                                           z = stats::qnorm(0.975),
                                           multiplier = 1000)
  # Point estimates: hand-standardize the per-cell weighted rates
  rate_cell <- function(a, s) {
    sel <- df$exposure == a & df$stratum == s
    sum(df$outcome[sel] * df$iptw[sel]) /
      sum(df$person_years[sel] * df$iptw[sel])
  }
  for (a in c("0", "1")) {
    ir_hand <- 1000 * sum(w_s * vapply(names(w_s), function(s)
      rate_cell(as.integer(a), s), numeric(1)))
    expect_equal(ctr$ir$ir[ctr$ir$exp_val == a], ir_hand, tolerance = 1e-8)
  }
  expect_true(all(is.finite(ctr$ir$lci)) && all(ctr$ir$lci < ctr$ir$ir) &&
                all(ctr$ir$ir < ctr$ir$uci))
  expect_true(ctr$irr$lci < ctr$irr$est && ctr$irr$est < ctr$irr$uci)
})

test_that("stratified engine survives a 0-event cell (CI stays defined)", {
  df <- make_toy_df(seed = 5)
  # Force a zero-event cell: stratum "c" exposed loses all its events
  df$outcome[df$stratum == "c" & df$exposure == 1L] <- 0L

  fitc <- suppressWarnings(rwetools:::fit_weighted_rate_cells(
    df, exp_var = "exposure", out_var = "outcome",
    py_var = "person_years", wt_var = "iptw", strat_var = "stratum"
  ))
  expect_false(is.null(fitc))
  expect_identical(nrow(fitc$cells), 6L)
  expect_true(any(fitc$cells$zero_events))

  pt_s <- vapply(split(df$person_years * df$iptw, df$stratum), sum, numeric(1))
  w_s  <- pt_s / sum(pt_s)
  ctr <- rwetools:::wtd_std_rate_contrasts(fitc, std_w = w_s,
                                           z = stats::qnorm(0.975),
                                           multiplier = 1000)

  # A zero-event CELL contributes ~nothing to the standardized estimate and
  # gradient (mirroring the unweighted standardized variance and Fay-Feuer),
  # so the standardized arm CI remains defined and ordered, and the IRD/IRR
  # contrasts stay finite.
  expect_true(is.finite(ctr$ird$se))
  a1 <- ctr$ir[ctr$ir$exp_val == "1", ]
  expect_true(is.finite(a1$lci) && is.finite(a1$uci))
  expect_true(a1$lci < a1$ir && a1$ir < a1$uci)
  expect_true(is.finite(ctr$irr$ln_se))
})

# direct-standardized IRR (unweighted r11-A) ============================

test_that("direct_std_irr matches hand calculation", {
  ev1 <- c(5, 9);  py1 <- c(120, 300)
  ev0 <- c(3, 14); py0 <- c(100, 420)
  w   <- c(0.4, 0.6)
  z   <- stats::qnorm(0.975)

  res <- rwetools:::direct_std_irr(ev1, py1, ev0, py0, w, z)

  ir1 <- sum(w * ev1 / py1); ir0 <- sum(w * ev0 / py0)
  v1  <- sum(w^2 * ev1 / py1^2); v0 <- sum(w^2 * ev0 / py0^2)
  ln_se <- sqrt(v1 / ir1^2 + v0 / ir0^2)

  expect_equal(res$est, ir1 / ir0, tolerance = 1e-12)
  expect_equal(res$ln_se, ln_se, tolerance = 1e-12)
  expect_equal(res$lci, exp(log(ir1 / ir0) - z * ln_se), tolerance = 1e-12)
})

test_that("direct_std_irr returns NA CIs when an arm has no events", {
  res <- rwetools:::direct_std_irr(c(0, 0), c(10, 20), c(2, 3), c(15, 25),
                                   c(0.5, 0.5), z = 1.96)
  expect_true(is.na(res$lci) && is.na(res$uci))
})

# pair-level resample ====================================================

test_that("pair_resample_index keeps matched sets intact and is seedable", {
  match_ids <- rep(seq_len(50), each = 2)   # 50 pairs of 2

  set.seed(99)
  idx <- rwetools:::pair_resample_index(match_ids)
  expect_identical(length(idx), 100L)
  # every sampled set contributes complete pairs
  expect_true(all(table(match_ids[idx]) %% 2 == 0))

  set.seed(99)
  idx2 <- rwetools:::pair_resample_index(match_ids)
  expect_identical(idx, idx2)

  # differs from a row resample of the same seed
  set.seed(99)
  row_idx <- sample.int(100, 100, replace = TRUE)
  expect_false(identical(sort(idx), sort(row_idx)))
})
