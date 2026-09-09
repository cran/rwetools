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
