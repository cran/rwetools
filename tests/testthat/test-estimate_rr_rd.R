# Tests for the v0.4.0 estimate_rr_rd (API v2: data blocks, always-on
# analytical CIs with cloglog Risk/CIF intervals, opt-in fixed-weight
# bootstrap; Fine-Gray moved to estimate_hr_ir).

fx <- readRDS(testthat::test_path("fixtures", "v030_fixtures.rds"))

test_that("estimate_rr_rd returns results for a crude block (KM)", {
  skip_on_cran()
  result <- estimate_rr_rd(
    in_df_crude          = small_df,
    exposure_var         = "exposure",
    outcome_var          = "outcome",
    followuptime_var     = "follow_up_days",
    time_unit            = "days",
    rr_rd_at_timepoint   = 365,
    risk_per_individuals = 1000,
    if_bootstrap_count   = 20,
    if_bootstrap_n_cores = 1,
    verbose              = FALSE
  )

  expect_type(result, "list")
  expect_true("Risk_Estimator" %in% names(result$estimates))
  expect_equal(unique(result$estimates$Risk_Estimator), "KM")
  expect_identical(result$estimates$Analysis, "Crude")
  # analytical columns are always present; bootstrap columns filled here
  expect_true(is.finite(result$estimates$RR_LCI_Analytical))
  expect_true(is.finite(result$estimates$RR_LCI_Boot))
  expect_true("Risk_SE" %in% names(result$cumulative_incidence))
  expect_false("Survival" %in% names(result$cumulative_incidence))
})

test_that("estimate_rr_rd runs a weighted block", {
  skip_on_cran()
  ps_df_rr <- estimate_ps(
    in_df = small_df, exposure_var = "exposure",
    class_vars = class_vars, cont_vars = cont_vars,
    ps_var = "ps", verbose = FALSE
  )
  wt_df_rr <- create_matching_weights(
    in_df = ps_df_rr, exposure_var = "exposure",
    ps_var = "ps", weight_var = "mw_wt", verbose = FALSE
  )

  result <- estimate_rr_rd(
    in_df_crude          = small_df,
    in_df_weight         = wt_df_rr,
    if_weight_weight_var = "mw_wt",
    exposure_var         = "exposure",
    outcome_var          = "outcome",
    followuptime_var     = "follow_up_days",
    time_unit            = "days",
    rr_rd_at_timepoint   = 365,
    risk_per_individuals = 1000,
    verbose              = FALSE
  )
  expect_identical(result$estimates$Analysis, c("Crude", "Weighted"))
})

test_that("fixture equivalence: RR/RD analytical reproduce; Risk CI is cloglog", {
  df <- read.csv(system.file("extdata", "sample_data.csv", package = "rwetools"))
  res <- estimate_rr_rd(
    in_df_crude = df, exposure_var = "exposure", outcome_var = "outcome",
    followuptime_var = "follow_up_days", time_unit = "days",
    rr_rd_at_timepoint = 365, risk_per_individuals = 1000, verbose = FALSE
  )
  fkm <- fx$rr_rd_km$estimates[1, ]   # old "Unweighted" row == new Crude block
  for (cn in c("RiskperN_Ref", "RiskperN_Exp", "RR_Analytical",
               "RR_LCI_Analytical", "RR_UCI_Analytical", "lnRR_SE_Analytical",
               "RDperN_Analytical", "RDperN_LCI_Analytical",
               "RDperN_UCI_Analytical", "RDperN_SE_Analytical")) {
    expect_equal(res$estimates[[cn]], fkm[[cn]], tolerance = 1e-10)
  }

  # Risk arm CI: cloglog == survfit(conf.type = "log-log"), != old log-scale
  fci <- fx$rr_rd_km$cumulative_incidence[1:2, ]
  expect_false(isTRUE(all.equal(res$cumulative_incidence$Risk_LCI,
                                fci$Risk_LCI, tolerance = 1e-6)))
  for (a in 0:1) {
    d_a <- df[df$exposure == a, ]
    sf <- survival::survfit(survival::Surv(follow_up_days, outcome) ~ 1,
                            data = d_a, conf.type = "log-log")
    ss <- summary(sf, times = 365, extend = TRUE)
    row <- res$cumulative_incidence[
      res$cumulative_incidence$Exposure_Value == as.character(a), ]
    expect_equal(row$Risk_LCI, 1 - ss$upper, tolerance = 1e-10)
    expect_equal(row$Risk_UCI, 1 - ss$lower, tolerance = 1e-10)
  }
})

test_that("stratified (standardized) estimates reproduce the fixture", {
  skip_on_cran()
  df <- read.csv(system.file("extdata", "sample_data.csv", package = "rwetools"))
  res <- estimate_rr_rd(
    in_df_crude = df, exposure_var = "exposure", outcome_var = "outcome",
    followuptime_var = "follow_up_days", time_unit = "days",
    rr_rd_at_timepoint = 365, risk_per_individuals = 1000,
    stratification_var = "cat1", verbose = FALSE
  )
  fkms <- fx$rr_rd_km_strat$estimates[1, ]
  expect_identical(res$estimates$Analysis, "Crude (Standardized)")
  for (cn in c("RR_Analytical", "RR_LCI_Analytical", "RDperN_Analytical",
               "RDperN_SE_Analytical")) {
    expect_equal(res$estimates[[cn]], fkms[[cn]], tolerance = 1e-10)
  }
  expect_equal(res$stratum_details$w_k, fx$rr_rd_km_strat$stratum_details$w_k[1:4],
               tolerance = 1e-12)
})

test_that("AJ estimator reproduces the fixture and enforces its argument rules", {
  skip_on_cran()
  df <- read.csv(system.file("extdata", "sample_data.csv", package = "rwetools"))
  res <- estimate_rr_rd(
    in_df_crude = df, exposure_var = "exposure", outcome_var = "outcome",
    followuptime_var = "follow_up_days", time_unit = "days",
    rr_rd_at_timepoint = 365, risk_per_individuals = 1000,
    risk_estimator = "AJ", if_aj_competing_event_var = "competing_event",
    verbose = FALSE
  )
  faj <- fx$rr_rd_aj_fg$estimates[1, ]
  expect_equal(res$estimates$RR_Analytical, faj$RR_Analytical, tolerance = 1e-10)
  expect_equal(res$estimates$RDperN_SE_Analytical, faj$RDperN_SE_Analytical,
               tolerance = 1e-10)
  expect_null(res$subdist_hazard)   # Fine-Gray lives in estimate_hr_ir now

  expect_error(
    estimate_rr_rd(in_df_crude = df, exposure_var = "exposure",
                   outcome_var = "outcome", followuptime_var = "follow_up_days",
                   risk_estimator = "AJ", verbose = FALSE),
    "if_aj_competing_event_var"
  )
  expect_message(
    estimate_rr_rd(in_df_crude = small_df, exposure_var = "exposure",
                   outcome_var = "outcome", followuptime_var = "follow_up_days",
                   risk_estimator = "KM",
                   if_aj_competing_event_var = "competing_event", verbose = TRUE),
    "ignored because risk_estimator"
  )
})

test_that("estimate_rr_rd errors when no data provided", {
  expect_error(
    estimate_rr_rd(
      outcome_var      = "outcome",
      followuptime_var = "follow_up_days",
      verbose          = FALSE
    ),
    "At least one of in_df_crude"
  )
})

test_that("matched block: pair-resample bootstrap differs from row resample", {
  skip_on_cran()
  set.seed(3)
  md <- small_df[sample.int(nrow(small_df), 300), ]
  md$pair_id <- rep(seq_len(150), each = 2)
  run_boot <- function(id) {
    suppressWarnings(suppressMessages(estimate_rr_rd(
      in_df_match = md, if_match_match_id = id,
      exposure_var = "exposure", outcome_var = "outcome",
      followuptime_var = "follow_up_days", rr_rd_at_timepoint = 365,
      if_bootstrap_count = 15, if_bootstrap_n_cores = 1,
      if_bootstrap_seed = 11, verbose = FALSE
    )))
  }
  b_pair <- run_boot("pair_id")
  b_row  <- run_boot(NULL)
  expect_true(is.finite(b_pair$estimates$RR_LCI_Boot))
  expect_false(identical(b_pair$estimates$RR_LCI_Boot,
                         b_row$estimates$RR_LCI_Boot))
  # seed reproducibility + analytical columns unaffected by the bootstrap
  b_pair2 <- run_boot("pair_id")
  expect_identical(b_pair$estimates$RR_LCI_Boot, b_pair2$estimates$RR_LCI_Boot)
  r0 <- estimate_rr_rd(in_df_match = md, if_match_match_id = "pair_id",
                       exposure_var = "exposure", outcome_var = "outcome",
                       followuptime_var = "follow_up_days",
                       rr_rd_at_timepoint = 365, verbose = FALSE)
  expect_equal(b_pair$estimates$RR_LCI_Analytical,
               r0$estimates$RR_LCI_Analytical, tolerance = 1e-12)
})

test_that("Risk on the 0 boundary yields NA cloglog CI plus a message", {
  # timepoint before any event: both arm risks are exactly 0
  expect_message(
    res <- estimate_rr_rd(
      in_df_crude = small_df, exposure_var = "exposure",
      outcome_var = "outcome", followuptime_var = "follow_up_days",
      rr_rd_at_timepoint = 0.5, verbose = FALSE
    ),
    "cloglog CI is undefined"
  )
  expect_true(all(res$cumulative_incidence$Risk == 0))
  expect_true(all(is.na(res$cumulative_incidence$Risk_LCI)))
  expect_true(all(is.na(res$cumulative_incidence$Risk_UCI)))
})

test_that("bootstrap is seed-reproducible on 2 cores as well as 1", {
  skip_on_cran()
  run2 <- function() {
    suppressWarnings(suppressMessages(estimate_rr_rd(
      in_df_crude = small_df, exposure_var = "exposure",
      outcome_var = "outcome", followuptime_var = "follow_up_days",
      rr_rd_at_timepoint = 365,
      if_bootstrap_count = 12, if_bootstrap_n_cores = 2,
      if_bootstrap_seed = 21, verbose = FALSE
    )))
  }
  b1 <- run2(); b2 <- run2()
  expect_identical(b1$estimates$RR_LCI_Boot, b2$estimates$RR_LCI_Boot)
  expect_identical(b1$estimates$RDperN_SE_Boot, b2$estimates$RDperN_SE_Boot)
})
