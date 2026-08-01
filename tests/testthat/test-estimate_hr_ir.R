# Tests for the v0.4.0 estimate_hr_ir (API v2: data blocks, fixed-rule
# analytical CIs, Fine-Gray move-in, opt-in fixed-weight bootstrap).
# Fixture v030_fixtures.rds holds the v0.3.0 outputs; method-unchanged
# quantities must reproduce exactly, intended changes are asserted as changed.

fixture_path <- testthat::test_path("fixtures", "v030_fixtures.rds")
fx <- readRDS(fixture_path)

sample_df <- function() {
  read.csv(system.file("extdata", "sample_data.csv", package = "rwetools"))
}
iptw_df <- function(df) {
  ps <- estimate_ps(in_df = df, exposure_var = "exposure",
                    class_vars = paste0("cat", 1:4),
                    cont_vars = paste0("cont", 1:7),
                    ps_var = "ps", verbose = FALSE)
  create_iptw(in_df = ps, exposure_var = "exposure", ps_var = "ps",
              weight_var = "iptw", verbose = FALSE)
}

test_that("crude block reproduces the v0.3.0 unweighted outputs exactly", {
  df <- sample_df()
  res <- estimate_hr_ir(
    in_df_crude = df, exposure_var = "exposure", outcome_var = "outcome",
    followuptime_var = "follow_up_days", time_unit = "days",
    ir_per_pyears = 1000, verbose = FALSE
  )
  fir <- fx$hr_ir$incidence_rates
  for (b in c("N_Subjects", "N_Events", "PersonYears", "IR_per_Npy",
              "IR_LCI", "IR_UCI", "IR_SE", "IRD_per_Npy", "IRD_LCI",
              "IRD_UCI", "IRD_SE")) {
    expect_equal(res$incidence_rates[[paste0(b, "_Crude")]],
                 fir[[paste0(b, "_Unwt")]], tolerance = 1e-12)
  }
  # HR: crude block uses the model SE (old unweighted "auto" rule)
  fhr <- fx$hr_ir$hazard_ratios
  expect_equal(res$hazard_ratios$HR, fhr$HR[1], tolerance = 1e-12)
  expect_equal(res$hazard_ratios$lnHR_SE, fhr$lnHR_SE[1], tolerance = 1e-12)
  expect_equal(res$hazard_ratios$HR_LCI, fhr$HR_LCI[1], tolerance = 1e-12)
})

test_that("crude HR equals a direct model-SE coxph fit (fixed rule)", {
  df <- sample_df()
  res <- estimate_hr_ir(
    in_df_crude = df, exposure_var = "exposure", outcome_var = "outcome",
    followuptime_var = "follow_up_days", verbose = FALSE
  )
  m <- survival::coxph(survival::Surv(follow_up_days, outcome) ~ exposure,
                       data = df, robust = FALSE)
  sm <- summary(m)
  expect_equal(res$hazard_ratios$HR, unname(sm$conf.int[, "exp(coef)"]))
  expect_equal(res$hazard_ratios$lnHR_SE, unname(sm$coefficients[, "se(coef)"]))
})

test_that("weighted block: robust HR rule and robust IR CIs", {
  df <- sample_df(); wt <- iptw_df(df)
  res <- estimate_hr_ir(
    in_df_crude = df, in_df_weight = wt, if_weight_weight_var = "iptw",
    exposure_var = "exposure", outcome_var = "outcome",
    followuptime_var = "follow_up_days", verbose = FALSE
  )
  fir <- fx$hr_ir$incidence_rates
  # points reproduce; robust arm CIs are WIDER than the old naive qgamma (P1)
  expect_equal(res$incidence_rates$IR_per_Npy_Weighted,
               fir$IR_per_Npy_Wted, tolerance = 1e-12)
  expect_true(all(
    (res$incidence_rates$IR_UCI_Weighted[1:2] - res$incidence_rates$IR_LCI_Weighted[1:2]) >
      (fir$IR_UCI_Wted[1:2] - fir$IR_LCI_Wted[1:2])
  ))
  # weighted HR: robust SE (== old "auto" weighted rule, fixture-exact)
  fhr <- fx$hr_ir$hazard_ratios
  expect_equal(res$hazard_ratios$HR[2], fhr$HR[2], tolerance = 1e-12)
  expect_equal(res$hazard_ratios$lnHR_SE[2], fhr$lnHR_SE[2], tolerance = 1e-12)
  expect_match(res$hazard_ratios$Analysis[2], "robust")
})

test_that("stratified: HRs reproduce; standardized arm CI is now Fay-Feuer", {
  df <- sample_df(); wt <- iptw_df(df)
  res <- estimate_hr_ir(
    in_df_crude = df, in_df_weight = wt, if_weight_weight_var = "iptw",
    exposure_var = "exposure", outcome_var = "outcome",
    followuptime_var = "follow_up_days", stratification_var = "cat1",
    verbose = FALSE
  )
  firs <- fx$hr_ir_strat$incidence_rates
  fhrs <- fx$hr_ir_strat$hazard_ratios
  expect_equal(res$hazard_ratios$HR, fhrs$HR, tolerance = 1e-12)
  expect_equal(res$hazard_ratios$lnHR_SE, fhrs$lnHR_SE, tolerance = 1e-12)
  expect_equal(res$incidence_rates$IR_per_Npy_Crude,
               firs$IR_per_Npy_Unwt, tolerance = 1e-12)
  expect_equal(res$incidence_rates$IRD_LCI_Crude,
               firs$IRD_LCI_Unwt, tolerance = 1e-12)  # IRD delta unchanged
  # arm CI intended change: FF != old log-normal
  expect_false(isTRUE(all.equal(res$incidence_rates$IR_LCI_Crude[1:2],
                                firs$IR_LCI_Unwt[1:2], tolerance = 1e-6)))
})

test_that("Fine-Gray moved in: sHR equals the old estimate_rr_rd FG output", {
  df <- sample_df(); wt <- iptw_df(df)
  res <- estimate_hr_ir(
    in_df_crude = df, in_df_weight = wt, if_weight_weight_var = "iptw",
    exposure_var = "exposure", outcome_var = "outcome",
    followuptime_var = "follow_up_days",
    hr_model = "Fine-Gray", if_fg_competing_event_var = "competing_event",
    verbose = FALSE
  )
  ffg <- fx$rr_rd_aj_fg$subdist_hazard
  expect_equal(res$subdist_hazard$SHR, ffg$SHR, tolerance = 1e-10)
  expect_equal(res$subdist_hazard$lnSHR_SE, ffg$lnSHR_SE, tolerance = 1e-10)
  expect_null(res$hazard_ratios)
  # Fine-Gray requires the competing-event variable
  expect_error(
    estimate_hr_ir(in_df_crude = df, exposure_var = "exposure",
                   outcome_var = "outcome", followuptime_var = "follow_up_days",
                   hr_model = "Fine-Gray", verbose = FALSE),
    "if_fg_competing_event_var"
  )
  # Cox + competing var -> ignore message
  expect_message(
    estimate_hr_ir(in_df_crude = df, exposure_var = "exposure",
                   outcome_var = "outcome", followuptime_var = "follow_up_days",
                   hr_model = "Cox", if_fg_competing_event_var = "competing_event",
                   verbose = TRUE),
    "ignored"
  )
})

test_that("matched block: robust + match-id cluster; pair rules enforced", {
  df <- sample_df()
  set.seed(1)
  md <- df[sample.int(nrow(df), 600), ]
  md$pair_id <- rep(seq_len(300), each = 2)
  res <- estimate_hr_ir(
    in_df_crude = df, in_df_match = md, if_match_match_id = "pair_id",
    exposure_var = "exposure", outcome_var = "outcome",
    followuptime_var = "follow_up_days", verbose = FALSE
  )
  expect_match(res$hazard_ratios$Analysis[2], "match-id cluster")
  expect_true(all(c("IR_per_Npy_Crude", "IR_per_Npy_Matched") %in%
                    names(res$incidence_rates)))
  # matched HR equals a direct robust cluster fit
  m <- survival::coxph(
    survival::Surv(follow_up_days, outcome) ~ exposure + cluster(pair_id),
    data = md, robust = TRUE
  )
  expect_equal(res$hazard_ratios$HR[2],
               unname(exp(stats::coef(m)["exposure"])), tolerance = 1e-12)
  expect_equal(res$hazard_ratios$lnHR_SE[2],
               unname(summary(m)$coefficients[, "robust se"]), tolerance = 1e-12)
})

test_that("block validation rules (XOR, crude-less, orphans)", {
  df <- sample_df(); wt <- iptw_df(df)
  md <- df[1:100, ]; md$pair_id <- rep(1:50, each = 2)
  expect_error(
    estimate_hr_ir(in_df_crude = df, in_df_weight = wt, in_df_match = md,
                   if_weight_weight_var = "iptw", exposure_var = "exposure",
                   outcome_var = "outcome", followuptime_var = "follow_up_days",
                   verbose = FALSE),
    "not both"
  )
  expect_error(
    estimate_hr_ir(in_df_weight = wt, exposure_var = "exposure",
                   outcome_var = "outcome", followuptime_var = "follow_up_days",
                   verbose = FALSE),
    "if_weight_weight_var"
  )
  expect_error(
    estimate_hr_ir(exposure_var = "exposure", outcome_var = "outcome",
                   followuptime_var = "follow_up_days", verbose = FALSE),
    "At least one"
  )
  # crude-less weight-only call is allowed
  res <- estimate_hr_ir(in_df_weight = wt, if_weight_weight_var = "iptw",
                        exposure_var = "exposure", outcome_var = "outcome",
                        followuptime_var = "follow_up_days", verbose = FALSE)
  expect_identical(nrow(res$hazard_ratios), 1L)
  # orphan if_* arguments -> ignore messages
  expect_message(
    estimate_hr_ir(in_df_crude = df, if_weight_weight_var = "iptw",
                   exposure_var = "exposure", outcome_var = "outcome",
                   followuptime_var = "follow_up_days", verbose = TRUE),
    "if_weight_weight_var is ignored"
  )
  expect_message(
    estimate_hr_ir(in_df_crude = df, if_match_match_id = "pair_id",
                   exposure_var = "exposure", outcome_var = "outcome",
                   followuptime_var = "follow_up_days", verbose = TRUE),
    "if_match_match_id is ignored"
  )
})

test_that("time_unit months (PY/12) matches the days pathway", {
  df <- sample_df()
  df_m <- df; df_m$fu_months <- df$follow_up_days / 365.25 * 12
  r_mo <- estimate_hr_ir(in_df_crude = df_m, exposure_var = "exposure",
                         outcome_var = "outcome", followuptime_var = "fu_months",
                         time_unit = "months", verbose = FALSE)
  r_dy <- estimate_hr_ir(in_df_crude = df, exposure_var = "exposure",
                         outcome_var = "outcome", followuptime_var = "follow_up_days",
                         time_unit = "days", verbose = FALSE)
  expect_equal(r_mo$incidence_rates$IR_per_Npy_Crude,
               r_dy$incidence_rates$IR_per_Npy_Crude, tolerance = 1e-10)
})

test_that("bootstrap: opt-in, fixed-weight message, seed-reproducible, analytical unchanged", {
  df <- sample_df()
  expect_message(
    estimate_hr_ir(in_df_crude = df, exposure_var = "exposure",
                   outcome_var = "outcome", followuptime_var = "follow_up_days",
                   if_bootstrap_seed = 1, verbose = TRUE),
    "if_bootstrap_seed is ignored"
  )
  run_boot <- function() {
    suppressMessages(estimate_hr_ir(
      in_df_crude = df, exposure_var = "exposure", outcome_var = "outcome",
      followuptime_var = "follow_up_days",
      if_bootstrap_count = 20, if_bootstrap_n_cores = 1,
      if_bootstrap_seed = 42, verbose = FALSE
    ))
  }
  expect_message(
    estimate_hr_ir(in_df_crude = df[1:400, ], exposure_var = "exposure",
                   outcome_var = "outcome", followuptime_var = "follow_up_days",
                   if_bootstrap_count = 5, if_bootstrap_n_cores = 1,
                   verbose = FALSE),
    "fixed-weight"
  )
  b1 <- run_boot(); b2 <- run_boot()
  expect_true(all(c("IR_BootLCI_Crude", "IRD_BootSE_Crude") %in%
                    names(b1$incidence_rates)))
  expect_identical(b1$hazard_ratios$HR_BootLCI, b2$hazard_ratios$HR_BootLCI)
  expect_identical(b1$incidence_rates$IR_BootLCI_Crude,
                   b2$incidence_rates$IR_BootLCI_Crude)
  # analytical columns unchanged by running the bootstrap
  r0 <- estimate_hr_ir(in_df_crude = df, exposure_var = "exposure",
                       outcome_var = "outcome", followuptime_var = "follow_up_days",
                       verbose = FALSE)
  expect_equal(b1$incidence_rates$IR_LCI_Crude, r0$incidence_rates$IR_LCI_Crude,
               tolerance = 1e-12)
  expect_equal(b1$hazard_ratios$HR, r0$hazard_ratios$HR, tolerance = 1e-12)
})

test_that("matched pair-resample differs from row resample (same seed)", {
  df <- sample_df()
  set.seed(2)
  md <- df[sample.int(nrow(df), 400), ]
  md$pair_id <- rep(seq_len(200), each = 2)
  # small resamples can produce degenerate Cox fits ("coefficient may be
  # infinite"); the assertion is about resample indices, not fit quality
  boot_pair <- suppressWarnings(suppressMessages(estimate_hr_ir(
    in_df_match = md, if_match_match_id = "pair_id",
    exposure_var = "exposure", outcome_var = "outcome",
    followuptime_var = "follow_up_days",
    if_bootstrap_count = 15, if_bootstrap_n_cores = 1,
    if_bootstrap_seed = 7, verbose = FALSE
  )))
  boot_row <- suppressWarnings(suppressMessages(estimate_hr_ir(
    in_df_match = md,
    exposure_var = "exposure", outcome_var = "outcome",
    followuptime_var = "follow_up_days",
    if_bootstrap_count = 15, if_bootstrap_n_cores = 1,
    if_bootstrap_seed = 7, verbose = FALSE
  )))
  expect_false(identical(boot_pair$incidence_rates$IR_BootLCI_Matched,
                         boot_row$incidence_rates$IR_BootLCI_Matched))
})

test_that("weighted IR robust CI width vs naive qgamma matches Phase C (~1.3x)", {
  df <- sample_df(); wt <- iptw_df(df)
  res <- estimate_hr_ir(
    in_df_crude = df, in_df_weight = wt, if_weight_weight_var = "iptw",
    exposure_var = "exposure", outcome_var = "outcome",
    followuptime_var = "follow_up_days", verbose = FALSE
  )
  fir <- fx$hr_ir$incidence_rates
  ratio <- (res$incidence_rates$IR_UCI_Weighted[1:2] -
              res$incidence_rates$IR_LCI_Weighted[1:2]) /
    (fir$IR_UCI_Wted[1:2] - fir$IR_LCI_Wted[1:2])
  # P1: the naive interval under-covered; robust arm intervals are wider by a
  # data-dependent factor (~1.26 and ~1.75 on the sample data). Generous
  # window to avoid platform flakiness while pinning the direction and order.
  expect_true(all(ratio > 1.1 & ratio < 2.5))
})
