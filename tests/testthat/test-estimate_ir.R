# Tests for estimate_ir(), split out of the v0.4.0 estimate_hr_ir() at 0.5.0.
# Fixture v030_fixtures.rds holds the v0.3.0 outputs; method-unchanged
# quantities must reproduce exactly, intended changes are asserted as changed.
# The split itself must be behaviour-frozen: every rate estimate below is also
# pinned against the v0.4.0 golden baseline by
# "rwetools audit/v040_baseline/compare_to_v040_baseline.R".

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

test_that("crude block reproduces the v0.3.0 unweighted rate outputs exactly", {
  df <- sample_df()
  res <- estimate_ir(
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
  # estimate_ir returns no hazard estimate (it moved to estimate_hr)
  expect_null(res$hazard_ratios)
  expect_null(res$subdist_hazard)
  expect_null(res$models)
})

test_that("estimate_ir has no stratification argument (0.5.0 contract)", {
  df <- sample_df()
  base <- list(in_df_crude = df, exposure_var = "exposure",
               outcome_var = "outcome",
               followuptime_var = "follow_up_days", verbose = FALSE)
  # neither the 0.4.0 name nor the new estimate_hr name is accepted
  expect_error(do.call(estimate_ir, c(base, list(stratification_var = "cat1"))),
               "unused argument")
  expect_error(do.call(estimate_ir, c(base, list(strata_var = "cat1"))),
               "unused argument")
  fm <- names(formals(estimate_ir))
  expect_false(any(c("stratification_var", "strata_var") %in% fm))
  expect_identical(fm[length(fm)], "verbose")
})

test_that("exp_value/ref_value are honoured on an already-0/1 column (Item 8)", {
  # Up to 0.4.0 the exposure recode was skipped when the column already looked
  # like 0/1, so exp_value = 0 / ref_value = 1 was silently ignored and the
  # WRONG arm was treated as exposed. Ground truth: flipping the two values
  # must invert the IRR and negate the IRD.
  df <- sample_df()
  call <- function(ev, rv, d = df) suppressMessages(estimate_ir(
    in_df_crude = d, exposure_var = "exposure", exp_value = ev, ref_value = rv,
    outcome_var = "outcome", followuptime_var = "follow_up_days",
    time_unit = "days", ir_per_pyears = 1000, verbose = FALSE))

  a <- call(1, 0)
  b <- call(0, 1)

  # hand arithmetic on the raw data
  ev1 <- sum(df$outcome[df$exposure == 1]); ev0 <- sum(df$outcome[df$exposure == 0])
  py1 <- sum(df$follow_up_days[df$exposure == 1]) / 365.25
  py0 <- sum(df$follow_up_days[df$exposure == 0]) / 365.25
  irr <- (ev1 / py1) / (ev0 / py0)
  ird <- 1000 * (ev1 / py1 - ev0 / py0)

  # IRR comes from a Poisson GLM, so it matches the raw ratio to ~1e-8, not
  # to machine precision; IRD is a raw sum and does match exactly.
  expect_equal(a$incidence_rate_ratios$IRR, irr, tolerance = 1e-6)
  expect_equal(a$incidence_rates$IRD_per_Npy_Crude[1], ird, tolerance = 1e-10)
  # the flip must actually flip
  expect_equal(b$incidence_rate_ratios$IRR, 1 / irr, tolerance = 1e-6)
  expect_equal(b$incidence_rates$IRD_per_Npy_Crude[1], -ird, tolerance = 1e-10)
  # per-arm rates swap between the two calls
  ir_exp <- function(r) r$incidence_rates$IR_per_Npy_Crude[
    r$incidence_rates$Exposure_Group == "Exposed"]
  ir_ref <- function(r) r$incidence_rates$IR_per_Npy_Crude[
    r$incidence_rates$Exposure_Group == "Reference"]
  expect_equal(ir_exp(b), ir_ref(a), tolerance = 1e-12)
  expect_equal(ir_ref(b), ir_exp(a), tolerance = 1e-12)

  # a non-0/1 labelled column must give exactly the same two answers
  dl <- df; dl$exposure <- ifelse(df$exposure == 1, "TRT", "REF")
  expect_equal(call("TRT", "REF", dl)$incidence_rates$IRD_per_Npy_Crude[1],
               ird, tolerance = 1e-12)
  expect_equal(call("REF", "TRT", dl)$incidence_rates$IRD_per_Npy_Crude[1],
               -ird, tolerance = 1e-12)
})

test_that("exposure recoding guards: equal values, unmatched rows, empty arm", {
  df <- sample_df()
  expect_error(
    estimate_ir(in_df_crude = df, exposure_var = "exposure",
                exp_value = 1, ref_value = 1, outcome_var = "outcome",
                followuptime_var = "follow_up_days", verbose = FALSE),
    "must differ"
  )
  # a third exposure level is dropped, warned about twice: once as unmatched
  # by exp_value/ref_value, once as an incomplete-case removal
  d3 <- df; d3$exposure[1:25] <- 2
  w <- testthat::capture_warnings(suppressMessages(estimate_ir(
    in_df_crude = d3, exposure_var = "exposure", outcome_var = "outcome",
    followuptime_var = "follow_up_days", verbose = FALSE)))
  expect_true(any(grepl("25 row\\(s\\).*match neither", w)))
  expect_true(any(grepl("25 observation\\(s\\) with missing values", w)))
  # an exp_value that matches nothing leaves no exposed arm -> stop
  expect_error(
    suppressWarnings(suppressMessages(estimate_ir(
      in_df_crude = df, exposure_var = "exposure", exp_value = 7,
      ref_value = 0, outcome_var = "outcome",
      followuptime_var = "follow_up_days", verbose = FALSE))),
    "no exposed .* rows"
  )
})

test_that("weighted block: robust IR CIs", {
  df <- sample_df(); wt <- iptw_df(df)
  res <- estimate_ir(
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
})

test_that("matched block reports its own IR columns", {
  df <- sample_df()
  md <- make_test_matched_pairs(300, seed = 1)
  res <- estimate_ir(
    in_df_crude = df, in_df_match = md, if_match_match_id = "pair_id",
    exposure_var = "exposure", outcome_var = "outcome",
    followuptime_var = "follow_up_days", verbose = FALSE
  )
  expect_true(all(c("IR_per_Npy_Crude", "IR_per_Npy_Matched") %in%
                    names(res$incidence_rates)))
})

test_that("block validation rules (XOR, crude-less, orphans)", {
  df <- sample_df(); wt <- iptw_df(df)
  md <- make_test_matched_pairs(50, seed = 3)
  expect_error(
    estimate_ir(in_df_crude = df, in_df_weight = wt, in_df_match = md,
                if_weight_weight_var = "iptw", exposure_var = "exposure",
                outcome_var = "outcome", followuptime_var = "follow_up_days",
                verbose = FALSE),
    "not both"
  )
  expect_error(
    estimate_ir(in_df_weight = wt, exposure_var = "exposure",
                outcome_var = "outcome", followuptime_var = "follow_up_days",
                verbose = FALSE),
    "if_weight_weight_var"
  )
  expect_error(
    estimate_ir(exposure_var = "exposure", outcome_var = "outcome",
                followuptime_var = "follow_up_days", verbose = FALSE),
    "At least one"
  )
  # crude-less weight-only call is allowed
  res <- estimate_ir(in_df_weight = wt, if_weight_weight_var = "iptw",
                     exposure_var = "exposure", outcome_var = "outcome",
                     followuptime_var = "follow_up_days", verbose = FALSE)
  expect_identical(nrow(res$incidence_rate_ratios), 1L)
  # orphan if_* arguments -> ignore messages
  expect_message(
    estimate_ir(in_df_crude = df, if_weight_weight_var = "iptw",
                exposure_var = "exposure", outcome_var = "outcome",
                followuptime_var = "follow_up_days", verbose = TRUE),
    "if_weight_weight_var is ignored"
  )
  expect_message(
    estimate_ir(in_df_crude = df, if_match_match_id = "pair_id",
                exposure_var = "exposure", outcome_var = "outcome",
                followuptime_var = "follow_up_days", verbose = TRUE),
    "if_match_match_id is ignored"
  )
})

test_that("ir_per_pyears is validated against the allowed multipliers", {
  df <- sample_df()
  expect_error(
    estimate_ir(in_df_crude = df, exposure_var = "exposure",
                outcome_var = "outcome", followuptime_var = "follow_up_days",
                ir_per_pyears = 500, verbose = FALSE),
    "ir_per_pyears must be one of"
  )
})

test_that("time_unit months (PY/12) matches the days pathway", {
  df <- sample_df()
  df_m <- df; df_m$fu_months <- df$follow_up_days / 365.25 * 12
  r_mo <- estimate_ir(in_df_crude = df_m, exposure_var = "exposure",
                      outcome_var = "outcome", followuptime_var = "fu_months",
                      time_unit = "months", verbose = FALSE)
  r_dy <- estimate_ir(in_df_crude = df, exposure_var = "exposure",
                      outcome_var = "outcome", followuptime_var = "follow_up_days",
                      time_unit = "days", verbose = FALSE)
  expect_equal(r_mo$incidence_rates$IR_per_Npy_Crude,
               r_dy$incidence_rates$IR_per_Npy_Crude, tolerance = 1e-10)
})

test_that("bootstrap: opt-in, fixed-weight message, seed-reproducible, analytical unchanged", {
  df <- sample_df()
  expect_message(
    estimate_ir(in_df_crude = df, exposure_var = "exposure",
                outcome_var = "outcome", followuptime_var = "follow_up_days",
                if_bootstrap_seed = 1, verbose = TRUE),
    "if_bootstrap_seed is ignored"
  )
  run_boot <- function() {
    suppressMessages(estimate_ir(
      in_df_crude = df, exposure_var = "exposure", outcome_var = "outcome",
      followuptime_var = "follow_up_days",
      if_bootstrap_count = 20, if_bootstrap_n_cores = 1,
      if_bootstrap_seed = 42, verbose = FALSE
    ))
  }
  expect_message(
    estimate_ir(in_df_crude = df[1:400, ], exposure_var = "exposure",
                outcome_var = "outcome", followuptime_var = "follow_up_days",
                if_bootstrap_count = 5, if_bootstrap_n_cores = 1,
                verbose = FALSE),
    "fixed-weight"
  )
  b1 <- run_boot(); b2 <- run_boot()
  expect_true(all(c("IR_BootLCI_Crude", "IRD_BootSE_Crude") %in%
                    names(b1$incidence_rates)))
  expect_identical(b1$incidence_rates$IR_BootLCI_Crude,
                   b2$incidence_rates$IR_BootLCI_Crude)
  expect_identical(colnames(b1$bootstrap$Crude),
                   c("IR1", "IR0", "IRD", "IRR"))
  # analytical columns unchanged by running the bootstrap
  r0 <- estimate_ir(in_df_crude = df, exposure_var = "exposure",
                    outcome_var = "outcome", followuptime_var = "follow_up_days",
                    verbose = FALSE)
  expect_equal(b1$incidence_rates$IR_LCI_Crude, r0$incidence_rates$IR_LCI_Crude,
               tolerance = 1e-12)
})

test_that("matched pair-resample differs from row resample (same seed)", {
  md <- make_test_matched_pairs(200, seed = 2)
  boot_pair <- suppressWarnings(suppressMessages(estimate_ir(
    in_df_match = md, if_match_match_id = "pair_id",
    exposure_var = "exposure", outcome_var = "outcome",
    followuptime_var = "follow_up_days",
    if_bootstrap_count = 15, if_bootstrap_n_cores = 1,
    if_bootstrap_seed = 7, verbose = FALSE
  )))
  boot_row <- suppressWarnings(suppressMessages(estimate_ir(
    in_df_match = md,
    exposure_var = "exposure", outcome_var = "outcome",
    followuptime_var = "follow_up_days",
    if_bootstrap_count = 15, if_bootstrap_n_cores = 1,
    if_bootstrap_seed = 7, verbose = FALSE
  )))
  expect_false(identical(boot_pair$incidence_rates$IR_BootLCI_Matched,
                         boot_row$incidence_rates$IR_BootLCI_Matched))
})
