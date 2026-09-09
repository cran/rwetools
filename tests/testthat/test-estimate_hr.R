# Tests for estimate_hr(), split out of the v0.4.0 estimate_hr_ir() at 0.5.0.
# Fixture v030_fixtures.rds holds the v0.3.0 outputs; method-unchanged
# quantities must reproduce exactly, intended changes are asserted as changed.
# The split itself must be behaviour-frozen: every hazard estimate below is
# also pinned against the v0.4.0 golden baseline by
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

test_that("crude block reproduces the v0.3.0 hazard output exactly", {
  df <- sample_df()
  res <- estimate_hr(
    in_df_crude = df, exposure_var = "exposure", outcome_var = "outcome",
    followuptime_var = "follow_up_days", verbose = FALSE
  )
  # crude block uses the model SE (old unweighted "auto" rule)
  fhr <- fx$hr_ir$hazard_ratios
  expect_equal(res$hazard_ratios$HR, fhr$HR[1], tolerance = 1e-12)
  expect_equal(res$hazard_ratios$lnHR_SE, fhr$lnHR_SE[1], tolerance = 1e-12)
  expect_equal(res$hazard_ratios$HR_LCI, fhr$HR_LCI[1], tolerance = 1e-12)
  # estimate_hr returns no rate measures (they moved to estimate_ir)
  expect_null(res$incidence_rates)
  expect_null(res$incidence_rate_ratios)
})

test_that("crude HR equals a direct model-SE coxph fit (fixed rule)", {
  df <- sample_df()
  res <- estimate_hr(
    in_df_crude = df, exposure_var = "exposure", outcome_var = "outcome",
    followuptime_var = "follow_up_days", verbose = FALSE
  )
  m <- survival::coxph(survival::Surv(follow_up_days, outcome) ~ exposure,
                       data = df, robust = FALSE)
  sm <- summary(m)
  expect_equal(res$hazard_ratios$HR, unname(sm$conf.int[, "exp(coef)"]))
  expect_equal(res$hazard_ratios$lnHR_SE, unname(sm$coefficients[, "se(coef)"]))
})

test_that("weighted block: robust HR rule", {
  df <- sample_df(); wt <- iptw_df(df)
  res <- estimate_hr(
    in_df_crude = df, in_df_weight = wt, if_weight_weight_var = "iptw",
    exposure_var = "exposure", outcome_var = "outcome",
    followuptime_var = "follow_up_days", verbose = FALSE
  )
  fhr <- fx$hr_ir$hazard_ratios
  expect_equal(res$hazard_ratios$HR[2], fhr$HR[2], tolerance = 1e-12)
  expect_equal(res$hazard_ratios$lnHR_SE[2], fhr$lnHR_SE[2], tolerance = 1e-12)
  expect_match(res$hazard_ratios$Analysis[2], "robust")
})

test_that("strata_var stratifies the hazard model (0.5.0 contract)", {
  df <- sample_df(); wt <- iptw_df(df)
  args <- list(
    in_df_crude = df, in_df_weight = wt, if_weight_weight_var = "iptw",
    exposure_var = "exposure", outcome_var = "outcome",
    followuptime_var = "follow_up_days", verbose = FALSE
  )
  res_s <- do.call(estimate_hr, c(args, list(strata_var = "cat1")))
  res_u <- do.call(estimate_hr, args)

  # strata() Cox HR reproduces the v0.3.0 stratified fixture
  fhrs <- fx$hr_ir_strat$hazard_ratios
  expect_equal(res_s$hazard_ratios$HR, fhrs$HR, tolerance = 1e-12)
  expect_equal(res_s$hazard_ratios$lnHR_SE, fhrs$lnHR_SE, tolerance = 1e-12)
  # and really does differ from the unstratified fit
  expect_false(isTRUE(all.equal(res_s$hazard_ratios$HR,
                                res_u$hazard_ratios$HR, tolerance = 1e-6)))
  expect_match(res_s$hazard_ratios$Analysis[1], "stratified by cat1")
  expect_identical(res_s$strata_var, "cat1")
  expect_null(res_u$strata_var)
})

test_that("Item 3': a strata level missing an arm is rejected, not absorbed", {
  df <- sample_df()
  # cohort-entry-year style stratifier where one period has no reference rows
  df$period <- cut(seq_len(nrow(df)), breaks = 4, labels = FALSE)
  bad <- df
  bad <- bad[!(bad$period == 3 & bad$exposure == 0), ]

  expect_error(
    suppressMessages(estimate_hr(
      in_df_crude = bad, exposure_var = "exposure", outcome_var = "outcome",
      followuptime_var = "follow_up_days", strata_var = "period",
      verbose = FALSE)),
    "no rows in one exposure arm"
  )
  # the message names the level and both counts so the fix is obvious
  err <- tryCatch(suppressMessages(estimate_hr(
    in_df_crude = bad, exposure_var = "exposure", outcome_var = "outcome",
    followuptime_var = "follow_up_days", strata_var = "period",
    verbose = FALSE)), error = function(e) conditionMessage(e))
  expect_match(err, "'3' \\(0 reference,")
  expect_match(err, "Collapse or exclude")

  # the full-support version runs, and every level really does have both arms
  res <- suppressMessages(estimate_hr(
    in_df_crude = df, exposure_var = "exposure", outcome_var = "outcome",
    followuptime_var = "follow_up_days", strata_var = "period",
    verbose = FALSE))
  expect_true(all(table(df$period, df$exposure) > 0))
  expect_identical(nrow(res$hazard_ratios), 1L)
})

test_that("exp_value/ref_value are honoured on an already-0/1 column (Item 8)", {
  # 0.4.0 skipped the exposure recode when the column already looked like 0/1,
  # so exp_value = 0 was ignored and the HR came out un-inverted.
  df <- sample_df()
  call <- function(ev, rv, d = df) suppressMessages(estimate_hr(
    in_df_crude = d, exposure_var = "exposure", exp_value = ev, ref_value = rv,
    outcome_var = "outcome", followuptime_var = "follow_up_days",
    verbose = FALSE))

  hr_a <- call(1, 0)$hazard_ratios$HR
  hr_b <- call(0, 1)$hazard_ratios$HR

  # reference: a direct coxph on the flipped indicator
  df$exp_flip <- 1L - df$exposure
  m <- survival::coxph(survival::Surv(follow_up_days, outcome) ~ exp_flip,
                       data = df, robust = FALSE)
  expect_equal(hr_b, unname(exp(stats::coef(m)["exp_flip"])), tolerance = 1e-10)
  # and it is the inverse of the default direction
  expect_equal(hr_b, 1 / hr_a, tolerance = 1e-10)

  # a non-0/1 labelled column gives the same two answers
  dl <- df; dl$exposure <- ifelse(df$exposure == 1, "TRT", "REF")
  expect_equal(call("TRT", "REF", dl)$hazard_ratios$HR, hr_a, tolerance = 1e-10)
  expect_equal(call("REF", "TRT", dl)$hazard_ratios$HR, hr_b, tolerance = 1e-10)
})

test_that("0.4.0 argument names are gone (0.5.0 breaking renames)", {
  df <- sample_df()
  base <- list(in_df_crude = df, exposure_var = "exposure",
               outcome_var = "outcome",
               followuptime_var = "follow_up_days", verbose = FALSE)
  # stratification_var -> strata_var
  expect_error(do.call(estimate_hr, c(base, list(stratification_var = "cat1"))),
               "unused argument")
  # time_unit dropped: inert on the hazard path
  expect_error(do.call(estimate_hr, c(base, list(time_unit = "days"))),
               "unused argument")
  # ir_per_pyears belongs to estimate_ir
  expect_error(do.call(estimate_hr, c(base, list(ir_per_pyears = 1000))),
               "unused argument")
  fm <- names(formals(estimate_hr))
  expect_true("strata_var" %in% fm)
  expect_false(any(c("stratification_var", "time_unit", "ir_per_pyears") %in% fm))
  expect_identical(fm[length(fm)], "verbose")
})

test_that("Fine-Gray: sHR equals the pre-0.4.0 estimate_rr_rd FG output", {
  df <- sample_df(); wt <- iptw_df(df)
  res <- estimate_hr(
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
    estimate_hr(in_df_crude = df, exposure_var = "exposure",
                outcome_var = "outcome", followuptime_var = "follow_up_days",
                hr_model = "Fine-Gray", verbose = FALSE),
    "if_fg_competing_event_var"
  )
  # Cox + competing var -> ignore message
  expect_message(
    estimate_hr(in_df_crude = df, exposure_var = "exposure",
                outcome_var = "outcome", followuptime_var = "follow_up_days",
                hr_model = "Cox", if_fg_competing_event_var = "competing_event",
                verbose = TRUE),
    "ignored"
  )
})

test_that("matched block: robust + match-id cluster", {
  df <- sample_df()
  md <- make_test_matched_pairs(300, seed = 1)
  res <- estimate_hr(
    in_df_crude = df, in_df_match = md, if_match_match_id = "pair_id",
    exposure_var = "exposure", outcome_var = "outcome",
    followuptime_var = "follow_up_days", verbose = FALSE
  )
  expect_match(res$hazard_ratios$Analysis[2], "match-id cluster")
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
  md <- make_test_matched_pairs(50, seed = 3)
  expect_error(
    estimate_hr(in_df_crude = df, in_df_weight = wt, in_df_match = md,
                if_weight_weight_var = "iptw", exposure_var = "exposure",
                outcome_var = "outcome", followuptime_var = "follow_up_days",
                verbose = FALSE),
    "not both"
  )
  expect_error(
    estimate_hr(in_df_weight = wt, exposure_var = "exposure",
                outcome_var = "outcome", followuptime_var = "follow_up_days",
                verbose = FALSE),
    "if_weight_weight_var"
  )
  expect_error(
    estimate_hr(exposure_var = "exposure", outcome_var = "outcome",
                followuptime_var = "follow_up_days", verbose = FALSE),
    "At least one"
  )
  # crude-less weight-only call is allowed
  res <- estimate_hr(in_df_weight = wt, if_weight_weight_var = "iptw",
                     exposure_var = "exposure", outcome_var = "outcome",
                     followuptime_var = "follow_up_days", verbose = FALSE)
  expect_identical(nrow(res$hazard_ratios), 1L)
  # orphan if_* arguments -> ignore messages
  expect_message(
    estimate_hr(in_df_crude = df, if_weight_weight_var = "iptw",
                exposure_var = "exposure", outcome_var = "outcome",
                followuptime_var = "follow_up_days", verbose = TRUE),
    "if_weight_weight_var is ignored"
  )
  expect_message(
    estimate_hr(in_df_crude = df, if_match_match_id = "pair_id",
                exposure_var = "exposure", outcome_var = "outcome",
                followuptime_var = "follow_up_days", verbose = TRUE),
    "if_match_match_id is ignored"
  )
})

test_that("bootstrap: opt-in, fixed-weight message, seed-reproducible, analytical unchanged", {
  df <- sample_df()
  expect_message(
    estimate_hr(in_df_crude = df, exposure_var = "exposure",
                outcome_var = "outcome", followuptime_var = "follow_up_days",
                if_bootstrap_seed = 1, verbose = TRUE),
    "if_bootstrap_seed is ignored"
  )
  run_boot <- function() {
    suppressMessages(estimate_hr(
      in_df_crude = df, exposure_var = "exposure", outcome_var = "outcome",
      followuptime_var = "follow_up_days",
      if_bootstrap_count = 20, if_bootstrap_n_cores = 1,
      if_bootstrap_seed = 42, verbose = FALSE
    ))
  }
  expect_message(
    estimate_hr(in_df_crude = df[1:400, ], exposure_var = "exposure",
                outcome_var = "outcome", followuptime_var = "follow_up_days",
                if_bootstrap_count = 5, if_bootstrap_n_cores = 1,
                verbose = FALSE),
    "fixed-weight"
  )
  b1 <- run_boot(); b2 <- run_boot()
  expect_true(all(c("HR_BootLCI", "HR_BootSE") %in% names(b1$hazard_ratios)))
  expect_identical(b1$hazard_ratios$HR_BootLCI, b2$hazard_ratios$HR_BootLCI)
  expect_identical(colnames(b1$bootstrap$Crude), "HR")
  # analytical columns unchanged by running the bootstrap
  r0 <- estimate_hr(in_df_crude = df, exposure_var = "exposure",
                    outcome_var = "outcome", followuptime_var = "follow_up_days",
                    verbose = FALSE)
  expect_equal(b1$hazard_ratios$HR, r0$hazard_ratios$HR, tolerance = 1e-12)
})
