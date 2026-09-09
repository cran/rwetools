# Guards on which matched cohorts the matched block will accept
# (backlog Items 5 and 6, both landing in validate_effect_blocks()).
#
# The matched block analyses rows UNWEIGHTED and clusters inference on the
# match id. That is only valid when (a) every row belongs to exactly one
# matched set and (b) the matching weights carry nothing the point estimate
# needs -- i.e. every value is exactly 1. Cohorts violating
# either are rejected with a message naming the correct alternative.

sample_df <- function() {
  read.csv(system.file("extdata", "sample_data.csv", package = "rwetools"))
}
ps_df <- function(df) {
  estimate_ps(in_df = df, exposure_var = "exposure",
              class_vars = paste0("cat", 1:4),
              cont_vars = paste0("cont", 1:7),
              ps_var = "ps", verbose = FALSE)
}
matched <- function(ps, ...) {
  suppressWarnings(suppressMessages(
    create_ps_matched_cohort(in_df = ps, exposure_var = "exposure",
                             ps_var = "ps", verbose = FALSE, ...)))
}
call_ir <- function(d, ...) suppressWarnings(suppressMessages(estimate_ir(
  in_df_match = d, exposure_var = "exposure", outcome_var = "outcome",
  followuptime_var = "follow_up_days", time_unit = "days",
  verbose = FALSE, ...)))

test_that("1:1 nearest matching is accepted (weights are all 1)", {
  m <- matched(ps_df(sample_df()), method = "nearest", ratio = 1, caliper = 0.2)
  expect_true(all(m$.match_weights == 1))
  expect_false(any(grepl(";", m$match_id, fixed = TRUE)))
  res <- call_ir(m, if_match_match_id = "match_id")
  expect_true("IR_per_Npy_Matched" %in% names(res$incidence_rates))
})

test_that("Item 6: with-replacement match ids are rejected", {
  m <- matched(ps_df(sample_df()), method = "nearest", ratio = 1,
               replace = TRUE, caliper = 0.2)
  # create_ps_matched_cohort joins the ids of every set a reused control
  # belongs to with ";"
  expect_true(any(grepl(";", m$match_id, fixed = TRUE)))
  expect_error(call_ir(m, if_match_match_id = "match_id"),
               "name more than one matched set")
  # dropping the match id does not sneak the cohort past: its non-unit weights
  # also trigger Item 5's guard (belt and suspenders)
  m2 <- m; m2$match_id <- NULL
  expect_error(call_ir(m2), "every \\.match_weights value must equal 1")
  # and the same guards apply in estimate_hr / estimate_risk
  expect_error(
    suppressMessages(estimate_hr(
      in_df_match = m, if_match_match_id = "match_id",
      exposure_var = "exposure", outcome_var = "outcome",
      followuptime_var = "follow_up_days", verbose = FALSE)),
    "name more than one matched set"
  )
  expect_error(
    suppressMessages(estimate_risk(
      in_df_match = m, if_match_match_id = "match_id",
      exposure_var = "exposure", outcome_var = "outcome",
      followuptime_var = "follow_up_days", verbose = FALSE)),
    "name more than one matched set"
  )
})

test_that("Item 5: non-unit matching weights are rejected", {
  ps <- ps_df(sample_df())
  for (nm in names(list(subclass = 1, variable_ratio = 1, ratio2 = 1))) {
    m <- switch(nm,
      subclass       = matched(ps, method = "subclass", subclass_n = 5),
      variable_ratio = matched(ps, method = "nearest", ratio = 2,
                               min_controls = 1, max_controls = 3,
                               caliper = 0.2),
      ratio2         = matched(ps, method = "nearest", ratio = 2, caliper = 0.2)
    )
    # the premise: weights really do vary inside the reference arm
    w0 <- m$.match_weights[m$exposure == 0]
    expect_gt(diff(range(w0)), 0, label = paste0(nm, ": within-arm weight range"))
    expect_error(call_ir(m, if_match_match_id = "match_id"),
                 "every \\.match_weights value must equal 1",
                 info = nm)
  }
})

test_that("Item 5: the documented alternative (weighted block) does run", {
  m <- matched(ps_df(sample_df()), method = "subclass", subclass_n = 5)
  res <- suppressWarnings(suppressMessages(estimate_ir(
    in_df_weight = m, if_weight_weight_var = ".match_weights",
    exposure_var = "exposure", outcome_var = "outcome",
    followuptime_var = "follow_up_days", time_unit = "days",
    verbose = FALSE)))
  expect_true("IR_per_Npy_Weighted" %in% names(res$incidence_rates))
  # and it is NOT the same answer as ignoring the weights would have given
  m_unw <- m; m_unw$.match_weights <- 1
  res_unw <- suppressWarnings(suppressMessages(estimate_ir(
    in_df_weight = m_unw, if_weight_weight_var = ".match_weights",
    exposure_var = "exposure", outcome_var = "outcome",
    followuptime_var = "follow_up_days", time_unit = "days",
    verbose = FALSE)))
  expect_false(isTRUE(all.equal(res$incidence_rate_ratios$IRR,
                               res_unw$incidence_rate_ratios$IRR,
                               tolerance = 1e-4)))
})

test_that("arm-constant but unequal matching weights are rejected", {
  # These constants cancel from arm-level rates/risks but not generally from
  # the Cox partial likelihood, so the shared matched-block contract rejects
  # the ambiguous input for all three effect functions.
  m <- matched(ps_df(sample_df()), method = "nearest", ratio = 1, caliper = 0.2)
  m$.match_weights <- ifelse(m$exposure == 1, 1, 0.5)
  expect_error(call_ir(m, if_match_match_id = "match_id"),
               "every \\.match_weights value must equal 1")
  expect_error(suppressMessages(estimate_hr(
    in_df_match = m, if_match_match_id = "match_id",
    exposure_var = "exposure", outcome_var = "outcome",
    followuptime_var = "follow_up_days", verbose = FALSE)),
    "every \\.match_weights value must equal 1")
})

test_that("all-1 weights do not let a 1:2 set enter the matched block", {
  d <- data.frame(
    exposure = c(1L, 0L, 0L, 1L, 0L, 0L),
    outcome = c(1L, 0L, 0L, 0L, 1L, 0L),
    follow_up_days = c(100, 110, 120, 130, 140, 150),
    match_id = rep(c("A", "B"), each = 3),
    .match_weights = 1
  )
  expect_error(call_ir(d, if_match_match_id = "match_id"),
               "strict 1:1 matching only.*2 reference, 1 exposed")
})

test_that("create_ps_matched_cohort warns at creation time for replace = TRUE", {
  expect_message(
    suppressWarnings(create_ps_matched_cohort(
      in_df = ps_df(sample_df()), exposure_var = "exposure", ps_var = "ps",
      method = "nearest", ratio = 1, replace = TRUE, caliper = 0.2,
      verbose = TRUE)),
    "cannot be analysed through the matched block"
  )
})
