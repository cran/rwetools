test_that("estimate_rr_rd returns results for unweighted analysis", {
  skip_on_cran()
  result <- estimate_rr_rd(
    in_df_unwt            = small_df,
    in_df_wted            = NULL,
    exposure_var          = "exposure",
    outcome_var           = "outcome",
    survival_time         = "follow_up_days",
    time_unit             = "days",
    rr_rd_at_timepoint    = 365,
    rr_rd_per_individuals = 1000,
    bootstrap_count       = 20,
    n_cores               = 1,
    verbose               = FALSE
  )

  expect_type(result, "list")
  has_df <- any(vapply(result, is.data.frame, logical(1)))
  expect_true(has_df)

  # Schema: Risk_Estimator present, Survival columns removed
  expect_true("Risk_Estimator" %in% names(result$estimates))
  expect_equal(unique(result$estimates$Risk_Estimator), "KM")
  expect_true("Risk_Estimator" %in% names(result$cumulative_incidence))
  expect_true("Risk_SE" %in% names(result$cumulative_incidence))
  expect_false("Survival" %in% names(result$cumulative_incidence))
  expect_false("Survival_SE" %in% names(result$cumulative_incidence))
})

test_that("estimate_rr_rd returns results for weighted analysis", {
  skip_on_cran()
  # Prepare weighted data
  ps_df_rr <- estimate_ps(
    in_df        = small_df,
    exposure_var = "exposure",
    class_vars   = class_vars,
    cont_vars    = cont_vars,
    ps_var       = "ps",
    verbose      = FALSE
  )
  wt_df_rr <- create_ps_weights(
    in_df         = ps_df_rr,
    exposure_var  = "exposure",
    ps_var        = "ps",
    weight_method = "mw",
    weight_var    = "mw_wt",
    verbose       = FALSE
  )

  result <- estimate_rr_rd(
    in_df_unwt            = small_df,
    in_df_wted            = wt_df_rr,
    exposure_var          = "exposure",
    outcome_var           = "outcome",
    survival_time         = "follow_up_days",
    weight_var            = "mw_wt",
    time_unit             = "days",
    rr_rd_at_timepoint    = 365,
    rr_rd_per_individuals = 1000,
    bootstrap_count       = 20,
    n_cores               = 1,
    verbose               = FALSE
  )

  expect_type(result, "list")
})

test_that("estimate_rr_rd errors when no data provided", {
  expect_error(
    estimate_rr_rd(
      in_df_unwt    = NULL,
      in_df_wted    = NULL,
      outcome_var   = "outcome",
      survival_time = "follow_up_days",
      verbose       = FALSE
    ),
    "At least one of in_df_unwt or in_df_wted"
  )
})
