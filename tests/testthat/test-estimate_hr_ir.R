test_that("estimate_hr_ir returns results for unweighted analysis", {
  result <- estimate_hr_ir(
    in_df_unwt    = small_df,
    in_df_wted    = NULL,
    exposure_var  = "exposure",
    outcome_var   = "outcome",
    survival_time = "follow_up_days",
    time_unit     = "days",
    ir_per_pyears = 1000,
    verbose       = FALSE
  )

  expect_type(result, "list")
  # Should contain at least one data frame of results
  has_df <- any(vapply(result, is.data.frame, logical(1)))
  expect_true(has_df)
})

test_that("estimate_hr_ir returns results for weighted analysis", {
  # Prepare weighted data
  ps_df_hr <- estimate_ps(
    in_df        = small_df,
    exposure_var = "exposure",
    class_vars   = class_vars,
    cont_vars    = cont_vars,
    ps_var       = "ps",
    verbose      = FALSE
  )
  wt_df_hr <- create_ps_weights(
    in_df         = ps_df_hr,
    exposure_var  = "exposure",
    ps_var        = "ps",
    weight_method = "mw",
    weight_var    = "mw_wt",
    verbose       = FALSE
  )

  result <- estimate_hr_ir(
    in_df_unwt    = small_df,
    in_df_wted    = wt_df_hr,
    exposure_var  = "exposure",
    outcome_var   = "outcome",
    survival_time = "follow_up_days",
    weight_var    = "mw_wt",
    time_unit     = "days",
    ir_per_pyears = 1000,
    verbose       = FALSE
  )

  expect_type(result, "list")
})

test_that("estimate_hr_ir errors when no data provided", {
  expect_error(
    estimate_hr_ir(
      in_df_unwt    = NULL,
      in_df_wted    = NULL,
      outcome_var   = "outcome",
      survival_time = "follow_up_days",
      verbose       = FALSE
    ),
    "At least one of in_df_unwt or in_df_wted"
  )
})
