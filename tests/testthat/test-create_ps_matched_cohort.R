# Estimate PS first
ps_df_match <- estimate_ps(
  in_df        = small_df,
  exposure_var = "exposure",
  class_vars   = class_vars,
  cont_vars    = cont_vars,
  ps_var       = "ps",
  verbose      = FALSE
)

test_that("create_ps_matched_cohort returns matched data frame", {
  skip_if_not_installed("MatchIt")

  result <- create_ps_matched_cohort(
    in_df        = ps_df_match,
    exposure_var = "exposure",
    ps_var       = "ps",
    ratio        = 1,
    caliper      = 0.2,
    verbose      = FALSE
  )

  expect_s3_class(result, "data.frame")
  # Matched cohort should have fewer or equal rows
  expect_lte(nrow(result), nrow(ps_df_match))
  # Should have both exposure groups
  expect_setequal(unique(result$exposure), c(0, 1))
})

test_that("create_ps_matched_cohort errors when no data provided", {
  expect_error(
    create_ps_matched_cohort(in_df = NULL, in_csvpath = NULL, verbose = FALSE),
    "Either in_df or in_csvpath"
  )
})

test_that("create_ps_matched_cohort rejects invalid ratio", {
  expect_error(
    create_ps_matched_cohort(
      in_df        = ps_df_match,
      exposure_var = "exposure",
      ps_var       = "ps",
      ratio        = -1,
      verbose      = FALSE
    ),
    "ratio must be a single number >= 1"
  )
})
