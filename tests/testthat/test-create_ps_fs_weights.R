# Estimate PS first
ps_df_fs <- estimate_ps(
  in_df        = small_df,
  exposure_var = "exposure",
  class_vars   = class_vars,
  cont_vars    = cont_vars,
  ps_var       = "ps",
  verbose      = FALSE
)

test_that("create_ps_fs_weights returns weighted data frame", {
  result <- create_ps_fs_weights(
    in_df                  = ps_df_fs,
    exposure_var           = "exposure",
    ps_var                 = "ps",
    weight_var             = "fs_wt",
    number_of_strata       = 10,
    stratification_method  = "exposure",
    estimand               = "ATT",
    trim_nonoverlap_region = TRUE,
    verbose                = FALSE
  )

  expect_s3_class(result, "data.frame")
  expect_true("fs_wt" %in% names(result))
  # Trimming may remove rows
  expect_lte(nrow(result), nrow(ps_df_fs))
  expect_true(all(result$fs_wt > 0, na.rm = TRUE))
})

test_that("create_ps_fs_weights errors when ps_var is NULL", {
  expect_error(
    create_ps_fs_weights(
      in_df        = ps_df_fs,
      exposure_var = "exposure",
      ps_var       = NULL,
      verbose      = FALSE
    ),
    "ps_var must be specified"
  )
})

test_that("create_ps_fs_weights errors when ps_var not in data", {
  expect_error(
    create_ps_fs_weights(
      in_df        = ps_df_fs,
      exposure_var = "exposure",
      ps_var       = "nonexistent_col",
      verbose      = FALSE
    ),
    "not found in dataset"
  )
})
