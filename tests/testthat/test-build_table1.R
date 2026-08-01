test_that("build_table1 returns data frame for unweighted analysis", {
  result <- build_table1(
    in_df        = small_df,
    exposure_var = "exposure",
    cont_vars    = cont_vars,
    binary_vars  = binary_vars,
    cat_vars     = class_vars,
    verbose      = FALSE
  )

  expect_s3_class(result, "data.frame")
  expect_true(nrow(result) > 0)
  expect_true("Variable" %in% names(result))
})

test_that("build_table1 returns data frame for weighted analysis", {
  ps_df_t1 <- estimate_ps(
    in_df        = small_df,
    exposure_var = "exposure",
    class_vars   = class_vars,
    cont_vars    = cont_vars,
    ps_var       = "ps",
    verbose      = FALSE
  )
  wt_df_t1 <- create_matching_weights(
    in_df         = ps_df_t1,
    exposure_var  = "exposure",
    ps_var        = "ps",
    weight_var    = "mw_wt",
    verbose       = FALSE
  )

  result <- build_table1(
    in_df        = wt_df_t1,
    exposure_var = "exposure",
    use_weights  = TRUE,
    weight_var   = "mw_wt",
    cont_vars    = cont_vars,
    binary_vars  = binary_vars,
    cat_vars     = class_vars,
    verbose      = FALSE
  )

  expect_s3_class(result, "data.frame")
  expect_true(nrow(result) > 0)
})

test_that("build_table1 writes Excel file when out_xlsxpath provided", {
  skip_if_not_installed("openxlsx")

  out_xlsx <- tempfile(fileext = ".xlsx")
  on.exit(unlink(out_xlsx))

  result <- build_table1(
    in_df        = small_df,
    out_xlsxpath = out_xlsx,
    exposure_var = "exposure",
    cont_vars    = cont_vars,
    binary_vars  = binary_vars,
    verbose      = FALSE
  )

  expect_true(file.exists(out_xlsx))
})

test_that("build_table1 errors when exposure_var not in data", {
  expect_error(
    build_table1(
      in_df        = small_df,
      exposure_var = "nonexistent_var",
      verbose      = FALSE
    )
  )
})
