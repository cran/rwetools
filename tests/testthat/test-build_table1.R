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

test_that("fractional weights drive every publication-facing weighted count", {
  d <- data.frame(
    exposure = c(0, 0, 1, 1),
    binary = c(0, 1, 0, 1),
    w = c(0.25, 0.75, 0.50, 0.50)
  )
  result <- build_table1(
    in_df = d, exposure_var = "exposure", use_weights = TRUE,
    weight_var = "w", binary_vars = "binary", n_decimal = 2,
    verbose = FALSE
  )

  n_row <- result[result$Variable == "n_of_patients", ]
  b_row <- result[result$Variable == "binary", ]
  expect_identical(n_row[c("Total", "Exp", "Ref")],
                   data.frame(Total = "2.00", Exp = "1.00", Ref = "1.00"))
  expect_identical(b_row$Total, "1.25 (62.5%)")
  expect_identical(b_row$Exp, "0.50 (50.0%)")
  expect_identical(b_row$Ref, "0.75 (75.0%)")
})

test_that("weighted Table 1 rejects missing or invalid weight columns", {
  base <- data.frame(exposure = c(0, 1), x = c(1, 2))
  expect_error(build_table1(base, exposure_var = "exposure", use_weights = TRUE,
                            weight_var = "absent", cont_vars = "x", verbose = FALSE),
               "not found")

  invalid <- list(
    character = c("1", "1"),
    missing = c(1, NA_real_),
    infinite = c(1, Inf),
    negative = c(1, -0.1)
  )
  for (w in invalid) {
    d <- base
    d$w <- w
    expect_error(build_table1(d, exposure_var = "exposure", use_weights = TRUE,
                              weight_var = "w", cont_vars = "x", verbose = FALSE))
  }
})

test_that("zero-weight rows are removed before survey statistics", {
  d <- data.frame(
    exposure = c(0, 0, 0, 1, 1, 1),
    binary = c(0, 1, 1, 0, 1, 1),
    value = c(10, 20, 1000, 30, 40, 2000),
    w = c(1, 1, 0, 1, 1, 0)
  )
  args <- list(exposure_var = "exposure", use_weights = TRUE,
               weight_var = "w", cont_vars = "value", binary_vars = "binary",
               verbose = FALSE)

  expect_message(with_zero <- do.call(build_table1, c(list(in_df = d), args)),
                 "Removed 2 zero-weight row")
  manually_filtered <- do.call(build_table1, c(list(in_df = d[d$w > 0, ]), args))
  expect_identical(with_zero, manually_filtered)
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
