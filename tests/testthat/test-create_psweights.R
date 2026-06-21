# First, estimate PS on small_df so we have a df with a ps column
ps_df <- estimate_ps(
  in_df        = small_df,
  exposure_var = "exposure",
  class_vars   = class_vars,
  cont_vars    = cont_vars,
  ps_var       = "ps",
  verbose      = FALSE
)

# --- create_ps_weights tests ---

test_that("create_ps_weights adds weight column (mw)", {
  result <- create_ps_weights(
    in_df         = ps_df,
    exposure_var  = "exposure",
    ps_var        = "ps",
    weight_method = "mw",
    weight_var    = "mw_wt",
    verbose       = FALSE
  )

  expect_s3_class(result, "data.frame")
  expect_true("mw_wt" %in% names(result))
  expect_equal(nrow(result), nrow(ps_df))
  expect_true(all(result$mw_wt >= 0, na.rm = TRUE))
})

test_that("create_ps_weights works for ow method", {
  result <- create_ps_weights(
    in_df         = ps_df,
    exposure_var  = "exposure",
    ps_var        = "ps",
    weight_method = "ow",
    weight_var    = "ow_wt",
    verbose       = FALSE
  )

  expect_true("ow_wt" %in% names(result))
  expect_true(all(result$ow_wt >= 0, na.rm = TRUE))
})

test_that("create_ps_weights works for iptw method", {
  result <- create_ps_weights(
    in_df         = ps_df,
    exposure_var  = "exposure",
    ps_var        = "ps",
    weight_method = "iptw",
    weight_var    = "iptw_wt",
    estimand      = "ATE",
    verbose       = FALSE
  )

  expect_true("iptw_wt" %in% names(result))
  expect_true(all(result$iptw_wt > 0, na.rm = TRUE))
})

test_that("create_ps_weights works for stabiptw method", {
  result <- create_ps_weights(
    in_df         = ps_df,
    exposure_var  = "exposure",
    ps_var        = "ps",
    weight_method = "stabiptw",
    weight_var    = "sipw_wt",
    estimand      = "ATT",
    verbose       = FALSE
  )

  expect_true("sipw_wt" %in% names(result))
  expect_true(all(result$sipw_wt > 0, na.rm = TRUE))
})

test_that("create_ps_weights rejects invalid weight_method", {
  expect_error(
    create_ps_weights(
      in_df         = ps_df,
      exposure_var  = "exposure",
      ps_var        = "ps",
      weight_method = "invalid_method",
      verbose       = FALSE
    ),
    "Invalid weight_method"
  )
})

test_that("create_ps_weights errors when no data provided", {
  expect_error(
    create_ps_weights(in_df = NULL, in_csvpath = NULL, verbose = FALSE),
    "Either in_df or in_csvpath"
  )
})

test_that("create_ps_weights emits FS-style distribution plots", {
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("openxlsx")
  d <- file.path(tempdir(), "rwetools_wt_plots")
  dir.create(d, showWarnings = FALSE)
  on.exit(unlink(d, recursive = TRUE), add = TRUE)
  suppressWarnings(create_ps_weights(
    in_df = ps_df, exposure_var = "exposure", ps_var = "ps", weight_method = "mw",
    weight_var = "mw_wt", out_dir_plots = d,
    out_xlsxpath_report = file.path(d, "rep.xlsx"), verbose = FALSE))
  # legacy faceted plots retained
  expect_true(file.exists(file.path(d, "ps_distribution_unweighted.png")))
  expect_true(file.exists(file.path(d, "ps_distribution_weighted.png")))
  # new unified FS-style set
  expect_true(file.exists(file.path(d, "rep_ps_distr_unwt.png")))
  expect_true(file.exists(file.path(d, "rep_ps_distr_wt.png")))
  expect_true(file.exists(file.path(d, "rep_ps_distr_density_wt.png")))
  expect_true(file.exists(file.path(d, "rep_boxplot.png")))
})
