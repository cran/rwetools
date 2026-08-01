# First, estimate PS on small_df so we have a df with a ps column
ps_df <- estimate_ps(
  in_df        = small_df,
  exposure_var = "exposure",
  class_vars   = class_vars,
  cont_vars    = cont_vars,
  ps_var       = "ps",
  verbose      = FALSE
)

# --- weighting family: create_iptw / create_matching_weights / create_overlap_weights ---

test_that("create_matching_weights adds a bounded weight column (ATM)", {
  result <- create_matching_weights(
    in_df        = ps_df,
    exposure_var = "exposure",
    ps_var       = "ps",
    weight_var   = "mw_wt",
    verbose      = FALSE
  )

  expect_s3_class(result, "data.frame")
  expect_true("mw_wt" %in% names(result))
  expect_equal(nrow(result), nrow(ps_df))
  # matching weights are bounded in [0, 1]
  expect_true(all(result$mw_wt >= 0 & result$mw_wt <= 1, na.rm = TRUE))
})

test_that("create_overlap_weights adds a bounded weight column (ATO)", {
  result <- create_overlap_weights(
    in_df        = ps_df,
    exposure_var = "exposure",
    ps_var       = "ps",
    weight_var   = "ow_wt",
    verbose      = FALSE
  )

  expect_true("ow_wt" %in% names(result))
  # overlap weights equal exp*(1-ps) + (1-exp)*ps, bounded in [0, 1]
  expect_true(all(result$ow_wt >= 0 & result$ow_wt <= 1, na.rm = TRUE))
  exp01 <- ifelse(ps_df$exposure == 1, 1, 0)
  expect_equal(result$ow_wt, exp01 * (1 - ps_df$ps) + (1 - exp01) * ps_df$ps)
})

test_that("create_iptw adds positive ATE weights", {
  result <- create_iptw(
    in_df        = ps_df,
    exposure_var = "exposure",
    ps_var       = "ps",
    weight_var   = "iptw_wt",
    verbose      = FALSE
  )

  expect_true("iptw_wt" %in% names(result))
  expect_true(all(result$iptw_wt > 0, na.rm = TRUE))
  # ATE IPTW: exposed get 1/ps, unexposed get 1/(1-ps)
  exp01 <- ifelse(ps_df$exposure == 1, 1, 0)
  expect_equal(result$iptw_wt, ifelse(exp01 == 1, 1 / ps_df$ps, 1 / (1 - ps_df$ps)))
})

test_that("create_iptw(stabilize = TRUE) returns stabilized weights", {
  result <- create_iptw(
    in_df        = ps_df,
    exposure_var = "exposure",
    ps_var       = "ps",
    weight_var   = "siptw_wt",
    stabilize    = TRUE,
    verbose      = FALSE
  )

  expect_true("siptw_wt" %in% names(result))
  expect_true(all(result$siptw_wt > 0, na.rm = TRUE))
  exp01 <- ifelse(ps_df$exposure == 1, 1, 0)
  p_bar <- mean(exp01)
  expect_equal(result$siptw_wt,
               ifelse(exp01 == 1, p_bar / ps_df$ps, (1 - p_bar) / (1 - ps_df$ps)))
})

test_that("weighting functions error when no data provided", {
  expect_error(
    create_iptw(in_df = NULL, in_csvpath = NULL, verbose = FALSE),
    "Either in_df or in_csvpath"
  )
  expect_error(
    create_matching_weights(in_df = NULL, in_csvpath = NULL, verbose = FALSE),
    "Either in_df or in_csvpath"
  )
})

test_that("the removed create_ps_weights() is gone", {
  expect_false(exists("create_ps_weights", mode = "function"))
})

test_that("create_iptw emits FS-style distribution plots", {
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("openxlsx")
  d <- file.path(tempdir(), "rwetools_wt_plots")
  dir.create(d, showWarnings = FALSE)
  on.exit(unlink(d, recursive = TRUE), add = TRUE)
  suppressWarnings(create_iptw(
    in_df = ps_df, exposure_var = "exposure", ps_var = "ps",
    weight_var = "iptw_wt", out_dir_plots = d,
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
