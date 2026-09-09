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

test_that("weighting stages honor labelled exposure and keep it row-aligned", {
  labelled <- ps_df
  labelled$exposure <- ifelse(labelled$exposure == 1, "GLP1", "DPP4")
  labelled$.check_id <- seq_len(nrow(labelled))

  # Row-aligned, not merely value-preserving: every returned row must carry the
  # label that ITS OWN id had on input. Comparing the two runs' exposure
  # columns to each other cannot detect a shared row-order error.
  aligned <- function(res) {
    identical(res$exposure,
              labelled$exposure[match(res$.check_id, labelled$.check_id)])
  }

  calls <- list(
    iptw = list(fun = create_iptw, weight = "iptw_wt"),
    mw = list(fun = create_matching_weights, weight = "mw_wt"),
    ow = list(fun = create_overlap_weights, weight = "ow_wt")
  )
  for (nm in names(calls)) {
    spec <- calls[[nm]]
    labelled_result <- spec$fun(
      in_df = labelled, exposure_var = "exposure",
      exp_value = "GLP1", ref_value = "DPP4", ps_var = "ps",
      weight_var = spec$weight, verbose = FALSE
    )
    numeric_result <- spec$fun(
      in_df = ps_df, exposure_var = "exposure",
      exp_value = 1, ref_value = 0, ps_var = "ps",
      weight_var = spec$weight, verbose = FALSE
    )

    expect_true(aligned(labelled_result), info = nm)
    expect_equal(labelled_result[[spec$weight]], numeric_result[[spec$weight]],
                 info = nm)
  }

  # trim_method != "none" is the only path that subsets the saved original
  # labels, so it must be exercised at least once.
  trimmed <- create_iptw(
    in_df = labelled, exposure_var = "exposure",
    exp_value = "GLP1", ref_value = "DPP4", ps_var = "ps",
    weight_var = "iptw_wt", trim_method = "crump", verbose = FALSE
  )
  expect_lt(nrow(trimmed), nrow(labelled))
  expect_true(aligned(trimmed))
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
