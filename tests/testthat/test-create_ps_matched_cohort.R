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

test_that("create_ps_matched_cohort supports caliper_scale options", {
  skip_if_not_installed("MatchIt")

  # default scale is logit_ps_sd (PS from estimate_ps is strictly within (0,1))
  res_logit <- create_ps_matched_cohort(
    in_df = ps_df_match, exposure_var = "exposure", ps_var = "ps",
    ratio = 1, verbose = FALSE)
  expect_s3_class(res_logit, "data.frame")
  expect_setequal(unique(res_logit$exposure), c(0, 1))

  # raw_ps_sd reproduces the <= 0.1.x behavior; raw applies a flat caliper
  res_rawsd <- create_ps_matched_cohort(
    in_df = ps_df_match, exposure_var = "exposure", ps_var = "ps",
    ratio = 1, caliper = 0.2, caliper_scale = "raw_ps_sd", verbose = FALSE)
  res_raw <- create_ps_matched_cohort(
    in_df = ps_df_match, exposure_var = "exposure", ps_var = "ps",
    ratio = 1, caliper = 0.05, caliper_scale = "raw", verbose = FALSE)
  expect_s3_class(res_rawsd, "data.frame")
  expect_s3_class(res_raw, "data.frame")
  expect_lte(nrow(res_raw), nrow(ps_df_match))
})

test_that("create_ps_matched_cohort rejects invalid caliper_scale", {
  expect_error(
    create_ps_matched_cohort(
      in_df = ps_df_match, exposure_var = "exposure", ps_var = "ps",
      caliper_scale = "bogus", verbose = FALSE),
    "should be one of"
  )
})

test_that("create_ps_matched_cohort: logit_ps_sd requires PS within (0,1)", {
  skip_if_not_installed("MatchIt")
  bad <- ps_df_match
  bad$ps[1:2] <- c(0, 1)
  expect_error(
    create_ps_matched_cohort(
      in_df = bad, exposure_var = "exposure", ps_var = "ps",
      caliper_scale = "logit_ps_sd", verbose = FALSE),
    "strictly within"
  )
})

test_that("create_ps_matched_cohort supports trimming before matching", {
  skip_if_not_installed("MatchIt")
  res <- create_ps_matched_cohort(
    in_df = ps_df_match, exposure_var = "exposure", ps_var = "ps", ratio = 1,
    trim_method = "crump", trim_crump_alpha = 0.1, verbose = FALSE)
  expect_s3_class(res, "data.frame")
  expect_true(all(res$ps >= 0.1 & res$ps <= 0.9))
})

test_that("create_ps_matched_cohort: subclass matching returns subclasses", {
  skip_if_not_installed("MatchIt")
  res <- create_ps_matched_cohort(
    in_df = ps_df_match, exposure_var = "exposure", ps_var = "ps",
    method = "subclass", subclass_n = 5, verbose = FALSE)
  expect_s3_class(res, "data.frame")
  expect_true(".match_weights" %in% names(res))
  # match_id carries the subclass label for subclass matching
  expect_true("match_id" %in% names(res))
  expect_setequal(unique(res$exposure), c(0, 1))
})

test_that("create_ps_matched_cohort rejects invalid method", {
  expect_error(
    create_ps_matched_cohort(
      in_df = ps_df_match, exposure_var = "exposure", ps_var = "ps",
      method = "bogus", verbose = FALSE),
    "should be one of")
})

test_that("create_ps_matched_cohort emits FS-style distribution plots", {
  skip_if_not_installed("MatchIt")
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("openxlsx")
  d <- file.path(tempdir(), "rwetools_match_plots")
  dir.create(d, showWarnings = FALSE)
  on.exit(unlink(d, recursive = TRUE), add = TRUE)
  suppressWarnings(create_ps_matched_cohort(
    in_df = ps_df_match, exposure_var = "exposure", ps_var = "ps", ratio = 1,
    out_dir_plots = d, out_xlsxpath_report = file.path(d, "rep.xlsx"), verbose = FALSE))
  # legacy comparison/mirror plots retained
  expect_true(file.exists(file.path(d, "ps_distribution_comparison.png")))
  # new unified FS-style set (pre/post-matching panels, no weight box plot)
  expect_true(file.exists(file.path(d, "rep_ps_distr_prematch.png")))
  expect_true(file.exists(file.path(d, "rep_ps_distr_postmatch.png")))
  expect_true(file.exists(file.path(d, "rep_ps_distr_density_prematch.png")))
  expect_true(file.exists(file.path(d, "rep_ps_distr_histog_postmatch.png")))
})
