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

test_that("create_ps_matched_cohort shares the common exposure canonicalization", {
  skip_if_not_installed("MatchIt")
  labelled <- ps_df_match
  labelled$exposure <- ifelse(labelled$exposure == 1, "GLP1", "DPP4")
  labelled$.check_id <- seq_len(nrow(labelled))

  result <- create_ps_matched_cohort(
    in_df = labelled, exposure_var = "exposure",
    exp_value = "GLP1", ref_value = "DPP4",
    ps_var = "ps", ratio = 1, caliper = 0.2, verbose = FALSE
  )
  # original labels, row-aligned to the id each row came in with
  expect_setequal(unique(result$exposure), c("GLP1", "DPP4"))
  expect_identical(
    result$exposure,
    labelled$exposure[match(result$.check_id, labelled$.check_id)]
  )
  # exp_value really selected the GLP1 arm as treated: the labelled call must
  # reproduce the numeric-coded call row-for-row, and naming the other arm as
  # exposed must NOT give the same cohort (ATT is arm-specific).
  numeric_in <- ps_df_match
  numeric_in$.check_id <- seq_len(nrow(numeric_in))
  same_direction <- create_ps_matched_cohort(
    in_df = numeric_in, exposure_var = "exposure",
    exp_value = 1, ref_value = 0,
    ps_var = "ps", ratio = 1, caliper = 0.2, verbose = FALSE
  )
  # naming the larger arm as exposed legitimately warns about control supply
  other_direction <- suppressWarnings(create_ps_matched_cohort(
    in_df = labelled, exposure_var = "exposure",
    exp_value = "DPP4", ref_value = "GLP1",
    ps_var = "ps", ratio = 1, caliper = 0.2, verbose = FALSE
  ))
  expect_identical(result$.check_id, same_direction$.check_id)
  expect_false(identical(result$.check_id, other_direction$.check_id))

  # the shared helper's value/type semantics, not a local variant
  expect_error(
    create_ps_matched_cohort(
      in_df = labelled, exposure_var = "exposure",
      exp_value = "GLP1", ref_value = "GLP1", ps_var = "ps", verbose = FALSE),
    "must differ"
  )
  expect_error(
    create_ps_matched_cohort(
      in_df = labelled, exposure_var = "exposure",
      exp_value = "SGLT2i", ref_value = "DPP4", ps_var = "ps", verbose = FALSE),
    "not found in exposure"
  )
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

test_that("built-in Table 1 honors non-unit nearest-matching weights", {
  skip_if_not_installed("MatchIt")
  skip_if_not_installed("openxlsx")

  report_path <- tempfile(fileext = ".xlsx")
  on.exit(unlink(report_path), add = TRUE)

  res <- suppressWarnings(create_ps_matched_cohort(
    in_df = ps_df_match,
    exposure_var = "exposure",
    ps_var = "ps",
    ratio = 2,
    min_controls = 1,
    max_controls = 3,
    make_crude_matched_table1 = TRUE,
    table1_cont_vars = "cont1",
    table1_binary_vars = character(),
    table1_cat_vars = character(),
    out_xlsxpath_report = report_path,
    verbose = FALSE
  ))

  expect_true(any(abs(res$.match_weights - 1) > 1e-8))

  observed <- openxlsx::read.xlsx(report_path, sheet = "Matched_Table1")
  expected <- build_table1(
    in_df = res,
    exposure_var = "exposure",
    cont_vars = "cont1",
    binary_vars = character(),
    cat_vars = character(),
    use_weights = TRUE,
    weight_var = ".match_weights",
    verbose = FALSE
  )
  compare_cols <- c("Total", "Exp", "Ref", "Crude_diff", "Std_diff")
  expect_equal(
    observed[observed$Variable == "cont1", compare_cols],
    expected[expected$Variable == "cont1", compare_cols],
    ignore_attr = TRUE
  )
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
