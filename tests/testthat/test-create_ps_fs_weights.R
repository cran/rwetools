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

test_that("create_ps_fs_weights keeps labelled exposure row-aligned and orients its diagnostics", {
  skip_if_not_installed("openxlsx")
  labelled <- ps_df_fs
  labelled$exposure <- ifelse(labelled$exposure == 1, "GLP1", "DPP4")
  labelled$.check_id <- seq_len(nrow(labelled))
  n_glp1 <- sum(labelled$exposure == "GLP1")
  n_dpp4 <- sum(labelled$exposure == "DPP4")
  # a swapped Exp/Ref column can only be caught if the two arms differ in size
  expect_false(n_glp1 == n_dpp4)

  report <- tempfile(fileext = ".xlsx")
  common <- list(
    ps_var = "ps", weight_var = "fs_wt", number_of_strata = 10,
    stratification_method = "exposure", estimand = "ATT",
    trim_nonoverlap_region = TRUE,
    make_unwt_wt_table1 = TRUE, table1_cont_vars = "cont1",
    table1_binary_vars = "binary1", verbose = FALSE
  )
  labelled_result <- do.call(create_ps_fs_weights, c(list(
    in_df = labelled, exposure_var = "exposure",
    exp_value = "GLP1", ref_value = "DPP4",
    out_xlsxpath_report = report
  ), common))
  numeric_result <- do.call(create_ps_fs_weights, c(list(
    in_df = ps_df_fs, exposure_var = "exposure",
    exp_value = 1, ref_value = 0
  ), common))

  # Trimming drops rows, so this asserts the row-order restore and not just the
  # label mapping: every returned row must carry the label ITS OWN id had.
  # Comparing the two runs' exposure columns to each other cannot do that --
  # both pass through the same restore line, so a shared error would cancel.
  expect_identical(
    labelled_result$exposure,
    labelled$exposure[match(labelled_result$.check_id, labelled$.check_id)]
  )
  expect_lt(nrow(labelled_result), nrow(ps_df_fs))
  expect_equal(labelled_result$fs_wt, numeric_result$fs_wt)

  # The internal diagnostics run on the canonical 1/0 coding. Pin the
  # direction: the crude Exp column must be the GLP1 arm, Ref the DPP4 arm.
  bal <- openxlsx::read.xlsx(report, sheet = "Balance Table", startRow = 4)
  n_row <- bal[bal$Variable == "n_of_patients", , drop = FALSE]
  expect_identical(nrow(n_row), 1L)
  expect_equal(as.numeric(n_row$Exp_Crude), n_glp1)
  expect_equal(as.numeric(n_row$Ref_Crude), n_dpp4)
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

test_that("create_ps_fs_weights emits FS-style distribution plots", {
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("openxlsx")
  d <- file.path(tempdir(), "rwetools_fs_plots")
  dir.create(d, showWarnings = FALSE)
  on.exit(unlink(d, recursive = TRUE), add = TRUE)
  suppressWarnings(create_ps_fs_weights(
    in_df = ps_df_fs, exposure_var = "exposure", ps_var = "ps", weight_var = "fs_wt",
    number_of_strata = 10, stratification_method = "exposure", estimand = "ATT",
    trim_nonoverlap_region = TRUE, out_dir_plots = d,
    out_xlsxpath_report = file.path(d, "rep.xlsx"), verbose = FALSE))
  expect_true(file.exists(file.path(d, "rep_ps_distr_unwt.png")))
  expect_true(file.exists(file.path(d, "rep_ps_distr_wt.png")))
  expect_true(file.exists(file.path(d, "rep_ps_distr_density_unwt.png")))
  expect_true(file.exists(file.path(d, "rep_boxplot.png")))
})
