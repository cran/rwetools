# Trimming / truncation / SMD-scale regression tests for the weighting family.
# Uses the shared small_df / class_vars / cont_vars from helper.R.

ps_df_feat <- estimate_ps(
  in_df        = small_df,
  exposure_var = "exposure",
  class_vars   = class_vars,
  cont_vars    = cont_vars,
  ps_var       = "ps",
  verbose      = FALSE
)

# ---- .trim_ps ----

test_that(".trim_ps('none') keeps everything", {
  keep <- .trim_ps(ps_df_feat$ps, ifelse(ps_df_feat$exposure == 1, 1, 0),
                   method = "none", verbose = FALSE)
  expect_true(all(keep))
  expect_length(keep, nrow(ps_df_feat))
})

test_that("Crump trimming keeps PS within [alpha, 1 - alpha]", {
  exp01 <- ifelse(ps_df_feat$exposure == 1, 1, 0)
  keep <- .trim_ps(ps_df_feat$ps, exp01, method = "crump",
                   crump_alpha = 0.1, verbose = FALSE)
  kept_ps <- ps_df_feat$ps[keep]
  expect_true(all(kept_ps >= 0.1 & kept_ps <= 0.9))
  # anything outside is dropped
  expect_equal(sum(!keep), sum(ps_df_feat$ps < 0.1 | ps_df_feat$ps > 0.9))
})

test_that("Sturmer trimming uses exposure-group percentile bounds", {
  exp01 <- ifelse(ps_df_feat$exposure == 1, 1, 0)
  p <- 0.05
  keep <- .trim_ps(ps_df_feat$ps, exp01, method = "sturmer",
                   sturmer_p = p, verbose = FALSE)
  lo <- as.numeric(stats::quantile(ps_df_feat$ps[exp01 == 1], p))
  hi <- as.numeric(stats::quantile(ps_df_feat$ps[exp01 == 0], 1 - p))
  kept_ps <- ps_df_feat$ps[keep]
  expect_true(all(kept_ps >= lo & kept_ps <= hi))
})

test_that(".trim_ps rejects out-of-range tuning parameters", {
  exp01 <- ifelse(ps_df_feat$exposure == 1, 1, 0)
  expect_error(.trim_ps(ps_df_feat$ps, exp01, method = "crump",
                        crump_alpha = 0.6, verbose = FALSE), "0, 0.5")
  expect_error(.trim_ps(ps_df_feat$ps, exp01, method = "sturmer",
                        sturmer_p = 0, verbose = FALSE), "0, 0.5")
})

test_that("create_iptw with Crump trimming removes rows", {
  full <- create_iptw(ps_df_feat, exposure_var = "exposure", ps_var = "ps",
                      verbose = FALSE)
  trimmed <- create_iptw(ps_df_feat, exposure_var = "exposure", ps_var = "ps",
                         trim_method = "crump", trim_crump_alpha = 0.1,
                         verbose = FALSE)
  expect_lte(nrow(trimmed), nrow(full))
  expect_true(all(trimmed$ps >= 0.1 & trimmed$ps <= 0.9))
})

# ---- .truncate_ps_weights ----

test_that("percentile truncation winsorizes to the quantile cut points", {
  w <- c(1, 2, 3, 50, 100)
  res <- .truncate_ps_weights(w, method = "percentile",
                              percentile = c(0, 0.8), verbose = FALSE)
  hi <- as.numeric(stats::quantile(w, 0.8))
  expect_true(max(res$w) <= hi + 1e-9)
  expect_equal(res$cut[2], hi)
})

test_that("cap truncation applies an absolute upper bound", {
  w <- c(0.5, 5, 20)
  res <- .truncate_ps_weights(w, method = "cap", cap = 10, verbose = FALSE)
  expect_equal(max(res$w), 10)
  expect_error(.truncate_ps_weights(w, method = "cap", cap = -1, verbose = FALSE),
               "positive number")
})

test_that("create_iptw applies weight truncation", {
  res <- create_iptw(ps_df_feat, exposure_var = "exposure", ps_var = "ps",
                     weight_var = "wt", truncate_method = "cap",
                     truncate_cap = 5, verbose = FALSE)
  expect_lte(max(res$wt, na.rm = TRUE), 5)
})

# ---- SMD scale regression (.is_balanced) ----

test_that(".is_balanced compares on the raw SMD scale (no x100 bug)", {
  # |SMD| = 0.45 must NOT be 'balanced' at threshold 0.1
  expect_false(.is_balanced(0.45, 0.1))
  expect_false(.is_balanced(-0.45, 0.1))
  expect_true(.is_balanced(0.05, 0.1))
  expect_true(is.na(.is_balanced(NA_real_, 0.1)))
})

# ---- Love plot reference line / label (raw scale) ----

test_that("create_love_plot draws its reference line on the raw SMD scale", {
  skip_if_not_installed("ggplot2")
  out <- tempfile(fileext = ".png")
  on.exit(unlink(out), add = TRUE)
  p <- create_love_plot(
    variable_names    = c("a", "b"),
    crude_std_diff    = c(0.45, 0.30),
    adjusted_std_diff = c(0.05, 0.02),
    crude_label = "Unweighted", adjusted_label = "Weighted",
    title = "t", output_path = out, std_diff_threshold = 0.1
  )
  expect_s3_class(p, "ggplot")
  # axis label must not advertise a percentage scale
  expect_false(grepl("%", p$labels$x))
  # the dashed reference line is at +/- 0.1, never at 10
  xint <- unlist(lapply(p$layers, function(l) l$data$xintercept))
  xint <- xint[is.finite(xint)]
  expect_true(any(abs(abs(xint) - 0.1) < 1e-9))
  expect_false(any(abs(xint - 10) < 1e-9))
})
