test_that("estimate_ps returns data.frame with PS column", {
  result <- estimate_ps(
    in_df        = small_df,
    exposure_var = "exposure",
    class_vars   = class_vars,
    cont_vars    = cont_vars,
    ps_var       = "ps",
    verbose      = FALSE
  )

  expect_s3_class(result, "data.frame")
  expect_true("ps" %in% names(result))
  expect_equal(nrow(result), nrow(small_df))
  expect_true(all(result$ps >= 0 & result$ps <= 1))
})

test_that("estimate_ps errors when no data provided", {
  expect_error(
    estimate_ps(in_df = NULL, in_csvpath = NULL, verbose = FALSE),
    "Either in_df or in_csvpath"
  )
})

test_that("estimate_ps warns when both in_df and in_csvpath given", {
  expect_warning(
    estimate_ps(
      in_df        = small_df,
      in_csvpath   = sample_csv,
      exposure_var = "exposure",
      class_vars   = class_vars,
      cont_vars    = cont_vars,
      verbose      = FALSE
    ),
    "Both in_df and in_csvpath"
  )
})

test_that("estimate_ps custom ps_var name works", {
  result <- estimate_ps(
    in_df        = small_df,
    exposure_var = "exposure",
    class_vars   = class_vars,
    cont_vars    = cont_vars,
    ps_var       = "my_ps",
    verbose      = FALSE
  )

  expect_true("my_ps" %in% names(result))
})

# Regression test for the Rule 1 (Sparse level) zero-count detection fix.
# A categorical level that appears in only one exposure group (n = 0 in the
# other) used to slip past the cnt < 5 screen because table(useNA = "no")
# silently dropped absent levels. The union-of-levels patch makes the 0-count
# cell visible to the check.
#
# Deterministic toy frame (no set.seed): 60 rows, exposure 30/30.
#   - cont1: well-separated continuous covariate
#   - good_cat: present in both groups, no sparse level
#   - race_cat: levels {0, 1} in exposure == 1; levels {0, 1, 2} in exposure == 0.
#     Level 2 has n = 0 in the exposed group -> must be detected as extreme.
.toy_zero_count_df <- function() {
  exposure <- c(rep(1L, 30), rep(0L, 30))
  cont1    <- c(seq(0.10, 3.00, length.out = 30),
                seq(0.05, 2.95, length.out = 30))
  good_cat <- c(rep(c(0L, 1L), times = 15),
                rep(c(0L, 1L), times = 15))
  race_cat <- c(rep(0L, 15), rep(1L, 15),                  # exp = 1: only {0, 1}
                rep(0L, 10), rep(1L, 10), rep(2L, 10))     # exp = 0: {0, 1, 2}
  data.frame(
    exposure = exposure,
    cont1    = cont1,
    good_cat = good_cat,
    race_cat = race_cat,
    stringsAsFactors = FALSE
  )
}

test_that("estimate_ps flags categorical level with zero count in one group", {
  toy_df <- .toy_zero_count_df()

  # With auto-exclusion on, race_cat must be detected and excluded; the call
  # then succeeds using only cont1 + good_cat. The verbose output emits the
  # "Auto-excluding" header and the variable name on separate message() lines,
  # so capture all messages and assert both substrings rather than requiring
  # them in a single message.
  msgs <- capture_messages(
    result <- estimate_ps(
      in_df                               = toy_df,
      exposure_var                        = "exposure",
      class_vars                          = c("good_cat", "race_cat"),
      cont_vars                           = "cont1",
      exclude_vars_w_extreme_distribution = TRUE,
      verbose                             = TRUE
    )
  )
  combined <- paste(msgs, collapse = "")
  expect_match(combined, "Auto-excluding")
  expect_match(combined, "race_cat")
  expect_s3_class(result, "data.frame")
  expect_true("ps" %in% names(result))
})

test_that("estimate_ps errors on zero-count level when auto-exclude is off", {
  toy_df <- .toy_zero_count_df()

  # With auto-exclusion off, the same zero-count level must trigger an error
  # rather than slip through silently (this is the bug the union-of-levels
  # patch fixes).
  expect_error(
    estimate_ps(
      in_df                               = toy_df,
      exposure_var                        = "exposure",
      class_vars                          = c("good_cat", "race_cat"),
      cont_vars                           = "cont1",
      exclude_vars_w_extreme_distribution = FALSE,
      verbose                             = FALSE
    ),
    "extreme distribution"
  )
})
