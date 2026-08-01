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

# ---------------------------------------------------------------------------
# NA handling (audit: "Fix estimate_ps() to handle NA explicitly").
# glm()'s previous default (na.omit) returned a predict() vector shorter than
# the data, so the PS column assignment crashed / mis-recycled. The fix uses
# na.action = na.exclude so predict() pads back to the full input length with
# NA in the dropped rows, preserving nrow and row alignment.
#
# Deterministic toy: 60 rows (exposure 30/30), well-overlapping cont1 and a
# balanced binary cat1 (no separation, no sparse level), with cont1 set to NA
# at rows 5 (exposed) and 40 (reference).
.toy_na_df <- function() {
  exposure <- c(rep(1L, 30), rep(0L, 30))
  cont1    <- c(seq(0.10, 3.00, length.out = 30),
                seq(0.05, 2.95, length.out = 30))
  cat1     <- rep(c(0L, 1L), times = 30)
  df <- data.frame(exposure = exposure, cont1 = cont1, cat1 = cat1,
                   stringsAsFactors = FALSE)
  df$cont1[c(5L, 40L)] <- NA  # inject missingness at known rows
  df
}

test_that("estimate_ps handles NA: nrow preserved and PS is NA on missing rows", {
  toy_df       <- .toy_na_df()
  missing_rows <- which(is.na(toy_df$cont1))   # c(5, 40)

  expect_warning(
    result <- estimate_ps(
      in_df        = toy_df,
      exposure_var = "exposure",
      class_vars   = "cat1",
      cont_vars    = "cont1",
      verbose      = FALSE
    ),
    "missing"
  )

  # Row count preserved
  expect_equal(nrow(result), nrow(toy_df))
  expect_length(result$ps, nrow(toy_df))
  # PS is NA exactly at the rows with a missing covariate
  expect_true(all(is.na(result$ps[missing_rows])))
  expect_false(any(is.na(result$ps[-missing_rows])))
  # All other PS values are finite probabilities
  expect_true(all(result$ps[-missing_rows] >= 0 & result$ps[-missing_rows] <= 1))
})

# ---------------------------------------------------------------------------
# Joint-design (quasi-)separation diagnostic
# (audit: "PS_estimation (a) joint-design separation check").
# The marginal Step-0 screen is blind to separation arising jointly across
# covariates. Cube design: the 8 corners of {-1,1}^3, exposure = 1 iff
# x1 + x2 + x3 > 0. Each covariate's marginal SMD is ~1.13 (< 1.5, so Step 0
# passes), but x1 + x2 + x3 separates exposure perfectly, so glm diverges and
# the post-fit diagnostic must flag it.
.toy_cube_separation_df <- function(reps = 6L) {
  corners <- expand.grid(x1 = c(-1, 1), x2 = c(-1, 1), x3 = c(-1, 1))
  s <- corners$x1 + corners$x2 + corners$x3
  corners$exposure <- ifelse(s > 0, 1L, 0L)
  df <- corners[rep(seq_len(nrow(corners)), each = reps),
                c("exposure", "x1", "x2", "x3")]
  rownames(df) <- NULL
  df
}

test_that("estimate_ps warns on joint separation (separation_action = 'warn')", {
  cube_df <- .toy_cube_separation_df()

  # Cube passes Step 0 (no extreme-distribution error) yet jointly separates.
  expect_warning(
    estimate_ps(
      in_df             = cube_df,
      exposure_var      = "exposure",
      cont_vars         = c("x1", "x2", "x3"),
      separation_action = "warn",
      verbose           = FALSE
    ),
    "separation"
  )
})

test_that("estimate_ps errors on joint separation (separation_action = 'error')", {
  cube_df <- .toy_cube_separation_df()

  expect_error(
    estimate_ps(
      in_df             = cube_df,
      exposure_var      = "exposure",
      cont_vars         = c("x1", "x2", "x3"),
      separation_action = "error",
      verbose           = FALSE
    ),
    "separation|near-positivity"
  )
})

test_that("estimate_ps stays silent on separation (separation_action = 'ignore')", {
  cube_df <- .toy_cube_separation_df()

  # 'ignore' must suppress BOTH the consolidated diagnostic and the underlying
  # base-R glm separation warnings (muffled via withCallingHandlers).
  expect_no_warning(
    estimate_ps(
      in_df             = cube_df,
      exposure_var      = "exposure",
      cont_vars         = c("x1", "x2", "x3"),
      separation_action = "ignore",
      verbose           = FALSE
    )
  )
})

test_that("estimate_ps does not warn on clean, non-separated data (backward compat)", {
  # Complete data with no separation must behave exactly as before: no missing
  # warning and no separation warning under the default action.
  expect_no_warning(
    estimate_ps(
      in_df        = small_df,
      exposure_var = "exposure",
      class_vars   = class_vars,
      cont_vars    = cont_vars,
      verbose      = FALSE
    )
  )
})
