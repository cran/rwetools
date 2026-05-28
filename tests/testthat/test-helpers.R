# --- fmt_num ---
test_that("fmt_num formats numbers correctly", {
  expect_equal(rwetools:::fmt_num(1.23456, digits = 3), "1.235")
  expect_equal(rwetools:::fmt_num(NA), "NA")
  expect_equal(rwetools:::fmt_num(-5.5, digits = 1, use_abs = TRUE), "5.5")
  expect_equal(rwetools:::fmt_num(0, digits = 2), "0.00")
})

# --- fmt_mean_sd ---
test_that("fmt_mean_sd formats mean (SD) correctly", {
  expect_equal(rwetools:::fmt_mean_sd(10.5, 2.3), "10.50 (2.30)")
  expect_equal(rwetools:::fmt_mean_sd(NA, 2.3), "NA")
  expect_equal(
    rwetools:::fmt_mean_sd(-3.1, 1.2, mean_digits = 1, sd_digits = 1, use_abs = TRUE),
    "3.1 (1.2)"
  )
})

# --- fmt_n_pct ---
test_that("fmt_n_pct formats count (percentage) correctly", {
  expect_equal(rwetools:::fmt_n_pct(25, 100), "25 (25.0%)")
  expect_equal(rwetools:::fmt_n_pct(NA, 100), "NA")
  expect_equal(rwetools:::fmt_n_pct(10, 0), "NA")
})

# --- fmt_pct ---
test_that("fmt_pct formats proportions as percentages", {
  expect_equal(rwetools:::fmt_pct(0.5), "50.0%")
  expect_equal(rwetools:::fmt_pct(NA), "--")
  expect_equal(rwetools:::fmt_pct(-0.25, digits = 2, use_abs = TRUE), "25.00%")
})

# --- canonize_levels ---
test_that("canonize_levels produces correct factors", {
  expect_equal(levels(rwetools:::canonize_levels(c(TRUE, FALSE, TRUE))), c("0", "1"))
  expect_equal(levels(rwetools:::canonize_levels(c(1, 2, 3))), c("1", "2", "3"))
  expect_equal(levels(rwetools:::canonize_levels(c("A", "B"))), c("A", "B"))
})

# --- levels_excl_zero ---
test_that("levels_excl_zero excludes baseline levels", {
  f <- factor(c("0", "1", "2"))
  expect_equal(rwetools:::levels_excl_zero(f), c("1", "2"))

  f2 <- factor(c("No", "Yes"))
  expect_equal(rwetools:::levels_excl_zero(f2), "Yes")
})

# --- get_var_types ---
test_that("get_var_types classifies variables from sample data", {
  vt <- rwetools:::get_var_types(small_df)

  # cont1 (many unique values) should be continuous
  expect_true("cont1" %in% vt$cont_vars)
  # binary1 (0/1) should be binary
  expect_true("binary1" %in% vt$binary_vars)
})

# --- summarize_ps_by_group ---
test_that("summarize_ps_by_group returns summary per group", {
  ps_vals <- runif(100)
  grp     <- rep(c(0, 1), each = 50)

  result <- rwetools:::summarize_ps_by_group(ps_vals, grp)

  expect_s3_class(result, "data.frame")
  expect_equal(nrow(result), 2)
  expect_true(all(c("mean_ps", "sd_ps", "min_ps", "max_ps") %in% names(result)))
})

# --- calc_c_statistic ---
test_that("calc_c_statistic returns c-stat and se", {
  set.seed(1)
  n <- 200
  exposure <- rep(c(0, 1), each = n / 2)
  ps <- ifelse(exposure == 1, rbeta(n / 2, 3, 2), rbeta(n / 2, 2, 3))

  result <- rwetools:::calc_c_statistic(exposure, ps, verbose = FALSE)

  expect_type(result, "list")
  expect_true("c_stat" %in% names(result))
  expect_true("se" %in% names(result))
  expect_true(result$c_stat >= 0 & result$c_stat <= 1)
})
