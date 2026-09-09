# Tests for the always-on marginal IRR in the v0.4.0 estimate_ir.
# Fixed methods (no user args): crude/matched = marginal Poisson rate model
# (model SE; matrix r10-A), weighted = svyglm quasipoisson robust (r12-A),
# stratified = direct-standardized (r11-A / r13-A).

fx <- readRDS(testthat::test_path("fixtures", "v030_fixtures.rds"))

sample_df <- function() {
  read.csv(system.file("extdata", "sample_data.csv", package = "rwetools"))
}
iptw_df <- function(df) {
  ps <- estimate_ps(in_df = df, exposure_var = "exposure",
                    class_vars = paste0("cat", 1:4),
                    cont_vars = paste0("cont", 1:7),
                    ps_var = "ps", verbose = FALSE)
  create_iptw(in_df = ps, exposure_var = "exposure", ps_var = "ps",
              weight_var = "iptw", verbose = FALSE)
}

test_that("IRR is always reported; crude = Poisson model SE (r10-A)", {
  df <- sample_df()
  res <- estimate_ir(in_df_crude = df, exposure_var = "exposure",
                        outcome_var = "outcome",
                        followuptime_var = "follow_up_days", verbose = FALSE)
  irr <- res$incidence_rate_ratios
  expect_s3_class(irr, "data.frame")
  expect_identical(irr$Analysis, "Crude")
  expect_identical(irr$Model, "poisson")

  # equals a direct marginal Poisson rate model
  df2 <- df; df2$person_years <- df2$follow_up_days / 365.25
  m <- stats::glm(outcome ~ exposure + offset(log(person_years)),
                  data = df2, family = stats::poisson())
  sm <- summary(m)
  expect_equal(irr$lnIRR,    unname(stats::coef(m)["exposure"]))
  expect_equal(irr$lnIRR_SE, unname(sm$coefficients["exposure", "Std. Error"]))

  # point == v0.3.0 fixture; SE ~= sqrt(1/D1 + 1/D0) (A == B up to IRLS tol)
  firr <- fx$hr_ir$incidence_rate_ratios
  fir  <- fx$hr_ir$incidence_rates
  expect_equal(irr$IRR, firr$IRR[1], tolerance = 1e-10)
  D1 <- fir$N_Events_Unwt[2]; D0 <- fir$N_Events_Unwt[1]
  expect_equal(irr$lnIRR_SE, sqrt(1 / D1 + 1 / D0), tolerance = 2e-4)
})
test_that("weighted IRR = svyglm quasipoisson robust, fixture-exact", {
  df <- sample_df(); wt <- iptw_df(df)
  res <- estimate_ir(in_df_crude = df, in_df_weight = wt,
                        if_weight_weight_var = "iptw",
                        exposure_var = "exposure", outcome_var = "outcome",
                        followuptime_var = "follow_up_days", verbose = FALSE)
  irr <- res$incidence_rate_ratios
  expect_identical(irr$Analysis, c("Crude", "Weighted"))
  expect_match(irr$SE_type[2], "svyglm")
  firr <- fx$hr_ir$incidence_rate_ratios
  expect_equal(irr$IRR[2],      firr$IRR[2],      tolerance = 1e-10)
  expect_equal(irr$lnIRR_SE[2], firr$lnIRR_SE[2], tolerance = 1e-10)
})
