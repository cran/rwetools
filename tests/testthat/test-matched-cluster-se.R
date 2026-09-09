# Backlog Item 2: matched-block IR / IRD / IRR analytical SEs are cluster-robust
# on the match id.
#
# The acceptance criteria deliberately do NOT assert a direction. A cluster SE
# is not necessarily smaller than the independence SE: it shrinks under
# positive within-set covariance and grows under negative covariance. What is
# asserted instead is equivalence -- to a directly specified
# survey::svydesign(ids = ~match_id) fit, and to hand-computed algebra -- plus
# preservation of the point estimates.

sample_df <- function() {
  read.csv(system.file("extdata", "sample_data.csv", package = "rwetools"))
}

# a 1:1 matched cohort with all-1 matching weights, so only the clustering
# distinguishes it from an unmatched analysis
paired <- function(n_pairs = 400, seed = 11, rho = c("pos", "neg", "none")) {
  rho <- match.arg(rho)
  set.seed(seed)
  df <- sample_df()
  d <- df[sample.int(nrow(df), 2 * n_pairs), ]
  d$pair_id <- rep(sprintf("P%04d", seq_len(n_pairs)), each = 2)
  d$exposure <- rep(c(0L, 1L), times = n_pairs)
  # induce a within-pair correlation in the outcome
  pair_u <- stats::rbinom(n_pairs, 1, 0.25)
  d$outcome <- switch(rho,
    pos  = rep(pair_u, each = 2),                       # both members agree
    neg  = as.integer(c(rbind(pair_u, 1L - pair_u))),   # members disagree
    none = stats::rbinom(2 * n_pairs, 1, 0.25)
  )
  d$follow_up_days <- pmax(30, d$follow_up_days)
  d
}

ir_matched <- function(d, with_id = TRUE) suppressMessages(estimate_ir(
  in_df_match = d,
  if_match_match_id = if (with_id) "pair_id" else NULL,
  exposure_var = "exposure", outcome_var = "outcome",
  followuptime_var = "follow_up_days", time_unit = "days",
  ir_per_pyears = 1000, verbose = FALSE))

test_that("matched IRR equals a direct svydesign(ids = ~match_id) svyglm fit", {
  d <- paired(rho = "pos")
  res <- ir_matched(d)
  irr <- res$incidence_rate_ratios

  d$person_years <- d$follow_up_days / 365.25
  # suppressWarnings: svydesign notes the absence of weights; unit weights
  # are exactly what the package's own unweighted design uses
  des <- suppressWarnings(survey::svydesign(ids = ~pair_id, data = d))
  fit <- survey::svyglm(outcome ~ exposure + offset(log(person_years)),
                        design = des, family = stats::quasipoisson())
  cf <- summary(fit)$coefficients["exposure", ]

  expect_equal(irr$lnIRR, unname(cf[1]), tolerance = 1e-10)
  expect_equal(irr$lnIRR_SE, unname(cf[2]), tolerance = 1e-10)
  expect_equal(irr$IRR, exp(unname(cf[1])), tolerance = 1e-10)
  expect_match(irr$SE_type, "match-id cluster")

  # ids = ~1 is the no-cluster declaration and must NOT reproduce it
  des0 <- suppressWarnings(survey::svydesign(ids = ~1, data = d))
  fit0 <- survey::svyglm(outcome ~ exposure + offset(log(person_years)),
                         design = des0, family = stats::quasipoisson())
  expect_false(isTRUE(all.equal(irr$lnIRR_SE,
                               unname(summary(fit0)$coefficients["exposure", 2]),
                               tolerance = 1e-6)))
})

test_that("matched IR / IRD equal the joint cluster-robust cell algebra", {
  d <- paired(rho = "pos")
  res <- ir_matched(d)$incidence_rates

  d$person_years <- d$follow_up_days / 365.25
  d$.cell <- factor(d$exposure)
  des <- suppressWarnings(survey::svydesign(ids = ~pair_id, data = d))
  fit <- survey::svyglm(outcome ~ 0 + .cell + offset(log(person_years)),
                        design = des, family = stats::quasipoisson())
  b <- stats::coef(fit); V <- stats::vcov(fit)
  r <- 1000 * exp(b)                     # per-1000 arm rates
  # arm CI: robust log-Wald
  z <- stats::qnorm(0.975)
  se_log <- sqrt(diag(V))
  lci <- r * exp(-z * se_log); uci <- r * exp(z * se_log)
  # IRD: joint delta on the identity scale
  g <- c(-r[[".cell0"]], r[[".cell1"]])   # d/dbeta of (r1 - r0)
  se_ird <- sqrt(drop(t(g) %*% V[c(".cell0", ".cell1"), c(".cell0", ".cell1")] %*% g))

  ex <- res$Exposure_Group
  expect_equal(res$IR_per_Npy_Matched[ex == "Reference"], unname(r[".cell0"]),
               tolerance = 1e-8)
  expect_equal(res$IR_per_Npy_Matched[ex == "Exposed"], unname(r[".cell1"]),
               tolerance = 1e-8)
  expect_equal(res$IR_LCI_Matched[ex == "Exposed"], unname(lci[".cell1"]),
               tolerance = 1e-8)
  expect_equal(res$IR_UCI_Matched[ex == "Reference"], unname(uci[".cell0"]),
               tolerance = 1e-8)
  expect_equal(unique(res$IRD_SE_Matched), se_ird, tolerance = 1e-8)
})

test_that("clustering moves only the intervals, never the point estimates", {
  for (rho in c("pos", "neg", "none")) {
    d <- paired(rho = rho)
    with_id <- ir_matched(d, with_id = TRUE)
    no_id   <- ir_matched(d, with_id = FALSE)
    # point estimates identical
    expect_equal(with_id$incidence_rates$IR_per_Npy_Matched,
                 no_id$incidence_rates$IR_per_Npy_Matched,
                 tolerance = 1e-10, info = rho)
    expect_equal(with_id$incidence_rates$IRD_per_Npy_Matched,
                 no_id$incidence_rates$IRD_per_Npy_Matched,
                 tolerance = 1e-10, info = rho)
    expect_equal(with_id$incidence_rate_ratios$IRR,
                 no_id$incidence_rate_ratios$IRR,
                 tolerance = 1e-6, info = rho)
    # intervals differ
    expect_false(isTRUE(all.equal(with_id$incidence_rate_ratios$lnIRR_SE,
                                  no_id$incidence_rate_ratios$lnIRR_SE,
                                  tolerance = 1e-6)), info = rho)
  }
})

test_that("the cluster SE is NOT required to be smaller (both signs occur)", {
  # This is the corrected acceptance criterion. The original plan asserted the
  # matched CI could only narrow; that is false. Under positive within-pair
  # covariance the IRD/IRR cluster SE shrinks, under negative covariance it
  # grows. Both directions are pinned here so neither is mistaken for a bug.
  se <- function(rho) {
    d <- paired(rho = rho)
    c(clustered = ir_matched(d, TRUE)$incidence_rate_ratios$lnIRR_SE,
      independent = ir_matched(d, FALSE)$incidence_rate_ratios$lnIRR_SE)
  }
  s_pos <- se("pos")
  s_neg <- se("neg")
  expect_lt(s_pos[["clustered"]], s_pos[["independent"]])
  expect_gt(s_neg[["clustered"]], s_neg[["independent"]])
})

test_that("the Total row gets a real CI (intercept-only cell fit, Item 9)", {
  # The Total row is fitted on a constant "cell" column. Up to 0.4.0 that went
  # through `~ 0 + .cell`, which errors on a one-level factor ("contrasts can
  # be applied only to factors with 2 or more levels"); the caller's tryCatch
  # swallowed it, so the WEIGHTED block's Total-row IR CI was silently NA in
  # every release. A single-level cell now uses a plain intercept.
  d <- paired(rho = "pos")
  res <- ir_matched(d)$incidence_rates
  tot <- res[res$Exposure_Group == "Total", ]
  expect_true(is.finite(tot$IR_LCI_Matched))
  expect_true(is.finite(tot$IR_UCI_Matched))
  expect_true(is.finite(tot$IR_SE_Matched))
  expect_lt(tot$IR_LCI_Matched, tot$IR_per_Npy_Matched)
  expect_gt(tot$IR_UCI_Matched, tot$IR_per_Npy_Matched)

  # and it equals a direct clustered intercept fit
  d$person_years <- d$follow_up_days / 365.25
  des <- suppressWarnings(survey::svydesign(ids = ~pair_id, weights = ~1, data = d))
  fit <- survey::svyglm(outcome ~ 1 + offset(log(person_years)),
                        design = des, family = stats::quasipoisson())
  # vcov() of a one-parameter fit is a 1x1 matrix; drop() so the comparison is
  # against a plain numeric
  b <- stats::coef(fit); v <- drop(stats::vcov(fit)); z <- stats::qnorm(0.975)
  expect_equal(tot$IR_per_Npy_Matched, unname(1000 * exp(b)), tolerance = 1e-8)
  expect_equal(tot$IR_LCI_Matched, unname(1000 * exp(b - z * sqrt(v))),
               tolerance = 1e-8)
  expect_equal(tot$IR_UCI_Matched, unname(1000 * exp(b + z * sqrt(v))),
               tolerance = 1e-8)
})

test_that("the weighted block's Total row is no longer NA (Item 9)", {
  df <- sample_df()
  ps <- estimate_ps(in_df = df, exposure_var = "exposure",
                    class_vars = paste0("cat", 1:4),
                    cont_vars = paste0("cont", 1:7),
                    ps_var = "ps", verbose = FALSE)
  wt <- create_iptw(in_df = ps, exposure_var = "exposure", ps_var = "ps",
                    weight_var = "iptw", verbose = FALSE)
  res <- suppressMessages(estimate_ir(
    in_df_weight = wt, if_weight_weight_var = "iptw",
    exposure_var = "exposure", outcome_var = "outcome",
    followuptime_var = "follow_up_days", time_unit = "days",
    ir_per_pyears = 1000, verbose = FALSE))$incidence_rates
  tot <- res[res$Exposure_Group == "Total", ]
  expect_true(is.finite(tot$IR_LCI_Weighted))
  expect_true(is.finite(tot$IR_UCI_Weighted))
  expect_true(is.finite(tot$IR_SE_Weighted))
})

test_that("matched block without a match id keeps the independence methods", {
  d <- paired(rho = "pos")
  res <- ir_matched(d, with_id = FALSE)
  expect_match(res$incidence_rate_ratios$SE_type, "model \\(Poisson")
  # Garwood arm CI: qgamma on the raw counts
  ex <- res$incidence_rates$Exposure_Group
  n_ev <- res$incidence_rates$N_Events_Matched[ex == "Exposed"]
  py   <- res$incidence_rates$PersonYears_Matched[ex == "Exposed"]
  expect_equal(res$incidence_rates$IR_LCI_Matched[ex == "Exposed"],
               1000 * stats::qgamma(0.025, n_ev) / py, tolerance = 1e-10)
})
