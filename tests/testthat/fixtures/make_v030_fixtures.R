# make_v030_fixtures.R -----------------------------------------------------
# Historical fixture-generation script written against the rwetools 0.3.0
# API. It is retained as provenance and is not sourced or run by the test suite.
#
# Step 0d: capture v0.3.0 numeric outputs (branch point of feature/v040-api-
# redesign, commit 1151194, tag v0.3.0) as .rds fixtures for the v0.4.0
# equivalence regression (EFFECT_CI_V2_PLAN.md Step 0 / Step 10).
#
# Captured (all analytical, fully deterministic -- no bootstrap anywhere):
#   hr_ir         : crude + IPTW-weighted IR/IRD, Cox HR (auto: model/robust),
#                   IRR quasipoisson (crude r10 + weighted r12)
#   hr_ir_strat   : stratified (direct-standardized IR/IRD, strata() Cox;
#                   IRR not allowed with stratify_by in 0.3.0)
#   rr_rd_km      : KM Risk/RR/RD analytical, crude + weighted
#   rr_rd_km_strat: KM stratified (standardized risk)
#   rr_rd_aj_fg   : AJ Risk/RR/RD analytical + Fine-Gray sHR, crude + weighted
#
# v0.4.0 must reproduce all point estimates exactly; the ONLY outputs allowed
# to differ are: weighted IR/IRD CI (-> robust), standardized IR arm CI
# (-> Fay-Feuer), Risk/CIF arm CI (-> cloglog).

pkg <- "D:/03_Dev/rwetools dev/rwetools"
devtools::load_all(pkg, quiet = TRUE)

df <- read.csv(file.path(pkg, "inst", "extdata", "sample_data.csv"),
               stringsAsFactors = FALSE)

cont_vars  <- paste0("cont", 1:7)
class_vars <- paste0("cat", 1:4)

# Deterministic PS + IPTW (ATE) weights
ps_df <- rwetools::estimate_ps(
  in_df = df, exposure_var = "exposure",
  class_vars = class_vars, cont_vars = cont_vars,
  ps_var = "ps", verbose = FALSE
)
wt_df <- rwetools::create_iptw(
  in_df = ps_df, exposure_var = "exposure",
  ps_var = "ps", weight_var = "iptw", verbose = FALSE
)

fixtures <- list()

# Full coxph objects are ~0.5 MB each (environments) -- too big to ship in the
# tarball. Keep coef + vcov instead; every formatted number lives in the
# output data frames anyway.
slim_models <- function(res) {
  for (m in c("unweighted_model", "weighted_model")) {
    if (!is.null(res[[m]])) {
      res[[paste0(m, "_coef")]] <- stats::coef(res[[m]])
      res[[paste0(m, "_vcov")]] <- stats::vcov(res[[m]])
      res[[m]] <- NULL
    }
  }
  res
}

fixtures$hr_ir <- rwetools::estimate_hr_ir(
  in_df_unwt = df, in_df_wted = wt_df,
  exposure_var = "exposure", outcome_var = "outcome",
  survival_time = "follow_up_days", weight_var = "iptw",
  time_unit = "days", ir_per_pyears = 1000,
  irr_method = "quasipoisson", verbose = FALSE
)
fixtures$hr_ir <- slim_models(fixtures$hr_ir)

fixtures$hr_ir_strat <- rwetools::estimate_hr_ir(
  in_df_unwt = df, in_df_wted = wt_df,
  exposure_var = "exposure", outcome_var = "outcome",
  survival_time = "follow_up_days", weight_var = "iptw",
  time_unit = "days", ir_per_pyears = 1000,
  stratify_by = "cat1", verbose = FALSE
)
fixtures$hr_ir_strat <- slim_models(fixtures$hr_ir_strat)

fixtures$rr_rd_km <- rwetools::estimate_rr_rd(
  in_df_unwt = df, in_df_wted = wt_df,
  exposure_var = "exposure", outcome_var = "outcome",
  survival_time = "follow_up_days", weight_var = "iptw",
  time_unit = "days", rr_rd_at_timepoint = 365, rr_rd_per_individuals = 1000,
  conf_int_method = "analytical", risk_estimator = "KM", verbose = FALSE
)

fixtures$rr_rd_km_strat <- rwetools::estimate_rr_rd(
  in_df_unwt = df, in_df_wted = wt_df,
  exposure_var = "exposure", outcome_var = "outcome",
  survival_time = "follow_up_days", weight_var = "iptw",
  time_unit = "days", rr_rd_at_timepoint = 365, rr_rd_per_individuals = 1000,
  conf_int_method = "analytical", risk_estimator = "KM",
  stratify_by = "cat1", verbose = FALSE
)

fixtures$rr_rd_aj_fg <- rwetools::estimate_rr_rd(
  in_df_unwt = df, in_df_wted = wt_df,
  exposure_var = "exposure", outcome_var = "outcome",
  survival_time = "follow_up_days", weight_var = "iptw",
  time_unit = "days", rr_rd_at_timepoint = 365, rr_rd_per_individuals = 1000,
  conf_int_method = "analytical", risk_estimator = "AJ",
  competing_event_var = "competing_event", fg_regression = TRUE,
  verbose = FALSE
)

fixtures$meta <- list(
  package_version = as.character(utils::packageVersion("rwetools")),
  source_commit   = "1151194 (tag v0.3.0)",
  r_version       = R.version.string,
  created         = format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"),
  weights         = "estimate_ps(cat1-4, cont1-7) + create_iptw(ATE, unstabilized, var 'iptw')",
  calls           = "see make_v030_fixtures.R alongside this file"
)

out_dir <- file.path(pkg, "tests", "testthat", "fixtures")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
out_rds <- file.path(out_dir, "v030_fixtures.rds")
saveRDS(fixtures, out_rds)  # regenerating this session's own artifact (slimmed)

cat("\n== fixture components ==\n")
print(names(fixtures))
for (nm in setdiff(names(fixtures), "meta")) {
  cat("\n--", nm, ":", paste(names(fixtures[[nm]]), collapse = ", "), "\n")
}
cat("\nSaved:", out_rds, "size:", file.size(out_rds), "bytes\n")
