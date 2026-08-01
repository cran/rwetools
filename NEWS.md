# rwetools 0.4.0

Complete redesign of the two effect-estimation functions,
`estimate_hr_ir()` and `estimate_rr_rd()` (API v2). This is a BREAKING
release for both functions; the propensity-score, weighting, matching, and
Table-1 functions are unchanged.

## Breaking API changes (both effect functions)

* Data blocks replace the old unweighted/weighted pair: `in_df_crude` plus
  either `in_df_weight` (with `if_weight_weight_var`) or `in_df_match`
  (optionally with `if_match_match_id`) — supplying both weighted and
  matched blocks is an error; weight- or match-only calls without a crude
  block are allowed. Output column suffixes are now `_Crude` / `_Weighted` /
  `_Matched` (previously `_Unwt` / `_Wted`).
* Renames: `survival_time` -> `followuptime_var`; `stratify_by` ->
  `stratification_var`; `rr_rd_per_individuals` -> `risk_per_individuals`;
  `weight_var` -> `if_weight_weight_var`.
* Analytical confidence intervals are ALWAYS computed and the CI-method
  selector arguments are gone: `hr_conf_int_method`, `cluster_var`,
  `conf_int_method`, `irr_method`, and `irr_covariates` have been removed.
  Each estimate uses one fixed, documented method (below).
* The percentile bootstrap is opt-in via `if_bootstrap_count` (plus
  `if_bootstrap_n_cores` / `if_bootstrap_seed`) and is always fixed-weight
  ("frozen"): analysis weights and standardization shares are held at their
  original values, so propensity-score-estimation uncertainty is not
  propagated (a message states this on every bootstrap run). Matched blocks
  resample matched sets by `if_match_match_id`. Only percentile intervals
  are provided (no BCa).
* `time_unit` gains `"months"` (person-years = months / 12) in
  `estimate_hr_ir()`.

## Fixed analytical methods (per estimate)

* Crude/matched IR: Garwood exact-Poisson (gamma) intervals (unchanged).
* Standardized IR arm CIs: Fay-Feuer (1997) gamma intervals replace the
  log-normal interval; a single stratum collapses exactly to the Garwood
  interval. The Total row keeps the crude Garwood interval; the IRD keeps
  the normal-Wald delta interval.
* Weighted IR / IRD: design-based robust (sandwich) variance from a joint
  weighted quasi-Poisson cell model (Lumley 2004), replacing the previous
  exact-Poisson interval on the rounded weighted event sum, which ignored
  the weighting variance and could severely under-cover under IPTW.
  Standardized weighted variants use the joint sandwich/delta method.
* HR (`estimate_hr_ir()`): crude block = Cox model SE; weighted block =
  robust sandwich SE; matched block = robust SE clustered on the match id
  (Lin-Wei 1989; Austin 2016). The old "auto" dispatch is now the only
  behavior.
* Fine-Gray subdistribution hazards MOVED from `estimate_rr_rd()`
  (`fg_regression` removed) into `estimate_hr_ir(hr_model = "Fine-Gray")`
  with `if_fg_competing_event_var`; weighted expansions keep
  fgwt = IPCW x PS weight, and matched blocks cluster the robust SE on the
  match id. `estimate_rr_rd()` keeps `risk_estimator = c("KM", "AJ")` with
  the competing-event indicator renamed to `if_aj_competing_event_var`
  (AJ-only).
* IRR is now always reported by `estimate_hr_ir()` as a marginal rate
  ratio: crude/matched blocks use a Poisson rate model (model SE, equal to
  sqrt(1/D1 + 1/D0)); the weighted block uses `survey::svyglm`
  quasi-Poisson (robust SE). The former quasi-Poisson option for
  unweighted fits is gone. Stratified IRRs, previously an error, are now
  direct-standardized.
* Risk / CIF arm CIs (`estimate_rr_rd()`): the log-scale interval is
  replaced by the complementary log-log (cloglog) interval applied to the
  final (standardized/weighted) estimate with its combined SE (Kalbfleisch
  & Prentice 2002). A crude Kaplan-Meier risk now reproduces
  `survfit(conf.type = "log-log")` exactly. A risk of exactly 0 or 1 has
  no defined cloglog interval: NA is returned and a message recommends the
  bootstrap. RD (normal-Wald) and RR (log-delta) intervals are unchanged.

## Other

* New internal analytical engines with unit tests: Fay-Feuer gamma CI,
  cloglog risk CI, joint-sandwich weighted rate cells with delta-method
  contrasts, direct-standardized IRR, and matched-set bootstrap
  resampling.
* v0.3.0 numeric fixtures ship in `tests/testthat/fixtures/` and guard the
  redesign: every method-unchanged output reproduces the 0.3.0 values
  exactly; only the weighted IR/IRD CIs (robust), standardized IR arm CIs
  (Fay-Feuer), Risk/CIF arm CIs (cloglog), and the crude IRR SE (Poisson)
  changed, each by design.
* The package now declares `Depends: R (>= 3.5.0)` (serialized test
  fixtures use serialization format version 3).

# rwetools 0.3.0

(Not submitted to CRAN; folded into 0.4.0.)

* BREAKING: `create_ps_weights()` has been removed. The single function with
  `weight_method` + `estimand` arguments is replaced by one function per
  estimand:
  - `create_iptw()` — inverse-probability weights for the ATE
    (`stabilize = TRUE` for stabilized IPTW).
  - `create_matching_weights()` — matching weights, now correctly labeled as
    targeting the ATM (Li & Greene 2013), not the ATT.
  - `create_overlap_weights()` — overlap weights for the ATO (Li et al. 2018).
  IPTW for the ATT (SMR weights) is not supported in this version. The broken
  stabilized-IPTW + ATO path has been removed.
* `create_iptw()`, `create_matching_weights()`, `create_overlap_weights()`, and
  `create_ps_matched_cohort()` gain optional propensity-score trimming
  (`trim_method`: `"crump"` symmetric or `"sturmer"` asymmetric). `create_iptw()`
  additionally gains weight truncation/winsorizing (`truncate_method`:
  `"percentile"` or `"cap"`). All default to off; both trimming and truncation
  change the analytic population and the interpretation of the estimand.
* `create_ps_matched_cohort()` gains a `method` argument
  (`"nearest"`/`"subclass"`), variable-ratio matching (`min_controls` /
  `max_controls` / `m_order`), and `subclass_n`. For `"subclass"` the
  returned cohort carries matching weights in `.match_weights`, which must be
  used for the matched Table 1 and any downstream effect estimation.
* Balance diagnostics and Love plots now use the raw (0-1) standardized-mean-
  difference scale consistently, fixing a bug where the "balanced" flag compared
  `|SMD|` against `threshold * 100` and the Love-plot reference line was drawn
  at 10 (off-screen) on raw-scale data.
* `estimate_hr_ir()` gained a `hr_conf_int_method` argument controlling how the
  hazard-ratio confidence interval (and the reported `lnHR_SE`) was computed:
  - `"auto"` (default): model-based SE for the unweighted fit and robust
    sandwich SE for the weighted fit (Austin 2016).
  - `"model"`: model-based (naive) Cox SE with a Wald CI for both fits.
  - `"robust"`: robust sandwich SE with a Wald CI for both fits
    (`"sandwich"` accepted as an alias).
  - `"bootstrap"`: percentile bootstrap CI with a bootstrap SE
    (`lnHR_SE = sd(log(HR))`); no p-value reported.
  (Removed again in 0.4.0: the "auto" dispatch is now the only behavior.)
* `estimate_hr_ir()` gained `cluster_var`, `bootstrap_count`, `n_cores`, and
  `seed` arguments for the robust/bootstrap methods (removed again in 0.4.0
  in favor of the fixed per-block methods and `if_bootstrap_*`).
* `estimate_hr_ir()` gained an optional incidence-rate-ratio rate model via
  `irr_method` (`"none"`/`"poisson"`/`"quasipoisson"`) with optional
  `irr_covariates`; weighted fits used `survey::svyglm(quasipoisson)`
  (reworked in 0.4.0: IRR is always-on and marginal).
* `estimate_rr_rd()` gained a Fine-Gray subdistribution-hazard regression via
  `fg_regression` + `competing_event_var` (moved to `estimate_hr_ir()` in
  0.4.0).
* `estimate_ps()` gained a separation diagnostic for the propensity-score
  logistic model and uses `na.exclude` so the returned PS column aligns with
  the input rows.

# rwetools 0.2.0

* `create_ps_weights()` and `create_ps_matched_cohort()` now also emit the
  unified PS-distribution figures (density, within-group histogram, and a
  histogram+density overlay for the unweighted/pre and weighted/post samples,
  plus a weight box plot for the weighting methods) that `create_ps_fs_weights()`
  already produced. The existing faceted / mirror plots are still written. A new
  internal helper `.plot_ps_distribution_set()` is shared by all three
  functions, so the figure style is now consistent across PS methods.
  `create_ps_fs_weights()` figure output is unchanged.
* `create_ps_matched_cohort()` gains a `caliper_scale` argument controlling the
  scale on which matching and the caliper are applied:
  `"logit_ps_sd"` matches on logit(PS) with a caliper of `caliper` x
  SD(logit(PS)) (Austin 2011); `"raw_ps_sd"` matches on PS with a caliper of
  `caliper` x SD(PS); `"raw"` applies a flat caliper on the raw PS scale
  (e.g. `caliper = 0.01`).
* BREAKING: the default matching caliper scale changed. rwetools <= 0.1.x
  matched on the raw PS with a standardized caliper (`caliper` x SD(PS)); the
  new default (`caliper_scale = "logit_ps_sd"`) matches on logit(PS)
  (`caliper` x SD(logit(PS))). Pass `caliper_scale = "raw_ps_sd"` to reproduce
  the previous results exactly. The verbose/report label that previously
  described the raw-PS caliper as "SD of logit(PS)" is now accurate for each
  scale.

# rwetools 0.1.2

* CRAN resubmission fixes:
  - DESCRIPTION: single-quoted `'Excel'` per CRAN software-name convention.
  - DESCRIPTION: added method references in the `Description` field
    (Rosenbaum and Rubin 1983 <doi:10.1093/biomet/70.1.41>;
    Austin 2011 <doi:10.1080/00273171.2011.568786>;
    Desai et al. 2017 <doi:10.1097/EDE.0000000000000595>).
  - `estimate_rr_rd()`: removed hardcoded RNG seed from the internal
    bootstrap helper. New `seed` argument (default `NULL`) lets the user
    opt in to reproducibility; when `NULL` the function does not touch
    the global RNG state. In the sequential branch, `.Random.seed` is
    saved and restored on exit when a seed is supplied.
* tests: fixed the zero-count-level regression test (test-only change).
  The 0.1.1 test used `expect_message(..., "Auto-excluding.*race_cat")`,
  which required both substrings in a single `message()` call, but the
  verbose path emits the header and the variable name on separate
  lines. Switched to `capture_messages()` and asserts both substrings
  individually.

# rwetools 0.1.1

* fix: `estimate_ps()` categorical extreme-distribution screen now detects
  zero-count levels via union of levels across exposure groups. Previously,
  a level present in only one exposure group was silently dropped by
  `table()` and could slip past the `cnt < 5` rule, triggering
  quasi-separation in `glm()`.

# rwetools 0.1.0

* Initial CRAN release.
