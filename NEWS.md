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
