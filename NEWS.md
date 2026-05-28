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
