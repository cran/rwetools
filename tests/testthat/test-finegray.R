# Fine-Gray subdistribution-hazard regression (SHR), now via
# estimate_hr_ir(hr_model = "Fine-Gray") after the v0.4.0 move-in.
# finegray(weights = psweight) -> fgwt already = IPCW x PS weight; the Cox
# fit uses weights = fgwt with a robust SE clustered on the original subject
# (crude/weighted) or on the match id (matched block).

# Helper: build the multi-state status the package uses internally.
.fg_status <- function(df) {
  factor(ifelse(df$outcome == 1, 1L, ifelse(df$competing_event == 1, 2L, 0L)),
         levels = c(0L, 1L, 2L), labels = c("censor", "event", "compete"))
}

test_that("crude Fine-Gray SHR matches a direct finegray + coxph fit", {
  res <- estimate_hr_ir(
    in_df_crude = small_df, exposure_var = "exposure",
    outcome_var = "outcome", followuptime_var = "follow_up_days",
    hr_model = "Fine-Gray", if_fg_competing_event_var = "competing_event",
    verbose = FALSE
  )
  shr <- res$subdist_hazard
  expect_s3_class(shr, "data.frame")
  expect_identical(shr$Analysis, "Crude")

  # Manual replication (same row order -> same .fg_id clustering)
  df <- small_df
  df$.ms    <- .fg_status(df)
  df$.fg_id <- seq_len(nrow(df))
  fg  <- survival::finegray(survival::Surv(follow_up_days, .ms) ~ exposure + .fg_id,
                            data = df, etype = "event")
  fit <- survival::coxph(survival::Surv(fgstart, fgstop, fgstatus) ~ exposure,
                         data = fg, weights = fgwt, cluster = .fg_id, robust = TRUE)
  sm  <- summary(fit, conf.int = 0.95)

  expect_equal(shr$SHR,      unname(sm$conf.int[, "exp(coef)"]))
  expect_equal(shr$SHR_LCI,  unname(sm$conf.int[, 3]))
  expect_equal(shr$SHR_UCI,  unname(sm$conf.int[, 4]))
  expect_equal(shr$lnSHR_SE, unname(sm$coefficients[, "robust se"]))
})

test_that("weighted Fine-Gray multiplies PS weight into fgwt exactly once", {
  ps_df <- estimate_ps(
    in_df = small_df, exposure_var = "exposure",
    class_vars = class_vars, cont_vars = cont_vars, ps_var = "ps", verbose = FALSE
  )
  wt_df <- create_matching_weights(
    in_df = ps_df, exposure_var = "exposure", ps_var = "ps",
    weight_var = "mw_wt", verbose = FALSE
  )

  res <- estimate_hr_ir(
    in_df_weight = wt_df, if_weight_weight_var = "mw_wt",
    exposure_var = "exposure", outcome_var = "outcome",
    followuptime_var = "follow_up_days",
    hr_model = "Fine-Gray", if_fg_competing_event_var = "competing_event",
    verbose = FALSE
  )
  shr <- res$subdist_hazard
  expect_identical(shr$Analysis, "Weighted")

  # Manual: PS weights to finegray (fgwt = IPCW x mw_wt), Cox uses weights = fgwt
  df <- wt_df
  df$.ms    <- .fg_status(df)
  df$.fg_id <- seq_len(nrow(df))
  fg  <- survival::finegray(survival::Surv(follow_up_days, .ms) ~ exposure + .fg_id,
                            data = df, weights = df$mw_wt, etype = "event")
  fit <- survival::coxph(survival::Surv(fgstart, fgstop, fgstatus) ~ exposure,
                         data = fg, weights = fgwt, cluster = .fg_id, robust = TRUE)
  sm  <- summary(fit, conf.int = 0.95)

  expect_equal(shr$SHR,      unname(sm$conf.int[, "exp(coef)"]))
  expect_equal(shr$lnSHR_SE, unname(sm$coefficients[, "robust se"]))
})

test_that("matched Fine-Gray clusters the robust SE on the match id", {
  set.seed(4)
  md <- small_df[sample.int(nrow(small_df), 300), ]
  md$pair_id <- rep(seq_len(150), each = 2)

  res <- estimate_hr_ir(
    in_df_match = md, if_match_match_id = "pair_id",
    exposure_var = "exposure", outcome_var = "outcome",
    followuptime_var = "follow_up_days",
    hr_model = "Fine-Gray", if_fg_competing_event_var = "competing_event",
    verbose = FALSE
  )
  shr <- res$subdist_hazard

  # Manual: carry pair_id through the expansion and cluster on it
  df <- md
  df$.ms    <- .fg_status(df)
  df$.fg_id <- seq_len(nrow(df))
  fg  <- survival::finegray(
    survival::Surv(follow_up_days, .ms) ~ exposure + .fg_id + pair_id,
    data = df, etype = "event"
  )
  fit <- survival::coxph(survival::Surv(fgstart, fgstop, fgstatus) ~ exposure,
                         data = fg, weights = fgwt, cluster = pair_id,
                         robust = TRUE)
  sm  <- summary(fit, conf.int = 0.95)

  expect_equal(shr$SHR,      unname(sm$conf.int[, "exp(coef)"]))
  expect_equal(shr$lnSHR_SE, unname(sm$coefficients[, "robust se"]))
})

test_that("Fine-Gray rejects rows that are both event and competing event", {
  bad <- small_df
  bad$competing_event[which(bad$outcome == 1)[1]] <- 1   # one ambiguous row
  expect_error(
    estimate_hr_ir(
      in_df_crude = bad, exposure_var = "exposure",
      outcome_var = "outcome", followuptime_var = "follow_up_days",
      hr_model = "Fine-Gray", if_fg_competing_event_var = "competing_event",
      verbose = FALSE
    ),
    "ambiguous"
  )
})
