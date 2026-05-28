# Suppress R CMD check NOTEs for ggplot2 non-standard evaluation
# These are column names used as bare aesthetics in ggplot2::aes() calls
utils::globalVariables(c(
  # rlang .data pronoun
  ".data",
  # create_psweights.R plot_data columns (PS distribution plots)
  "ps", "exposure", "weight", "group",
  # create_psweights.R (PS check / fs-weight plots)
  "exposure_factor", "exposure_group",
  # helpers_ps_psweights.R Love plot columns
  "Std_Diff", "Variable", "Type",
  # build_table1.R dplyr grouping variable
  ".grp",
  # create_ps_fs_weights.R / create_ps_matched_cohort.R histogram aesthetics
  "count", "density"
))
