# PS HELPERS (used by multiple PS functions) ##########
#### Summarize PS distribution by exposure group) ##########
#' Summarize propensity-score distribution by exposure group
#'
#' @param ps_vals Numeric vector of propensity scores
#' @param group_vals Vector of group indicators (0/1)
#' @param group_labels Named character vector, e.g. c("0"="Control","1"="Treated"). 
#'   If NULL, uses numeric group values.
#' @return data.frame with one row per group containing n, mean, SD, min,
#'   quartiles, and max of the propensity score.
#' @keywords internal
#' @noRd
summarize_ps_by_group <- function(ps_vals, group_vals, group_labels = NULL) {
  groups <- sort(unique(group_vals[!is.na(group_vals)]))
  result <- do.call(rbind, lapply(groups, function(g) {
    ps_g <- ps_vals[group_vals == g & !is.na(group_vals)]
    data.frame(
      group = if (!is.null(group_labels)) group_labels[as.character(g)] else g,
      n = length(ps_g),
      mean_ps = mean(ps_g, na.rm = TRUE),
      sd_ps = stats::sd(ps_g, na.rm = TRUE),
      min_ps = min(ps_g, na.rm = TRUE),
      q25_ps = stats::quantile(ps_g, 0.25, na.rm = TRUE),
      median_ps = stats::median(ps_g, na.rm = TRUE),
      q75_ps = stats::quantile(ps_g, 0.75, na.rm = TRUE),
      max_ps = max(ps_g, na.rm = TRUE),
      stringsAsFactors = FALSE
    )
  }))
  rownames(result) <- NULL
  result
}

#### Calculate C-statistic (concordance) for PS model discrimination ####
#' Calculate C-statistic (concordance) for PS model discrimination
#'
#' Uses \code{survival::concordance} which is O(n log n) and supports sample
#' weights natively.
#'
#' @param exposure_vec Binary exposure vector (0/1)
#' @param ps_vec Propensity score vector
#' @param weights_vec Optional numeric weights vector (for post-weighting / post-matching c-stat)
#' @param label Character label printed in console (default "PS Model")
#' @param verbose Logical. Print progress messages (default TRUE).
#' @return A list with elements \code{c_stat} (numeric) and \code{se} (numeric),
#'         or NULL if survival is not available.
#' @keywords internal
#' @noRd
calc_c_statistic <- function(exposure_vec, ps_vec, weights_vec = NULL, label = "PS Model",
                             verbose = TRUE) {
  if (!requireNamespace("survival", quietly = TRUE)) {
    if (verbose) message("\nNote: Package 'survival' is required to calculate C-statistic")
    return(NULL)
  }
  tryCatch({
    df_tmp <- data.frame(y = as.numeric(exposure_vec), x = as.numeric(ps_vec))
    if (!is.null(weights_vec)) {
      conc <- survival::concordance(y ~ x, data = df_tmp,
                                    weights = as.numeric(weights_vec))
    } else {
      conc <- survival::concordance(y ~ x, data = df_tmp)
    }
    c_stat <- as.numeric(conc$concordance)
    c_se   <- sqrt(as.numeric(conc$var))
    if (verbose) message(sprintf("\n%s C-statistic: %.3f (SE: %.4f)", label, c_stat, c_se))
    list(c_stat = c_stat, se = c_se)
  }, error = function(e) {
    if (verbose) message(sprintf("\nNote: Could not compute %s C-statistic: %s", label, conditionMessage(e)))
    NULL
  })
}

#### Add a README sheet to an openxlsx workbook ####
#' Add a README sheet to an openxlsx workbook
#'
#' @param wb An openxlsx workbook object
#' @param readme_text Character string (may contain newlines)
#' @param verbose Logical. Print progress messages (default TRUE).
#' @return Called for its side effect (adds a sheet to \code{wb}). Returns
#'   \code{NULL} invisibly.
#' @keywords internal
#' @noRd
add_readme_sheet <- function(wb, readme_text,
                             verbose = TRUE) {
  openxlsx::addWorksheet(wb, "README")
  readme_lines <- data.frame(V1 = strsplit(readme_text, "\n")[[1]])
  openxlsx::writeData(wb, "README", readme_lines, colNames = FALSE, startCol = 1, startRow = 1)
  openxlsx::setColWidths(wb, "README", cols = 1, widths = 80)
  if (verbose) message("  Added README sheet")
}

#### Create a balance/Love plot comparing standardized differences #####
#' Create a balance (Love) plot comparing standardized differences
#'
#' @param variable_names Character vector of variable names
#' @param crude_std_diff Numeric vector of crude/pre standardized differences
#' @param adjusted_std_diff Numeric vector of adjusted/post standardized differences
#' @param crude_label Label for crude/pre group (e.g. "Unweighted", "Pre-matching", "Crude")
#' @param adjusted_label Label for adjusted/post group (e.g. "Weighted", "Post-matching")
#' @param title Plot title
#' @param output_path Full file path for saved plot
#' @param colors Named character vector of length 2 with crude_label and adjusted_label as names
#' @param shapes Optional named integer vector of length 2 with crude_label and
#'   adjusted_label as names giving ggplot2 shape codes (e.g., 16 = filled circle,
#'   17 = filled triangle). If NULL, ggplot2 default shapes are used.
#' @param use_absolute Logical. Plot absolute values of std diff? (default FALSE)
#' @param std_diff_threshold Numeric. Balance threshold on the raw standardized
#'   difference scale (default 0.1) at which the dashed reference line(s) are
#'   drawn. Inputs are expected on the raw (0-1) SMD scale produced by
#'   \code{\link{build_table1}}, not a percentage.
#' @return Invisibly returns a \code{ggplot} object, or \code{invisible(NULL)}
#'   if \pkg{ggplot2} is unavailable or no data remain after removing NAs.
#' @section Side Effects:
#' Saves the plot to \code{output_path} via \code{ggplot2::ggsave}.
#' @keywords internal
#' @noRd
create_love_plot <- function(variable_names, crude_std_diff, adjusted_std_diff,
                             crude_label, adjusted_label, title, output_path,
                             colors = NULL, shapes = NULL, use_absolute = FALSE,
                             std_diff_threshold = 0.1) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    warning("Package 'ggplot2' is required for Love plot. Skipping.")
    return(invisible(NULL))
  }
  
  if (is.null(colors)) {
    colors <- stats::setNames(c("#377eb8", "#e41a1c"), c(crude_label, adjusted_label))
  }
  if (length(variable_names) != length(crude_std_diff) ||
      length(crude_std_diff) != length(adjusted_std_diff)) {
    stop("variable_names, crude_std_diff, and adjusted_std_diff must have the same length")
  }
  
  # Build long-format data
  balance_long <- rbind(
    data.frame(Variable = variable_names, Std_Diff = crude_std_diff,
               Type = crude_label, stringsAsFactors = FALSE),
    data.frame(Variable = variable_names, Std_Diff = adjusted_std_diff,
               Type = adjusted_label, stringsAsFactors = FALSE)
  )
  balance_long <- balance_long[!is.na(balance_long$Std_Diff), , drop = FALSE]
  
  if (nrow(balance_long) == 0) return(invisible(NULL))
  
  # Order by absolute crude std diff
  var_order <- variable_names[order(abs(crude_std_diff), na.last = TRUE)]
  balance_long$Variable <- factor(balance_long$Variable, levels = var_order)
  balance_long$Type <- factor(balance_long$Type, levels = c(crude_label, adjusted_label))
  
  n_vars <- length(unique(balance_long$Variable))
  plot_height <- max(8, min(20, n_vars * 0.25))
  
  if (use_absolute) {
    balance_long$Std_Diff <- abs(balance_long$Std_Diff)
    x_lab <- "Absolute Standardized Difference"
    vline_x <- std_diff_threshold
  } else {
    x_lab <- "Standardized Difference"
    vline_x <- c(-std_diff_threshold, std_diff_threshold)
  }
  
  p <- ggplot2::ggplot(balance_long,
                       ggplot2::aes(x = Std_Diff, y = Variable, color = Type, shape = Type)) +
    ggplot2::geom_point(size = 2, alpha = 0.8) +
    ggplot2::geom_vline(xintercept = vline_x, linetype = "dashed", color = "red", alpha = 0.5) +
    ggplot2::geom_vline(xintercept = 0, linetype = "solid", color = "gray50") +
    ggplot2::labs(title = title, x = x_lab, y = "Variable",
                  color = "Sample", shape = "Sample") +
    ggplot2::theme_minimal() +
    ggplot2::theme(legend.position = "bottom",
                   axis.text.y = ggplot2::element_text(size = max(6, min(10, 300 / n_vars)))) +
    ggplot2::scale_color_manual(values = colors)

  if (!is.null(shapes)) {
    p <- p + ggplot2::scale_shape_manual(values = shapes)
  }

  ggplot2::ggsave(output_path, plot = p, width = 10, height = plot_height, dpi = 150)
  invisible(p)
}

#### Emit the unified FS-style PS-distribution plot set #####
#' Emit the unified propensity-score distribution diagnostic plot set
#'
#' Internal helper that renders the standard rwetools propensity-score (PS)
#' distribution figures - a density plot, a within-group histogram, and a
#' histogram+density overlay - for a "panel 1" (unweighted / pre) sample and a
#' "panel 2" (weighted / post) sample, plus an optional weight box plot. Shared
#' by the weighting functions ([create_iptw()], [create_matching_weights()],
#' [create_overlap_weights()]), [create_ps_matched_cohort()], and
#' [create_ps_fs_weights()] so that all emit a consistent figure style.
#'
#' Exposure groups are coloured red (\code{"Reference"}, exposure 0) and blue
#' (\code{"Exposure"}, exposure 1).
#'
#' @param crude_df data.frame for panel 1 (unweighted / pre). Must contain
#'   \code{ps_var} and a 0/1 \code{exposure_var} column.
#' @param weighted_df data.frame for panel 2 (weighted / post). Must contain the
#'   same columns and, when \code{weight_var} is supplied, the weight column.
#' @param ps_var Character. Name of the propensity-score column.
#' @param exposure_var Character. Name of the 0/1 exposure column.
#' @param weight_var Character or \code{NULL}. Weight column in
#'   \code{weighted_df}. When \code{NULL}, the panel-2 plots are drawn as
#'   unweighted counts of \code{weighted_df} (e.g. a matched subset, where each
#'   retained subject carries weight 1) and the box plot is skipped.
#' @param out_dir_plots Character. Output directory (created if needed).
#' @param plot_prefix Character. File-name prefix for the saved PNGs.
#' @param unwt_title,wt_title Character. Title stems for panel 1 / panel 2. The
#'   density and histogram variants append \code{" - Density"} /
#'   \code{" - Histogram"}.
#' @param box_title Character. Title for the weight box plot.
#' @param panel1_suffix,panel2_suffix Character. File-name tokens for the two
#'   panels (default \code{"unwt"} / \code{"wt"}).
#' @param make_boxplot Logical. Draw the weight box plot? (ignored when
#'   \code{weight_var} is \code{NULL}).
#' @param verbose Logical. Print progress messages (default \code{FALSE}).
#' @return Invisibly \code{NULL}. Called for the PNG files written under
#'   \code{out_dir_plots}.
#' @keywords internal
#' @noRd
.plot_ps_distribution_set <- function(crude_df, weighted_df, ps_var, exposure_var,
                                      weight_var = NULL, out_dir_plots, plot_prefix,
                                      unwt_title = "Unweighted PS Distribution",
                                      wt_title   = "Weighted PS Distribution",
                                      box_title  = "Distribution of Propensity Score Weights",
                                      panel1_suffix = "unwt", panel2_suffix = "wt",
                                      make_boxplot = TRUE, verbose = FALSE) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    warning("Package 'ggplot2' is required for plots. Skipping plot generation.")
    return(invisible(NULL))
  }
  if (!dir.exists(out_dir_plots)) dir.create(out_dir_plots, recursive = TRUE)

  # Colour palette
  fill_vals  <- c("0" = "#e41a1c", "1" = "#377eb8")
  fill_labs  <- c("0" = "Reference", "1" = "Exposure")

  # ---- Panel 1: unweighted / pre ----
  if (verbose) message("  Creating panel-1 PS distribution plots...")
  unwt_data <- crude_df
  unwt_data$exposure_factor <- factor(unwt_data[[exposure_var]], levels = c(0, 1), labels = c("0", "1"))

  # density
  p <- ggplot2::ggplot(unwt_data, ggplot2::aes(x = .data[[ps_var]], fill = exposure_factor, color = exposure_factor)) +
    ggplot2::geom_density(alpha = 0.3, adjust = 1.5) +
    ggplot2::scale_fill_manual(values = fill_vals, labels = fill_labs) +
    ggplot2::scale_color_manual(values = fill_vals, labels = fill_labs) +
    ggplot2::labs(title = paste0(unwt_title, " - Density"), x = "PS", y = "Density") +
    ggplot2::theme_minimal() + ggplot2::theme(plot.title = ggplot2::element_text(size = 14, face = "bold", hjust = 0.5),
                                              legend.position = "top") + ggplot2::xlim(0, 1)
  ggplot2::ggsave(file.path(out_dir_plots, paste0(plot_prefix, "_ps_distr_density_", panel1_suffix, ".png")),
                  plot = p, width = 8, height = 6, dpi = 150)

  # histogram
  p <- ggplot2::ggplot(unwt_data, ggplot2::aes(x = .data[[ps_var]],
                                               y = ggplot2::after_stat(count / tapply(count, group, sum)[group] * 100),
                                               fill = exposure_factor)) +
    ggplot2::geom_histogram(alpha = 0.5, position = "identity", bins = 100) +
    ggplot2::scale_fill_manual(values = fill_vals, labels = fill_labs) +
    ggplot2::labs(title = paste0(unwt_title, " - Histogram"), x = "PS", y = "% patients (within group)") +
    ggplot2::theme_minimal() + ggplot2::theme(plot.title = ggplot2::element_text(size = 14, face = "bold", hjust = 0.5),
                                              legend.position = "top") + ggplot2::xlim(0, 1)
  ggplot2::ggsave(file.path(out_dir_plots, paste0(plot_prefix, "_ps_distr_histog_", panel1_suffix, ".png")),
                  plot = p, width = 8, height = 6, dpi = 150)

  # both
  p <- ggplot2::ggplot(unwt_data, ggplot2::aes(x = .data[[ps_var]])) +
    ggplot2::geom_histogram(ggplot2::aes(y = ggplot2::after_stat(density), fill = exposure_factor),
                            alpha = 0.5, position = "identity", bins = 100) +
    ggplot2::geom_density(ggplot2::aes(y = ggplot2::after_stat(density), color = exposure_factor), linewidth = 1.2, adjust = 1.5) +
    ggplot2::scale_fill_manual(values = fill_vals, labels = fill_labs) +
    ggplot2::scale_color_manual(values = fill_vals, labels = fill_labs, guide = "none") +
    ggplot2::labs(title = unwt_title, x = "PS", y = "Density") +
    ggplot2::theme_minimal() + ggplot2::theme(plot.title = ggplot2::element_text(size = 14, face = "bold", hjust = 0.5),
                                              legend.position = "top") + ggplot2::xlim(0, 1)
  ggplot2::ggsave(file.path(out_dir_plots, paste0(plot_prefix, "_ps_distr_", panel1_suffix, ".png")),
                  plot = p, width = 8, height = 6, dpi = 150)

  # ---- Panel 2: weighted / post ----
  if (verbose) message("  Creating panel-2 PS distribution plots...")
  wt_data <- weighted_df
  wt_data$exposure_factor <- factor(wt_data[[exposure_var]], levels = c(0, 1), labels = c("0", "1"))
  has_w <- !is.null(weight_var)

  # density
  aes_density <- if (has_w) {
    ggplot2::aes(x = .data[[ps_var]], weight = .data[[weight_var]], fill = exposure_factor, color = exposure_factor)
  } else {
    ggplot2::aes(x = .data[[ps_var]], fill = exposure_factor, color = exposure_factor)
  }
  p <- ggplot2::ggplot(wt_data, aes_density) +
    ggplot2::geom_density(alpha = 0.3, adjust = 1.5) +
    ggplot2::scale_fill_manual(values = fill_vals, labels = fill_labs) +
    ggplot2::scale_color_manual(values = fill_vals, labels = fill_labs) +
    ggplot2::labs(title = paste0(wt_title, " - Density"), x = "PS", y = "Density") +
    ggplot2::theme_minimal() + ggplot2::theme(plot.title = ggplot2::element_text(size = 14, face = "bold", hjust = 0.5),
                                              legend.position = "top") + ggplot2::xlim(0, 1)
  ggplot2::ggsave(file.path(out_dir_plots, paste0(plot_prefix, "_ps_distr_density_", panel2_suffix, ".png")),
                  plot = p, width = 8, height = 6, dpi = 150)

  # histogram
  aes_histog <- if (has_w) {
    ggplot2::aes(x = .data[[ps_var]], weight = .data[[weight_var]], fill = exposure_factor)
  } else {
    ggplot2::aes(x = .data[[ps_var]], fill = exposure_factor)
  }
  p <- ggplot2::ggplot(wt_data, aes_histog) +
    ggplot2::geom_histogram(ggplot2::aes(y = ggplot2::after_stat(count / tapply(count, group, sum)[group] * 100)),
                            alpha = 0.5, position = "identity", bins = 100) +
    ggplot2::scale_fill_manual(values = fill_vals, labels = fill_labs) +
    ggplot2::labs(title = paste0(wt_title, " - Histogram"), x = "PS", y = "% patients (within group)") +
    ggplot2::theme_minimal() + ggplot2::theme(plot.title = ggplot2::element_text(size = 14, face = "bold", hjust = 0.5),
                                              legend.position = "top") + ggplot2::xlim(0, 1)
  ggplot2::ggsave(file.path(out_dir_plots, paste0(plot_prefix, "_ps_distr_histog_", panel2_suffix, ".png")),
                  plot = p, width = 8, height = 6, dpi = 150)

  # both
  aes_both <- if (has_w) {
    ggplot2::aes(x = .data[[ps_var]], weight = .data[[weight_var]])
  } else {
    ggplot2::aes(x = .data[[ps_var]])
  }
  p <- ggplot2::ggplot(wt_data, aes_both) +
    ggplot2::geom_histogram(ggplot2::aes(y = ggplot2::after_stat(density), fill = exposure_factor),
                            alpha = 0.5, position = "identity", bins = 100) +
    ggplot2::geom_density(ggplot2::aes(y = ggplot2::after_stat(density), color = exposure_factor), linewidth = 1.2, adjust = 1.5) +
    ggplot2::scale_fill_manual(values = fill_vals, labels = fill_labs) +
    ggplot2::scale_color_manual(values = fill_vals, labels = fill_labs, guide = "none") +
    ggplot2::labs(title = wt_title, x = "PS", y = "Density") +
    ggplot2::theme_minimal() + ggplot2::theme(plot.title = ggplot2::element_text(size = 14, face = "bold", hjust = 0.5),
                                              legend.position = "top") + ggplot2::xlim(0, 1)
  ggplot2::ggsave(file.path(out_dir_plots, paste0(plot_prefix, "_ps_distr_", panel2_suffix, ".png")),
                  plot = p, width = 8, height = 6, dpi = 150)

  # ---- Weight box plot (panel 2 weights only) ----
  if (make_boxplot && has_w) {
    if (verbose) message("  Creating weight distribution box plot...")
    box_data <- weighted_df
    box_data$exposure_group <- factor(box_data[[exposure_var]], levels = c(0, 1),
                                      labels = c("Reference", "Exposure"))
    y_limit <- stats::quantile(weighted_df[[weight_var]], 0.99, na.rm = TRUE) * 1.1

    p <- ggplot2::ggplot(box_data, ggplot2::aes(x = exposure_group, y = .data[[weight_var]], fill = exposure_group)) +
      ggplot2::geom_boxplot(alpha = 0.7, outlier.alpha = 0.3) +
      ggplot2::geom_hline(yintercept = 1, linetype = "dashed", color = "gray50", linewidth = 1) +
      ggplot2::labs(title = box_title, x = "Group", y = "Weight") +
      ggplot2::theme_minimal() + ggplot2::theme(plot.title = ggplot2::element_text(size = 14, face = "bold", hjust = 0.5),
                                                legend.position = "none") +
      ggplot2::scale_fill_manual(values = c("Exposure" = "#377eb8", "Reference" = "#e41a1c")) +
      ggplot2::coord_cartesian(ylim = c(0, y_limit))
    ggplot2::ggsave(file.path(out_dir_plots, paste0(plot_prefix, "_boxplot.png")),
                    plot = p, width = 8, height = 6, dpi = 150)
  }

  invisible(NULL)
}

#### Balance verdict on the raw SMD scale ####
#' Is a covariate balanced at the given threshold?
#'
#' Compares the absolute standardized mean difference against the threshold on
#' the \emph{raw} (0-1) SMD scale produced by \code{\link{build_table1}}. This is
#' the single source of truth for the "Balanced" verdict used by the weighting,
#' matching, and fine-stratification balance tables (it replaces an earlier
#' \code{abs(SMD) < threshold * 100} comparison that silently passed every
#' covariate on raw-scale input).
#'
#' @param std_diff Numeric (scalar or vector) standardized mean difference, raw scale.
#' @param threshold Numeric balance threshold on the raw scale (e.g. 0.1).
#' @return Logical of the same length as \code{std_diff}; \code{NA} propagates.
#' @keywords internal
#' @noRd
.is_balanced <- function(std_diff, threshold) {
  abs(std_diff) < threshold
}

#### PS trimming (Crump symmetric / Sturmer asymmetric) ####
#' Trim a propensity-score distribution (point removal)
#'
#' Returns a logical \code{keep} mask flagging which observations fall inside the
#' retained propensity-score region. Two conventions are supported:
#' \itemize{
#'   \item \code{"crump"} - symmetric Crump et al. (2009): keep
#'     \code{ps in [alpha, 1 - alpha]}.
#'   \item \code{"sturmer"} - asymmetric Sturmer et al. (2010/2021): keep
#'     \code{ps in [Q(ps | exposed, p), Q(ps | reference, 1 - p)]}, where the
#'     lower bound is the \code{p}-th percentile of the PS among the exposed
#'     (\code{exp == 1}) and the upper bound is the \code{(1 - p)}-th percentile
#'     among the reference group (\code{exp == 0}). This assumes \code{exp == 1}
#'     denotes the treated / new-treatment group.
#' }
#' This is \emph{point removal} and is distinct from the common-support range
#' trimming in \code{\link{create_ps_fs_weights}}
#' (\code{trim_nonoverlap_region}). Trimming changes the analytic population:
#' the resulting estimand refers to the trimmed population, not the original
#' cohort.
#'
#' @param ps Numeric propensity-score vector.
#' @param exp Numeric 0/1 exposure vector (1 = exposed / treated).
#' @param method One of \code{"none"}, \code{"crump"}, \code{"sturmer"}.
#' @param crump_alpha Numeric in (0, 0.5). Symmetric Crump bound (default 0.1).
#' @param sturmer_p Numeric in (0, 0.5). Sturmer tail percentile (default 0.05).
#' @param verbose Logical. Print the number removed and the interpretation note.
#' @return Logical vector (length of \code{ps}); \code{TRUE} = keep. \code{NA}
#'   propensity scores are dropped (\code{keep = FALSE}).
#' @keywords internal
#' @noRd
.trim_ps <- function(ps, exp, method = c("none", "crump", "sturmer"),
                     crump_alpha = 0.1, sturmer_p = 0.05, verbose = TRUE) {
  method <- match.arg(method)
  n <- length(ps)
  if (method == "none") {
    return(rep(TRUE, n))
  }

  if (method == "crump") {
    if (!is.numeric(crump_alpha) || length(crump_alpha) != 1 ||
        crump_alpha <= 0 || crump_alpha >= 0.5) {
      stop("trim_crump_alpha must be a single number in (0, 0.5)")
    }
    lower <- crump_alpha
    upper <- 1 - crump_alpha
    desc  <- sprintf("Crump symmetric: keep PS in [%.3f, %.3f]", lower, upper)
  } else {  # sturmer
    if (!is.numeric(sturmer_p) || length(sturmer_p) != 1 ||
        sturmer_p <= 0 || sturmer_p >= 0.5) {
      stop("trim_sturmer_p must be a single number in (0, 0.5)")
    }
    ps_exp <- ps[exp == 1 & !is.na(exp)]
    ps_ref <- ps[exp == 0 & !is.na(exp)]
    if (length(ps_exp) == 0 || length(ps_ref) == 0) {
      stop("Sturmer trimming requires both exposure groups to be non-empty")
    }
    lower <- as.numeric(stats::quantile(ps_exp, sturmer_p, na.rm = TRUE))
    upper <- as.numeric(stats::quantile(ps_ref, 1 - sturmer_p, na.rm = TRUE))
    desc  <- sprintf("Sturmer asymmetric: keep PS in [%.4f, %.4f] (Q exposed %.2f, Q reference %.2f)",
                     lower, upper, sturmer_p, 1 - sturmer_p)
  }

  keep <- !is.na(ps) & ps >= lower & ps <= upper
  n_removed <- sum(!keep)
  if (verbose) {
    message(sprintf("PS trimming (%s)", desc))
    message(sprintf("  Removed %d of %d observations (%.1f%%)",
                    n_removed, n, 100 * n_removed / n))
    message("  Note: trimming changes the analytic population; estimated effects ",
            "refer to the trimmed population, not the original cohort.")
  }
  keep
}

#### Weight truncation / winsorizing (IPTW-family weights) ####
#' Truncate or winsorize a weight vector
#'
#' Bounds extreme weights for inverse-probability-style weights. Matching and
#' overlap weights are bounded by construction and do not use this. Truncation
#' changes the effective analytic population / estimand interpretation and
#' should be reported as a sensitivity analysis.
#'
#' @param weights Numeric weight vector.
#' @param method One of \code{"none"}, \code{"percentile"}, \code{"cap"}.
#' @param percentile Length-2 numeric \code{c(lower, upper)} in `[0, 1]` used when
#'   \code{method = "percentile"} (winsorize to those pooled quantiles). Upper-only
#'   truncation is expressed as \code{c(0, 0.99)}.
#' @param cap Single positive number used when \code{method = "cap"} (absolute
#'   upper cap on the weights).
#' @param verbose Logical. Print the cut points and number affected.
#' @return A list with \code{w} (truncated weights) and \code{cut} (the applied
#'   cut points: \code{c(lower, upper)} for percentile, or \code{c(NA, cap)} for cap).
#' @keywords internal
#' @noRd
.truncate_ps_weights <- function(weights, method = c("none", "percentile", "cap"),
                                 percentile = c(0.01, 0.99), cap = NULL,
                                 verbose = TRUE) {
  method <- match.arg(method)
  if (method == "none") {
    return(list(w = weights, cut = c(NA_real_, NA_real_)))
  }

  if (method == "percentile") {
    if (!is.numeric(percentile) || length(percentile) != 2 ||
        any(percentile < 0) || any(percentile > 1) || percentile[1] >= percentile[2]) {
      stop("truncate_percentile must be c(lower, upper) with 0 <= lower < upper <= 1")
    }
    lo <- as.numeric(stats::quantile(weights, percentile[1], na.rm = TRUE))
    hi <- as.numeric(stats::quantile(weights, percentile[2], na.rm = TRUE))
    n_low  <- sum(weights < lo, na.rm = TRUE)
    n_high <- sum(weights > hi, na.rm = TRUE)
    w <- pmin(pmax(weights, lo), hi)
    if (verbose) {
      message(sprintf("Weight truncation (percentile [%.3f, %.3f]): cut at [%.4f, %.4f]",
                      percentile[1], percentile[2], lo, hi))
      message(sprintf("  Winsorized %d low and %d high weights", n_low, n_high))
    }
    return(list(w = w, cut = c(lo, hi)))
  }

  # method == "cap"
  if (!is.numeric(cap) || length(cap) != 1 || cap <= 0) {
    stop("truncate_cap must be a single positive number when truncate_method = 'cap'")
  }
  n_capped <- sum(weights > cap, na.rm = TRUE)
  w <- pmin(weights, cap)
  if (verbose) {
    message(sprintf("Weight truncation (absolute cap = %.4f): capped %d weights", cap, n_capped))
  }
  list(w = w, cut = c(NA_real_, cap))
}
