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
#' @return Invisibly returns a \code{ggplot} object, or \code{invisible(NULL)}
#'   if \pkg{ggplot2} is unavailable or no data remain after removing NAs.
#' @section Side Effects:
#' Saves the plot to \code{output_path} via \code{ggplot2::ggsave}.
#' @keywords internal
create_love_plot <- function(variable_names, crude_std_diff, adjusted_std_diff,
                             crude_label, adjusted_label, title, output_path,
                             colors = NULL, shapes = NULL, use_absolute = FALSE) {
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
    x_lab <- "Absolute Standardized Difference (%)"
    vline_x <- 10
  } else {
    x_lab <- "Standardized Difference (%)"
    vline_x <- c(-10, 10)
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