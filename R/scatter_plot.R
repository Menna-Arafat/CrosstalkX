#' Create a Scatter Plot Between Two Crosstalk Metrics With Correlation Annotation
#'
#' Generates a scatter plot for two selected columns from a crosstalk results data frame,
#' optionally applies log transformations, overlays a smoothing curve, and annotates the
#' plot with a correlation coefficient and p-value. The plot is also saved as a PNG file.
#'
#' @param result A data.frame containing crosstalk results (e.g., output of
#'   \code{\link{crosstalk}} or \code{\link{compute_crosstalk}}). Must contain the columns
#'   specified by `x` and `y`.
#'
#' @param x Character. Column name to use for the x-axis. Default is `"Jaccard"`.
#'
#' @param y Character. Column name to use for the y-axis. Default is `"pval_deg_norm"`.
#'
#' @param log_transform_x Logical. If `TRUE`, transforms x-values as \code{log10(x + 0.001)}.
#'   Default is `TRUE`.
#'
#' @param log_transform_y Logical. If `TRUE`, transforms y-values as \code{-log10(y + 0.001)}.
#'   Default is `TRUE`.
#'
#' @param model.fit Character. Correlation method passed to \code{stats::cor.test()}.
#'   Typical options include `"spearman"`, `"pearson"`, and `"kendall"`. Default is `"spearman"`.
#'
#' @param model.visual Character. Smoother method passed to \code{ggplot2::geom_smooth()}.
#'   Default is `"loess"`.
#' @param prefix name of the output
#'
#' @details
#' The function creates two temporary columns `x_val` and `y_val` after optional
#' transformation. Correlation is computed using \code{stats::cor.test()} on these
#' transformed values. A PNG file is written to the working directory named
#' \code{"scatter_plot_<x>_vs_<y>.png"}.
#'
#' Note: The annotation label currently reports "Spearman r" regardless of `model.fit`.
#' If you use `"pearson"` or `"kendall"`, consider updating the label accordingly.
#'
#' @return A \pkg{ggplot2} object (the scatter plot). The plot is also saved to disk
#'   as a PNG file.
#'
#' @importFrom ggplot2 ggplot aes geom_point geom_smooth theme_minimal labs annotate ggsave
#' @importFrom dplyr mutate
#' @importFrom stats cor.test
#' @importFrom rlang .data
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # result is a data.frame with columns "Jaccard" and "pval_deg_norm"
#' p <- scatter_plot(result, x = "Jaccard", y = "pval_deg_norm",
#'                   model.fit = "spearman", model.visual = "loess")
#' p
#' }


scatter_plot <- function(result, x = "Jaccard", y = "pval_deg_norm",
                         log_transform_x = TRUE, log_transform_y = TRUE,
                         model.fit = "spearman", model.visual = "loess", prefix = "res") {

        result <- result |>
          dplyr::mutate(
            x_val = if (log_transform_x) log10(.data[[x]] + .001) else .data[[x]],
            y_val = if (log_transform_y) -log10(.data[[y]] + .001) else .data[[y]]
          )

        y_lab <- if (log_transform_y) paste0("log10(", y, ")") else y
        x_lab <- if (log_transform_x) paste0("-log10(", x, ")") else x

        # Compute correlation
        cor_test <- stats::cor.test(result$x_val, result$y_val, method = model.fit)
        cor_coeff <- round(cor_test$estimate, 2)
        p_val <- signif(cor_test$p.value, 2)

        # Create plot
        plot <- ggplot2::ggplot(result, ggplot2::aes(x = .data$x_val, y = .data$y_val)) +
          ggplot2::geom_point(alpha = 0.7, color = "blue4") +
          ggplot2::geom_smooth(method = model.visual, se = TRUE, color = "red2", linetype = "dashed") +
          ggplot2::theme_minimal() +
          ggplot2::labs(
            x = x,
            y = paste0("-log10(", y, ")"),
            title = paste0("Correlation between ", x, " & -log10(", y, ")")
          ) +
          ggplot2::annotate(
            "text",
            x = max(result$x_val, na.rm = TRUE),
            y = max(result$y_val, na.rm = TRUE),
            label = paste0(" r = ", cor_coeff, "\npvalue = ", p_val),
            hjust = 1.1,
            vjust = 1.1,
            size = 5
          )

        ggplot2::ggsave(
          filename = paste0(prefix, "_scatter_plot_", x, "_vs_", y, ".png"),
          plot = plot, width = 10, height = 7, dpi = 600
        )

        plot
}
