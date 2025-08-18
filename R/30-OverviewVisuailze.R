# ----  Visualization ----

# ---- Screened cell fraction(+/-/N)-sample/source stacked graph ----

#' @title Visualization of Cell Screening Fractions
#'
#' @description
#' Generates stacked bar plots showing the proportion of cells classified as
#' Positive/Negative/Neutral by single-cell screening algorithms (Scissor, scPAS,
#' scPP, or scAB) across different sample groups.
#'
#' @usage
#' ScreenFractionPlot(
#'   screened_seurat,
#'   group_by = "Source",
#'   screen_type = c("scissor", "scPAS", "scPP", "scAB"),
#'   show_null = FALSE,
#'   plot_color = NULL,
#'   show_plot = TRUE,
#'   plot_title = "Screen Fraction",
#'   stack_width = 0.85,
#'   x_text_angle = 45,
#'   axis_linewidth = 0.8,
#'   legend_position = "right"
#' )
#'
#' @param screened_seurat A Seurat object containing screening results in metadata.
#'        Must contain columns corresponding to `screen_type`.
#' @param group_by Metadata column name for grouping samples (default: "Source").
#' @param screen_type Screening algorithm used. Must match a metadata column
#'        (case-sensitive, e.g., "scissor" for Scissor results).
#' @param show_null Logical whether to show groups with zero cells (default: FALSE).
#' @param plot_color Custom color palette (named vector format):
#'        - Required names: "Positive", "Negative", "Neutral"
#'        - Default: c("Neutral"="#CECECE", "Positive"="#ff3333", "Negative"="#386c9b")
#' @param show_plot Logical to immediately display plot (default: TRUE).
#' @param plot_title Plot title (default: "Screen Fraction").
#' @param stack_width Bar width (default: 0.85).
#' @param x_text_angle X-axis label angle (default: 45).
#' @param axis_linewidth Axis line thickness (default: 0.8).
#' @param legend_position Legend position (default: "right").
#' @param x_lab X-axis label (default: NULL).
#' @param y_lab Y-axis label (default: "Status Fraction")
#'
#' @return A list containing:
#' \itemize{
#'   \item{stats: A data frame with screening statistics including:
#'     \itemize{
#'       \item Grouping variable counts
#'       \item Raw cell counts
#'       \item Percentage fractions
#'     }
#'   }
#'   \item{plot: A ggplot2 object of the stacked bar plot}
#' }
#'
#' @section Visualization Details:
#' - Bars are ordered by descending Positive fraction
#' - Y-axis shows percentage (0-100%)
#' - Zero-fraction groups are automatically hidden unless `show_null=TRUE`
#'
#' @examples
#' \dontrun{
#' # Basic usage with Scissor results
#' res <- ScreenFractionPlot(
#'   screened_seurat = scissor_result,
#'   group_by = "PatientID",
#'   screen_type = "scissor"
#' )
#'
#' # Customized plot
#' custom_plot <- ScreenFractionPlot(
#'   screened_seurat = scAB_result,
#'   group_by = "TissueType",
#'   screen_type = "scAB",
#'   plot_color = c("Positive"="red", "Negative"="blue", "Neutral"="grey80"),
#'   plot_title = "scAB Screening Results",
#'   x_text_angle = 90
#' )
#' }
#'
#' @export
#' @importFrom dplyr count group_by mutate ungroup pull filter arrange
#' @importFrom tidyr complete
#' @importFrom ggplot2 ggplot aes geom_col scale_y_continuous scale_fill_manual
#' @importFrom ggplot2 theme_classic labs theme element_text element_line
#' @importFrom scales percent_format
#' @importFrom glue glue
#'
#'
ScreenFractionPlot = function(
    screened_seurat,
    group_by = "Source",
    screen_type = c("scissor", "scPAS", "scPP", "scAB"),
    show_null = FALSE,
    plot_color = NULL,
    show_plot = TRUE,
    plot_title = "Screen Fraction",
    stack_width = 0.85,
    x_text_angle = 45,
    axis_linewidth = 0.8,
    legend_position = "right",
    x_lab = NULL,
    y_lab = "Fraction of Status"
) {
    library(dplyr)

    if (!inherits(screened_seurat, "Seurat")) {
        cli::cli_abort(
            c("x" = "{.var screened_seurat} must be a Seurat object"),
            class = "TypeError"
        )
    }

    if (length(screen_type) != 1) {
        cli::cli_abort(c(
            "x" = "Please refer one screen algorithm type({.var screen_type}).",
            "i" = "Available screen types: ",
            grep(
                "scissor$|scPAS$|scPP$|scAB.*$",
                screened_seurat$scRNA_data@meta.data %>% names(),
                value = TRUE
            )
        ))
    }

    plot_color <- plot_color %||%
        c("Neutral" = "#CECECE", "Positive" = "#ff3333", "Negative" = "#386c9b")

    result <- list() # return

    stats_df <- screened_seurat@meta.data %>%
        dplyr::count(!!sym(group_by), !!sym(screen_type)) %>%
        tidyr::complete(
            !!sym(group_by),
            !!sym(screen_type),
            fill = list(n = 0)
        ) %>%
        dplyr::group_by(!!sym(group_by)) %>%
        dplyr::mutate(Total = sum(n)) %>% # total number of each group
        dplyr::ungroup() %>%
        dplyr::mutate(
            Fraction = ifelse(Total == 0, 0, n / Total)
        ) %>%
        {
            result$stats <<- .
            .
        } %>%
        dplyr::select(!!sym(group_by), !!sym(screen_type), Fraction)

    plot_order <- stats_df %>%
        dplyr::filter(!!sym(screen_type) == "Positive") %>%
        dplyr::arrange(desc(Fraction)) %>%
        dplyr::pull(!!sym(group_by))

    # for plot labels
    label_type = screened_seurat@misc[[grep(
        screen_type,
        names(screened_seurat@misc),
        value = TRUE
    )]]

    # filter null records
    if (!show_null) {
        stats_df <- stats_df %>% dplyr::filter(Fraction > 0)
    }

    result$plot <- ggplot2::ggplot(
        stats_df,
        ggplot2::aes(
            x = factor(!!sym(group_by), levels = plot_order),
            y = `Fraction`,
            fill = !!sym(screen_type),
            title = plot_title
        )
    ) +
        ggplot2::geom_col(position = "stack", width = stack_width) +
        ggplot2::scale_y_continuous(
            labels = scales::percent_format(accuracy = 1),
            expand = c(0, 0),
            breaks = seq(0, 1, 0.1)
        ) +
        ggplot2::scale_fill_manual(
            values = plot_color
        ) +
        ggplot2::theme_classic(base_size = 14) +
        ggplot2::labs(
            x = x_lab,
            y = y_lab,
            fill = label_type
        ) +
        ggplot2::theme(
            axis.text.x = ggplot2::element_text(
                angle = x_text_angle,
                hjust = 1,
                vjust = 1
            ),
            axis.text = ggplot2::element_text(color = "black"),
            legend.position = legend_position,
            axis.line = ggplot2::element_line(linewidth = axis_linewidth)
        )

    if (show_plot) {
        print(result$plot)
    }

    return(result)
}
