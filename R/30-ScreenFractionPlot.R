# ---- Visualization ----

# ---- Screened cell fraction(+/-/N)-sample/source stacked graph ----

#' @title Visualization of Cell Screening Fractions
#'
#' @description
#' Generates stacked bar plots showing the proportion of cells classified as
#' Positive/Negative/Neutral by single-cell screening algorithms (Scissor, scPAS,
#' scPP, or scAB) across different sample groups. Supports both single and multiple
#' screen types visualization.
#'
#' @usage
#' ScreenFractionPlot(
#'   screened_seurat,
#'   group_by = "Source", # x-axis label
#'   screen_type = NULL, # default: search all available screens
#'   show_null = FALSE,
#'   plot_color = NULL, # stack bar color of each status
#'   show_plot = FALSE,
#'   plot_title = "Screen Fraction",
#'   stack_width = 0.85,
#'   x_text_angle = 45,
#'   axis_linewidth = 0.8,
#'   legend_position = "right",
#'   x_lab = NULL,
#'   y_lab = "Fraction of Status",
#'   ncol = 2, # number of columns for facet wrap
#'   nrow = NULL, # number of rows for facet wrap
#'   ... # unused
#' )
#'
#' @param screened_seurat A Seurat object containing screening results in metadata.
#'        Must contain columns corresponding to `screen_type`.
#' @param group_by Metadata column name for grouping samples (default: "Source").
#' @param screen_type Screening algorithm(s) used. Can be a single value or vector.
#'        Must match metadata column(s) (case-sensitive, e.g., "scissor" for Scissor results).
#' @param show_null Logical whether to show groups with zero cells (default: FALSE).
#' @param plot_color Custom color palette (named vector format):
#'        - Required names: "Positive", "Negative", "Neutral"
#'        - Default: c("Neutral"="#CECECE", "Other"="#CECECE", Positive"="#ff3333", "Negative"="#386c9b")
#' @param show_plot Logical to immediately display plot (default: TRUE).
#' @param plot_title Plot title (default: "Screen Fraction"). When multiple screen types,
#'        can be a vector of titles or single title (will append screen type).
#' @param stack_width Bar width (default: 0.85).
#' @param x_text_angle X-axis label angle (default: 45).
#' @param axis_linewidth Axis line thickness (default: 0.8).
#' @param legend_position Legend position (default: "right").
#' @param x_lab X-axis label (default: NULL).
#' @param y_lab Y-axis label (default: "Fraction of Status").
#' @param ncol Number of columns for facet wrap when multiple screen types (default: 2).
#' @param nrow Number of rows for facet wrap when multiple screen types (default: NULL).
#' @param.... Useless
#'
#' @return A list containing:
#' \itemize{
#'   \item{stats: A data frame (single screen) or list of data frames (multiple screens)
#'     with screening statistics including:
#'     \itemize{
#'       \item Grouping variable counts
#'       \item Raw cell counts
#'       \item Percentage fractions
#'     }
#'   }
#'   \item{plot: A ggplot2 object (single screen) or list of ggplot2 objects (multiple screens)}
#'   \item{combined_plot: A combined plot using patchwork (only for multiple screens)}
#' }
#'
#' @section Visualization Details:
#' - Bars are ordered by descending Positive fraction
#' - Y-axis shows percentage (0-100%)
#' - Zero-fraction groups are automatically hidden unless `show_null=TRUE`
#' - For multiple screen types, plots can be combined using patchwork
#'
#' @examples
#' \dontrun{
#' # Single screen type usage
#' res <- ScreenFractionPlot(
#'   screened_seurat = scissor_result,
#'   group_by = "PatientID",
#'   screen_type = "scissor"
#' )
#'
#' # Multiple screen types at once
#' multi_res <- ScreenFractionPlot(
#'   screened_seurat = multi_screened_result,
#'   group_by = "TissueType",
#'   screen_type = c("scissor", "scPAS", "scPP", "scAB", "DEGAS"),
#'   ncol = 2
#' )
#' }
#'
#' @export
#' @importFrom dplyr count group_by mutate ungroup pull filter arrange
#' @importFrom tidyr complete
#' @importFrom ggplot2 ggplot aes geom_col scale_y_continuous scale_fill_manual
#' @importFrom ggplot2 theme_classic labs theme element_text element_line
#' @importFrom glue glue
#' @importFrom patchwork wrap_plots plot_annotation
#'
#' @family visualization_function
#'
ScreenFractionPlot <- function(
    screened_seurat,
    group_by = "Source",
    screen_type = NULL,
    show_null = FALSE,
    plot_color = NULL,
    show_plot = FALSE,
    plot_title = "Screen Fraction",
    stack_width = 0.85,
    x_text_angle = 45L,
    axis_linewidth = 0.8,
    legend_position = "right",
    x_lab = NULL,
    y_lab = "Fraction of Status",
    ncol = 2L,
    nrow = NULL,
    ...
) {
    chk::chk_is(screened_seurat, "Seurat")
    chk::chk_character(group_by)
    chk::chk_length(group_by, 1)
    chk::chk_flag(show_null)
    if (!is.null(plot_color)) {
        chk::chk_vector(plot_color)
    }

    meta_data <- screened_seurat[[]]
    all_screen_types <- colnames(meta_data)
    available_screens <- grep(
        "sc[a-zA-Z_]+$|DEGAS$",
        names(meta_data),
        value = TRUE
    )
    if (!group_by %chin% all_screen_types) {
        cli::cli_abort(c(
            "x" = "Grouping variable not found in metadata.",
            ">" = "Available grouping variables: {.val {all_screen_types}}"
        ))
    }
    # Check available screen types in the Seurat object
    if (is.null(screen_type)) {
        screen_type <- available_screens
    } else if (
        !all(purrr::map_vec(
            screen_type,
            ~ . %chin% all_screen_types # We don't use `available_screens` for the sake of compatibility to other groups
        ))
    ) {
        cli::cli_abort(c(
            "x" = "Screen type(s) not found in metadata.",
            ">" = "Available screen types: {.val {available_screens}}"
        ))
    }

    plot_color <- plot_color %||%
        c(
            "Neutral" = "#CECECE",
            "Other" = "#CECECE",
            "Positive" = "#ff3333",
            "Negative" = "#386c9b"
        )

    # Function to create plot for single screen type
    SinglePlot <- function(single_screen_type, title_suffix = "") {
        stats_df <- meta_data %>%
            dplyr::count(!!sym(group_by), !!sym(single_screen_type)) %>%
            tidyr::complete(
                !!sym(group_by),
                !!sym(single_screen_type),
                fill = list(`n` = 0)
            ) %>%
            dplyr::group_by(!!sym(group_by)) %>%
            dplyr::mutate(Total = sum(`n`)) %>%
            dplyr::ungroup() %>%
            dplyr::mutate(
                Fraction = ifelse(Total == 0, 0, `n` / Total)
            )

        plot_order <- stats_df %>%
            dplyr::filter(!!sym(single_screen_type) == "Positive") %>%
            dplyr::arrange(dplyr::desc(Fraction)) %>%
            dplyr::pull(!!sym(group_by))

        # Get label type from misc
        label_type <- screened_seurat@misc[[grep(
            glue::glue("{single_screen_type}_type"),
            names(screened_seurat@misc),
            value = TRUE
        )[1]]]

        if (is.null(label_type)) {
            label_type <- single_screen_type
        }

        # Filter null records
        plot_df <- if (!show_null) {
            dplyr::filter(stats_df, Fraction > 0)
        } else {
            stats_df
        }

        # Create plot title
        current_title <- if (
            length(plot_title) > 1 && length(plot_title) == length(screen_type)
        ) {
            plot_title[which(screen_type == single_screen_type)]
        } else if (title_suffix != "") {
            glue::glue("{plot_title} - {title_suffix}")
        } else {
            plot_title
        }

        plot <- ggplot2::ggplot(
            plot_df,
            ggplot2::aes(
                x = factor(!!sym(group_by), levels = plot_order),
                y = `Fraction`,
                fill = !!sym(single_screen_type)
            )
        ) +
            ggplot2::geom_col(position = "stack", width = stack_width) +
            ggplot2::scale_y_continuous(
                labels = function(x) paste0(round(x * 100L, 0L), "%"),
                expand = c(0, 0),
                breaks = seq(0, 1, 0.1)
            ) +
            ggplot2::scale_fill_manual(
                values = plot_color
            ) +
            ggplot2::theme_classic(base_size = 14L) +
            ggplot2::labs(
                x = x_lab,
                y = y_lab,
                fill = label_type,
                title = current_title
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

        return(list(stats = stats_df, plot = plot))
    }

    # Process single or multiple screen types
    if (length(screen_type) == 1) {
        # Single screen type - maintain backward compatibility
        result <- SinglePlot(screen_type)

        if (show_plot) {
            print(result$plot)
        }

        return(result)
    } else {
        # Multiple screen types
        cli::cli_inform(
            "Creating plots for {.val {length(screen_type)}} screen types..."
        )

        # Create plots for all screen types in parallel if possible
        plot_results <- lapply(screen_type, function(st) {
            SinglePlot(st, title_suffix = st)
        })

        # Extract stats and plots
        stats_list <- lapply(plot_results, function(x) x$stats)
        names(stats_list) <- screen_type

        plots_list <- lapply(plot_results, function(x) x$plot)
        names(plots_list) <- screen_type

        # Combine plots using patchwork
        combined_plot <- patchwork::wrap_plots(
            plots_list,
            ncol = ncol,
            nrow = nrow
        )

        if (show_plot) {
            print(combined_plot)
        }

        return(list(
            stats = stats_list,
            plot = plots_list,
            combined_plot = combined_plot
        ))
    }
}
