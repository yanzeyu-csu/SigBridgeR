#' @title Perform PIPET Screening Analysis
#' @description
#' Predicts cell subpopulations in single-cell data by matching expression profiles
#' to predefined marker gene templates using various distance/similarity metrics.
#' This function implements a template-based classification approach with permutation
#' testing for significance assessment.
#'
#'
#' @param matched_bulk Normalized bulk expression matrix (features × samples).
#'        Column names must match `phenotype` identifiers.
#' @param sc_data Seurat object containing single-cell RNA-seq data.
#' @param phenotype Clinical outcome data. Can be:
#'        - Vector: named with sample IDs
#'        - Data frame: with row names matching bulk columns
#' @param phenotype_class Analysis mode:
#'        - `"binary"`: Case-control design (e.g., responder/non-responder)
#'        - `"continuous"`: Continuous outcome (e.g.,)
#'        - `"survival"`: Time-to-event analysis data.frame
#' @param label_type Character specifying phenotype label type (e.g., "SBS1", "time"), stored in `scRNA_data@misc`
#' @param lg2FC In the DESeq differential expression analysis results, the cutoff value of lg2FC. The default value is `1L`.
#' @param p_adjust In the DESeq differential expression analysis results, the cutoff value of adjust P. The default value is `0.05`.
#' @param show_lg2FC Select whether to show log2 fold changes. The default value is `TRUE`.
#' @param group A character, name of one metadata column to group cells by (for example, orig.ident). The default value is `NULL`
#' @param rm_NA Select Whether to remove NA values. The default value is `TRUE`.
#' @param freq_counts An integer, keep genes expressed in more than a certain number of cells. The default value is `NULL`, which means no filtering.
#' @param normalize Select whether to perform normalization of count data. The default value is `TRUE`.
#' @param scale Select whether to scale and center features in the dataset. The default value is `TRUE`.
#' @param nPerm An integer, number of permutations to do. The default value is `1000L`.
#' @param distance A character, the distance algorithm must be included in "cosine", "pearson", "spearman", "kendall","euclidean","maximum". default value is `NULL`, which means `"cosine"`.
#' @param ... Additional arguments to be passed to `PIPET.optimized`.
#' - seed: Random seed for reproducibility
#' - verbose: Whether to show progress messages
#' - parallel: Whether to use parallel processing
#' - parallel.type: Type of parallel backend
#' - workers: Number of parallel workers
#'
#' @family screen_method
#' @family PIPET
#' @export
DoPIPET <- function(
    matched_bulk,
    sc_data,
    phenotype,
    phenotype_class = c("binary", "continuous"),
    label_type = "PIPET",
    lg2FC = 1,
    p_adjust = 0.05,
    show_lg2FC = TRUE,
    group = NULL,
    rm_NA = TRUE,
    freq_counts = NULL,
    normalize = TRUE,
    scale = TRUE,
    nPerm = 1000L,
    distance = NULL,
    ...
) {
    # * Input validation
    chk::chk_is(matched_bulk, c("matrix", "data.frame"))
    chk::chk_is(sc_data, "Seurat")
    chk::chk_character(label_type)
    phenotype_class <- SigBridgeRUtils::MatchArg(
        phenotype_class,
        c("binary", "continuous"),
        NULL
    )
    purrr::walk(c(show_lg2FC, rm_NA, normalize, scale), chk::chk_flag)
    purrr::walk(
        c(group, distance),
        ~ {
            if (!is.null(.x)) chk::chk_character(.x)
        }
    )
    chk::chk_integer(nPerm)
    purrr::walk(c(lg2FC, p_adjust), chk::chk_double)

    # * Default params
    ## * general params
    dots <- rlang::list2(...)
    verbose <- dots$verbose %||% getFuncOption("verbose")
    seed <- dots$seed %||% getFuncOption("seed")
    parallel <- dots$parallel %||% getFuncOption("parallel")
    workers <- dots$workers %||% getFuncOption("workers")
    ## * PIPET params
    distance <- SigBridgeRUtils::MatchArg(
        distance,
        c(
            "cosine",
            "pearson",
            "spearman",
            "kendall",
            "euclidean",
            "maximum"
        )
    )

    if (verbose) {
        ts_cli$cli_alert_info(cli::col_green("Starting PIPET screen"))
        ts_cli$cli_alert_info("Creating marker genes from bulk data...")
    }
    markers <- PIPET::Create_Markers(
        bulk_data = matched_bulk,
        colData = phenotype_df,
        class_col = "class",
        lg2FC = lg2FC,
        p.adjust = p_adjust,
        show_lg2FC = show_lg2FC,
        verbose = verbose,
        seed = seed,
        parallel = parallel
    )

    # Run PIPET core algorithm
    if (verbose) {
        ts_cli$cli_alert_info("Running PIPET correlation analysis...")
    }

    pipet_result <- PIPET::PIPET.optimized(
        sc_data = sc_data,
        markers = markers,
        group = group,
        rm_NA = rm_NA,
        freq_counts = freq_counts,
        normalize = normalize,
        scale = scale,
        nPerm = nPerm,
        distance = distance,
        verbose = verbose,
        seed = seed,
        parallel = parallel,
        workers = workers
    )

    # Add results to Seurat object if applicable

    sc_data_final <- Seurat::AddMetaData(
        sc_data,
        metadata = pipet_result
    ) %>%
        AddMisc(
            PIPET_type = label_type,
            PIPET_params = list(
                lg2FC = lg2FC,
                p.adjust = p_adjust,
                distance = distance,
                nPerm = nPerm
            ),
            cover = FALSE
        )

    if (verbose) {
        ts_cli$cli_alert_success(cli::col_green("PIPET screening completed"))
    }

    list(
        scRNA_data = sc_data_final,
        pipet_result = pipet_result,
        markers = markers
    )
}
