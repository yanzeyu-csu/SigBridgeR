DoPIPET <- function(
    matched_bulk,
    sc_data,
    phenotype,
    label_type = "PIPET",
    lg2FC = 1,
    p_adjust = 0.05,
    show_lg2FC = TRUE,
    group = NULL,
    rm_NA = TRUE,
    freq_counts = NULL,
    normalize = TRUE,
    scale = TRUE,
    nPerm = 1000,
    distance = "cosine",
    ...
) {
    # Input validation
    chk::chk_is(matched_bulk, c("matrix", "data.frame"))
    chk::chk_is(sc_data, "Seurat")
    chk::chk_character(label_type)

    # * Default params
    dots <- rlang::list2(...)
    verbose <- dots$verbose %||% getFuncOption("verbose")
    seed <- dots$seed %||% getFuncOption("seed")
    parallel <- dots$parallel %||% getFuncOption("parallel")
    workers <- dots$workers %||% getFuncOption("workers")

    if (verbose) {
        ts_cli$cli_alert_info(cli::col_green("Starting PIPET screen"))
    }

    # Create markers from bulk data
    if (verbose) {
        ts_cli$cli_alert_info("Creating marker genes from bulk data...")
    }
    markers <- PIPET::PIPET_CreateMarkers(
        bulk_data = matched_bulk,
        colData = phenotype_df,
        class_col = "class",
        lg2FC = lg2FC,
        p.adjust = p_adjust,
        show_lg2FC = show_lg2FC
    )

    # Run PIPET core algorithm
    ts_cli$cli_alert_info("Running PIPET correlation analysis...")
    pipet_result <- PIPET::PIPET.optimized(
        SC_data = sc_data,
        markers = markers,
        group = group,
        rm_NA = rm_NA,
        freq_counts = freq_counts,
        normalize = normalize,
        scale = scale,
        nPerm = nPerm,
        distance = distance,
        nCores = nCores
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
