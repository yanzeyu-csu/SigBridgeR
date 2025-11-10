#' @title Perform scPP screening analysis
#'
#' @description
#' This function performs scPP screening on single-cell data using matched bulk data and phenotype information.
#' It supports binary, continuous, and survival phenotype types.
#'
#' @usage
#' DoscPP(
#'   matched_bulk,
#'   sc_data,
#'   phenotype,
#'   label_type = "scPP",
#'   phenotype_class = c("Binary", "Continuous", "Survival"),
#'   ref_group = 0,
#'   Log2FC_cutoff = 0.585,
#'   estimate_cutoff = 0.2,
#'   probs = c(0.2, NULL),
#'   ...
#' )
#'
#' @param matched_bulk Bulk expression data (genes × samples) where:
#'        - Column names must match `phenotype` row names
#' @param sc_data Seurat object containing preprocessed single-cell data:
#'        - Normalized counts in `RNA` assay
#' @param phenotype Data frame or tibble or named vector with:
#'        - Rownames matching `matched_bulk` columns
#'        - For survival: must contain time and status columns
#' @param label_type Character specifying phenotype label type (e.g., "SBS1"), stored in `scRNA_data@misc`
#' @param phenotype_class Analysis type (case-sensitive):
#'        - `"Binary"`: Case-control studies (e.g., tumor/normal)
#'        - `"Continuous"`: Quantitative traits (e.g., drug response)
#'        - `"Survival"`: Time-to-event data (requires time/status columns)
#' @param ref_group Reference group or baseline for **binary** comparisons, e.g. "Normal" for Tumor/Normal studies and 0 for 0/1 case-control studies. (default: 0)
#' @param Log2FC_cutoff Minimum log2 fold-change for binary markers (default: 0.585)
#' @param estimate_cutoff Effect size threshold for **continuous** traits (default: 0.2)
#' @param probs A numeric value indicating the quantile cutoff for cell classification. This parameter can also be a numeric vector, in which case an optimal threshold will be selected based on the AUC and enrichment score.(default: 0.2)
#' @param ... Additional arguments. Currently supports:
#'    - `verbose`: Logical indicating whether to print progress messages. Defaults to `TRUE`.
#'    - `seed`: For reproducibility, default is `123L`
#'    - `parallel`: Logical. When `probs` is a numeric vector or NULL, whether to enable parallel mode to search for the optimal `probs` value. Default: `FALSE`.
#'    - `workers`: Number of workers to use in parallel mode. Default: `NULL` (use all 4 cores).
#'
#' @return A list containing:
#' \describe{
#'   \item{scRNA_data}{Seurat object with added metadata:
#'     \describe{
#'       \item{ScPP}{"Positive"/"Negative"/"Neutral" classification}
#'     }
#'   }
#'   \item{gene_list}{List of genes used for screening}
#'   \item{AUC}{A data.frame with area under the ROC curve:
#'     \describe{
#'         \item{scPP_AUCup}{AUC for positive}
#'         \item{scPP_AUCdown}{AUC for negative}
#'     }
#'   }
#' }
#'
#' @section Algorithm Steps:
#' 1. Data Validation: Checks sample alignment between bulk and phenotype data
#' 2. Marker Selection: Identifies phenotype-associated genes from bulk data
#' 3. Single-cell Screening: Projects bulk markers onto single-cell data
#' 4. Cell Classification: Categorizes cells based on phenotype association
#'
#' @section Reference:
#' WangX-Lab/ScPP \[Internet\]. \[cited 2025 Aug 31\]. Available from: https://github.com/WangX-Lab/ScPP
#'
#' @examples
#' \dontrun{
#' # Binary phenotype analysis
#' res <- DoscPP(
#'   matched_bulk = bulk_data,
#'   sc_data = seurat_obj,
#'   phenotype = ms_data,
#'   label_type = "SBS1",
#'   phenotype_class = "Binary"
#' )
#'
#' # Survival analysis
#' surv_res <- DoscPP(
#'   sc_data = seurat_obj,
#'   matched_bulk = bulk_data,
#'   phenotype = surv_df,
#'   label_type = "OS_status",
#'   phenotype_class = "Survival"
#' )
#' }
#'
#'
#' @family screen_method
#' @family scPP
#'
DoscPP <- function(
    matched_bulk,
    sc_data,
    phenotype,
    label_type = "scPP",
    phenotype_class = c("Binary", "Continuous", "Survival"),
    ref_group = 0,
    Log2FC_cutoff = 0.585,
    estimate_cutoff = 0.2,
    probs = c(0.2, NULL),
    ...
) {
    chk::chk_is(sc_data, "Seurat")
    chk::chk_character(label_type)
    phenotype_class <- MatchArg(
        phenotype_class,
        c("Binary", "Continuous", "Survival"),
        NULL
    )
    chk::chk_number(ref_group)
    chk::chk_range(Log2FC_cutoff)
    chk::chk_range(estimate_cutoff)
    if (!is.null(probs)) {
        chk::chk_range(probs)
    }
    # scPP can't tolerate NA
    chk::chk_not_any_na(matched_bulk)
    chk::chk_not_any_na(phenotype)

    # robust, scPP is more strict than scissor and scPAS
    if (phenotype_class == "Survival") {
        if (!all(rownames(phenotype) == colnames(matched_bulk))) {
            cli::cli_abort(c(
                "x" = "Please check the rownames of {.var phenotype} and colnames of {.var bulk_dataset}, they should be the same."
            ))
        }
    } else {
        if (!all(names(phenotype) == colnames(matched_bulk))) {
            cli::cli_abort(c(
                "x" = "Please check the names of {.var phenotype} and colnames of {.var bulk_dataset}, they should be the same."
            ))
        }
    }

    dots <- rlang::list2(...)
    verbose <- dots$verbose %||% getFuncOption("verbose")
    parallel <- dots$parallel %||% getFuncOption("parallel")
    workers <- dots$workers %||% getFuncOption("workers")
    seed <- dots$seed %||% getFuncOption("seed")

    if (verbose) {
        ts_cli$cli_alert_info(cli::col_green("Start scPP screening."))
        ts_cli$cli_alert_info("Finding overall markers...")
    }

    matched_bulk <- as.data.frame(matched_bulk)
    # decide which type of phenotype data is used
    if (is.vector(phenotype)) {
        phenotype <- as.data.frame(phenotype) %>%
            Rownames2Col("Sample") %>%
            dplyr::rename("Feature" = 2) %>%
            dplyr::mutate(Feature = as.numeric(`Feature`))
    }
    gene_list <- if (tolower(phenotype_class) == "binary") {
        Check0VarRows(matched_bulk)
        ScPP::marker_Binary(
            bulk_data = matched_bulk,
            features = phenotype,
            ref_group = ref_group,
            Log2FC_cutoff = Log2FC_cutoff
        )
    } else if (tolower(phenotype_class) == "continuous") {
        ScPP::marker_Continuous(
            bulk_data = matched_bulk,
            features = phenotype$Feature,
            estimate_cutoff = estimate_cutoff
        )
    } else {
        ScPP::marker_Survival(
            bulk_data = matched_bulk,
            survival_data = phenotype
        )
    }

    l <- lapply(gene_list, length)
    pos_null <- FALSE
    neg_null <- FALSE
    if ("gene_pos" %chin% names(l)) {
        # Cannot combine the conditions due to the feature of `gene_list`
        if (l[["gene_pos"]] == 0 && verbose) {
            ts_cli$cli_alert_info("No significant positive genes found")
            pos_null <- TRUE
        }
    }
    if ("gene_neg" %chin% names(l)) {
        if (l[["gene_neg"]] == 0 && verbose) {
            ts_cli$cli_alert_info("No significant negative genes found")
            neg_null <- TRUE
        }
    }
    if (pos_null && neg_null) {
        cli::cli_warn(c(
            "scPP is not applicable to the current data. Returning {.val NULL}",
            "scPP screening exits 1."
        ))
        return(NULL)
    }

    if (verbose) {
        ts_cli$cli_alert_info("Screening...")
    }

    # *Start screen
    scPP_result <- ScPP.optimized(
        sc_dataset = sc_data,
        geneList = gene_list,
        probs = probs,
        verbose = verbose,
        parallel = parallel,
        workers = workers,
        seed = seed
    )
    sc_data[[]] <- scPP_result$metadata
    sc_data <- AddMisc(sc_data, scPP_type = label_type, cover = FALSE)

    if (verbose) {
        ts_cli$cli_alert_success(cli::col_green("scPP screening done."))
    }

    list(
        scRNA_data = sc_data,
        gene_list = list(
            genes_pos = scPP_result$Genes_pos,
            genes_neg = scPP_result$Genes_neg
        ),
        AUC = dplyr::select(
            scPP_result$metadata,
            "scPP_AUCup",
            "scPP_AUCdown"
        )
    )
}


#' @title scPP Screening with Optimal Threshold Detection
#'
#' @description
#' Performs single-cell phenotype screening using gene set enrichment analysis.
#' Can either use a fixed probability threshold or automatically find the optimal
#' threshold by testing multiple values and maximizing NES difference.
#'
#' @param sc_dataset A Seurat object containing single-cell RNA-seq data.
#'   Must have RNA assay with normalized data.
#' @param geneList A named list containing gene sets:
#'   \itemize{
#'     \item `gene_pos` - Genes associated with positive phenotype (required)
#'     \item `gene_neg` - Genes associated with negative phenotype (required)
#'     \item `genes_sort` - Named numeric vector of ranked genes (required for optimization mode)
#'   }
#' @param probs Numeric value or vector of probability thresholds:
#'   \itemize{
#'     \item Single value (e.g., 0.2): Performs phenotype profiling with fixed threshold
#'     \item Multiple values (e.g., seq(0.2, 0.45, by = 0.05)): Finds optimal threshold
#'     \item NULL (default): Automatically searches optimal threshold using seq(0.2, 0.45, by = 0.05)
#'   }
#' @param verbose Logical, whether to print progress messages.
#'
#' @return
#' A list with three components:
#' \itemize{
#'   \item `metadata` - Data frame with cell metadata including scPP_AUCup, scPP_AUCdown, scPP
#'   \item `Genes_pos` - Genes upregulated in Positive vs Negative
#'   \item `Genes_neg` - Genes upregulated in Negative vs Positive
#' }
#'
#' @details
#' This function operates in two modes based on the `probs` parameter:
#'
#' ## Fixed Threshold Mode (length(probs) == 1):
#' 1. Computes AUCell scores for positive and negative gene sets
#' 2. Classifies cells based on the specified threshold
#' 3. Identifies differential markers between phenotype groups
#' 4. Returns complete results with metadata and marker genes
#'
#' ## Optimization Mode (length(probs) > 1 or NULL):
#' 1. Tests multiple probability thresholds
#' 2. For each threshold, classifies cells and finds markers
#' 3. Runs GSEA to calculate NES for marker sets
#' 4. Returns threshold with maximum NES difference (Positive - Negative)
#' 5. Requires `genes_sort` in geneList for GSEA analysis
#'
#' @note
#' - Fixed threshold mode: Faster, returns detailed results
#' - Optimization mode: Slower, requires genes_sort, but robust
#'
#' @examples
#' \dontrun{
#' # Fixed threshold mode
#' result <- ScPP.optimized(
#'   sc_dataset = seurat_obj,
#'   geneList = list(
#'     gene_pos = c("CD4", "IL7R"),
#'     gene_neg = c("CD8A", "CD8B")
#'   ),
#'   probs = 0.2
#' )
#'
#' # Optimization mode
#' result <- ScPP.optimized(
#'   sc_dataset = seurat_obj,
#'   geneList = list(
#'     gene_pos = c("CD4", "IL7R"),
#'     gene_neg = c("CD8A", "CD8B"),
#'     genes_sort = ranked_genes
#'   ),
#'   probs = NULL  # or seq(0.2, 0.45, by = 0.05)
#' )
#' }
#'
#' @seealso
#' [AUCell::AUCell_calcAUC()], [Seurat::FindMarkers()], [fgsea::fgsea()]
#'
#' @keywords internal
#' @family scPP
#'
ScPP.optimized <- function(
    sc_dataset,
    geneList,
    probs = c(0.2, NULL),
    verbose = TRUE,
    parallel = FALSE,
    workers = NULL,
    seed = 123L
) {
    # Set default probs if NULL, serach for optimal probs if vector
    probs <- probs %||% seq(0.2, 0.45, by = 0.05)

    # Validate probs
    if (any(probs <= 0 | probs >= 0.5)) {
        cli::cli_abort("{.arg probs} must be numeric values between 0 and 0.5")
    }

    # Determine mode based on probs length
    n_probs <- length(probs)
    is_optimization_mode <- n_probs > 1

    if (is_optimization_mode && "genes_sort" %chin% names(geneList)) {
        # ============================================================
        # OPTIMIZATION MODE: Find optimal probs
        # ============================================================
        probs <- OptimizationMode(
            sc_dataset = sc_dataset,
            geneList = geneList,
            probs = probs,
            parallel = parallel,
            workers = workers,
            seed = seed
        )
    }
    # ============================================================
    # FIXED THRESHOLD MODE: Single prob profiling
    # ============================================================
    FixedProbMode(
        sc_dataset = sc_dataset,
        geneList = geneList,
        probs = probs,
        verbose = verbose,
        seed = seed
    )
}


#' @keywords internal
#' @family scPP
#' @seealso [ScPP.optimized()]
FixedProbMode <- function(
    sc_dataset,
    geneList,
    probs = 0.2,
    verbose = TRUE,
    seed = 123L,
    ...
) {
    chk::chk_length(probs)
    set.seed(seed)

    # Extract gene sets
    geneList_AUC <- geneList[names(geneList) %chin% c("gene_pos", "gene_neg")]

    if (length(geneList_AUC) != 2) {
        cli::cli_abort("geneList must contain 'gene_pos' and 'gene_neg'")
    }
    if (verbose) {
        ts_cli$cli_alert_info(
            "Running fixed threshold mode with prob = {.val {probs}}"
        )
    }

    # Get RNA data
    rna_data <- SeuratObject::LayerData(sc_dataset)
    if (verbose) {
        ts_cli$cli_alert_info("Computing AUC scores...")
    }

    # Compute AUCell scores
    cellrankings <- AUCell::AUCell_buildRankings(rna_data, plotStats = FALSE)
    cellAUC <- AUCell::AUCell_calcAUC(geneList_AUC, cellrankings)

    auc_matrix <- AUCell::getAUC(cellAUC)
    auc_up <- as.numeric(auc_matrix["gene_pos", ])
    auc_down <- as.numeric(auc_matrix["gene_neg", ])

    # Create metadata table
    metadata_dt <- data.table::as.data.table(
        sc_dataset[[]],
        keep.rownames = "cell_id"
    )
    metadata_dt[, `:=`(
        scPP_AUCup = auc_up,
        scPP_AUCdown = auc_down
    )]

    # Calculate quantiles
    up_quantiles <- colQuantiles(
        matrix(c(auc_up, auc_down), ncol = 2),
        probs = c(probs, 1 - probs)
    )

    up_q1 <- up_quantiles[1, 1]
    up_q2 <- up_quantiles[1, 2]
    down_q1 <- up_quantiles[2, 1]
    down_q2 <- up_quantiles[2, 2]

    # Identify phenotype cells
    downcells1 <- metadata_dt[scPP_AUCup <= up_q1, cell_id]
    upcells1 <- metadata_dt[scPP_AUCup >= up_q2, cell_id]
    downcells2 <- metadata_dt[scPP_AUCdown >= down_q2, cell_id]
    upcells2 <- metadata_dt[scPP_AUCdown <= down_q1, cell_id]

    scPP_neg <- purrr::reduce(list(downcells1, downcells2), intersect)
    scPP_pos <- purrr::reduce(list(upcells1, upcells2), intersect)

    # Classify cells
    metadata_dt[, scPP := "Neutral"]
    metadata_dt[cell_id %chin% scPP_pos, scPP := "Positive"]
    metadata_dt[cell_id %chin% scPP_neg, scPP := "Negative"]
    if (verbose) {
        ts_cli$cli_alert_success(
            "Classified {.val {length(scPP_pos)}} Positive, {.val {length(scPP_neg)}} Negative cells"
        )
    }

    # Update Seurat object
    sc_dataset$scPP <- metadata_dt$scPP
    Seurat::Idents(sc_dataset) <- "scPP"
    if (verbose) {
        ts_cli$cli_alert_info(
            "Finding markers between `Positive` and `Negative` group..."
        )
    }

    # Find markers
    markers <- Seurat::FindMarkers(
        sc_dataset,
        ident.1 = "Positive",
        ident.2 = "Negative",
        verbose = FALSE
    )

    # Filter markers
    markers_mat <- as.matrix(markers[, c("avg_log2FC", "p_val_adj")])

    pos_mask <- markers_mat[, "avg_log2FC"] > 1 &
        markers_mat[, "p_val_adj"] < 0.05
    neg_mask <- markers_mat[, "avg_log2FC"] < -1 &
        markers_mat[, "p_val_adj"] < 0.05

    genes_pos <- rownames(markers)[pos_mask]
    genes_neg <- rownames(markers)[neg_mask]

    # Warnings for empty marker sets
    CheckGenes <- purrr::safely(function(genes, msg) {
        if (length(genes) == 0) cli::cli_warn(msg)
    })

    CheckGenes(
        genes_pos,
        "There are no genes significantly upregulated in `Positive` compared to `Negative`."
    )
    CheckGenes(
        genes_neg,
        "There are no genes significantly upregulated in `Negative` compared to `Positive`."
    )
    if (verbose) {
        ts_cli$cli_alert_success(
            "Found {.val {length(genes_pos)}} positive markers, {.val {length(genes_neg)}} negative markers"
        )
    }

    # Return results
    list(
        metadata = as.data.frame(metadata_dt) %>%
            Col2Rownames("cell_id"),
        Genes_pos = genes_pos,
        Genes_neg = genes_neg
    )
}


#' @keywords internal
#' @family scPP
#' @family scPP_optimal_param
#' @seealso [ScPP.optimized()]
OptimizationMode <- function(
    sc_dataset,
    geneList,
    probs,
    verbose = TRUE,
    parallel = FALSE,
    workers = NULL,
    seed = 123L
) {
    set.seed(seed)

    # Extract gene sets
    geneList_AUC <- geneList[names(geneList) %in% c("gene_pos", "gene_neg")]
    genes_sort <- geneList$genes_sort

    if (length(geneList_AUC) != 2) {
        cli::cli_abort("geneList must contain 'gene_pos' and 'gene_neg'")
    }

    if (is.null(genes_sort)) {
        cli::cli_abort(
            "Optimization mode requires 'genes_sort' in {.arg geneList} for GSEA analysis"
        )
    }

    n_probs <- length(probs)

    if (verbose) {
        ts_cli$cli_alert_info(
            "Running optimization mode: testing {.val {n_probs}} threshold{?s}"
        )
    }

    # Get RNA data
    rna_data <- SeuratObject::LayerData(sc_dataset)
    if (verbose) {
        ts_cli$cli_alert_info("Computing AUC scores...")
    }

    # Compute AUCell scores once
    cellrankings <- AUCell::AUCell_buildRankings(rna_data, plotStats = FALSE)
    cellAUC <- AUCell::AUCell_calcAUC(geneList_AUC, cellrankings)

    auc_matrix <- AUCell::getAUC(cellAUC)
    auc_up <- as.numeric(auc_matrix["gene_pos", ])
    auc_down <- as.numeric(auc_matrix["gene_neg", ])

    # Create metadata table
    metadata_dt <- data.table::as.data.table(
        sc_dataset[[]],
        keep.rownames = "cell_id"
    )
    metadata_dt[, `:=`(
        AUCup = auc_up,
        AUCdown = auc_down
    )]

    # Pre-compute all quantiles at once
    all_probs <- c(probs, 1 - probs)
    quantiles_up <- stats::quantile(auc_up, probs = all_probs)
    quantiles_down <- stats::quantile(auc_down, probs = all_probs)

    # Pre-allocate results matrix
    NES_dif_res <- matrix(NA_real_, nrow = n_probs, ncol = 2)
    colnames(NES_dif_res) <- c("prob", "NES_dif")
    NES_dif_res[, "prob"] <- probs

    if (verbose) {
        ts_cli$cli_alert_info(
            "Testing thresholds and computing NES differences..."
        )
    }

    # Iterate through probability thresholds
    if (parallel) {
        workers <- workers %||% 4L
        plan("multisession", workers = workers)

        if (verbose) {
            ts_cli$cli_alert_info(sprintf(
                "Using parallel processing with %d workers",
                workers
            ))
        }
        ProcessSingleProb <- SigBridgeR:::ProcessSingleProb

        ProcessAllProb <- function(i) {
            ProcessSingleProb(
                i = i,
                metadata_dt = metadata_dt,
                sc_dataset = sc_dataset,
                probs = probs,
                quantiles_up = quantiles_up,
                quantiles_down = quantiles_down,
                genes_sort = genes_sort
            )
        }

        results <- future_map(
            .x = seq_len(n_probs),
            .f = ProcessAllProb,
            .progress = verbose,
            .options = furrr_options(
                seed = seed,
                packages = c("Seurat", "data.table", "fgsea", "rlang"),
                globals = c(
                    "metadata_dt",
                    "sc_dataset",
                    "probs",
                    "quantiles_up",
                    "quantiles_down",
                    "genes_sort",
                    "ProcessSingleProb"
                )
            )
        )

        plan("sequential")
    } else {
        results <- purrr::map(
            seq_len(n_probs),
            ~ ProcessSingleProb(
                i = .x,
                metadata_dt = metadata_dt,
                sc_dataset = sc_dataset,
                probs = probs,
                quantiles_up = quantiles_up,
                quantiles_down = quantiles_down,
                genes_sort = genes_sort
            ),
            .progress = verbose
        )
    }

    for (result in results) {
        if (!is.na(result$NES_dif)) {
            NES_dif_res[result$index, "NES_dif"] <- result$NES_dif
        }
    }

    # Convert to data.frame and find optimal
    NES_dif_res <- as.data.frame(NES_dif_res)

    # Check if we have any valid results
    valid_results <- !is.na(NES_dif_res$NES_dif)
    if (!any(valid_results)) {
        cli::cli_abort(c(
            "x" = "No valid NES differences calculated. Try different probability thresholds.",
            ">" = "Current probs: {.val {probs}}"
        ))
    }

    # Find optimal probability
    opt_idx <- which.max(NES_dif_res$NES_dif)
    opt_prob <- NES_dif_res$prob[opt_idx]
    opt_nes <- NES_dif_res$NES_dif[opt_idx]

    if (verbose) {
        ts_cli$cli_alert_success(
            "Optimal threshold: {.val {opt_prob}} (NES difference: {.val {round(opt_nes, 3)}})"
        )
        # Show summary of valid results
        valid_summary <- NES_dif_res[valid_results, ]
        ts_cli$cli_alert_info(
            "Valid results: {.val {sum(valid_results)}}/{.val {n_probs}} thresholds"
        )
    }

    opt_prob
}

#' @keywords internal
#' @family scPP_optimal_param
ProcessSingleProb <- function(
    i,
    metadata_dt,
    sc_dataset,
    probs,
    quantiles_up,
    quantiles_down,
    genes_sort
) {
    prob_i <- probs[i]
    n_probs <- length(probs)

    # Get quantile thresholds
    up_low <- quantiles_up[i]
    up_high <- quantiles_up[n_probs + i]
    down_low <- quantiles_down[i]
    down_high <- quantiles_down[n_probs + i]

    # Identify cells with high AUCup and low AUCdown (Positive group)
    pos_mask <- metadata_dt$AUCup >= up_high &
        metadata_dt$AUCdown <= down_low
    # Identify cells with low AUCup and high AUCdown (Negative group)
    neg_mask <- metadata_dt$AUCup <= up_low &
        metadata_dt$AUCdown >= down_high

    scPP_pos <- metadata_dt$cell_id[pos_mask]
    scPP_neg <- metadata_dt$cell_id[neg_mask]

    # Skip if too few cells in either group
    if (length(scPP_pos) < 3 || length(scPP_neg) < 3) {
        cli::cli_warn(
            "Skipping prob {.val {prob_i}} because too few cells in either group"
        )
        return(list(index = i, NES_dif = NA_real_))
    }

    # * Copy one to avoid modifying
    # Classify cells
    metadata_dt_copy <- data.table::copy(metadata_dt)
    metadata_dt_copy[, scPP := "Neutral"]
    metadata_dt_copy[pos_mask, scPP := "Positive"]
    metadata_dt_copy[neg_mask, scPP := "Negative"]

    # Update Seurat object
    sc_dataset_copy <- sc_dataset
    sc_dataset_copy$scPP <- metadata_dt_copy$scPP
    Seurat::Idents(sc_dataset_copy) <- "scPP"

    # Find markers
    markers <- Seurat::FindMarkers(
        object = sc_dataset_copy,
        ident.1 = "Positive",
        ident.2 = "Negative",
        verbose = FALSE
    )

    if (is.null(markers) || nrow(markers) == 0) {
        cli::cli_warn(
            "Skipping prob {.val {prob_i}} because no markers found"
        )
        return(list(index = i, NES_dif = NA_real_))
    }

    # Filter markers
    markers_mat <- as.matrix(markers[, c("avg_log2FC", "p_val_adj")])

    pos_mask_genes <- markers_mat[, "avg_log2FC"] > 1 &
        markers_mat[, "p_val_adj"] < 0.05
    neg_mask_genes <- markers_mat[, "avg_log2FC"] < -1 &
        markers_mat[, "p_val_adj"] < 0.05

    genes_pos <- rownames(markers)[pos_mask_genes]
    genes_neg <- rownames(markers)[neg_mask_genes]

    # Check if we have markers in both directions
    if (length(genes_pos) == 0 || length(genes_neg) == 0) {
        cli::cli_warn(
            "Skipping prob {.val {prob_i}} because no markers found in both directions"
        )
        return(list(index = i, NES_dif = NA_real_))
    }

    # Run GSEA
    res <- list(
        Genes_pos = genes_pos,
        Genes_neg = genes_neg
    )

    fgsea_res <- rlang::try_fetch(
        fgsea::fgsea(
            pathways = res,
            stats = genes_sort
        ),
        error = function(e) NULL
    )

    if (is.null(fgsea_res)) {
        cli::cli_warn(
            "Skipping prob {.val {prob_i}} because GSEA failed"
        )
        return(list(index = i, NES_dif = NA_real_))
    }

    # Calculate NES difference
    nes_pos <- fgsea_res$NES[fgsea_res$pathway == "Genes_pos"]
    nes_neg <- fgsea_res$NES[fgsea_res$pathway == "Genes_neg"]

    if (length(nes_pos) > 0 && length(nes_neg) > 0) {
        return(list(index = i, NES_dif = nes_pos - nes_neg))
    }

    cli::cli_warn(
        "Skipping prob {.val {prob_i}} because no relevant NES found"
    )

    list(index = i, NES_dif = NA_real_)
}
