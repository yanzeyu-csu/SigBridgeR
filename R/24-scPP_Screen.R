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
#'   probs = 0.2
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
#' @param probs Quantile cutoff for cell classification (default: 0.2)
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
#'         \item{AUCdown}{AUC for negative}
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
#' @importFrom ScPP marker_Binary marker_Continuous marker_Survival ScPP
#' @importFrom Seurat AddMetaData
#' @importFrom dplyr %>% rename
#' @importFrom glue glue
#' @importFrom cli cli_abort cli_alert_info cli_alert_success col_green col_yellow
#' @importFrom tibble rownames_to_column
#' @importFrom stats quantile
#'
#' @family screen_method
#' @family scPP
#' @keywords internal
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
    probs = 0.2
) {
    chk::chk_is(matched_bulk, c("matrix", "data.frame"))
    chk::chk_is(sc_data, "Seurat")
    chk::chk_character(label_type)
    phenotype_class <- match.arg(phenotype_class)
    chk::chk_number(ref_group)
    chk::chk_range(Log2FC_cutoff)
    chk::chk_range(estimate_cutoff)
    chk::chk_number(probs)
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

    ts_cli$cli_alert_info(cli::col_green("Start scPP screening."))

    ts_cli$cli_alert_info("Finding markers...")

    matched_bulk <- as.data.frame(matched_bulk)
    # decide which type of phenotype data is used
    if (is.vector(phenotype)) {
        phenotype <- as.data.frame(phenotype) %>%
            tibble::rownames_to_column("Sample") %>%
            dplyr::rename("Feature" = 2) %>%
            dplyr::mutate(Feature = as.numeric(`Feature`))
    }
    if (tolower(phenotype_class) == "binary") {
        Check0VarRows(matched_bulk)
        gene_list <- ScPP::marker_Binary(
            bulk_data = matched_bulk,
            features = phenotype,
            ref_group = ref_group,
            Log2FC_cutoff = Log2FC_cutoff
        )
    } else if (tolower(phenotype_class) == "continuous") {
        gene_list <- ScPP::marker_Continuous(
            bulk_data = matched_bulk,
            features = phenotype$Feature,
            estimate_cutoff = estimate_cutoff
        )
    } else if (tolower(phenotype_class) == "survival") {
        gene_list <- ScPP::marker_Survival(
            bulk_data = matched_bulk,
            survival_data = phenotype
        )
    }

    l <- lapply(gene_list, length)
    pos_null <- FALSE
    neg_null <- FALSE
    if ("gene_pos" %chin% names(l)) {
        # Cannot combine the conditions due to the feature of `gene_list`
        if (l[["gene_pos"]] == 0) {
            ts_cli$cli_alert_info(" No significant positive genes found")
            pos_null <- TRUE
        }
    }
    if ("gene_neg" %chin% names(l)) {
        if (l[["gene_neg"]] == 0) {
            ts_cli$cli_alert_info(" No significant negative genes found")
            neg_null <- TRUE
        }
    }
    if (pos_null && neg_null) {
        ts_cli$cli_alert_info(
            cli::col_yellow("scPP screening exits 1.")
        )
        return(list(
            scRNA_data = "`scPP` is not applicable to the current data."
        ))
    }

    ts_cli$cli_alert_info("Screening...")

    # *Start screen
    tryCatch(
        scPP_result <- ScPP.optimized(sc_data, gene_list, probs = probs),
        error = function(e) {
            cli::cli_alert_danger(e$message)

            cli::cli_abort("scPP screening exits 1.")
        }
    )
    sc_data[[]] <- scPP_result$metadata
    sc_data <- AddMisc(sc_data, scPP_type = label_type, cover = FALSE)

    ts_cli$cli_alert_success(cli::col_green("scPP screening done."))

    return(
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
    )
}


#' @title Check for Rows with Zero Variance in a Matrix
#'
#' @description
#' This function checks each row of a matrix (including sparse matrices of class `dgCMatrix`)
#' for zero variance. Rows with zero variance or only `NA` values are identified, and an error
#' is thrown listing the names of these rows. This is useful for preprocessing data where
#' constant rows may cause issues in analyses (e.g., PCA, regression).
#'
#' @param mat A numeric matrix or a sparse matrix of class `dgCMatrix` (from the `Matrix` package).
#'   Rows represent features (e.g., genes), and columns represent observations.
#' @param call The environment from which the function was called, used for error reporting.
#'   Defaults to rlang::caller_env(). Most users can ignore this parameter.
#'
#' @return Invisibly returns a numeric vector of row variances. If zero-variance rows are found,
#'   the function throws an error with a message listing the problematic row names.
#'
#' @details
#' For dense matrices, variance is computed using an optimized `RowVars` function that
#' efficiently calculates row variances with proper NA handling. For sparse matrices of
#' class `dgCMatrix`, variance is computed using a mathematical identity that avoids
#' creating large intermediate matrices. Rows with fewer than 2 non-zero observations
#' are treated as zero-variance.
#'
#' This implementation is memory-efficient and handles large matrices better than
#' the original version.
#'
#' @examples
#' \dontrun{
#' # Dense matrix example
#' set.seed(123)
#' mat_dense <- matrix(rnorm(100), nrow = 10)
#' rownames(mat_dense) <- paste0("Gene", 1:10)
#' Check0VarRows(mat_dense) # No error if all rows have variance
#'
#' # Introduce zero variance
#' mat_dense[1, ] <- rep(5, 10) # First row is constant
#' Check0VarRows(mat_dense)     # Throws error listing "Gene1"
#'
#' # Sparse matrix example
#' library(Matrix)
#' mat_sparse <- as(matrix(rpois(100, 0.5), nrow = 10), "dgCMatrix")
#' rownames(mat_sparse) <- paste0("Gene", 1:10)
#' Check0VarRows(mat_sparse)
#' }
#'
#' @importFrom Matrix rowSums
#' @importFrom rlang caller_env
#' @importFrom cli cli_abort
#'
#' @keywords internal
#' @family scPP
#'
Check0VarRows <- function(mat, call = rlang::caller_env()) {
    if (inherits(mat, "dgCMatrix")) {
        n <- ncol(mat)
        row_sums <- rowSums(mat)
        row_sums_sq <- rowSums(mat^2)

        # Var = (Σx² - (Σx)²/n) / (n-1)
        vars <- (row_sums_sq - (row_sums^2) / n) / (n - 1)

        zero_counts <- rowSums(mat != 0)
        vars[zero_counts <= 1] <- 0
        vars[is.nan(vars)] <- 0
    } else {
        vars <- RowVars(mat, na.rm = TRUE)
    }

    zero_idx <- which(vars == 0 | is.na(vars))
    if (length(zero_idx)) {
        bad_genes <- rownames(mat)[zero_idx]
        cli::cli_abort(
            c(
                "Detected {.val {length(bad_genes)}} gene(s) with zero variance:",
                "i" = "{.val {bad_genes}}"
            ),
            class = "ZeroVarianceError",
            call = call
        )
    }
    invisible(vars)
}


#' @title Single-Cell Phenotype Profiling (Optimized)
#'
#' @description
#' Performs optimized single-cell phenotype profiling using gene set enrichment
#' analysis and differential expression to identify phenotype-associated cells
#' and marker genes. This implementation uses AUCell for efficient gene set
#' scoring and provides comprehensive phenotype characterization.
#'
#' @param sc_dataset A Seurat object containing single-cell RNA-seq data.
#'   Must have RNA assay with normalized data.
#' @param geneList A named list of length 2 containing gene sets for phenotype
#'   characterization. The list should contain:
#'   \itemize{
#'     \item `gene_pos` - Genes associated with positive phenotype
#'     \item `gene_neg` - Genes associated with negative phenotype
#'   }
#' @param probs Numeric value between 0 and 1 specifying the quantile threshold
#'   for phenotype classification. Default: 0.2 (20th and 80th percentiles).
#'
#' @return
#' A list with three components:
#' \itemize{
#'   \item `metadata` - Data frame containing cell metadata with added columns:
#'     \itemize{
#'       \item `scPP_AUC` - AUCell scores for positive gene set
#'       \item `scPP_AUCwn` - AUCell scores for negative gene set
#'       \item `ScPP` - Phenotype classification: "Phenotype+", "Phenotype-", or "Background"
#'     }
#'   \item `Genes_pos` - Character vector of genes significantly upregulated
#'         in Phenotype+ cells compared to Phenotype- cells
#'   \item `Genes_neg` - Character vector of genes significantly upregulated
#'         in Phenotype- cells compared to Phenotype+ cells
#' }
#'
#' @details
#' This function implements an optimized workflow for single-cell phenotype
#' profiling that combines gene set enrichment analysis with differential
#' expression to identify and characterize cell phenotypes:
#'
#' ## Method Overview:
#' 1. **Gene Set Scoring**: Uses AUCell to calculate enrichment scores for
#'    both positive and negative gene sets in each cell
#' 2. **Phenotype Classification**: Identifies phenotype-positive and
#'    phenotype-negative cells based on quantile thresholds of AUC scores
#' 3. **Differential Expression**: Finds marker genes distinguishing the
#'    two phenotype groups using Seurat's FindMarkers
#'
#' ## Phenotype Classification Criteria:
#' - **Phenotype+**: High positive gene set score AND low negative gene set score
#' - **Phenotype-**: Low positive gene set score AND high negative gene set score
#' - **Background**: All other cells not meeting above criteria
#'
#' ## Differential Expression Thresholds:
#' - Absolute average log2 fold change > 1
#' - Adjusted p-value < 0.05
#'
#' @note
#' The function requires the Seurat object to have normalized data in the RNA
#' assay. For Seurat version compatibility, it handles both v4 and v5 data
#' structures automatically.
#'
#' @examples
#' \dontrun{
#' # Example using a Seurat object and gene lists
#' gene_list <- list(
#'   gene_pos = c("CD4", "IL7R", "CCR7"),
#'   gene_neg = c("CD8A", "CD8B", "GZMB")
#' )
#'
#' result <- ScPP.optimized(
#'   sc_dataset = seurat_obj,
#'   geneList = gene_list,
#'   probs = 0.2
#' )
#'
#' # Access results
#' head(result$metadata)
#' result$Genes_pos
#' result$Genes_neg
#' }
#'
#' @seealso
#' [AUCell::AUCell_buildRankings()], [AUCell::AUCell_calcAUC()],
#' [Seurat::FindMarkers()]
#'
#' @references
#' \url{https://github.com/WangX-Lab/ScPP/}, function ScPP()
#'
#' @keywords internal
#' @family scPP
#
ScPP.optimized <- function(sc_dataset, geneList, probs = 0.2) {
    chk::chk_length(geneList, 2)
    chk::chk_is(sc_dataset, "Seurat")
    chk::chk_range(probs)

    rna_data <- if (utils::packageVersion("Seurat") >= "5.0.0") {
        sc_dataset@assays$RNA$data
    } else {
        sc_dataset@assays$RNA@data
    }

    ts_cli$cli_alert_info("Computing AUC scores...")

    cellrankings <- AUCell::AUCell_buildRankings(rna_data, plotStats = FALSE)
    cellAUC <- AUCell::AUCell_calcAUC(geneList, cellrankings)

    auc_matrix <- AUCell::getAUC(cellAUC)
    auc_up <- as.numeric(auc_matrix["gene_pos", ])
    auc_down <- as.numeric(auc_matrix["gene_neg", ])

    metadata_dt <- data.table::as.data.table(
        sc_dataset[[]],
        keep.rownames = "cell_id"
    )
    metadata_dt[, `:=`(
        scPP_AUCup = auc_up,
        scPP_AUCdown = auc_down
    )]

    up_quantiles <- matrixStats::colQuantiles(
        matrix(c(auc_up, auc_down), ncol = 2),
        probs = c(probs, 1 - probs)
    )

    up_q1 <- up_quantiles[1, 1]
    up_q2 <- up_quantiles[1, 2]
    down_q1 <- up_quantiles[2, 1]
    down_q2 <- up_quantiles[2, 2]

    downcells1 <- metadata_dt[scPP_AUCup <= up_q1, cell_id]
    upcells1 <- metadata_dt[scPP_AUCup >= up_q2, cell_id]
    downcells2 <- metadata_dt[scPP_AUCdown >= down_q2, cell_id]
    upcells2 <- metadata_dt[scPP_AUCdown <= down_q1, cell_id]

    scPP_neg <- purrr::reduce(list(downcells1, downcells2), intersect)
    scPP_pos <- purrr::reduce(list(upcells1, upcells2), intersect)

    metadata_dt[, scPP := "Neutral"]
    metadata_dt[cell_id %chin% scPP_pos, scPP := "Positive"]
    metadata_dt[cell_id %chin% scPP_neg, scPP := "Negative"]

    sc_dataset$scPP <- metadata_dt$scPP
    Seurat::Idents(sc_dataset) <- "scPP"

    ts_cli$cli_alert_info("Finding markers...")

    markers <- Seurat::FindMarkers(
        sc_dataset,
        ident.1 = "Positive",
        ident.2 = "Negative"
    )

    markers_mat <- as.matrix(markers[, c("avg_log2FC", "p_val_adj")])

    pos_mask <- markers_mat[, "avg_log2FC"] > 1 &
        markers_mat[, "p_val_adj"] < 0.05
    neg_mask <- markers_mat[, "avg_log2FC"] < -1 &
        markers_mat[, "p_val_adj"] < 0.05

    genes_pos <- rownames(markers)[pos_mask]
    genes_neg <- rownames(markers)[neg_mask]

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

    return(list(
        metadata = as.data.frame(metadata_dt) %>%
            tibble::column_to_rownames("cell_id"),
        Genes_pos = genes_pos,
        Genes_neg = genes_neg
    ))
}
