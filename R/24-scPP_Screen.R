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
#'       \item{label_type}{Outcome label used}
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
#' @importFrom crayon green
#' @importFrom cli cli_abort cli_alert_info cli_alert_success
#' @importFrom tibble rownames_to_column
#' @importFrom data.table as.data.table fifelse
#' @importFrom stats quantile
#'
#' @family screen_method
#'
#' @keywords internal
#' @export
#'
DoscPP = function(
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
    chk::chk_subset(phenotype_class, c("Binary", "Continuous", "Survival"))
    chk::chk_length(phenotype_class, 1)
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

    cli::cli_alert_info(c(
        "[{TimeStamp()}]",
        crayon::green(" Start scPP screening.")
    ))

    cli::cli_alert_info(c(
        "[{TimeStamp()}]",
        " Finding markers..."
    ))

    matched_bulk = as.data.frame(matched_bulk)
    # decide which type of phenotype data is used
    if (is.vector(phenotype)) {
        phenotype = as.data.frame(phenotype) %>%
            tibble::rownames_to_column("Sample") %>%
            dplyr::rename("Feature" = 2) %>%
            dplyr::mutate(Feature = as.numeric(Feature))
    }
    if (tolower(phenotype_class) == "binary") {
        Check0VarRows(matched_bulk)
        gene_list = ScPP::marker_Binary(
            bulk_data = matched_bulk,
            features = phenotype,
            ref_group = ref_group,
            Log2FC_cutoff = Log2FC_cutoff
        )
    } else if (tolower(phenotype_class) == "continuous") {
        gene_list = ScPP::marker_Continuous(
            bulk_data = matched_bulk,
            features = phenotype$Feature,
            estimate_cutoff = estimate_cutoff
        )
    } else if (tolower(phenotype_class) == "survival") {
        gene_list = ScPP::marker_Survival(
            bulk_data = matched_bulk,
            survival_data = phenotype
        )
    }

    l = lapply(gene_list, length)
    pos_null = FALSE
    neg_null = FALSE
    if ("gene_pos" %in% names(l)) {
        # Cannot combine the conditions due to the feature of `gene_list`
        if (l[["gene_pos"]] == 0) {
            cli::cli_alert_info(c(
                "[{TimeStamp()}]",
                " No significant positive genes found"
            ))
            pos_null = TRUE
        }
    }
    if ("gene_neg" %in% names(l)) {
        if (l[["gene_neg"]] == 0) {
            cli::cli_alert_info(c(
                "[{TimeStamp()}]",
                " No significant negative genes found"
            ))
            neg_null = TRUE
        }
    }
    if (pos_null & neg_null) {
        cli::cli_alert_info(c(
            "[{TimeStamp()}]",
            crayon::yellow(" scPP screening exit.")
        ))
        return(list(scRNA_data = NULL))
    }

    cli::cli_alert_info(c(
        "[{TimeStamp()}]",
        " Screening..."
    ))

    # *Start screen
    err_flag = FALSE
    tryCatch(
        scPP_result <- ScPP::ScPP(sc_data, gene_list, probs = probs),
        error = function(e) {
            cli::cli_alert_danger(c("[{TimeStamp()}] ", e$message))

            cli::cli_alert_info(c(
                "[{TimeStamp()}]",
                crayon::yellow(" scPP screening exit.")
            ))

            err_flag = TRUE
        }
    )
    if (err_flag) {
        return(list(
            scRNA_data = "`scPP` is not applicable to the current data."
        ))
    }

    sc_data@meta.data[, "scPP"] <- data.table::as.data.table(
        scPP_result$metadata
    )[,
        .(
            scPP = data.table::fifelse(
                ScPP == "Phenotype+",
                "Positive",
                data.table::fifelse(
                    ScPP == "Phenotype-",
                    "Negative",
                    "Neutral"
                )
            )
        )
    ]$scPP

    sc_data <- sc_data %>%
        AddMisc(scPP_type = label_type, cover = FALSE)

    cli::cli_alert_success(c(
        "[{TimeStamp()}]",
        crayon::green(" scPP screening done.")
    ))

    return(list(scRNA_data = sc_data, gene_list = gene_list))
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
#' For dense matrices, variance is computed using an optimized `rowVars` function that
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
#' @export
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
        vars <- rowVars(mat, na.rm = TRUE)
    }

    zero_idx <- which(vars == 0 | is.na(vars))
    if (length(zero_idx)) {
        bad_genes <- rownames(mat)[zero_idx]
        cli::cli_abort(
            c(
                "Detected {.val {length(bad_genes)}} gene(s) with zero variance:",
                "i" = "{.val {bad_genes}}"
            ),
            class = "zero_variance_error",
            call = call
        )
    }
    invisible(vars)
}
