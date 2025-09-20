# ---- 2. Do scPAS ----

#' @title Perform scPAS Screening Analysis
#' @description
#' This function performs scPAS screening analysis by integrating bulk and single-cell RNA-seq data.
#' It includes data filtering steps and wraps the core scPAS::scPAS function.
#'
#' @param matched_bulk Bulk RNA-seq data (genes x samples)
#' @param sc_data Single-cell RNA-seq data (Seurat object and preprocessed)
#' @param phenotype Phenotype data frame with sample annotations
#' @param label_type Character specifying phenotype label type (e.g., "SBS1", "time"), stored in `scRNA_data@misc`.  (default: "scPAS")
#' @param assay Assay to use from sc_data (default: 'RNA')
#' @param imputation Logical, whether to perform imputation (default: FALSE)
#' @param imputation_method Character. Name of alternative method for imputation. (options: "KNN", "ALRA")
#' @param nfeature Number of features to select (default: 3000, indicating that the top 3000 highly variable genes are selected for model training
#' @param alpha Numeric. Significance threshold. Parameter used to balance the effect of the l1 norm and the network-based penalties. It can be a number or a searching vector. If alpha = NULL, a default searching vector is used. The range of alpha is in `[0,1]`. A larger alpha lays more emphasis on the l1 norm. (default: 0.01)
#' @param network_class Network class to use (default: 'SC', indicating gene-gene similarity networks derived from single-cell data. The other one is 'bulk'.)
#' @param scPAS_family Model family for analysis (options: "cox", "gaussian", "binomial")
#' @param permutation_times Number of permutations to perform (default: 2000)
#' @param FDR_threshold Numeric. FDR value threshold for identifying phenotype-associated cells (default: 0.05)
#' @param independent Logical. The background distribution of risk scores is constructed independently of each cell. (default: TRUE)
#' @param ... Additional arguments passed to `DoscPAS` functions
#'
#' @return A Seurat object from scPAS analysis
#'
#' @references
#' Xie A, Wang H, Zhao J, Wang Z, Xu J, Xu Y. scPAS: single-cell phenotype-associated subpopulation identifier. Briefings in Bioinformatics. 2024 Nov 22;26(1):bbae655.
#'
#' @section LICENSE:
#' Licensed under the GNU General Public License version 3 (GPL-3.0).
#' A copy of the license is available at <https://www.gnu.org/licenses/gpl-3.0.en.html>.
#'
#'
#' @importFrom cli cli_alert_info cli_alert_success
#' @importFrom Matrix rowSums as.matrix
#' @import dplyr
#'
#' @family screen_method
#'
#' @keywords internal
#' @export
#'
DoscPAS = function(
    matched_bulk,
    sc_data,
    phenotype,
    label_type = "scPAS",
    assay = 'RNA',
    imputation = F,
    imputation_method = c("KNN", "ALRA"),
    nfeature = 3000,
    alpha = c(0.01, NULL),
    network_class = c("SC", "bulk"),
    scPAS_family = c("cox", "gaussian", "binomial"),
    permutation_times = 2000,
    FDR_threshold = 0.05,
    independent = TRUE,
    ...
) {
    chk::chk_is(matched_bulk, c("matrix", "data.frame"))
    chk::chk_is(sc_data, "Seurat")
    chk::chk_character(label_type)
    chk::chk_flag(imputation)
    if (imputation) {
        chk::chk_length(imputation_method, 1)
        chk::chk_subset(imputation_method, c("KNN", "ALRA"))
    }
    chk::chk_number(nfeature)
    chk::chk_null_or(alpha, chk::chk_range)
    chk::chk_subset(network_class, c("SC", "bulk"))
    if (length(network_class) > 1) {
        network_class <- network_class[1]
    }
    chk::chk_subset(scPAS_family, c("cox", "gaussian", "binomial"))
    chk::chk_length(scPAS_family, 1)
    chk::chk_number(permutation_times)
    chk::chk_number(FDR_threshold)
    chk::chk_flag(independent)

    # robust
    if (scPAS_family == "cox") {
        if (is.null(intersect(colnames(matched_bulk), rownames(phenotype)))) {
            cli::cli_abort(c(
                "x" = "No intersection between the rownames of {.var phenotype} and colnames of {.var matched_bulk}."
            ))
        }
    } else {
        if (is.null(intersect(colnames(matched_bulk), names(phenotype)))) {
            cli::cli_abort(c(
                "x" = "No intersection between the names of {.var phenotype} and colnames of {.var matched_bulk}."
            ))
        }
    }

    cli::cli_alert_info(c(
        "[{TimeStamp()}]",
        crayon::green(" Start scPAS screening.")
    ))

    scPAS_result <- scPAS::scPAS(
        bulk_dataset = as.matrix(matched_bulk),
        sc_dataset = sc_data,
        assay = 'RNA',
        tag = label_type,
        phenotype = phenotype,
        imputation = imputation,
        imputation_method = imputation_method,
        nfeature = nfeature,
        alpha = alpha,
        network_class = network_class,
        family = scPAS_family,
        independent = independent,
        permutation_times = permutation_times,
        FDR.threshold = FDR_threshold,
        ...
    ) %>%
        AddMisc(scPAS_type = label_type, cover = FALSE)

    # *rename level
    scPAS_result@meta.data = dplyr::mutate(
        scPAS_result@meta.data,
        `scPAS` = dplyr::case_when(
            scPAS == "scPAS+" ~ "Positive",
            scPAS == "scPAS-" ~ "Negative",
            TRUE ~ "Neutral"
        )
    )

    detailed_info <- dplyr::select(
        scPAS_result@meta.data,
        dplyr::contains("scPAS_")
    )

    cli::cli_alert_success(c(
        "[{TimeStamp()}]",
        crayon::green(" scPAS screening done.")
    ))

    return(list(scRNA_data = scPAS_result, stats = detailed_info))
}
