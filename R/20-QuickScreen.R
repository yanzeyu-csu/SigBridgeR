#' @title Single-Cell Data Screening
#'
#' @description
#' Integrates matched bulk expression data and phenotype information to identify
#' phenotype-associated cell populations in single-cell RNA-seq data using one of
#' four computational methods. Ensures consistency between bulk and phenotype data
#' before analysis.
#'
#' @param matched_bulk Matrix or data frame of preprocessed bulk RNA-seq expression
#'        data (genes x samples). Column names must match names/IDs in `phenotype`.
#' @param sc_data A Seurat object containing scRNA-seq data to be screened.
#' @param phenotype Phenotype data, either:
#'        - Named vector (names match `matched_bulk` columns), or
#'        - Data frame with row names matching `matched_bulk` columns
#' @param label_type Character specifying phenotype label type (e.g., "SBS1", "time")
#' @param phenotype_class Type of phenotypic outcome (must be consistent with input data):
#'        - `"binary"`: Binary traits (e.g., case/control)
#'        - `"continuous"`: Continuous measurements (only for `Scissor`, `scPAS`, `scPP`)
#'        - `"survival"`: Survival objects
#' @param screen_method Screening algorithm to use, there are four options:
#'        - `"Scissor"`: see also `DoScissor()`
#'        - `"scPP"`: see also `DoscPP()`
#'        - `"scPAS"`: see also `DoscPAS()`
#'        - `"scAB"`: see also `DoscAB()`, no continuous support
#' @param ... Additional method-specific parameters:
#' \describe{
#'   \item{Scissor}{\describe{
#'     \item{scissor_alpha}{(numeric or NULL) Significance threshold. When NULL, alpha will keep increasing iteratively until the corresponding cells are screened out, default 0.05}
#'     \item{scissor_cutoff}{(numeric) A threshold for terminating the iteration of alpha, only work when `scissor_alpha` is NULL, default 0.2}
#'     \item{path2load_scissor_cache}{(character) default `NULL`}
#'     \item{path2save_scissor_inputs}{(character) A path to save the intermediary data. By using `path2load_scissor_cache`,  the intermediary data can be loaded from the specified path. default `"Scissor_inputs.RData"`}
#'     \item{nfold}{(integer) Cross-validation folds for reliability test, default 10}
#'     \item{reliability_test}{(logical) Whether to perform reliability test, default FALSE}
#'   }}
#'   \item{scPP}{\describe{
#'     \item{ref_group}{(integer or character) Reference group or baseline for **binary** comparisons, e.g. "Normal" for Tumor/Normal studies and 0 for 0/1 case-control studies. default: 0}
#'     \item{Log2FC_cutoff}{(numeric) Minimum log2 fold-change for binary markers, default 0.585}
#'     \item{estimate_cutoff}{(numeric) Effect size threshold for **continuous** traits, default 0.2}
#'     \item{probs}{(numeric) Quantile cutoff for cell classification, default 0.2}
#'   }}
#'   \item{scPAS}{\describe{
#'     \item{assay}{(character) Assay to use from sc_data, default "RNA"}
#'     \item{imputation}{(logical) Whether to perform imputation, default FALSE}
#'     \item{nfeature}{(integer) Number of features to select, default 3000}
#'     \item{alpha}{(numeric or NULL) Significance threshold, When NULL, alpha will keep increasing iteratively until the corresponding cells are screened out, default 0.01}
#'     \item{independent}{(logical) The background distribution of risk scores is constructed independently of each cell. default: TRUE}
#'     \item{network_class}{(character) Network class to use. default: 'SC', indicating gene-gene similarity networks derived from single-cell data. The other one is 'bulk'.}
#'     \item{permutation_times}{(integer) Number of permutations, default 2000}
#'     \item{FDR_threshold}{(numeric) FDR value threshold for identifying phenotype-associated cells default 0.05}
#'   }}
#'   \item{scAB}{\describe{
#'     \item{alpha}{(numeric) Coefficient of phenotype regularization ,default 0.005}
#'     \item{alpha_2}{(numeric) Coefficent of cell-cell similarity regularization, default 5e-05}
#'     \item{maxiter}{(integer) NMF optimization iterations, default 2000}
#'     \item{tred}{(integer) Z-score threshold, default 2}
#'   }}
#' }
#'
#' @return A list containing:
#' \describe{
#'   \item{scRNA_data}{Filtered Seurat object with phenotype-associated cells}
#'   \item{Some screen_result}{Important information about the screened result related to the selected method}
#' }
#'
#'
#' @section Data Matching Requirements:
#' - `matched_bulk` column names and `phenotype` names/rownames must be identical
#' - Phenotype values must correspond to bulk samples (not directly to single cells)
#' - Mismatches will trigger an error before analysis begins, and there is a built-in pre-run check.
#'
#' @section Method Compatibility:
#'
#' | **Method** | **Supported Phenotypes**      | **Additional Parameters**      |
#' |------------|-------------------------------|---------------------------------|
#' | `Scissor`  | All three types               | `scissor_alpha`, `scissor_cutoff`, `path2load_scissor_cache`, `path2save_scissor_inputs`, `nfold`, `reliability_test`, `reliability_test_n` |           |
#' | `scPP`     | All three types               | `ref_group`, `Log2FC_cutoff`, `estimate_cutoff`, `probs`                |
#' | `scPAS`    | All three types               | `n_components` ,`assay`, `imputation`,`nfeature`, `alpha`,`network_class`,`permutation_times`,`FDR_threshold`,`independent`               |
#' | `scAB`     | Binary/Survival               | `alpha`, `alpha_2`, `maxiter`, `tred`            |
#'
#'
#'
#' @seealso Associated functions:
#' \itemize{
#'   \item \code{\link{DoScissor}}
#'   \item \code{\link{DoscPP}}
#'   \item \code{\link{DoscPAS}}
#'   \item \code{\link{DoscAB}}
#' }
#'
#'
#' @export
#' @import dplyr
#' @importFrom glue glue
#'
Screen <- function(
    matched_bulk,
    sc_data,
    phenotype,
    label_type = NULL,
    phenotype_class = c("binary", "survival", "continuous"),
    screen_method = c("Scissor", "scPP", "scPAS", "scAB"),
    ...
) {
    chk::chk_subset(phenotype_class, c("binary", "survival", "continuous"))
    chk::chk_subset(screen_method, c("Scissor", "scPP", "scPAS", "scAB"))
    chk::chk_is(sc_data, "Seurat")

    if (is.null(label_type) || length(label_type) != 1) {
        cli::cli_alert_info(c(
            "{.var label_type} not specified or not of length {.val 1}, using {.val {screen_method}}"
        ))
        label_type = screen_method
    }

    screened_result =
        switch(
            screen_method,
            "Scissor" = {
                family = switch(
                    phenotype_class,
                    "binary" = "binomial",
                    "survival" = "cox",
                    "continuous" = "gaussian"
                )

                DoScissor(
                    sc_data = sc_data,
                    matched_bulk = matched_bulk,
                    phenotype = phenotype,
                    label_type = label_type,
                    scissor_family = family, # "gaussian", "binomial", "cox"
                    ...
                )
            },
            "scPAS" = {
                family = switch(
                    phenotype_class,
                    "binary" = "binomial",
                    "survival" = "cox",
                    "continuous" = "gaussian",
                )

                DoscPAS(
                    sc_data = sc_data,
                    matched_bulk = matched_bulk,
                    phenotype = phenotype,
                    label_type = label_type,
                    scPAS_family = family, # "gaussian", "binomial", "cox"
                    ...
                )
            },
            "scPP" = {
                phenotype_class = glue::glue(
                    toupper(substr(phenotype_class, 1, 1)),
                    tolower(substr(phenotype_class, 2, nchar(phenotype_class)))
                )

                DoscPP(
                    sc_data = sc_data,
                    matched_bulk = matched_bulk,
                    phenotype = phenotype,
                    label_type = label_type,
                    phenotype_class = phenotype_class, # "Binary", "Continuous", "Survival"
                    ...
                )
            },
            "scAB" = {
                if (phenotype_class == "continuous") {
                    cli::cli_abort(c(
                        "x" = "{.strong scAB} does not support continuous phenotype."
                    ))
                }

                DoscAB(
                    sc_data = sc_data,
                    matched_bulk = matched_bulk,
                    phenotype = phenotype,
                    label_type = label_type,
                    phenotype_class = phenotype_class, # "Binary", "Survival"
                    ...
                )
            },
            TRUE ~
                cli::cli_abort(c(
                    "x" = "Screen method not found.",
                    "i" = "Available methods: {.val Scissor}, {.val scPP}, {.val scPAS}, {.val scAB}"
                ))
        )

    return(screened_result)
}
