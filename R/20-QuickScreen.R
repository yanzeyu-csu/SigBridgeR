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
#'     \item{scissor_alpha}{(numeric) default 0.05}
#'     \item{scissor_cutoff}{(numeric) default 0.2}
#'     \item{path2load_scissor_cache}{(character) default `NULL`}
#'     \item{path2save_scissor_inputs}{(character) default `"Scissor_inputs.RData"`}
#'     \item{nfold}{(integer) Cross-validation folds for reliability test, default 10}
#'     \item{reliability_test}{(logical) default FALSE}
#'   }}
#'   \item{scPP}{\describe{
#'     \item{ref_group}{(integer) Reference group for binary comparisons, default 1}
#'     \item{Log2FC_cutoff}{(numeric) Minimum log2 fold-change for binary markers, default 0.585}
#'     \item{estimate_cutoff}{(numeric) Effect size threshold for continuous traits, default 0.2}
#'     \item{probs}{(numeric) Quantile cutoff for cell classification, default 0.2}
#'   }}
#'   \item{scPAS}{\describe{
#'     \item{assay}{(character) Assay to use from sc_data, default "RNA"}
#'     \item{imputation}{(logical) Whether to perform imputation, default FALSE}
#'     \item{nfeature}{(integer) Number of features to select, default 3000}
#'     \item{alpha}{(numeric) Significance threshold, default 0.01}
#'     \item{extra_filter}{(logical) Whether to perform extra filtering based on RNA count, default FALSE}
#'     \item{gene_RNAcount_filter}{(numeric) Only used when `extra_filter=TRUE`, default 20}
#'     \item{bulk_0_filter_thresh}{(numeric) Only used when `extra_filter=TRUE`,default 0.25}
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
#'   \item{matched_samples}{Vector of samples used in the analysis}
#'   \item{method_output}{Method-specific output objects}
#' }
#'
#'
#' @section Data Matching Requirements:
#' - `matched_bulk` column names and `phenotype` names/rownames must be identical
#' - Phenotype values must correspond to bulk samples (not directly to single cells)
#' - Mismatches will trigger an error before analysis begins
#'
#' @section Method Compatibility:
#'
#' | **Method** | **Supported Phenotypes**      | **Additional Parameters**      |
#' |------------|-------------------------------|---------------------------------|
#' | `Scissor`  | All three types               | `alpha`, `lambda`               |
#' | `scPP`     | All three types               | `embedding_type`                |
#' | `scPAS`    | All three types               | `n_components`                  |
#' | `scAB`     | Binary/Survival               | `bandwidth`                     |
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
    library(dplyr)
    if (length(screen_method) != 1) {
        cli::cli_abort(c("x" = "Only one {.arg screen_method} is allowed."))
    }

    phenotype_class = tolower(phenotype_class)
    if (length(phenotype_class) != 1) {
        cli::cli_abort(c("x" = "Only one {.arg phenotype_class} is allowed."))
    } else if (!phenotype_class %in% c("binary", "survival", "continuous")) {
        cli::cli_abort(c(
            "x" = "Invalid {.arg phenotype_class = {phenotype_class}}.",
            "i" = " Must be one of {.val binary}, {.val survival}, or {.val continuous}."
        ))
    }

    if (is.null(label_type) || length(label_type) != 1) {
        cli::cli_warn(c(
            "{.var label_type} not specified or not of length {.val 1}, using {.val {screen_method}}"
        ))
        label_type = screen_method
    }

    screened_result = tolower(screen_method) %>%
        switch(
            "scissor" = {
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
            "scpas" = {
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
            "scpp" = {
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
            "scab" = {
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
            TRUE ~ stop("Screen method not found.")
        )

    return(screened_result)
}
