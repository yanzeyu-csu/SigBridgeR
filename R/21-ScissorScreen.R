# ? ---- 2. DO SCISSOR ----

#' @title Perform Scissor Screening Analysis
#' @description
#' Identifies phenotype-associated cell subpopulations in single-cell data using
#' regularized regression on matched bulk expression profiles.
#'
#' @usage
#' DoScissor(
#'    path2load_scissor_cache = NULL,
#'    path2save_scissor_inputs = "Scissor_inputs.RData",
#'    matched_bulk,
#'    sc_data,
#'    phenotype,
#'    label_type = "scissor",
#'    alpha = c(0.05,NULL),
#'    cutoff = 0.2,
#'    scissor_family = c("gaussian", "binomial", "cox"),
#'    reliability_test = FALSE,
#'    reliability_test.n = 10,
#'    reliability_test.nfold = 10,
#'    cell_evaluation = FALSE,
#'    cell_evaluation.benchmark_data = "path_to_file.RData",
#'    cell_evaluation.FDR = 0.05,
#'    cell_evaluation.bootstrap_n = 100,
#'    ...
#' )
#'
#' @param path2load_scissor_cache Path to precomputed Scissor inputs (RData file).
#'        If provided, skips recomputation (default: NULL).
#' @param path2save_scissor_inputs Path to save intermediate files (default: "Scissor_inputs.RData").
#' @param matched_bulk Normalized bulk expression matrix (features × samples).
#'        Column names must match `phenotype` identifiers.
#' @param sc_data Seurat object containing single-cell RNA-seq data.
#' @param phenotype Clinical outcome data. Can be:
#'        - Vector: named with sample IDs
#'        - Data frame: with row names matching bulk columns
#' @param label_type Character specifying phenotype label type (e.g., "SBS1", "time"), stored in `scRNA_data@misc`
#' @param alpha Parameter used to balance the effect of the l1 norm and the network-based penalties. It can be a number or a searching vector. If alpha = NULL, a default searching vector is used. The range of alpha is between 0 and 1. A larger alpha lays more emphasis on the l1 norm.
#' @param cutoff  (default: 0.2). When `alpha=NULL`, the cutoff is used to determine the optimal alpha.
#'        Higher values increase specificity.
#' @param scissor_family Model family for outcome type:
#'        - "gaussian": Continuous outcomes
#'        - "binomial": Binary outcomes (default)
#'        - "cox": Survival outcomes
#' @param reliability_test Logical to perform stability assessment when `scissor_alpha` is specified as a value between 0 and 1.(default: TRUE).
#' @param reliability_test.n Number of CV folds for reliability test (default: 10).
#'
#' @param reliability_test.nfold Cross-validation folds for reliability test (default: 10).
#' @param cell_evaluation Logical to perform cell evaluation (default: FALSE).
#' @param cell_evaluation.benchmark_data Path to benchmark data (RData file).
#' @param cell_evaluation.FDR FDR threshold for cell evaluation (default: 0.05).
#' @param cell_evaluation.bootstrap_n Number of bootstrap iterations for cell evaluation (default: 100).
#' @param ... Additional arguments. Currently supports:
#'    - `verbose`: Logical indicating whether to print progress messages. Defaults to `TRUE`.
#'    - `seed`: For reproducibility, default is `123L`
#'
#' @return A list containing:
#' \describe{
#'   \item{scRNA_data}{A Seurat object with screened cells containing metadata:
#'     \describe{
#'       \item{scissor}{"Positive"/"Negative"/"Neutral" classification}
#'       \item{label_type}{Outcome label used}
#'     }
#'   }
#'   \item{scissor_result}{Raw Scissor results}
#'   \item{reliability_result}{If reliability_test=TRUE, contains:
#'     \describe{
#'       \item{statistic}{A value between 0 and 1}
#'       \item{p}{p-value of the test statistic}
#'       \item{AUC_test_real}{10 values of AUC for real data}
#'       \item{AUC_test_back}{A list of AUC for background data}
#'     }
#'   }
#'   \item{cell_evaluation}{If cell_evaluation=TRUE, contains:
#'     \describe{
#'       \item{evaluation_res}{A data.frame with some supporting information for each Scissor selected cell}
#'     }
#'   }
#' }
#'
#' @references
#' Sun D, Guan X, Moran AE, Wu LY, Qian DZ, Schedin P, et al. Identifying phenotype-associated subpopulations by integrating bulk and single-cell sequencing data. Nat Biotechnol. 2022 Apr;40(4):527–38.
#'
#' @section LICENSE:
#' Licensed under the GNU General Public License version 3 (GPL-3.0).
#' A copy of the license is available at <https://www.gnu.org/licenses/gpl-3.0.en.html>.
#'
#' @examples
#' \dontrun{
#' # Binary outcome example
#' res <- DoScissor(
#'   matched_bulk = bulk_matrix,
#'   sc_data = seurat_obj,
#'   phenotype = a_named_vector,
#'   scissor_family = "binomial"
#' )
#' }
#'
#' @export
#' @family screen_method
#' @family scissor
#'
DoScissor <- function(
    path2load_scissor_cache = NULL,
    path2save_scissor_inputs = "Scissor_inputs.RData",
    matched_bulk,
    sc_data,
    phenotype,
    label_type = "scissor",
    alpha = c(0.05, NULL),
    cutoff = 0.2,
    scissor_family = c("gaussian", "binomial", "cox"),
    reliability_test = FALSE,
    reliability_test.n = 10,
    reliability_test.nfold = 10,
    cell_evaluation = FALSE,
    cell_evaluation.benchmark_data = "path_to_file.RData",
    cell_evaluation.FDR = 0.05,
    cell_evaluation.bootstrap_n = 100,
    ...
) {
    # Input validation
    chk::chk_is(matched_bulk, c("matrix", "data.frame"))
    chk::chk_is(sc_data, "Seurat")
    chk::chk_character(label_type)
    chk::chk_range(cutoff)
    scissor_family <- SigBridgeRUtils::MatchArg(
        scissor_family,
        c("gaussian", "binomial", "cox"),
        NULL
    )
    chk::chk_flag(reliability_test)
    chk::chk_flag(cell_evaluation)

    # get from dots
    dots <- rlang::list2(...)
    verbose <- dots$verbose %||% SigBridgeRUtils::getFuncOption("verbose")
    seed <- dots$seed %||% SigBridgeRUtils::getFuncOption("seed")

    if (scissor_family %chin% c("binomial", "cox")) {
        label_type_scissor <- c(
            glue::glue("{label_type}_Negative"),
            glue::glue("{label_type}_Positive")
        )
    } else if (scissor_family == "gaussian") {
        n <- length(table(phenotype))
        label_type_scissor <- glue::glue("{label_type}_{seq_len(n)}")
    }

    if (!is.null(path2save_scissor_inputs)) {
        path <- dirname(path2save_scissor_inputs)
        if (!dir.exists(path)) {
            dir.create(path, recursive = TRUE)
        }
    }

    infos1 <- Scissor::Scissor.v5.optimized(
        bulk_dataset = matched_bulk,
        sc_dataset = sc_data,
        phenotype = phenotype,
        tag = label_type_scissor,
        alpha = alpha,
        cutoff = cutoff,
        family = scissor_family,
        Save_file = path2save_scissor_inputs,
        Load_file = path2load_scissor_cache,
        verbose = verbose,
        seed = seed
    )

    # meta.data to add
    sc_meta <- data.frame(
        scissor = rep("Neutral", ncol(sc_data)),
        row.names = colnames(sc_data)
    )
    sc_meta$scissor[rownames(sc_meta) %chin% infos1$Scissor_pos] <- "Positive"
    sc_meta$scissor[rownames(sc_meta) %chin% infos1$Scissor_neg] <- "Negative"
    sc_data <- Seurat::AddMetaData(sc_data, metadata = sc_meta) %>%
        SigBridgeRUtils::AddMisc(
            scissor_type = label_type,
            scissor_para = infos1$para,
            cover = FALSE
        )

    # *reliability test
    if (reliability_test) {
        purrr::walk(
            c(reliability_test.n, reliability_test.nfold, alpha),
            chk::chk_number
        )

        # indicate that Y has only two levels, both Pos and Neg cells exist
        if (length(table(infos1$Y)) < 2) {
            cli::cli_abort(c(
                "x" = "Error in reliability test:",
                ">" = "one of the Pos or Neg cells doesn't exist"
            ))
        }
        if (verbose) {
            ts_cli$cli_alert_info(
                cli::col_green("Start reliability test")
            )
        }

        reliability_result <- Scissor::reliability.test(
            infos1$X,
            infos1$Y,
            infos1$network,
            alpha = alpha,
            family = scissor_family,
            cell_num = length(infos1$Scissor_pos) +
                length(infos1$Scissor_neg),
            n = reliability_test.n,
            nfold = reliability_test.nfold,
            verbose = verbose
        )
        if (verbose) {
            ts_cli$cli_alert_info(
                cli::col_green("reliability test: Done")
            )
        }
    } else {
        reliability_result <- NULL
    }

    # * cell_evaluation
    if (cell_evaluation) {
        chk::chk_file(cell_evaluation.benchmark_data)
        chk::chk_number(cell_evaluation.bootstrap_n)
        chk::chk_range(cell_evaluation.FDR)
        if (verbose) {
            ts_cli$cli_alert_info(
                cli::col_green("Start cell evalutaion")
            )
        }

        evaluate_res <- Scissor::evaluate.cell(
            Load_file = cell_evaluation.benchmark_data,
            Scissor_result = infos1,
            FDR_cutoff = cell_evaluation.FDR,
            bootstrap_n = cell_evaluation.bootstrap_n
        )
        if (verbose) {
            ts_cli$cli_alert_info(
                cli::col_green("Cell evalutaion: Done")
            )
        }
    } else {
        evaluate_res <- NULL
    }

    list(
        scRNA_data = sc_data,
        scissor_result = infos1, # parameters included
        reliability_result = reliability_result,
        cell_evaluation = evaluate_res
    )
}
