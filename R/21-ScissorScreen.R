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
#' @param ... Additional arguments passed to `Scissor.v5.optimized`.
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
#' @importFrom dplyr %>%
#' @importFrom Seurat AddMetaData
#' @importFrom Scissor reliability.test
#' @importFrom stats cor
#' @importFrom glue glue
#' @importFrom cli cli_abort cli_alert_info cli_alert_success cli_alert_danger col_green
#'
#' @keywords internal
#' @family screen_method
#'
DoScissor = function(
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
    scissor_family <- match.arg(scissor_family)
    chk::chk_character(path2save_scissor_inputs)
    chk::chk_flag(reliability_test)
    chk::chk_flag(cell_evaluation)

    if (scissor_family %chin% c("binomial", "cox")) {
        label_type_scissor = c(
            glue::glue("{label_type}_Negative"),
            glue::glue("{label_type}_Positive")
        )
    } else if (scissor_family == "gaussian") {
        n = length(table(phenotype))
        label_type_scissor = glue::glue("{label_type}_{seq_len(n)}")
    } else if (length(scissor_family) != 1) {
        cli::cli_abort(
            c(
                "x" = "Please choose one scissor family, use parameter {.var scissor_family}."
            ),
            class = "FamilyError"
        )
    }
    if (!is.null(path2save_scissor_inputs)) {
        path = dirname(path2save_scissor_inputs)
        if (!dir.exists(path)) {
            dir.create(path, recursive = TRUE)
        }
    }

    infos1 <- Scissor.v5.optimized(
        bulk_dataset = matched_bulk,
        sc_dataset = sc_data,
        phenotype = phenotype,
        tag = label_type_scissor,
        alpha = alpha,
        cutoff = cutoff,
        family = scissor_family,
        Save_file = path2save_scissor_inputs,
        Load_file = path2load_scissor_cache,
        ...
    )

    # meta.data to add
    sc_meta <- data.frame(
        scissor = rep("Neutral", ncol(sc_data)),
        row.names = colnames(sc_data)
    )
    sc_meta$scissor[rownames(sc_meta) %chin% infos1$Scissor_pos] <- "Positive"
    sc_meta$scissor[rownames(sc_meta) %chin% infos1$Scissor_neg] <- "Negative"
    sc_data <- Seurat::AddMetaData(sc_data, metadata = sc_meta) %>%
        AddMisc(scissor_type = label_type, cover = FALSE)

    # *reliability test
    if (reliability_test) {
        chk::chk_number(reliability_test.n)
        chk::chk_number(reliability_test.nfold)
        chk::chk_number(alpha)

        # indicate that Y has only two levels, both Pos and Neg cells exist
        if (!length(table(infos1$Y)) < 2) {
            ts_cli$cli_alert_info(
                cli::col_green("Start reliability test")
            )

            reliability_result <- Scissor::reliability.test(
                infos1$X,
                infos1$Y,
                infos1$network,
                alpha = alpha,
                family = scissor_family,
                cell_num = length(infos1$Scissor_pos) +
                    length(infos1$Scissor_neg),
                n = reliability_test.n,
                nfold = reliability_test.nfold
            )

            cli::cli_alert_success(
                "reliability test: Done"
            )
        } else {
            cli::cli_abort(c(
                "x" = "Error in reliability test:",
                ">" = "one of the Pos or Neg cells doesn't exist"
            ))
        }
    } else {
        reliability_result <- NULL
    }

    # * cell_evaluation
    if (cell_evaluation) {
        chk::chk_file(cell_evaluation.benchmark_data)
        chk::chk_number(cell_evaluation.bootstrap_n)
        chk::chk_range(cell_evaluation.FDR)

        ts_cli$cli_alert_info(
            cli::col_green("Start cell evalutaion")
        )

        evaluate_res <- Scissor::evaluate.cell(
            Load_file = cell_evaluation.benchmark_data,
            Scissor_result = infos1,
            FDR_cutoff = cell_evaluation.FDR,
            bootstrap_n = cell_evaluation.bootstrap_n
        )

        cli::cli_alert_success(
            "cell evalutaion: Done"
        )
    } else {
        evaluate_res <- NULL
    }

    return(list(
        scRNA_data = sc_data,
        scissor_result = infos1, # parameters
        reliability_result = reliability_result,
        cell_evaluation = evaluate_res
    ))
}

#' @title Optimized Scissor Algorithm for Seurat ver5
#' @description
#' Scissor.v5 from `https://doi.org/10.1038/s41587-021-01091-3`and `https://github.com/sunduanchen/Scissor/issues/59`
#' Another version of Scissor.v5() to optimize memory usage and execution speed in preprocess.
#'
#' @references
#' Sun D, Guan X, Moran AE, Wu LY, Qian DZ, Schedin P, et al. Identifying phenotype-associated subpopulations by integrating bulk and single-cell sequencing data. Nat Biotechnol. 2022 Apr;40(4):527–38.
#'
#' @section LICENSE:
#' Licensed under the GNU General Public License version 3 (GPL-3.0).
#' A copy of the license is available at <https://www.gnu.org/licenses/gpl-3.0.en.html>.
#'
#'
#' @family screen_method
#'
#' @keywords internal
#'
Scissor.v5.optimized <- function(
    bulk_dataset,
    sc_dataset,
    phenotype,
    tag = NULL,
    alpha = NULL,
    cutoff = 0.2,
    family = c("gaussian", "binomial", "cox"),
    Save_file = "Scissor_inputs.RData",
    Load_file = NULL,
    workers = 4,
    ...
) {
    ts_cli$cli_alert_info(
        cli::col_green("Scissor start...")
    )

    if (is.null(Load_file)) {
        ts_cli$cli_alert_info("Start from raw data...")
        common = intersect(
            rownames(bulk_dataset),
            rownames(sc_dataset)
        )
        if (length(common) == 0) {
            cli::cli_abort(c(
                "x" = "There is no common genes between the given single-cell and bulk samples. Please check Scissor inputs."
            ))
        }

        if (inherits(sc_dataset, "Seurat")) {
            sc_exprs <- as.matrix(sc_dataset@assays$RNA$data)
            if ("RNA_snn" %chin% names(sc_dataset@graphs)) {
                network <- as.matrix(sc_dataset@graphs$RNA_snn)
                cli::cli_alert_info(
                    "Using {.val RNA_snn} graph for network."
                )
            } else if ("integrated_snn" %chin% names(sc_dataset@graphs)) {
                network <- as.matrix(sc_dataset@graphs$integrated_snn)
                cli::cli_alert_info(
                    "Using {.val integrated_snn} graph for network."
                )
            } else {
                cli::cli_abort(c(
                    "x" = "No `RNA_snn` or `integrated_snn` graph in the given Seurat object. Please check Scissor inputs."
                ))
            }
        } else {
            sc_exprs <- as.matrix(sc_dataset)
            Seurat_tmp <- SCPreProcess(
                sc_dataset,
                quality_control.pattern = c("^MT-"),
                verbose = FALSE
            )
            network <- as.matrix(Seurat_tmp@graphs$RNA_snn)
        }
        diag(network) <- 0
        network <- (network != 0) * 1

        bulk_mat <- as.matrix(bulk_dataset[common, ])
        sc_mat <- as.matrix(sc_exprs[common, ])
        dataset0 <- cbind(bulk_mat, sc_mat)

        ts_cli$cli_alert_info(
            "Normalizing quantiles of data..."
        )

        dataset1 <- preprocessCore::normalize.quantiles(as.matrix(dataset0))
        rownames(dataset1) <- common
        colnames(dataset1) <- c(colnames(bulk_mat), colnames(sc_mat))

        ts_cli$cli_alert_info(
            "Subsetting data..."
        )

        n_bulk <- ncol(bulk_mat)
        # gene-sample
        Expression_bulk <- dataset1[, 1:n_bulk, drop = FALSE]
        # gene-cell
        Expression_cell <- dataset1[, (n_bulk + 1):ncol(dataset1), drop = FALSE]

        gc(verbose = FALSE)

        ts_cli$cli_alert_info(
            "Calculating correlation..."
        )

        X <- cor(Expression_bulk, Expression_cell)
        quality_check <- matrixStats::colQuantiles(X, probs = seq(0, 1, 0.25))

        cat(strrep("-", floor(getOption("width") / 2)), "\n", sep = "")
        cli::cli_text("Five-number summary of correlations:")
        print(quality_check %>% asplit(2) %>% purrr::map_dbl(mean))
        cat(strrep("-", floor(getOption("width") / 2)), "\n", sep = "")
        # median
        if (quality_check[3] < 0.01) {
            cli::cli_warn(
                "The median correlation between the single-cell and bulk samples is relatively low."
            )
        }

        FamilyProcessor <- list(
            binomial = function() {
                Y <- as.numeric(phenotype)
                z <- table(Y)
                if (length(z) != length(tag)) {
                    cli::cli_abort(
                        "x" = "The length differs between tags and phenotypes. Please check Scissor inputs and selected regression type."
                    )
                }
                cli::cli_alert_info(
                    "Current phenotype contains {.val {z[1]}} {tag[1]} and {.val {z[2]}} {tag[2]} samples."
                )
                ts_cli$cli_alert_info(
                    "Perform logistic regression on the given phenotypes..."
                )
                Y
            },
            gaussian = function() {
                Y <- as.numeric(phenotype)
                z <- table(Y)
                if (length(z) != length(tag)) {
                    cli::cli_abort(
                        "x" = "The length differs between tags and phenotypes. Please check Scissor inputs and selected regression type.",
                        "i" = "length of tags: {.val {length(tag)}}",
                        "i" = "length of phenotypes: {.val {length(z)}}"
                    )
                }
                tmp <- paste(z, tag)
                cli::cli_alert_info(
                    "Current phenotype contains: {.val {length(tmp)}} samples."
                )
                cli::cli_text("Sample examples:")
                cli::cli_bullets(c(
                    " " = "{.val {head(tmp, 5)}}",
                    " " = "... ({length(tmp)-6} more samples)",
                    " " = "{.val {tail(tmp, 1)}}"
                ))
                ts_cli$cli_alert_info(
                    "Perform linear regression on the given phenotypes..."
                )
                Y
            },
            cox = function() {
                Y <- as.matrix(phenotype)
                if (ncol(Y) != 2) {
                    cli::cli_abort(
                        "x" = "The size of survival data is wrong. Please check Scissor inputs and selected regression type."
                    )
                }
                ts_cli$cli_alert_info(
                    "Perform cox regression on the given clinical outcomes..."
                )
                Y
            }
        )

        family <- match.arg(family)
        Y <- FamilyProcessor[[family]]()

        if (!is.null(Save_file)) {
            save(
                X,
                Y,
                network,
                Expression_bulk,
                Expression_cell,
                file = Save_file
            )
            ts_cli$cli_alert_success(glue::glue(
                "Statistics data saved to `{Save_file}`."
            ))
        }
    } else {
        # Load data from previous work
        ts_cli$cli_alert_info(
            glue::glue("Loading data from `{Load_file}`...")
        )
        load(Load_file)
    }

    # garbage collection
    rm(
        Expression_bulk,
        Expression_cell,
        sc_dataset,
        bulk_dataset,
        phenotype
    )

    ts_cli$cli_alert_info("Screening...")

    alpha <- alpha %||%
        c(0.005, 0.01, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9)

    for (i in seq_along(alpha)) {
        tryCatch(
            {
                set.seed(123)

                fit0 <- Scissor::APML1(
                    X,
                    Y,
                    family = family,
                    penalty = "Net",
                    alpha = alpha[i],
                    Omega = network,
                    nlambda = 100,
                    nfolds = min(10, nrow(X))
                )
                fit1 <- Scissor::APML1(
                    X,
                    Y,
                    family = family,
                    penalty = "Net",
                    alpha = alpha[i],
                    Omega = network,
                    lambda = fit0$lambda.min
                )
                if (family == "binomial") {
                    Coefs <- as.numeric(fit1$Beta[2:(ncol(X) + 1)])
                } else {
                    Coefs <- as.numeric(fit1$Beta)
                }
                pos_mask <- Coefs > 0
                neg_mask <- Coefs < 0
                Cell1 <- colnames(X)[pos_mask]
                Cell2 <- colnames(X)[neg_mask]
                percentage <- (length(Cell1) + length(Cell2)) / ncol(X)
                cli::cli_h2("At alpha = {.val {alpha[i]}}")
                cli::cli_text(sprintf(
                    "Scissor identified {.val {%d}} Scissor+ cells and {.val {%d}} Scissor- cells.",
                    length(Cell1),
                    length(Cell2)
                ))
                cli::cli_text(sprintf(
                    "The percentage of selected cell is: {.val {%s}}%%",
                    round(percentage * 100, digits = 3)
                ))
                if (percentage < cutoff) {
                    ts_cli$cli_alert_info(
                        cli::col_green("Scissor Ended.")
                    )
                    break
                }
                cat("\n")
            },
            error = function(e) {
                cli::cli_alert_danger(e$message)

                ts_cli$cli_alert_info(
                    cli::col_yellow("Scissor screening exit 1.")
                )
            }
        )
    }
    return(list(
        para = list(
            alpha = alpha[i],
            lambda = fit0$lambda.min,
            family = family
        ),
        Coefs = Coefs,
        Scissor_pos = Cell1,
        Scissor_neg = Cell2,
        X = X,
        Y = Y,
        network = network
    ))
}
