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
#'
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
    scissor_family <- MatchArg(
        scissor_family,
        c("gaussian", "binomial", "cox"),
        NULL
    )
    chk::chk_character(path2save_scissor_inputs)
    chk::chk_flag(reliability_test)
    chk::chk_flag(cell_evaluation)

    # get from dots
    dots <- rlang::list2(...)
    verbose <- dots$verbose %||% getFuncOption("verbose")
    seed <- dots$seed %||% getFuncOption("seed")

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
        AddMisc(
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
#' @family scissor
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
    verbose = getFuncOption("verbose"),
    seed = getFuncOption("seed"),
    ...
) {
    if (verbose) {
        ts_cli$cli_alert_info(
            cli::col_green("Scissor start...")
        )
    }

    if (is.null(Load_file)) {
        if (verbose) {
            ts_cli$cli_alert_info("Start from raw data...")
        }
        common <- intersect(
            rownames(bulk_dataset),
            rownames(sc_dataset)
        )
        if (length(common) == 0) {
            cli::cli_abort(c(
                "x" = "There is no common genes between the given single-cell and bulk samples. Please check Scissor inputs."
            ))
        }

        if (inherits(sc_dataset, "Seurat")) {
            # sc_exprs <- as.matrix(sc_dataset@assays$RNA$data)
            sc_exprs <- SeuratObject::LayerData(sc_dataset, layer = "data")

            if ("RNA_snn" %chin% names(sc_dataset@graphs)) {
                # network <- as.matrix(sc_dataset@graphs$RNA_snn)
                network <- SeuratObject::Graphs(
                    sc_dataset,
                    slot = "RNA_snn"
                )

                if (verbose) {
                    cli::cli_alert_info(
                        "Using {.val RNA_snn} graph for network."
                    )
                }
            } else if ("integrated_snn" %chin% names(sc_dataset@graphs)) {
                # network <- as.matrix(sc_dataset@graphs$integrated_snn)
                network <- SeuratObject::Graphs(
                    sc_dataset,
                    slot = "integrated_snn"
                )

                if (verbose) {
                    cli::cli_alert_info(
                        "Using {.val integrated_snn} graph for network."
                    )
                }
            } else {
                cli::cli_abort(c(
                    "x" = "No `RNA_snn` or `integrated_snn` graph in the given Seurat object. Please check Scissor inputs."
                ))
            }
        } else {
            sc_exprs <- Matrix::Matrix(as.matrix(sc_dataset))
            Seurat_tmp <- SCPreProcess(
                sc_dataset,
                quality_control.pattern = c("^MT-"),
                verbose = FALSE
            )
            network <- SeuratObject::Graphs(
                Seurat_tmp,
                slot = "RNA_snn"
            )
            rm(Seurat_tmp)
        }
        Matrix::diag(network) <- 0
        network <- as.matrix((network != 0) * 1)

        # bulk_mat <- as.matrix(bulk_dataset[common, ])
        bulk_mat <- Matrix::Matrix(as.matrix(bulk_dataset[common, ]))

        # sc_mat <- as.matrix(sc_exprs[common, ])
        sc_mat <- sc_exprs[common, ] # A dgCMatrix object

        dataset0 <- Matrix::cbind2(bulk_mat, sc_mat) # much smaller
        if (verbose) {
            ts_cli$cli_alert_info(
                "Normalizing quantiles of data"
            )
        }

        dataset1 <- normalize.quantiles(as.matrix(dataset0), keep.names = TRUE)
        # rownames(dataset1) <- common
        # colnames(dataset1) <- c(colnames(bulk_mat), colnames(sc_mat))
        if (verbose) {
            ts_cli$cli_alert_info(
                "Subsetting data"
            )
        }

        n_bulk <- ncol(bulk_mat)
        # gene-sample
        Expression_bulk <- dataset1[, seq_len(n_bulk), drop = FALSE]
        # gene-cell
        Expression_cell <- dataset1[, (n_bulk + 1):ncol(dataset1), drop = FALSE]

        gc(verbose = FALSE)
        if (verbose) {
            ts_cli$cli_alert_info(
                "Calculating correlation"
            )
        }

        X <- stats::cor(Expression_bulk, Expression_cell)

        quality_check <- colQuantiles(X, probs = seq(0, 1, 0.25))
        if (verbose) {
            cli::cli_text(
                strrep("-", floor(getOption("width") / 2)),
                "\n",
                sep = ""
            )
            cli::cli_text("Five-number summary of correlations:")
            quality_check %>%
                asplit(2) %>%
                purrr::map_dbl(mean) %>%
                round(digits = 6) %>%
                paste(sep = " ", collapse = " ") %>%
                cli::cli_text()
            cli::cli_text(
                strrep("-", floor(getOption("width") / 2)),
                "\n",
                sep = ""
            )
        }
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
                if (verbose) {
                    cli::cli_alert_info(
                        "Current phenotype contains {.val {z[1]}} {tag[1]} and {.val {z[2]}} {tag[2]} samples."
                    )
                    ts_cli$cli_alert_info(
                        "Perform logistic regression on the given phenotypes..."
                    )
                }
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
                if (verbose) {
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
                }
                Y
            },
            cox = function() {
                Y <- as.matrix(phenotype)
                if (ncol(Y) != 2) {
                    cli::cli_abort(
                        "x" = "The size of survival data is wrong. Please check Scissor inputs and selected regression type."
                    )
                }
                if (verbose) {
                    ts_cli$cli_alert_info(
                        "Perform cox regression on the given clinical outcomes..."
                    )
                }
                Y
            }
        )

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
            if (verbose) {
                ts_cli$cli_alert_success(
                    "Statistics data saved to {.file {Save_file}}."
                )
            }
        }
    } else {
        # Load data from previous work
        if (verbose) {
            ts_cli$cli_alert_info(
                "Loading data from {.file {Load_file}}"
            )
        }
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
    gc(verbose = FALSE)
    if (verbose) {
        ts_cli$cli_alert_info("Screening...")
    }

    alpha <- alpha %||%
        c(0.005, 0.01, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9)

    results <- list()

    for (i in seq_along(alpha)) {
        set.seed(seed)

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
        cells <- colnames(X)
        Cell1 <- cells[pos_mask]
        Cell2 <- cells[neg_mask]
        percentage <- (length(Cell1) + length(Cell2)) / length(cells)

        if (verbose) {
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
        }

        if (percentage < cutoff) {
            ts_cli$cli_alert_info(cli::col_green("Scissor Ended."))
            break
        }
    }

    list(
        para = list(
            alpha = alpha,
            lambda = fit0$lambda.min,
            family = family,
            Coefs = Coefs # for miscellaneous informationss
        ),
        Coefs = Coefs, # for cell evaluation
        Scissor_pos = Cell1,
        Scissor_neg = Cell2,
        X = X,
        Y = Y,
        network = network
    )
}

# #' @title Reliability Significance Test
# #' @description
# #' Performs the reliability significance test to determine whether the inferred phenotype-to-cell
# #' associations are reliable (statistical p-value less than 0.05) or are false positives.
# #'
# #' This statistical test can determine whether the inferred phenotype-to-cell associations are
# #' reliable (statistical p-value less than 0.05) or are false positives.
# #'
# #' The motivation for the reliability significance test is: if the chosen single-cell and bulk data are not suitable for
# #' the phenotype-to-cell associations, the correlations would be less informative and not well associated with the phenotype
# #' labels. Thus, the corresponding prediction performance would be poor and not be significantly distinguishable from the
# #' randomly permutated labels.
# #'
# #' The evaluation measurements used in the reliability significance test are the mean squared error (MSE) for linear regression,
# #' the area under the ROC curve (AUC) for classification, and the concordance index (c-index) for Cox regression.
# #'
# #' @param X Scissor calculated correlation matrix for each pair of cells and bulk samples.
# #' @param Y Scissor calculated response variable in the regression model.
# #' @param network Scissor calculated cell-cell similarity network.
# #' @param alpha Same parameter alpha used in Scissor.
# #' @param family Same parameter family used in Scissor.
# #' @param cell_num The number of the Scissor selected cells.
# #' @param n Permutation times.
# #' @param nfold The fold number in cross-validation.
# #'
# #' @return A list containing the following components:
# #'   \item{statistic}{The test statistic.}
# #'   \item{p}{The test p-value.}
# #'   \item{Measurement_test_real}{The evaluation measurement on each fold calculated using the real label.}
# #'   \item{Measurement_test_back}{A list with each component contains the evaluation measurements calculated using the permutated labels.}
# #'
# #' @keywords internal
# #' @family scissor
# #
# reliability.test2 <- function(
#     X,
#     Y,
#     network,
#     alpha,
#     family = c("gaussian", "binomial", "cox"),
#     cell_num,
#     n = 100,
#     nfold = 10,
#     verbose = TRUE
# ) {
#     library(progress)
#     # library(Matrix)
#     switch(
#         family,
#         'gaussian' = {
#             Scissor::test_lm(
#                 X,
#                 Y,
#                 network,
#                 alpha,
#                 cell_num,
#                 n,
#                 nfold,
#                 verbose = verbose
#             )
#         },
#         'binomial' = {
#             library(pROC)
#             Scissor::test_logit(X, Y, network, alpha, cell_num, n, nfold)
#         },
#         'cox' = {
#             library(survival)
#             Scissor::test_cox(X, Y, network, alpha, cell_num, n, nfold)
#         }
#     )
# }

# #' @keywords internal
# test_lm2 <- function(
#     X,
#     Y,
#     network,
#     alpha,
#     cell_num,
#     n = 100,
#     nfold = 10,
#     ...
# ) {
#     dots <- rlang::list2(...)
#     seed <- dots$seed %||% getFuncOption("seed")
#     verbose <- dots$verbose %||% getFuncOption("verbose")
#     parallel <- dots$parallel %||% getFuncOption("parallel")
#     workers <- dots$workers %||% getFuncOption("workers")
#     parallel_type <- dots$parallel.type %||% getFuncOption("parallel.type")

#     set.seed(seed)
#     m <- nrow(X)
#     index0 <- sample(cut(seq_len(m), breaks = nfold, labels = FALSE))

#     if (verbose) {
#         ts_cli$cli_alert_info(
#             "Performing {nfold}-fold cross-validation on X with true labels"
#         )
#     }
#     X <- Matrix::Matrix(X)
#     MSE_test_real <- numeric(nfold)

#     if (verbose) {
#         cli::cli_progress_bar("CV with true labels", total = nfold)
#     }
#     for (j in seq_len(nfold)) {
#         MSE_test_real[j] <- ComputeFold(
#             j,
#             X,
#             Y,
#             index0,
#             network,
#             alpha,
#             cell_num
#         )
#         if (verbose) {
#             cli::cli_progress_update()
#         }
#         Sys.sleep(1 / 100)
#     }
#     if (verbose) {
#         cli::cli_progress_done()
#     }

#     if (verbose) {
#         ts_cli$cli_alert_info(
#             "Perform cross-validation on X with permutated label"
#         )
#     }

#     PermutationValidate <- function(
#         i,
#         X,
#         Y,
#         index0,
#         network,
#         alpha,
#         cell_num,
#         nfold,
#         m,
#         seed
#     ) {
#         set.seed(i + seed)
#         Y2 <- Y[sample(m)]

#         mse_vec <- numeric(nfold)
#         names(mse_vec) <- paste0("nfold_", seq_len(nfold))
#         for (j in seq_len(nfold)) {
#             mse_vec[j] <- ComputeFold(
#                 j,
#                 X,
#                 Y2,
#                 index0,
#                 network,
#                 alpha,
#                 cell_num,
#                 seed
#             )
#         }

#         mse_vec
#     }

#     if (verbose) {
#         ts_cli$cli_alert_info("Permutation test")
#     }

#     MSE_test_back <- if (parallel) {
#         plan(parallel_type, workers = workers)
#         on.exit(plan("sequential"), add = TRUE)

#         future_map(
#             seq_len(n),
#             ~ PermutationValidate(
#                 .x, # i
#                 X,
#                 Y,
#                 index0,
#                 network,
#                 alpha,
#                 cell_num,
#                 nfold,
#                 m,
#                 seed
#             ),
#             .options = furrr_options(seed = TRUE),
#             .progress = verbose
#         )
#     } else {
#         purrr::map(
#             seq_len(n),
#             ~ PermutationValidate(
#                 .x, # i
#                 X,
#                 Y,
#                 index0,
#                 network,
#                 alpha,
#                 cell_num,
#                 nfold,
#                 m,
#                 seed
#             ),
#             .progress = verbose
#         )
#     }
#     names(MSE_test_back) <- paste0("n_", seq_len(n))
#     background <- purrr::map_dbl(MSE_test_back, mean)

#     statistic <- mean(background, na.rm = TRUE)
#     p <- sum(background < statistic, na.rm = TRUE) / sum(!is.na(background))

#     if (verbose) {
#         cli::cli_h3("Results")
#         cli::cli_alert_success(
#             "Test statistic (mean MSE) = {.val {round(statistic, 3)}}"
#         )
#         cli::cli_alert_success(
#             "Reliability test p-value = {.val {round(p, 3)}}"
#         )
#     }

#     if (p >= 0.05) {
#         cli::cli_warn(
#             "Model is NOT significantly better than random (p >= 0.05)"
#         )
#     }

#     list(
#         statistic = statistic,
#         p = p,
#         MSE_test_real = MSE_test_real,
#         MSE_test_back = MSE_test_back
#     )
# }

# #' @keywords internal
# ComputeFold <- function(j, X, Y, index0, network, alpha, cell_num, seed) {
#     c_index <- which(index0 == j)
#     X_train <- X[-c_index, , drop = FALSE]
#     Y_train <- Y[-c_index]
#     X_test <- X[c_index, , drop = FALSE]
#     Y_test <- Y[c_index]

#     # 拟合模型
#     fit <- NULL
#     while (is.null(fit$fit)) {
#         set.seed(seed)
#         fit <- Scissor::APML1(
#             X_train,
#             Y_train,
#             family = "gaussian",
#             penalty = "Net",
#             alpha = alpha,
#             Omega = network,
#             nlambda = 100
#         )
#     }

#     index <- which.min(abs(fit$fit$nzero - cell_num))
#     Coefs <- as.numeric(fit$Beta[, index])

#     mean((Y_test - X_test %*% Coefs)^2)
# }
