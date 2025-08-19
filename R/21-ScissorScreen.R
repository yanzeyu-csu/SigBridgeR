# ---- 2. DO SCISSOR ----

#' @title Perform Scissor Screening Analysis
#' @description
#' Identifies phenotype-associated cell subpopulations in single-cell data using
#' regularized regression on matched bulk expression profiles.
#'
#' @usage
#' DoScissor(
#'    path2load_scissor_cache = NULL,
#'    matched_bulk,
#'    sc_data,
#'    phenotype,
#'    label_type = "scissor",
#'    scissor_alpha = 0.05,
#'    scissor_cutoff = 0.2,
#'    scissor_family = c("gaussian", "binomial", "cox"),
#'    reliability_test = TRUE,
#'    path2save_scissor_inputs = "Scissor_inputs.RData",
#'    reliability_test_n = 10,
#'    nfold = 10,
#'    ...
#' )
#'
#' @param path2load_scissor_cache Path to precomputed Scissor inputs (RData file).
#'        If provided, skips recomputation (default: NULL).
#' @param matched_bulk Normalized bulk expression matrix (features × samples).
#'        Column names must match `phenotype` identifiers.
#' @param sc_data Seurat object containing single-cell RNA-seq data.
#' @param phenotype Clinical outcome data. Can be:
#'        - Vector: named with sample IDs
#'        - Data frame: with row names matching bulk columns
#' @param label_type Character specifying phenotype label type (e.g., "SBS1", "time"), stored in `scRNA_data@misc`
#' @param scissor_alpha Parameter used to balance the effect of the l1 norm and the network-based penalties. It can be a number or a searching vector. If alpha = NULL, a default searching vector is used. The range of alpha is between 0 and 1. A larger alpha lays more emphasis on the l1 norm.(default: 0.05).
#' @param scissor_cutoff  (default: 0.2).
#'        Higher values increase specificity.
#' @param scissor_family Model family for outcome type:
#'        - "gaussian": Continuous outcomes
#'        - "binomial": Binary outcomes (default)
#'        - "cox": Survival outcomes
#' @param reliability_test Logical to perform stability assessment when `scissor_alpha` is specified as a value between 0 and 1.(default: TRUE).
#' @param reliability_test_n Number of CV folds for reliability test (default: 10).
#'
#' @param path2save_scissor_inputs Path to save intermediate files (default: "Scissor_inputs.RData").
#' @param nfold Cross-validation folds for reliability test (default: 10).
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
#' }
#'
#' @section Reference:
#' Sun S et al. (2022). "Scissor identifies phenotype-associated cell subsets
#' in single-cell genomics." Nat Methods 19(5):600-608. \doi{10.1038/s41592-022-01452-1}
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
#' @importFrom cli cli_abort cli_alert_info cli_alert_success cli_alert_danger
#' @importFrom crayon green red
#'
#' @keywords SigBridgeR_internal
#' @export
#'
DoScissor = function(
    path2load_scissor_cache = NULL,
    matched_bulk,
    sc_data,
    phenotype,
    label_type = "scissor",
    scissor_alpha = 0.05,
    scissor_cutoff = 0.2,
    scissor_family = c("gaussian", "binomial", "cox"),
    reliability_test = TRUE,
    path2save_scissor_inputs = "Scissor_inputs.RData",
    reliability_test_n = 10,
    nfold = 10,
    ...
) {
    # Input validation
    scissor_family <- match.arg(scissor_family)

    if (scissor_family %in% c("binomial", "cox")) {
        label_type = c(
            glue::glue("{label_type}_Negative"),
            glue::glue("{label_type}_Positive")
        )
    } else if (scissor_family == "gaussian") {
        n = length(table(phenotype))
        label_type = glue::glue("{label_type}_{seq_len(n)}")
    } else if (length(scissor_family) != 1) {
        cli::cli_abort(c(
            "x" = "Please choose one scissor family, use parameter {.var scissor_family}.",
            class = "FamilyError"
        ))
    }
    path = dirname(path2save_scissor_inputs)
    if (!dir.exists(path)) {
        dir.create(path, recursive = TRUE)
    }

    infos1 <- Scissor.v5.optimized(
        bulk_dataset = matched_bulk,
        sc_dataset = sc_data,
        phenotype = phenotype,
        tag = label_type,
        alpha = scissor_alpha,
        cutoff = scissor_cutoff,
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
    sc_meta$scissor[rownames(sc_meta) %in% infos1$Scissor_pos] <- "Positive"
    sc_meta$scissor[rownames(sc_meta) %in% infos1$Scissor_neg] <- "Negative"
    sc_data <- Seurat::AddMetaData(sc_data, metadata = sc_meta) %>%
        AddMisc(scissor_type = label_type, cover = FALSE)

    # reliability test
    if (reliability_test) {
        ifelse(
            # indicate that Y has only two levels, both Pos and Neg cells exist
            !length(table(infos1$Y)) < 2,
            {
                cli::cli_alert_info(c(
                    "[{TimeStamp()}]",
                    crayon::green(" Start reliability test")
                ))
                reliability_result <- Scissor::reliability.test(
                    infos1$X,
                    infos1$Y,
                    infos1$network,
                    alpha = scissor_alpha,
                    family = scissor_family,
                    cell_num = length(infos1$Scissor_pos) +
                        length(infos1$Scissor_neg),
                    n = reliability_test_n,
                    nfold = nfold
                )
                cli::cli_alert_success(
                    "[{TimeStamp()}] reliability test: Done"
                )
            },
            {
                cli::cli_alert_danger(c(
                    "{crayon::red('Error in reliability test')}: one of the Pos or Neg cells doesn't exist"
                ))
                reliability_result <- NULL
            }
        )
    } else {
        reliability_result <- NULL
    }

    return(list(
        scRNA_data = sc_data,
        scissor_result = infos1,
        reliability_result = reliability_result
    ))
}

#' @title Optimized Scissor Algorithm for Seurat ver5
#' @description
#' Scissor.v5 from `https://doi.org/10.1038/s41587-021-01091-3`and `https://github.com/sunduanchen/Scissor/issues/59`
#' Another version of Scissor.v5() to optimize memory usage and execution speed in preprocess.
#'
#' @family screen method
#'
#' @importFrom Matrix Matrix
#'
#' @keywords SigBridgeR_internal
#' @noRd
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
    workers = 32,
    ...
) {
    # Scissor needs Matrix()

    cl <- parallel::makeCluster(min(workers, parallel::detectCores() - 1))
    doParallel::registerDoParallel(cl)

    cli::cli_alert_info(
        c("[{TimeStamp()}]", crayon::green(" Scissor start..."))
    )

    ifelse(
        is.null(Load_file),
        {
            cli::cli_alert_info(c(
                "[{TimeStamp()}]",
                crayon::bold(" Start from raw data...")
            ))
            common = intersect(
                rownames(bulk_dataset),
                rownames(sc_dataset)
            )
            if (length(common) == 0) {
                stop(
                    "There is",
                    crayon::bold(" no common genes "),
                    "between the given single-cell and bulk samples. Please check Scissor inputs."
                )
            }

            if (inherits(sc_dataset, "Seurat")) {
                sc_exprs <- as.matrix(sc_dataset@assays$RNA$data)
                if ("RNA_snn" %in% names(sc_dataset@graphs)) {
                    network <- as.matrix(sc_dataset@graphs$RNA_snn)
                    cli::cli_alert_info(
                        "Using {.val RNA_snn} graph for network."
                    )
                } else if ("integrated_snn" %in% names(sc_dataset@graphs)) {
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
                Seurat_tmp <- Seurat::CreateSeuratObject(
                    counts = sc_dataset,
                    verbose = FALSE
                ) %>%
                    Seurat::FindVariableFeatures(
                        selection.method = "vst",
                        nfeatures = 2000,
                        verbose = FALSE
                    ) %>%
                    Seurat::ScaleData(verbose = FALSE) %>%
                    Seurat::RunPCA(
                        features = Seurat::VariableFeatures(.),
                        verbose = FALSE
                    ) %>%
                    Seurat::FindNeighbors(dims = 1:10, verbose = FALSE)
                network <- as.matrix(Seurat_tmp@graphs$RNA_snn)
            }
            diag(network) <- 0
            network[which(network != 0)] <- 1
            dataset0 <- cbind(bulk_dataset[common, ], sc_exprs[common, ])

            cli::cli_alert_info(
                "[{TimeStamp()}] Normalizing quantiles of data..."
            )

            dataset1 <- preprocessCore::normalize.quantiles(as.matrix(dataset0))
            rownames(dataset1) <- rownames(dataset0)
            colnames(dataset1) <- colnames(dataset0)

            cli::cli_alert_info(
                "[{TimeStamp()}] Subsetting data..."
            )
            # gene-sample
            Expression_bulk <- dataset1[, 1:ncol(bulk_dataset)]
            # gene-cell
            Expression_cell <- dataset1[,
                (ncol(bulk_dataset) + 1):ncol(dataset1)
            ]
            gc(verbose = FALSE)

            cli::cli_alert_info(
                "[{TimeStamp()}] Calculating correlation..."
            )
            X <- cor(Expression_bulk, Expression_cell)
            quality_check <- stats::quantile(X)

            cat(strrep("-", floor(getOption("width") / 2)), "\n", sep = "")
            message(crayon::bold("Five-number summary of correlations:\n"))
            print(quality_check)
            cat(strrep("-", floor(getOption("width") / 2)), "\n", sep = "")
            # median
            if (quality_check[3] < 0.01) {
                cli::cli_alert_warning(crayon::yellow(
                    "The median correlation between the single-cell and bulk samples is relatively low."
                ))
            }

            switch(
                family,
                "binomial" = {
                    Y <- as.numeric(phenotype)
                    z <- table(Y)
                    if (length(z) != length(tag)) {
                        stop(
                            "The length differs between tags and phenotypes. Please check Scissor inputs and selected regression type.",
                            .call = FALSE
                        )
                    } else {
                        cli::cli_alert_info(
                            "Current phenotype contains {crayon::bold(z[1])} {tag[1]} and {crayon::bold(z[2])} {tag[2]} samples."
                        )
                        cli::cli_alert_info(
                            "[{TimeStamp()}] Perform logistic regression on the given phenotypes(It may take long):"
                        )
                    }
                },
                "gaussian" = {
                    Y <- as.numeric(phenotype)
                    z <- table(Y)
                    if (length(z) != length(tag)) {
                        cli::cli_abort(c(
                            "x" = "The length differs between tags and phenotypes. Please check Scissor inputs and selected regression type.",
                            "i" = "length of tags: {.val {length(tag)}}",
                            "i" = "length of phenotypes: {.val {length(z)}}"
                        ))
                    } else {
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
                        cli::cli_alert_info(
                            "[{TimeStamp()}] Perform linear regression on the given phenotypes:"
                        )
                    }
                },
                "cox" = {
                    Y <- as.matrix(phenotype)
                    if (ncol(Y) != 2) {
                        stop(
                            "The size of survival data is wrong. Please check Scissor inputs and selected regression type.",
                            .call = FALSE
                        )
                    } else {
                        cli::cli_alert_info(
                            "[{TimeStamp()}] Perform cox regression on the given clinical outcomes:"
                        )
                    }
                }
            )
            if (!is.null(Save_file)) {
                save(
                    X,
                    Y,
                    network,
                    Expression_bulk,
                    Expression_cell,
                    file = Save_file
                )
                cli::cli_alert_success(
                    "Statistics data saved to `{Save_file}`."
                )
            }
        },
        {
            # Load data from previous work
            cli::cli_alert_info(c(
                "[{TimeStamp()}]",
                crayon::bold(" Loading data from `{Load_file}`...")
            ))
            load(Load_file)
        }
    )
    # garbage collection
    rm(
        Expression_bulk,
        Expression_cell,
        sc_dataset,
        bulk_dataset,
        phenotype
    )

    cli::cli_alert_info(c(
        "[{TimeStamp()}]",
        crayon::bold(" Screening...")
    ))

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
                Cell1 <- colnames(X)[which(Coefs > 0)]
                Cell2 <- colnames(X)[which(Coefs < 0)]
                percentage <- (length(Cell1) + length(Cell2)) / ncol(X)
                print(c(sprintf("alpha = %s", alpha[i])))
                print(sprintf(
                    "Scissor identified %d Scissor+ cells and %d Scissor- cells.",
                    length(Cell1),
                    length(Cell2)
                ))
                print(sprintf(
                    "The percentage of selected cell is: %s%%",
                    formatC(percentage * 100, format = "f", digits = 3)
                ))
                if (percentage < cutoff) {
                    cli::cli_alert_info(
                        c("[{TimeStamp()}]", crayon::green(" Scissor Ended."))
                    )
                    break
                }
                cat("\n")
            },
            error = function(e) {
                cli::cli_alert_danger(
                    "[{TimeStamp()}]",
                    crayon::red("Error at alpha={alpha}:"),
                    e$message
                )
                cli::cli_alert_danger(conditionMessage(e))
            }
        )
    }
    cat(strrep("-", getOption("width")), "\n", sep = "")

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
