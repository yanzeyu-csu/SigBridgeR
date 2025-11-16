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
        # The reason why using data.frame instead of vector is to
        # keep the same input and output format with scPP
        phenotype <- as.data.frame(phenotype) %>%
            Rownames2Col("Sample") %>%
            dplyr::rename("Feature" = 2) %>%
            dplyr::mutate(Feature = as.numeric(`Feature`))
    }
    gene_list <- if (tolower(phenotype_class) == "binary") {
        marker_Binary.optimized(
            bulk_data = matched_bulk,
            features = phenotype,
            ref_group = ref_group,
            Log2FC_cutoff = Log2FC_cutoff
        )
    } else if (tolower(phenotype_class) == "continuous") {
        marker_Continuous.optimized(
            bulk_data = matched_bulk,
            features = phenotype$Feature,
            method = "spearman",
            estimate_cutoff = estimate_cutoff,
        )
    } else {
        marker_Survival2(
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

#' @title Find Binary Feature Markers
#'
#' @param bulk_data Log2-normalized bulk expression data with genes in row and samples in column.
#' @param features Feature data of bulk samples, column1 are sample names (colname is "Sample") and column2 are feature (colname is "Feature") labels of each sample.
#' @param ref_group A character to indicate which feature is the control group.
#' @param Log2FC_cutoff Absolute cutoff value of fold change, default is 0.585.
#'
#' @return A gene list of feature markers with two binary groups and ranked genes based on t-test's statistic.
#' @keywords internal
#' @family scPP
#'
marker_Binary.optimized <- function(
    bulk_data,
    features,
    ref_group,
    Log2FC_cutoff = 0.585
) {
    # # 输入验证
    # if (missing(ref_group)) {
    #     cli::cli_abort("{.arg ref_group }is missing or incorrect.")
    # }
    # if (missing(bulk_data) || !inherits(bulk_data, c("matrix", "data.frame"))) {
    #     cli::cli_abort("{.arg bulk_data} is missing or incorrect.")
    # }
    # if (missing(features) || !inherits(features, c("matrix", "data.frame"))) {
    #     cli::cli_abort("{.arg features} is missing or incorrect.")
    # }

    # 确保bulk_data是矩阵格式
    if (!is.matrix(bulk_data)) {
        bulk_data <- as.matrix(bulk_data)
    }

    # 使用data.table进行高效的样本筛选
    features_dt <- data.table::as.data.table(features)
    ref <- features_dt[Feature == ref_group, Sample]
    tes <- features_dt[Feature != ref_group, Sample]

    # 向量化查找位置
    ref_pos <- which(colnames(bulk_data) %chin% ref)
    tes_pos <- which(colnames(bulk_data) %chin% tes)

    log2FCs <- rowMeans(bulk_data[, tes_pos, drop = FALSE]) -
        rowMeans(bulk_data[, ref_pos, drop = FALSE])

    # 使用matrixTests包进行批量t检验(比apply快很多)
    # 如果没有安装,回退到向量化的apply
    if (requireNamespace("matrixTests", quietly = TRUE)) {
        row_t_welch <- getExportedValue("matrixTests", "row_t_welch")
        t_results <- row_t_welch(
            bulk_data[, tes_pos, drop = FALSE],
            bulk_data[, ref_pos, drop = FALSE]
        )
        pvalues <- t_results$pvalue
        statistics <- t_results$statistic
        names(pvalues) <- rownames(bulk_data)
        names(statistics) <- rownames(bulk_data)
    } else {
        # 优化的apply方法:预分配结果,使用vapply
        n_genes <- nrow(bulk_data)
        pvalues <- numeric(n_genes)
        statistics <- numeric(n_genes)
        names(pvalues) <- rownames(bulk_data)
        names(statistics) <- rownames(bulk_data)

        for (i in seq_len(n_genes)) {
            test_result <- rlang::try_fetch(
                t.test(bulk_data[i, tes_pos], bulk_data[i, ref_pos]),
                error = function(e) list(p.value = NA, statistic = NA)
            )
            pvalues[i] <- test_result$p.value
            statistics[i] <- test_result$statistic
        }
    }

    # 排序统计量
    genes_sort <- sort(statistics[!is.na(statistics)], decreasing = TRUE)

    # 使用data.table构建结果
    res <- data.table::data.table(
        gene = names(pvalues),
        pvalue = pvalues,
        log2FC = log2FCs
    )

    # 计算FDR
    res[, fdr := stats::p.adjust(pvalue, method = "fdr")]

    # 高效筛选基因
    gene_pos <- res[pvalue < 0.05 & log2FC > Log2FC_cutoff, gene]
    gene_neg <- res[pvalue < 0.05 & log2FC < -Log2FC_cutoff, gene]

    # 构建返回列表
    geneList <- list(
        gene_pos = gene_pos,
        gene_neg = gene_neg,
        genes_sort = genes_sort
    )

    # 优化的条件判断
    has_pos <- length(gene_pos) > 0
    has_neg <- length(gene_neg) > 0

    if (has_pos && has_neg) {
        return(geneList)
    } else if (!has_pos) {
        cli::cli_warn(
            "There are no genes positively correlated with the given feature in this bulk dataset."
        )
        return(list(gene_neg = gene_neg))
    } else {
        cli::cli_warn(
            "There are no genes negatively correlated with the given feature in this bulk dataset."
        )
        return(list(gene_pos = gene_pos))
    }
}

#' @title Find Continuous Feature Markers
#'
#' @param bulk_data Log2-normalized bulk expression data with genes in row and samples in column.
#' @param features Feature data of bulk samples, such as TMB or CNA values of each sample.
#' @param estimate_cutoff Absolute cutoff value of correlation coefficient, default is 0.2.
#' @param method Method uses for cor.test, default is "spearman", another choice is "pearson".
#'
#' @return A gene list of feature markers and ranked genes based on correlation coefficients.
#' @keywords internal
#' @family scPP
#'
marker_Continuous.optimized <- function(
    bulk_data,
    features,
    method = "spearman",
    estimate_cutoff = 0.2
) {
    # # 输入验证 (上游已检查)
    # if (missing(bulk_data) || !inherits(bulk_data, c("matrix", "data.frame"))) {
    #     cli::cli_abort("{.arg bulk_data} is missing or incorrect.")
    # }
    # if (missing(features) || !is.numeric(features)) {
    #     cli::cli_abort("{.arg features} is missing or incorrect.")
    # }

    # 确保bulk_data是矩阵格式(性能更好)
    if (!is.matrix(bulk_data)) {
        bulk_data <- as.matrix(bulk_data)
    }

    # 预处理features: log2转换
    features_log <- log2(as.numeric(features) + 1)
    n_genes <- nrow(bulk_data)

    # 预分配结果向量
    pvalues <- numeric(n_genes)
    estimates <- numeric(n_genes)
    names(pvalues) <- rownames(bulk_data)
    names(estimates) <- rownames(bulk_data)

    # 向量化相关性检验
    for (i in seq_len(n_genes)) {
        cor_result <- rlang::try_fetch(
            stats::cor.test(
                bulk_data[i, ],
                features_log,
                method = method
            ),
            error = function(e) list(p.value = NA, estimate = NA)
        )
        pvalues[i] <- cor_result$p.value
        estimates[i] <- cor_result$estimate
    }

    # 排序估计值(去除NA)
    genes_sort <- sort(estimates[!is.na(estimates)], decreasing = TRUE)

    # 使用data.table构建结果
    res <- data.table::data.table(
        gene = names(pvalues),
        pvalue = pvalues,
        estimate = estimates
    )

    # 按pvalue排序并计算FDR
    data.table::setorder(res, pvalue)
    res[, fdr := stats::p.adjust(pvalue, method = "fdr")]

    # 高效筛选基因
    gene_pos <- res[fdr < 0.05 & estimate > estimate_cutoff, gene]
    gene_neg <- res[fdr < 0.05 & estimate < -estimate_cutoff, gene]

    # 构建返回列表
    geneList <- list(
        gene_pos = gene_pos,
        gene_neg = gene_neg,
        genes_sort = genes_sort
    )

    # 优化的条件判断
    has_pos <- length(gene_pos) > 0
    has_neg <- length(gene_neg) > 0

    if (has_pos && has_neg) {
        return(geneList)
    } else if (!has_pos) {
        cli::cli_warn(
            "There are no genes positively correlated with the given feature in this bulk dataset."
        )
        return(list(gene_neg = gene_neg))
    } else {
        cli::cli_warn(
            "There are no genes negatively correlated with the given feature in this bulk dataset."
        )
        return(list(gene_pos = gene_pos))
    }
}

#' @title Find Survival-Associated Markers
#'
#' @param bulk_data Log2-normalized bulk expression data with genes in row and samples in column.
#' @param survival_data Survival data with time in column1 and status in column2. Rownames are sample name.
#'
#' @return A gene list of survival-associated markers and ranked genes based on Cox hazard ratios.
#' @keywords internal
#' @family scPP
#'
marker_Survival2 <- function(bulk_data, survival_data) {
    # # 输入验证 (上游已检查)
    # if (missing(bulk_data) || !inherits(bulk_data, c("matrix", "data.frame"))) {
    #     cli::cli_abort("{.arg bulk_data} is missing or incorrect.")
    # }
    # if (missing(survival_data) || !inherits(survival_data, c("matrix", "data.frame"))) {
    #     cli::cli_abort("{.arg survival_data} is missing or incorrect.")
    # }

    SurvivalData <- data.frame(cbind(survival_data, Matrix::t(bulk_data)))
    colnames(SurvivalData) = make.names(colnames(SurvivalData))
    var <- make.names(rownames(bulk_data))

    Model_Formula <- sapply(var, function(x) {
        stats::as.formula(paste("survival::Surv(time, status) ~", x))
    })

    Model_all <- lapply(Model_Formula, function(x) {
        survival::coxph(x, data = SurvivalData)
    })

    res <- lapply(seq_along(Model_all), function(i) {
        coef_summary <- Matrix::summary(Model_all[[i]])$coefficients
        data.frame(
            variable = var[i],
            pvalue = coef_summary[, 5],
            coef = coef_summary[, 2]
        )
    }) %>%
        dplyr::bind_rows()

    genes_sort <- res %>%
        dplyr::arrange(dplyr::desc(coef)) %>%
        dplyr::pull(coef, name = variable)

    res <- res[order(res$pvalue), ]
    res$fdr <- stats::p.adjust(res$pvalue, method = "fdr")

    gene_pos <- res %>%
        dplyr::filter(fdr < 0.05, coef > 1) %>%
        dplyr::pull(variable) # correalted with worse survival
    gene_neg <- res %>%
        dplyr::filter(fdr < 0.05, coef < 1) %>%
        dplyr::pull(variable) # correlated with better survival

    geneList <- list(
        gene_pos = gene_pos,
        gene_neg = gene_neg,
        genes_sort = genes_sort
    )

    # ? It's confusing but it is
    if (length(gene_pos) > 0 & length(gene_neg) > 0) {
        return(geneList)
    } else if (length(gene_pos) == 0) {
        cli::cli_warn(
            "There are no genes negatively correlated with patients' prognosis in this bulk dataset."
        )
        return(list(gene_pos = gene_neg))
    } else if (length(gene_neg) == 0) {
        cli::cli_warn(
            "There are no genes positively correlated with patients' prognosis in this bulk dataset."
        )
        return(list(gene_neg = gene_pos))
    }
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
    verbose = getFuncOption("verbose"),
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
    verbose = getFuncOption("verbose"),
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
    verbose = getFuncOption("verbose"),
    parallel = getFuncOption("parallel"),
    workers = getFuncOption("workers"),
    seed = getFuncOption("seed")
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
        plan(getFuncOption("parallel.type"), workers = workers)

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
