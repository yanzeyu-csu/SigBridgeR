#' @title Perform LP-SGL Screening Analysis
#' @description
#' Identifies phenotype-associated cell subpopulations using Lasso-Penalized
#' Sparse Group Lasso with Leiden community detection.
#'
#' @param matched_bulk Bulk expression matrix (features × samples)
#' @param sc_data Single-cell RNA-seq data (matrix or Seurat object)
#' @param phenotype Binary phenotype vector for bulk samples
#' @param label_type Character specifying phenotype label type (default: "LP_SGL")
#' @param resolution Resolution parameter for Leiden clustering (default: 0.6)
#' @param alpha Alpha parameter for SGL (default: 0.5)
#' @param nfold Number of folds for cross-validation (default: 5)
#' @param logFC_threshold Log fold change threshold for DEG filtering (default: 1)
#' @param pval_threshold P-value threshold for DEG filtering (default: 0.05)
#' @param ... Additional arguments passed to preprocessing functions
#'
#' @return A list containing screened results
#'
#' @export
DoLP_SGL <- function(
    matched_bulk,
    sc_data,
    phenotype,
    label_type = "LP_SGL",
    resolution = 0.6,
    alpha = 0.5,
    nfold = 5,
    logFC_threshold = 1,
    pval_threshold = 0.05,
    ...
) {
    # Input validation
    chk::chk_is(matched_bulk, c("matrix", "data.frame"))
    chk::chk_is(sc_data, c("matrix", "data.frame", "Seurat"))

    if (length(unique(phenotype)) != 2) {
        cli::cli_abort("Phenotype must contain exactly two groups")
    }

    ts_cli$cli_alert_info("Starting LP-SGL screening analysis")

    # Call core algorithm
    lpsgl_results <- LP_SGL.optimized(
        bulk_dataset = matched_bulk,
        sc_dataset = sc_data,
        phenotype = phenotype,
        resolution = resolution,
        alpha = alpha,
        nfold = nfold,
        logFC_threshold = logFC_threshold,
        pval_threshold = pval_threshold,
        ...
    )

    # Add results to Seurat object
    lpsgl_meta <- data.frame(
        LP_SGL = rep("Neutral", ncol(lpsgl_results$sc_data)),
        LP_SGL_coefficient = lpsgl_results$coefficients,
        row.names = colnames(lpsgl_results$sc_data)
    )
    lpsgl_meta$LP_SGL[lpsgl_results$coefficients > 0] <- "Positive"
    lpsgl_meta$LP_SGL[lpsgl_results$coefficients < 0] <- "Negative"

    sc_data_final <- Seurat::AddMetaData(
        lpsgl_results$sc_data,
        metadata = lpsgl_meta
    ) %>%
        AddMisc(LP_SGL_type = label_type, cover = FALSE)

    ts_cli$cli_alert_success("LP-SGL screening completed")

    return(list(
        scRNA_data = sc_data_final,
        lpsgl_result = lpsgl_results$model,
        coefficients = lpsgl_results$coefficients
    ))
}

#' @title Optimized LP-SGL Algorithm
#' @keywords internal
LP_SGL.optimized <- function(
    bulk_dataset,
    sc_dataset,
    phenotype,
    resolution = 0.6,
    alpha = 0.5,
    nfold = 5,
    logFC_threshold = 1,
    pval_threshold = 0.05,
    ...
) {
    # Preprocess bulk data
    processed_bulk <- BulkPreProcess(
        data = bulk_dataset,
        sample_info = data.frame(
            sample = colnames(bulk_dataset),
            condition = phenotype,
            stringsAsFactors = FALSE
        ),
        check = FALSE,
        gene_symbol_conversion = FALSE,
        verbose = FALSE,
        ...
    )

    # Update phenotype
    valid_samples <- colnames(processed_bulk)
    phenotype <- phenotype[names(phenotype) %in% valid_samples]

    # Preprocess single-cell data
    if (!inherits(sc_dataset, "Seurat")) {
        sc_processed <- SCPreProcess(sc = sc_dataset, verbose = FALSE, ...)
    } else {
        sc_processed <- sc_dataset
        if (!"RNA_snn" %chin% names(sc_processed@graphs)) {
            sc_processed <- SCPreProcess(
                sc = sc_processed,
                verbose = FALSE,
                ...
            )
        }
    }

    # DEG analysis
    phenotype1 <- which(phenotype == unique(phenotype)[1])
    phenotype2 <- which(phenotype == unique(phenotype)[2])
    group <- c(rep(c("a", "b"), c(length(phenotype1), length(phenotype2))))
    design <- model.matrix(~ 0 + factor(group))
    colnames(design) <- levels(factor(group))
    rownames(design) <- colnames(processed_bulk)
    contrast.matrix <- limma::makeContrasts("a-b", levels = design)
    fit <- limma::lmFit(processed_bulk, design)
    fit2 <- limma::contrasts.fit(fit, contrast.matrix)
    fit2 <- limma::eBayes(fit2)
    DEG_ot <- limma::topTable(fit2, coef = 1, n = Inf)
    DEGs <- subset(
        DEG_ot,
        abs(DEG_ot$logFC) > logFC_threshold & DEG_ot$P.Value < pval_threshold
    )

    # Correlation matrix
    shared_genes <- intersect(rownames(processed_bulk), rownames(sc_processed))
    sc_exprs <- as.matrix(SeuratObject::LayerData(sc_processed))
    correlation_matrix <- cor(
        processed_bulk[shared_genes, ],
        sc_exprs[shared_genes, ]
    )

    # Leiden clustering
    network <- as.matrix(sc_processed@graphs$RNA_snn)
    snn <- Matrix::as(network, "dgCMatrix")
    snn_df <- as.data.frame(Matrix::summary(snn))
    colnames(snn_df) <- c("from", "to", "weight")
    g <- igraph::graph_from_data_frame(snn_df, directed = FALSE)
    set.seed(0)
    leiden <- leidenAlg::leiden.community(g, resolution = resolution)
    membership <- leiden[["membership"]]
    index <- match(seq_len(ncol(sc_processed)), as.numeric(names(membership)))
    leidengroup <- as.numeric(membership[index])

    # SGL analysis
    data <- list(x = correlation_matrix, y = phenotype)
    fit_lpsgl <- SGL::SGL(data, leidengroup, type = "logit", alpha = alpha)
    lam <- fit_lpsgl[['lambdas']]
    cvfit <- SGL::cvSGL(
        data,
        leidengroup,
        type = 'logit',
        nfold = nfold,
        alpha = alpha,
        lambdas = lam
    )
    error <- cvfit$lldiff
    m <- min(error)
    h <- which(error == m)
    a <- fit_lpsgl[["beta"]]
    b <- a[, h]
    names(b) <- colnames(sc_processed)

    return(list(
        sc_data = sc_processed,
        coefficients = b,
        model = list(fit = fit_lpsgl, cvfit = cvfit),
        degs = DEGs
    ))
}
