# * ---- Preprocess bulk expression data ----

#' @title Convert Ensembles Version IDs & TCGA Version IDs to Genes in Bulk Expression Data
#'
#' @description
#' Preprocess bulk expression data: convert Ensembles version IDs and TCGA version IDs to genes. NA values are replaced with `unknown_k` format (k stands for the position of the NA value in the row).
#'
#' @param data bulk expression data (matrix or data.frame)
#'
#' @export
#'
SymbolConvert <- function(data) {
    row_names <- gsub("\\..*$", "", rownames(data))
    if (is.null(row_names)) {
        cli::cli_abort(c("x" = "Row names are missing in the data"))
    }
    options(
        IDConverter.datapath = system.file("extdata", package = "IDConverter")
    )
    gene_symbols <- IDConverter::convert_hm_genes(row_names)

    na_count <- sum(is.na(gene_symbols))
    if (na_count > 0) {
        cli::cli_warn(c(
            "Found {.val {na_count}} NA values in gene symbols during conversion."
        ))

        na_indices <- which(is.na(gene_symbols))
        gene_symbols[na_indices] <- glue::glue("unknown_{na_indices}")
        cli::cli_warn(
            "Replaced {.val {na_count}} NA values with `unknown_k` format."
        )
    }

    rownames(data) <- gene_symbols
    return(data)
}


#' @title Bulk RNA-seq Data Preprocessing and Quality Control Function
#'
#' @description
#' This function performs comprehensive preprocessing and quality control analysis
#' for bulk RNA-seq data, including data validation, filtering, batch effect detection,
#' principal component analysis, and visualization.
#'
#' @details
#' The function performs the following operations:
#' \enumerate{
#'   \item Data validation and format conversion
#'   \item Basic statistics calculation (missing values, read depth, gene detection)
#'   \item Optional detailed quality checks including:
#'     \itemize{
#'       \item Sample correlation analysis
#'       \item Principal Component Analysis (PCA)
#'       \item Outlier detection using Mahalanobis distance
#'       \item Batch effect detection using ANOVA
#'     }
#'   \item Data filtering based on count thresholds and quality metrics
#'   \item Optional gene symbol conversion
#'   \item Visualization generation (PCA plots)
#' }
#'
#' @param data Expression matrix with genes as rows and samples as columns, or a list containing count_matrix and sample_info
#' @param sample_info Sample information data frame (optional), ignored if data is a list. A qualified `sample_info` should contain both `sample` and `condition` columns (case-sensitive), and there are no specific requirements for the data type stored in the `condition` column.
#' @param gene_symbol_conversion Whether to convert Ensembles version IDs and TCGA version IDs to genes with IDConverter, default: FALSE
#' @param check Whether to perform detailed quality checks, default: TRUE
#' @param min_count_threshold Minimum count threshold for gene filtering, default: 10
#' @param min_gene_expressed Minimum number of samples a gene must be expressed in, default: 3
#' @param min_total_reads Minimum total reads per sample, default: 1e6
#' @param min_genes_detected Minimum number of genes detected per sample, default: 10000
#' @param min_correlation Minimum correlation threshold between samples, default: 0.8
#' @param n_top_genes Number of top variable genes for PCA analysis, default: 500
#' @param show_plot_results Whether to generate visualization plots, default: TRUE
#' @param verbose Whether to output detailed information, default: TRUE
#'
#'
#' @section Quality Metrics:
#' The function calculates and reports several quality metrics:
#' \describe{
#'   \item{Data Integrity}{Number of missing values}
#'   \item{Gene Count}{Total number of genes after filtering}
#'   \item{Sample Read Depth}{Total reads per sample}
#'   \item{Gene Detection Rate}{Number of genes detected per sample}
#'   \item{Sample Correlation}{Pearson correlation between samples}
#'   \item{PCA Variance}{Variance explained by first two principal components}
#'   \item{Batch Effects}{Proportion of genes significantly affected by batch}
#' }
#'
#'
#' @section Sample Information Format:
#' The sample_info data frame should contain:
#' \describe{
#'   \item{sample}{Character vector of unique sample identifiers}
#'   \item{condition}{Character vector of experimental conditions}
#'   \item{batch}{Optional character vector of batch identifiers}
#' }
#'
#' @section Filtering Criteria:
#' Genes are retained if:
#' \itemize{
#'   \item They have counts >= min_count_threshold in >= min_gene_expressed samples
#' }
#'
#' Samples are retained if:
#' \itemize{
#'   \item Total reads >= min_total_reads
#'   \item Detected genes >= min_genes_detected
#'   \item Mean correlation >= min_correlation (if check=TRUE)
#' }
#'
#'
#' @seealso
#' \code{\link[edgeR]{cpm}} for counts per million calculation,
#' \code{\link[stats]{prcomp}} for PCA analysis,
#' \code{\link[stats]{cor}} for correlation analysis,
#' \code{\link[SigBridgeR]{SymbolConvert}} for gene symbol conversion
#'
#' @return Filtered count matrix
#' @export
#'
BulkPreProcess <- function(
    data,
    sample_info = NULL,
    gene_symbol_conversion = FALSE,
    check = TRUE,
    min_count_threshold = 10,
    min_gene_expressed = 3,
    min_total_reads = 1e6,
    min_genes_detected = 10000,
    min_correlation = 0.8,
    n_top_genes = 500,
    show_plot_results = TRUE,
    verbose = TRUE
) {
    purrr::walk(
        list(
            min_count_threshold,
            min_gene_expressed,
            min_genes_detected,
            min_total_reads,
            min_correlation,
            n_top_genes
        ),
        ~ chk::chk_numeric
    )
    purrr::walk(list(check, show_plot_results, verbose), ~ chk::chk_flag)

    if (verbose) {
        ts_cli$cli_alert_info(
            cli::col_green("Starting data preprocessing...")
        )
    }

    # Handle input data format
    if (is.list(data) && !is.data.frame(data)) {
        if ("count_matrix" %chin% names(data)) {
            counts_matrix <- as.matrix(data$count_matrix)
            if ("sample_info" %chin% names(data) && is.null(sample_info)) {
                sample_info <- data$sample_info
            }
        } else {
            cli::cli_abort(c(
                "x" = "If data is a list, it must contain 'count_matrix' element"
            ))
        }
    } else {
        counts_matrix <- as.matrix(data)
    }

    # Data validation
    if (!is.numeric(counts_matrix)) {
        cli::cli_warn(
            "Expression matrix must be numeric, converted to numeric now."
        )
        # Attempt to convert to numeric
        counts_matrix <- rlang::try_fetch(
            {
                apply(counts_matrix, 2, as.numeric)
            },
            error = function(e) {
                cli::cli_abort(c(
                    "x" = "Failed to convert expression matrix to numeric format:",
                    e$message
                ))
            },
            warning = function(w) {
                if (verbose) {
                    cli::cli_warn(
                        "Warnings during numeric conversion: {w$message}"
                    )
                }
            }
        )

        # Check if conversion introduced NAs
        new_nas <- sum(is.na(counts_matrix))
        if (new_nas > 0) {
            if (verbose) {
                cli::cli_warn(
                    "Numeric conversion introduced {.val {new_nas}} NA values (will be imputed)"
                )
            }
        }
    }

    if (any(counts_matrix < 0, na.rm = TRUE)) {
        cli::cli_abort(c(
            "x" = "Expression matrix cannot contain negative values"
        ))
    }

    n_genes <- nrow(counts_matrix)
    n_samples <- ncol(counts_matrix)

    if (is.null(sample_info)) {
        # Create default sample information
        if (verbose) {
            ts_cli$cli_alert_info(
                "No sample info provided, using default settings."
            )
        }
        sample_info <- data.frame(
            sample = colnames(counts_matrix) %||%
                glue::glue("Sample_{seq_len(n_samples)}"),
            condition = rep("unknown", n_samples),
            stringsAsFactors = FALSE
        )
    } else {
        # Validate matching
        if (nrow(sample_info) != n_samples) {
            cli::cli_abort(c(
                "x" = "Number of rows in `sample_info` does not match number of columns in count matrix"
            ))
        }
    }

    if (verbose) {
        ts_cli$cli_alert_success(
            "Data loaded: {.val {n_genes}} genes * {.val {n_samples}} samples"
        )
    }

    # Basic Statistics
    if (verbose) {
        ts_cli$cli_alert_info("Computing basic statistics...")
    }

    missing_data <- sum(is.na(counts_matrix))
    total_reads_per_sample <- colSums(counts_matrix, na.rm = TRUE)
    genes_detected_per_sample <- colSums(counts_matrix > 0, na.rm = TRUE)
    genes_with_reads <- rowSums(counts_matrix > 0, na.rm = TRUE)
    genes_expressed_multiple <- sum(genes_with_reads >= min_gene_expressed)

    if (verbose) {
        ts_cli$cli_alert_success(
            "{.val {genes_expressed_multiple}} genes pass expression filter"
        )
    }

    # Detailed Quality Checks (Conditional Execution)
    if (check) {
        if (verbose) {
            ts_cli$cli_alert_info(
                "Starting detailed quality checks..."
            )
        }

        # Calculate CPM (compute once, reuse)
        cpm_values <- edgeR::cpm(counts_matrix, log = TRUE, prior.count = 1)

        # Sample correlation
        if (verbose) {
            ts_cli$cli_alert_info(
                "Computing sample correlations..."
            )
        }
        sample_cor <- stats::cor(
            cpm_values,
            method = "pearson",
            use = "pairwise.complete.obs"
        )

        # Quick identification of low correlations
        cor_lower_tri <- sample_cor[lower.tri(sample_cor)]
        min_cor <- min(cor_lower_tri, na.rm = TRUE)
        mean_cor <- mean(cor_lower_tri, na.rm = TRUE)
        low_cor_count <- sum(cor_lower_tri < min_correlation, na.rm = TRUE)

        if (verbose) {
            ts_cli$cli_alert_success(
                "Correlation analysis completed: min correlation = {.val {round(min_cor, 3)}}"
            )
            if (low_cor_count > 0) {
                cli::cli_warn(
                    "Found {.val {low_cor_count}} sample pairs with correlation < {.val min_correlation = {min_correlation}}"
                )
            }
        }

        # PCA Analysis
        if (verbose) {
            ts_cli$cli_alert_info(
                "Performing principal component analysis..."
            )
        }

        # Select highly variable genes
        n_genes_for_pca <- min(n_top_genes, nrow(cpm_values))
        gene_vars <- rowVars(cpm_values, na.rm = TRUE)
        top_var_idx <- order(gene_vars, decreasing = TRUE)[seq_len(
            n_genes_for_pca
        )]

        # PCA computation
        pca_data <- Matrix::t(cpm_values[top_var_idx, , drop = FALSE])
        pca_result <- stats::prcomp(pca_data, scale. = TRUE, center = TRUE)
        pca_summary <- summary(pca_result)

        # Outlier detection (use only first 2 PCs)
        if (ncol(pca_result$x) >= 2) {
            pc12 <- pca_result$x[, c(1, 2)]
            mahal_dist <- stats::mahalanobis(
                pc12,
                center = colMeans(pc12),
                cov = stats::cov(pc12)
            )
            outlier_threshold <- stats::qchisq(0.95, df = 2)
            outliers <- which(mahal_dist > outlier_threshold)
        } else {
            outliers <- integer(0)
        }

        if (verbose) {
            var_pc1 <- round(pca_summary$importance[2, 1] * 100, 2)
            var_pc2 <- round(pca_summary$importance[2, 2] * 100, 2)
            ts_cli$cli_alert_info(
                "PCA completed: PC1({.val {var_pc1}}%) PC2({.val {var_pc2}}%), {.val {length(outliers)}} outlier samples"
            )
            if (length(outliers) > 0) {
                outlier_names <- sample_info$sample[outliers]
                cli::cli_warn(sprintf(
                    "Outlier samples detected: {.emph {%s}}",
                    glue::glue_collapse(outlier_names, sep = ', ')
                ))
            }
        }

        # Batch effect detection (ANOVA)
        if ("batch" %chin% colnames(sample_info)) {
            batch_levels <- unique(sample_info$batch)
            if (length(batch_levels) > 1) {
                if (verbose) {
                    ts_cli$cli_alert_info(
                        "Detecting batch effects..."
                    )
                }

                # Fast batch effect detection (random sampling of genes)
                n_genes_batch <- min(5000, nrow(cpm_values))
                batch_test_idx <- sample(nrow(cpm_values), n_genes_batch)

                perform_anova_test <- function(gene_expression, batch_info) {
                    if (all(is.na(gene_expression))) {
                        return(1)
                    }
                    rlang::try_fetch(
                        {
                            aov_result <- stats::aov(
                                gene_expression ~ batch_info
                            )
                            summary(aov_result)[[1]]["Pr(>F)"][1, 1]
                        },
                        error = function(e) 1
                    )
                }

                batch_pvalues <- cpm_values[batch_test_idx, , drop = FALSE] %>%
                    asplit(1) %>% # by row
                    purrr::map_dbl(~ perform_anova_test(.x, sample_info$batch))

                batch_sig_prop <- mean(batch_pvalues < 0.05, na.rm = TRUE)

                if (verbose) {
                    ts_cli$cli_alert_success(sprintf(
                        "Batch effect detection completed: {.val {%s}}% genes affected",
                        round(batch_sig_prop * 100, 2)
                    ))
                    if (batch_sig_prop > 0.3) {
                        cli::cli_warn(
                            "Strong batch effects detected, consider batch correction."
                        )
                    }
                }
            }
        }

        # Generate key plots
        if (show_plot_results) {
            if (verbose) {
                ts_cli$cli_alert_info(
                    "Generating visualization plots..."
                )
            }
            # PCA graphics::plot
            if (exists("pca_result")) {
                pca_df <- data.frame(
                    PC1 = pca_result$x[, 1],
                    PC2 = pca_result$x[, 2],
                    sample = sample_info$sample,
                    condition = sample_info$condition,
                    stringsAsFactors = FALSE
                )

                if ("batch" %chin% colnames(sample_info)) {
                    pca_df$batch <- sample_info$batch
                }

                var_labels <- pca_summary$importance[2, c(1, 2)] * 100

                p_pca <- ggplot2::ggplot(
                    pca_df,
                    ggplot2::aes(x = `PC1`, y = `PC2`, color = `condition`)
                ) +
                    ggplot2::geom_point(size = 3, alpha = 0.8) +
                    ggplot2::labs(
                        title = "Principal Component Analysis (PCA)",
                        x = paste0("PC1 (", round(var_labels[1], 2), "%)"),
                        y = paste0("PC2 (", round(var_labels[2], 2), "%)")
                    ) +
                    ggplot2::theme_minimal() +
                    ggplot2::theme(legend.position = "right") +
                    ggforce::geom_mark_ellipse(
                        ggplot2::aes(fill = condition, group = condition),
                        alpha = 0.1,
                        expand = ggplot2::unit(3, "mm"),
                        show.legend = FALSE
                    )

                if ("batch" %chin% colnames(pca_df)) {
                    p_pca <- p_pca +
                        ggplot2::aes(shape = `batch`) +
                        ggplot2::scale_shape_manual(
                            values = c(
                                16,
                                17,
                                18,
                                19,
                                20
                            )[seq_len(length(unique(pca_df$batch)))]
                        )
                }

                methods::show(p_pca)
            }

            if (verbose) {
                ts_cli$cli_alert_success(
                    "Plot generation completed"
                )
            }
        }

        rm(cpm_values, pca_data, pca_result)
        if (exists("pc12")) {
            rm(pc12)
        }
        gc(verbose = FALSE)
    }

    # Data Filtering
    if (verbose) {
        ts_cli$cli_alert_info("Performing data filtering...")
    }

    # Gene filtering
    genes_pass_count <- rowSums(counts_matrix >= min_count_threshold) >=
        min_gene_expressed

    # Sample filtering
    samples_pass_reads <- total_reads_per_sample >= min_total_reads
    samples_pass_genes <- genes_detected_per_sample >= min_genes_detected
    samples_to_keep <- samples_pass_reads & samples_pass_genes

    # Correlation filtering (if correlation was computed)
    if (check && exists("sample_cor")) {
        sample_mean_cor <- colMeans(sample_cor, na.rm = TRUE)
        samples_pass_cor <- sample_mean_cor >= min_correlation
        samples_to_keep <- samples_to_keep & samples_pass_cor
    }

    filtered_counts <- counts_matrix[
        genes_pass_count,
        samples_to_keep,
        drop = FALSE
    ]

    n_genes_filtered <- sum(genes_pass_count)
    n_samples_filtered <- sum(samples_to_keep)

    if (verbose) {
        genes_removed <- n_genes - n_genes_filtered
        samples_removed <- n_samples - n_samples_filtered
        ts_cli$cli_alert_success("Data filtering completed:")
        cli::cli_alert_info(
            "  Genes: {.val {n_genes}} -> {.val {n_genes_filtered}} (removed {.val {genes_removed}})"
        )
        cli::cli_alert_info(
            "  Samples: {.val {n_samples}} -> {.val {n_samples_filtered}} (removed {.val {samples_removed}})"
        )
    }
    # Quality Report Generation
    if (verbose) {
        cli::cli_h2("Quality Assessment Summary:")

        # Data integrity
        if (missing_data == 0) {
            cli::cli_alert_success("Data integrity: No missing values")
        } else {
            cli::cli_warn(
                "Data integrity: {.val {missing_data}} missing values detected"
            )
        }

        # Sample quality
        if (all(samples_pass_reads[samples_to_keep])) {
            cli::cli_alert_success(
                "Sample read depth: All samples pass threshold (>={.val {min_total_reads}})"
            )
        } else {
            n_low_reads <- sum(!samples_pass_reads[samples_to_keep])
            cli::cli_warn(
                "Sample read depth: {.val {n_low_reads}} samples below threshold"
            )
        }

        if (all(samples_pass_genes[samples_to_keep])) {
            cli::cli_alert_success(
                "Gene detection: All samples pass threshold (>={.val {min_genes_detected}})"
            )
        } else {
            n_low_genes <- sum(!samples_pass_genes[samples_to_keep])
            cli::cli_warn(
                "Gene detection: {.val {n_low_genes}} samples below threshold"
            )
        }

        # Correlation and outliers
        if (check) {
            if (exists("min_cor") && min_cor >= min_correlation) {
                cli::cli_alert_success(
                    "Sample correlation: Good (minimum = {.val {round(min_cor, 3)}})"
                )
            } else if (exists("min_cor")) {
                cli::cli_warn(
                    "Sample correlation: Low (minimum = {.val {round(min_cor, 3)}})"
                )
            }

            if (exists("outliers") && length(outliers) == 0) {
                cli::cli_alert_success("Sample outliers: None detected")
            } else if (exists("outliers") && length(outliers) > 0) {
                cli::cli_warn(
                    "Sample outliers: {.val {length(outliers)}} sample(s) detected"
                )
            }
        }
    }

    if (gene_symbol_conversion) {
        filtered_counts <- SymbolConvert(filtered_counts)
        if (verbose) {
            ts_cli$cli_alert_success("Gene symbol conversion done")
        }
    }

    if (verbose) {
        ts_cli$cli_alert_success(
            cli::col_green("BulkPreProcess completed")
        )
    }

    filtered_counts
}
