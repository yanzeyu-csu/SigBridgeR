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
#' @param cutoff Numeric. Cutoff value for selecting the optimal alpha value when alpha = NULL. (default: 0.2)
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
#'
#' @family screen_method
#' @family scPAS
#'
DoscPAS <- function(
    matched_bulk,
    sc_data,
    phenotype,
    label_type = "scPAS",
    assay = 'RNA',
    imputation = FALSE,
    imputation_method = c("KNN", "ALRA"),
    nfeature = 3000L,
    alpha = c(0.01, NULL),
    cutoff = 0.2,
    network_class = c("SC", "bulk"),
    scPAS_family = c("cox", "gaussian", "binomial"),
    permutation_times = 2000L,
    FDR_threshold = 0.05,
    independent = TRUE,
    ...
) {
    chk::chk_is(matched_bulk, c("matrix", "data.frame"))
    chk::chk_is(sc_data, "Seurat")
    chk::chk_character(label_type)
    chk::chk_flag(imputation)
    if (imputation) {
        imputation_method %<>% MatchArg(c("KNN", "ALRA"))
    }
    if (!is.null(alpha)) {
        chk::chk_range(alpha)
    }
    network_class %<>% MatchArg(c("SC", "bulk"))
    scPAS_family %<>% MatchArg(c("cox", "gaussian", "binomial"), NULL)
    purrr::walk(
        list(nfeature, permutation_times, FDR_threshold),
        ~ chk::chk_number
    )
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

    ts_cli$cli_alert_info(cli::col_green("Start scPAS screening."))

    scPAS_result <- scPAS.optimized(
        bulk_dataset = as.matrix(matched_bulk),
        sc_dataset = sc_data,
        assay = 'RNA',
        tag = label_type,
        phenotype = phenotype,
        imputation = imputation,
        imputation_method = imputation_method,
        nfeature = nfeature,
        alpha = alpha,
        cutoff = cutoff,
        network_class = network_class,
        family = scPAS_family,
        independent = independent,
        permutation_times = permutation_times,
        FDR.threshold = FDR_threshold,
        ...
    ) %>%
        AddMisc(scPAS_type = label_type, cover = FALSE)

    detailed_info <- dplyr::select(
        scPAS_result[[]],
        dplyr::contains("scPAS_")
    )

    ts_cli$cli_alert_success(
        cli::col_green("scPAS screening done.")
    )

    return(list(
        scRNA_data = scPAS_result,
        stats = detailed_info,
        para = scPAS_result@misc$scPAS_para
    ))
}


#' @title scPAS: Single-Cell Phenotype-Associated Subpopulations (Optimized)
#'
#' @description
#' An optimized implementation of scPAS for identifying phenotype-associated
#' cell subpopulations from single-cell RNA-seq data by integrating bulk
#' transcriptomic data. This version includes performance optimizations,
#' memory-efficient matrix operations, and enhanced statistical testing.
#'
#' @param bulk_dataset A matrix or data frame containing bulk expression data.
#'   Each row represents a gene and each column represents a sample. Expression
#'   values should be continuous
#' @param sc_dataset A Seurat object or matrix containing single-cell RNA-seq
#'   expression data. If a matrix is provided, it will be automatically
#'   processed using Seurat's default pipeline.
#' @param phenotype Phenotype annotation for bulk samples. The format depends
#'   on the regression family:
#'   \itemize{
#'     \item For `family = "gaussian"`: A continuous numeric vector
#'     \item For `family = "binomial"`: A binary group indicator vector (0/1
#'           encoded) or factor with two levels
#'     \item For `family = "cox"`: A two-column matrix with columns named
#'           'time' and 'status' (1 = event, 0 = censored)
#'   }
#' @param assay Character string specifying the assay name in the Seurat object
#'   to use for analysis. Default: 'RNA'.
#' @param tag Optional character vector of length 2 specifying names for each
#'   phenotypic group. Used only for logistic regression (`family = "binomial"`).
#' @param nfeature Numeric value or character vector specifying the number of
#'   variable features to select, or a custom set of feature names. If `NULL`,
#'   all common genes between bulk and single-cell data are used.
#' @param imputation Logical indicating whether to perform imputation on
#'   single-cell data. Default: `TRUE`.
#' @param imputation_method Character string specifying the imputation method.
#'   One of: 'KNN', 'ALRA'. Default: 'KNN'.
#' @param alpha Numeric value or vector specifying the regularization parameter
#'   balancing L1 and network-based penalties. If `NULL`, a default sequence
#'   from 0.001 to 0.9 is used.
#' @param cutoff When `alpha = NULL`, the threshold for selecting the optimal
#'   alpha value. Default: 0.2
#' @param network_class Character string specifying the source for constructing
#'   the gene-gene similarity network. One of: 'SC' (single-cell data),
#'   'bulk' (bulk data). Default: 'SC'.
#' @param independent Logical indicating whether to construct background
#'   distributions independently for each cell. Default: `TRUE`.
#' @param family Character string specifying the regression family. One of:
#'   "gaussian" (linear regression), "binomial" (logistic regression),
#'   "cox" (Cox regression). Default: "gaussian".
#' @param permutation_times Numeric value specifying the number of permutations
#'   for statistical testing. Default: 2000.
#' @param FDR.threshold Numeric value specifying the false discovery rate
#'   threshold for identifying phenotype-associated cells. Default: 0.05.
#'
#' @return
#' Returns the input Seurat object with the following additions:
#' \itemize{
#'   \item **Metadata columns**:
#'     \itemize{
#'       \item `scPAS_RS` - Raw risk scores for each cell
#'       \item `scPAS_NRS` - Normalized risk scores (Z-statistics)
#'       \item `scPAS_Pvalue` - P-values from permutation testing
#'       \item `scPAS_FDR` - False discovery rate adjusted p-values
#'       \item `scPAS` - Cell classification labels: "Positive", "Negative", or "Neutral"
#'     }
#'   \item **Miscellaneous slot** (`sc_dataset@misc$scPAS_para`):
#'     \itemize{
#'       \item `alpha` - Alpha values used in model optimization
#'       \item `lambda` - Lambda values used in model optimization
#'       \item `family` - Regression family used
#'       \item `Coefs` - Final model coefficients for each gene
#'       \item `bulk` - Processed bulk expression matrix
#'       \item `phenotype` - Processed phenotype vector
#'       \item `Network` - Gene-gene similarity network used
#'     }
#' }
#'
#' @details
#' This optimized implementation of scPAS integrates bulk and single-cell
#' transcriptomic data to identify phenotype-associated cell subpopulations
#' through a comprehensive analytical workflow:
#'
#' ## Workflow Overview:
#'
#' 1. **Data Preprocessing**:
#'    - Identifies common genes between bulk and single-cell datasets
#'    - Filters ribosomal and mitochondrial genes
#'    - Performs quantile normalization on bulk data
#'    - Optionally imputes single-cell data using specified methods
#'
#' 2. **Network Construction**:
#'    - Builds gene-gene similarity networks from either single-cell or bulk data
#'    - Uses correlation-based similarity measures
#'    - Applies sparse neighborhood network (SNN) construction
#'
#' 3. **Regularized Regression**:
#'    - Implements network-regularized sparse regression (APML0)
#'    - Optimizes alpha and lambda parameters through cross-validation
#'    - Supports multiple regression families (gaussian, binomial, cox)
#'
#' 4. **Risk Score Calculation**:
#'    - Computes phenotype-associated risk scores for each cell
#'    - Uses matrix optimizations for efficient computation
#'
#' 5. **Statistical Validation**:
#'    - Performs permutation testing to assess significance
#'    - Calculates Z-statistics and false discovery rates
#'    - Classifies cells based on statistical thresholds
#'
#'
#' @note
#' The function requires both bulk and single-cell data from related biological
#' conditions. For survival analysis (`family = "cox"`), the phenotype must be
#' a properly formatted survival object or matrix with 'time' and 'status'
#' columns.
#'
#' @examples
#' \dontrun{
#' # Example with continuous phenotype (linear regression)
#' result <- scPAS.optimized(
#'   bulk_dataset = bulk_expr_matrix,
#'   sc_dataset = seurat_obj,
#'   phenotype = continuous_phenotype,
#'   family = "gaussian"
#' )
#'
#' # Example with binary phenotype (logistic regression)
#' result <- scPAS.optimized(
#'   bulk_dataset = bulk_expr_matrix,
#'   sc_dataset = seurat_obj,
#'   phenotype = binary_groups,
#'   family = "binomial",
#'   tag = c("Control", "Disease")
#' )
#'
#' # Example with custom parameters
#' result <- scPAS.optimized(
#'   bulk_dataset = bulk_expr_matrix,
#'   sc_dataset = seurat_obj,
#'   phenotype = survival_data,
#'   family = "cox",
#'   nfeature = 2000,
#'   permutation_times = 5000,
#'   FDR.threshold = 0.01
#' )
#' }
#'
#' @references Xie A, Wang H, Zhao J, Wang Z, Xu J, Xu Y. scPAS: single-cell phenotype-associated subpopulation identifier. Briefings in Bioinformatics. 2024 Nov 22;26(1):bbae655.
#'
#' @importFrom methods as
#'
#' @keywords internal
#' @family scPAS
#'
scPAS.optimized <- function(
    bulk_dataset,
    sc_dataset,
    phenotype,
    assay = 'RNA',
    tag = NULL,
    nfeature = NULL,
    imputation = TRUE,
    imputation_method = c('KNN', 'ALRA'),
    alpha = NULL,
    cutoff = 0.2,
    network_class = c('SC', 'bulk'),
    independent = TRUE,
    family = c("gaussian", "binomial", "cox"),
    permutation_times = 2000,
    FDR.threshold = 0.05
) {
    # Set default assay
    Seurat::DefaultAssay(sc_dataset) <- assay

    # Step 0: Common gene identification with optimized filtering
    common_genes <- if (inherits(sc_dataset, "Seurat")) {
        if (is.null(nfeature)) {
            intersect(rownames(bulk_dataset), rownames(sc_dataset))
        } else if (is.numeric(nfeature) && length(nfeature) == 1) {
            sc_dataset <- Seurat::FindVariableFeatures(
                sc_dataset,
                selection.method = "vst",
                verbose = FALSE,
                nfeatures = nfeature
            )
            var_features <- Seurat::VariableFeatures(sc_dataset)
            intersect(rownames(bulk_dataset), var_features)
        } else if (is.character(nfeature) && length(nfeature) > 1) {
            intersect(rownames(bulk_dataset), nfeature)
        }
    } else {
        ts_cli$cli_alert_info(
            "The single-cell data is not a Seurat object, running default Seurat pipeline."
        )
        sc_dataset <- SCPreProcess(
            sc_dataset,
            quality_control.pattern = c("^MT-"),
            verbose = FALSE
        )
        if (is.null(nfeature)) {
            intersect(rownames(bulk_dataset), rownames(sc_dataset))
        } else if (is.numeric(nfeature)) {
            sc_dataset <- Seurat::FindVariableFeatures(
                sc_dataset,
                selection.method = "vst",
                verbose = FALSE,
                nfeatures = nfeature
            )
            var_features <- Seurat::VariableFeatures(sc_dataset)
            intersect(rownames(bulk_dataset), var_features)
        } else {
            intersect(rownames(bulk_dataset), nfeature)
        }
    }

    # Filter out ribosomal and mitochondrial genes
    gene_patterns <- c("^RP[LS]", "^MT-")
    common_genes <- purrr::reduce(
        gene_patterns,
        function(genes, pattern) {
            genes[!grepl(pattern, genes)]
        },
        .init = common_genes
    )

    if (length(common_genes) == 0) {
        cli::cli_abort(c(
            "x" = "There is no common genes between the given single-cell and bulk samples."
        ))
    }

    # Step 1: Quantile normalization with matrix optimization
    ts_cli$cli_alert_info("Quantile normalization of bulk data.")
    Expression_bulk <- preprocessCore::normalize.quantiles(as.matrix(bulk_dataset[
        common_genes,
    ]))
    rownames(Expression_bulk) <- common_genes
    colnames(Expression_bulk) <- colnames(bulk_dataset)

    # Step 2: Single-cell expression processing with data.table efficiency
    if (imputation) {
        sc_dataset <- scPAS::imputation(
            sc_dataset,
            assay = assay,
            method = imputation_method
        )
        assay <- Seurat::DefaultAssay(sc_dataset)
    }

    ts_cli$cli_alert_info(
        "Extracting single-cell expression profiles..."
    )
    sc_exprs <- SeuratObject::LayerData(sc_dataset)
    Expression_cell <- sc_exprs[common_genes, ]

    # Clean up memory using data.table approach
    rm_vars <- c("sc_exprs", "bulk_dataset")
    rm_vars <- rm_vars[vapply(
        rm_vars,
        exists,
        inherits = TRUE,
        FUN.VALUE = logical(1)
    )]
    if (length(rm_vars) > 0) {
        rm(list = rm_vars)
    }

    # Prepare X matrix
    x <- Matrix::t(Expression_bulk)

    # Step 3: Network construction with matrix optimizations
    if (network_class == 'bulk') {
        ts_cli$cli_alert_info(
            "Constructing a gene-gene similarity by bulk data..."
        )
        cor.m <- cor(x)
    } else {
        ts_cli$cli_alert_info(
            "Constructing a gene-gene similarity by single cell data..."
        )
        # Use matrix operations for efficient correlation
        cor.m <- scPAS::sparse.cor(Matrix::t(Expression_cell))
    }

    # Network construction
    cor.m[cor.m < 0] <- 0
    SNN <- Seurat::FindNeighbors(1 - cor.m, distance.matrix = TRUE)
    Network <- as.matrix(SNN$snn)
    diag(Network) <- 0
    Network <- (Network > 0.2) * 1 # binarization

    # Clean up
    rm(cor.m, SNN)
    gc(verbose = FALSE)

    # Step 4: Model optimization with purrr functional programming
    ts_cli$cli_alert_info(
        "Optimizing the network-regularized sparse regression model..."
    )

    # Prepare Y based on family using purrr pattern matching
    family_processor <- list(
        binomial = function() {
            y <- as.numeric(phenotype)
            z <- table(y)
            ts_cli$cli_alert_info(
                "Current phenotype contains {.val {z[1]}} {tag[1]} and {.val {z[2]}} {tag[2]} samples."
            )
            ts_cli$cli_alert_info(
                "Perform {.strong logistic} regression on the given phenotypes..."
            )
            y
        },
        gaussian = function() {
            y <- as.numeric(phenotype)
            ts_cli$cli_alert_info(
                "Perform linear regression on the given phenotypes..."
            )
            y
        },
        cox = function() {
            y <- as.matrix(phenotype)
            if (ncol(y) != 2) {
                cli::cli_abort(
                    c(
                        "x" = "The size of survival data is wrong. Please check inputs and selected regression type."
                    ),
                    class = "IncorrectNumberOfColumns"
                )
            }
            ts_cli$cli_alert_info(
                "Perform cox regression on the given phenotypes..."
            )
            y
        }
    )

    y <- family_processor[[family]]()

    alpha <- alpha %||%
        c(0.001, 0.005, 0.01, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9)
    lambda <- c()
    for (i in seq_along(alpha)) {
        set.seed(123)

        fit0 <- scPAS::APML0(
            x = x,
            y = y,
            family = family,
            penalty = 'Net',
            Omega = Network,
            alpha = alpha[i],
            nlambda = 100,
            nfolds = min(10, nrow(x))
        )

        fit1 <- scPAS::APML0(
            x = x,
            y = y,
            family = family,
            penalty = 'Net',
            Omega = Network,
            alpha = alpha[i],
            lambda = fit0$lambda.min
        )
        lambda <- c(lambda, fit0$lambda.min)
        # Extract coefficients using matrix indexing
        Coefs <- if (family == "binomial") {
            as.numeric(fit1$Beta[2:(ncol(x) + 1)])
        } else {
            as.numeric(fit1$Beta)
        }

        names(Coefs) <- colnames(x)

        # Fast feature counting using logical indexing
        pos_features <- colnames(x)[Coefs > 0]
        neg_features <- colnames(x)[Coefs < 0]
        percentage <- (length(pos_features) + length(neg_features)) / ncol(x)

        cli::cli_h2("At alpha = {.val {alpha[i]}}")
        cli::cli_text("lambda = {.val {fit0$lambda.min}}")
        cli::cli_text(
            "scPAS identified {.val {length(pos_features)}} risk+ features and {.val {length(neg_features)}} risk- features."
        )
        percentage_show <- round(percentage * 100, digits = 3)
        cli::cli_text(
            "The percentage of selected feature is: {.val {percentage_show}}%"
        )
        if (percentage < cutoff) {
            break
        }
    }

    # Step 5: Risk score calculation with matrix optimizations
    ts_cli$cli_alert_info("Calculating quantified risk scores...")

    # Fast sparse matrix scaling and risk calculation
    scaled_exp <- Seurat:::FastSparseRowScale(
        Expression_cell,
        display_progress = FALSE
    )
    scaled_exp[is.na(scaled_exp)] <- 0
    scaled_exp <- methods::as(scaled_exp, "dgCMatrix")

    # Fast matrix multiplication for risk scores
    risk_score <- Matrix::crossprod(scaled_exp, Coefs)

    # Step 6: Permutation test with purrr and matrix optimizations
    ts_cli$cli_alert_info(
        "Qualitative identification by permutation test program with {.val {permutation_times}} times random perturbations..."
    )

    set.seed(12345)

    randomPermutation_list <- lapply(
        seq_len(permutation_times),
        function(i) {
            set.seed(1234 + i)
            sample(Coefs, length(Coefs), replace = FALSE)
        }
    )

    randomPermutation <- matrix(
        unlist(randomPermutation_list),
        nrow = length(Coefs),
        ncol = permutation_times,
        byrow = FALSE
    )
    randomPermutation <- as(randomPermutation, "dgCMatrix")
    rownames(randomPermutation) <- names(Coefs)

    # Matrix multiplication for background scores
    risk_score.background <- Matrix::crossprod(scaled_exp, randomPermutation)
    rm(randomPermutation_list, randomPermutation)
    # Calculate background statistics
    if (independent) {
        risk_bg_matrix <- as.matrix(risk_score.background)
        mean.background <- rowMeans(risk_bg_matrix)
        sd.background <- rowSds(risk_bg_matrix)
        rm(risk_bg_matrix)
    } else {
        risk_bg_vector <- as.vector(risk_score.background)
        mean.background <- mean(risk_bg_vector)
        sd.background <- stats::sd(risk_bg_vector)
        rm(risk_bg_vector)
    }
    gc(verbose = FALSE)

    # Z-score calculation using vectorized operations with numerical stability
    # Add small epsilon to avoid division by zero
    sd.background[sd.background == 0] <- .Machine$double.eps
    # Z-score calculation using vectorized operations
    Z <- (risk_score[, 1] - mean.background) / sd.background

    # Fast p-value and FDR calculation
    p.value <- stats::pnorm(abs(Z), mean = 0, sd = 1, lower.tail = FALSE)
    q.value <- stats::p.adjust(p.value, method = 'BH')

    # Use data.table for fast data frame creation
    risk_score_data.frame <- data.table::data.table(
        cell = colnames(Expression_cell),
        raw_score = risk_score[, 1],
        Z.statistics = Z,
        p.value = p.value,
        FDR = q.value
    )

    # Fast conditional labeling using data.table
    risk_score_data.frame[,
        cell_label := data.table::fcase(
            Z > 0 & q.value <= FDR.threshold ,
            "Positive"                       ,
            Z < 0 & q.value <= FDR.threshold ,
            "Negative"                       ,
            default = "Neutral"
        )
    ]

    sc_dataset <- AddMisc(
        sc_dataset,
        scPAS_para = list(
            alpha = alpha,
            lambda = lambda,
            family = family,
            Coefs = Coefs
            # ,bulk = x,
            # phenotype = y,
            # Network = Network
        ),
        cover = FALSE
    )

    sc_dataset$scPAS_RS <- risk_score_data.frame$raw_score
    sc_dataset$scPAS_NRS <- risk_score_data.frame$Z.statistics
    sc_dataset$scPAS_Pvalue <- risk_score_data.frame$p.value
    sc_dataset$scPAS_FDR <- risk_score_data.frame$FDR
    sc_dataset$scPAS <- risk_score_data.frame$cell_label

    return(sc_dataset)
}
