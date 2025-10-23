#' @title Perform scAB Screening Analysis
#'
#' @description
#' Implements the scAB algorithm to identify phenotype-associated cell subpopulations
#' in single-cell RNA-seq data by integrating matched bulk expression and phenotype
#' information. Uses non-negative matrix factorization (NMF) with dual regularization
#' for phenotype association and cell-cell similarity.
#'
#' @param matched_bulk Normalized bulk expression matrix (genes × samples) where:
#'        - Columns match `phenotype` row names
#'        - Genes match features in `sc_data`
#' @param sc_data Seurat object containing preprocessed single-cell data:
#' @param phenotype Data frame with clinical annotations where:
#'        - Rows correspond to `matched_bulk` columns
#'        - For survival: contains `time` and `status` columns
#' @param label_type Character specifying phenotype label type (e.g., "SBS1", "time"), stored in `scRNA_data@misc`
#' @param phenotype_class Analysis mode:
#'        - `"binary"`: Case-control design (e.g., responder/non-responder)
#'        - `"survival"`: Time-to-event analysis data.frame
#' @param alpha Coefficient of phenotype regularization (default=0.005).
#' @param alpha_2 Coefficent of cell-cell similarity regularization (default=5e-05).
#' @param maxiter NMF optimization iterations (default=2000).
#' @param tred Z-score threshold (default=2).
#'
#' @return A list containing:
#' \describe{
#'   \item{scRNA_data}{Filtered Seurat object with selected cells}
#'   \item{scAB_result}{scAB screening result}
#' }
#'
#' @references
#' Zhang Q, Jin S, Zou X. scAB detects multiresolution cell states with clinical significance by integrating single-cell genomics and bulk sequencing data. Nucleic Acids Research. 2022 Nov 28;50(21):12112–30.
#'
#' @section LICENSE:
#' Licensed under the GNU General Public License version 3 (GPL-3.0).
#' A copy of the license is available at <https://www.gnu.org/licenses/gpl-3.0.en.html>.
#'
#' @examples
#' \dontrun{
#' # Binary phenotype example
#' result <- DoscAB(
#'   matched_bulk = bulk_matrix,
#'   sc_data = seurat_obj,
#'   phenotype = clinical_df,
#'   label_type = "disease_status",
#'   phenotype_class = "binary",
#'   alpha = 0.005,
#'   alpha_2 = 5e-05,
#'   maxiter = 2000,
#'   tred = 2
#' )
#' }
#'
#' @importFrom scAB create_scAB select_K scAB findSubset
#' @importFrom cli cli_alert_info col_green
#'
#' @family screen_method
#' @family scAB
#'
DoscAB <- function(
    matched_bulk,
    sc_data,
    phenotype,
    label_type = "scAB",
    phenotype_class = c("binary", "survival"),
    alpha = 0.005,
    alpha_2 = 5e-05,
    maxiter = 2000,
    tred = 2
) {
    chk::chk_is(matched_bulk, c("matrix", "data.frame"))
    chk::chk_is(sc_data, "Seurat")
    chk::chk_character(label_type)
    phenotype_class %<>% MatchArg(c("binary", "survival"), NULL)
    chk::chk_range(alpha)
    chk::chk_range(alpha_2)
    chk::chk_number(maxiter)
    chk::chk_number(tred)
    # scAB can't tolerate NA
    chk::chk_not_any_na(matched_bulk)
    chk::chk_not_any_na(phenotype)

    # scAB is more strict than Scissor and scPAS
    if (phenotype_class == "survival") {
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

    ts_cli$cli_alert_info(cli::col_green("Start scAB screening."))

    scAB_obj <- create_scAB.v5(
        Object = sc_data,
        bulk_dataset = matched_bulk,
        phenotype = phenotype,
        method = phenotype_class
    )

    ts_cli$cli_alert_info("Selecting K...")

    k <- scAB::select_K(
        scAB_obj,
        K_max = 20L,
        repeat_times = 10L,
        maxiter = 2000L,
        seed = 0L,
        verbose = FALSE
    )

    ts_cli$cli_alert_info(
        "Run NMF with phenotype and cell-cell similarity regularization at",
        "K = {.val {k}}."
    )

    scAB_result <- scAB::scAB(
        Object = scAB_obj,
        K = k,
        alpha = alpha,
        alpha_2 = alpha_2,
        maxiter = maxiter
    )

    ts_cli$cli_alert_info("Screening cells...")

    sc_data <- scAB::findSubset(
        sc_data,
        scAB_Object = scAB_result,
        tred = tred
    ) %>%
        AddMisc(
            scAB_type = label_type,
            scAB_para = list(
                iter = scAB_result$iter,
                loss = scAB_result$loss,
                method = scAB_result$method,
                tred = tred
            ),
            cover = FALSE
        )

    sc_data[[]] <- dplyr::rename(sc_data[[]], scAB = `scAB_select`) %>%
        dplyr::mutate(
            scAB = dplyr::case_when(
                scAB == "Other cells" ~ "Other",
                scAB == "scAB+ cells" ~ "Positive",
                TRUE ~ "NULL"
            )
        )

    ts_cli$cli_alert_info(
        cli::col_green("scAB screening done.")
    )

    return(list(scRNA_data = sc_data, scAB_result = scAB_result))
}


#' @title scAB_data preprocess
#' @description
#' preprocess the single-cell data, bulk data, and phenotype data.
#'
#' @param Object Seurat object
#' @param bulk_dataset matrix of bulk data
#' @param phenotype Phenotype data, a matrix with two columns "time" and "state", or a vector
#' @param method method "survival" or "binary"
#'
#' @return a scAB_data
#' @importFrom preprocessCore normalize.quantiles
#'
#' @family scAB
#'
#' @keywords internal
#'
create_scAB.v5 <- function(
    Object,
    bulk_dataset,
    phenotype,
    method = c("survival", "binary")
) {
    # cell neighbors
    if ("RNA_snn" %chin% names(Object@graphs)) {
        A <- as.matrix(SeuratObject::Graphs(object = Object, slot = "RNA_snn"))
        cli::cli_alert_info(
            " Using {.val RNA_snn} graph for network."
        )
    } else if ("integrated_snn" %chin% names(Object@graphs)) {
        A <- as.matrix(SeuratObject::Graphs(
            object = Object,
            slot = "integrated_snn"
        ))
        cli::cli_alert_info(
            "Using {.val integrated_snn} graph for network."
        )
    } else {
        cli::cli_abort(c(
            "x" = "No `RNA_snn` or `integrated_snn` graph in the given Seurat object. Please check `Object@graphs`."
        ))
    }
    diag(A) <- 0
    A[which(A != 0)] <- 1
    degrees <- matrixStats::rowSums2(A)
    D <- diag(degrees)
    eps <- 2.2204e-256
    D12 <- diag(1 / sqrt(pmax(degrees, eps)))

    L <- D12 %*% (D - A) %*% D12 # Normalized Graph Laplacian
    Dhat <- D12 %*% (D) %*% D12
    Ahat <- D12 %*% (A) %*% D12

    # similarity matrix
    sc_exprs <- as.data.frame(SeuratObject::LayerData(Object))
    common <- intersect(rownames(bulk_dataset), rownames(sc_exprs))
    dataset0 <- cbind(bulk_dataset[common, ], sc_exprs[common, ]) # Dataset before quantile normalization.
    dataset1 <- preprocessCore::normalize.quantiles(as.matrix(dataset0)) # Dataset after quantile normalization.
    rownames(dataset1) <- common
    colnames(dataset1) <- colnames(dataset0)
    Expression_bulk <- dataset1[, seq_len(ncol(bulk_dataset))]
    Expression_cell <- dataset1[, (ncol(bulk_dataset) + 1):ncol(dataset1)]
    X <- stats::cor(Expression_bulk, Expression_cell)
    X <- X / norm(X, "F")

    # phenotype ranking
    if (method == "survival") {
        ss <- scAB::guanrank(phenotype[, c("time", "status")])
        S <- diag(1 - ss[rownames(phenotype), 3])
    } else {
        S <- diag(1 - phenotype)
    }
    # return
    obj <- list(
        X = X,
        S = S,
        L = L,
        D = Dhat,
        A = Ahat,
        phenotype = phenotype,
        method = method
    )
    class(obj) <- "scAB_data"
    return(obj)
}

#' @title Selection of parameter K
#'
#' @param Object a scAB_data object
#' @param k_max  the maximum value of the rank in the matrix factorization
#' @param repeat_times  the number of repetitions
#' @param seed random seed
#' @param verbose Logical, whether to print output message
#'
#' @return A integer value of K
#'
#' @family scAB
#'
#' @keywords internal
#'
select_K.optimized <- function(
    Object,
    K_max = 20L,
    repeat_times = 10L,
    maxiter = 2000L,
    seed = 0L,
    verbose = FALSE
) {
    X <- Object$X
    set.seed(seed)

    K_all <- 2:K_max
    n_K <- length(K_all)
    dist_K <- matrix(NA_real_, nrow = n_K, ncol = repeat_times)

    dist_all <- norm(X, "F")

    eii <- numeric(n_K)
    row_means <- numeric(n_K)

    for (Ki_idx in seq_along(K_all)) {
        Ki <- K_all[Ki_idx]

        for (Kj in seq_len(repeat_times)) {
            res_ij <- scAB::NMF(X = X, K = Ki, maxiter = maxiter)
            diff_matrix <- X - res_ij$W %*% res_ij$H
            dist_K[Ki_idx, Kj] <- matrixStats::sum2(diff_matrix^2) # equivalent to `norm(., "F")^2`
        }
        if (verbose) {
            message(sprintf(
                "loss of %d: %.6f",
                Ki,
                matrixStats::mean2(dist_K[Ki_idx, ])
            ))
        }

        if (Ki_idx == 1) {
            next
        } # Ki == 2
        row_means <- matrixStats::rowMeans2(dist_K)
        numerator <- row_means[Ki_idx - 1] - row_means[Ki_idx]
        denominator <- row_means[1] - row_means[Ki_idx]
        eii[Ki_idx] <- numerator / denominator

        if (numerator <= 0) {
            break
        }

        if (eii[Ki_idx] < 0.05) {
            break
        }
    }
    K <- Ki - 1
    return(K)
}
