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
#' @param verbose Logical indicating whether to print progress messages. Default: `TRUE`.
#' @param ... For future update.
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
    tred = 2,
    verbose = TRUE,
    ...
) {
    chk::chk_is(matched_bulk, c("matrix", "data.frame"))
    chk::chk_is(sc_data, "Seurat")
    chk::chk_character(label_type)
    phenotype_class <- MatchArg(phenotype_class, c("binary", "survival"), NULL)
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

    if (verbose) {
        ts_cli$cli_alert_info(cli::col_green("Start scAB screening."))
    }

    scAB_obj <- create_scAB.v5(
        Object = sc_data,
        bulk_dataset = matched_bulk,
        phenotype = phenotype,
        method = phenotype_class,
        verbose = verbose
    )

    if (verbose) {
        ts_cli$cli_alert_info("Selecting K...")
    }

    k <- select_K.optimized(
        scAB_obj,
        K_max = 20L,
        repeat_times = 10L,
        maxiter = 2000L, # default in scAB
        seed = 0L,
        verbose = verbose
    )

    if (verbose) {
        ts_cli$cli_alert_info(
            "Run NMF with phenotype and cell-cell similarity regularization at K = {.val {k}}"
        )
    }

    scAB_result <- scAB.optimized(
        Object = scAB_obj,
        K = k,
        alpha = alpha,
        alpha_2 = alpha_2,
        maxiter = maxiter
    )

    if (verbose) {
        ts_cli$cli_alert_info("Screening cells...")
    }

    sc_data <- findSubset.optimized(
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

    if (verbose) {
        ts_cli$cli_alert_info(
            cli::col_green("scAB screening done.")
        )
    }

    list(scRNA_data = sc_data, scAB_result = scAB_result)
}


#' @title scAB_data preprocess
#' @description
#' preprocess the single-cell data, bulk data, and phenotype data.
#'
#' @param Object Seurat object
#' @param bulk_dataset matrix of bulk data
#' @param phenotype Phenotype data, a matrix with two columns "time" and "state", or a vector
#' @param method method "survival" or "binary"
#' @param verbose Logical, whether to print messages.
#' @param ... For future updates.
#'
#' @return a scAB_data
#'
#' @family scAB
#'
#' @keywords internal
#'
create_scAB.v5 <- function(
    Object,
    bulk_dataset,
    phenotype,
    method = c("survival", "binary"),
    verbose = TRUE,
    ...
) {
    # cell neighbors
    if ("RNA_snn" %chin% names(Object@graphs)) {
        A <- Matrix::Matrix(SeuratObject::Graphs(
            object = Object,
            slot = "RNA_snn"
        ))
        if (verbose) {
            cli::cli_alert_info(
                " Using {.val RNA_snn} graph for network."
            )
        }
    } else if ("integrated_snn" %chin% names(Object@graphs)) {
        A <- Matrix::Matrix(SeuratObject::Graphs(
            object = Object,
            slot = "integrated_snn"
        ))

        if (verbose) {
            cli::cli_alert_info(
                "Using {.val integrated_snn} graph for network."
            )
        }
    } else {
        cli::cli_abort(c(
            "x" = "No `RNA_snn` or `integrated_snn` graph in the given Seurat object. Please check `Object@graphs`."
        ))
    }
    Matrix::diag(A) <- 0
    A@x[which(A@x != 0)] <- 1
    degrees <- Matrix::rowSums(A)
    D <- Matrix::diag(degrees)
    eps <- 2.2204e-256 # In the original implementation of scAB, a custom-defined `eps` is used instead of `.Machine$double.eps`.
    D12 <- Matrix::diag(1 / sqrt(pmax(degrees, eps)))

    L <- D12 %*% (D - A) %*% D12 # Normalized Graph Laplacian
    Dhat <- D12 %*% (D) %*% D12
    Ahat <- D12 %*% (A) %*% D12

    # similarity matrix
    sc_exprs <- Matrix::Matrix(SeuratObject::LayerData(Object))
    common <- intersect(rownames(bulk_dataset), rownames(sc_exprs))
    dataset0 <- cbind(bulk_dataset[common, ], sc_exprs[common, ]) # Dataset before quantile normalization.
    dataset1 <- normalize.quantiles(as.matrix(dataset0)) # Dataset after quantile normalization.
    dataset1 <- Matrix::Matrix(dataset1)
    rownames(dataset1) <- common
    colnames(dataset1) <- colnames(dataset0)

    ncol_bulk <- ncol(bulk_dataset)
    Expression_bulk <- as.matrix(dataset1[, seq_len(ncol_bulk)])
    Expression_cell <- as.matrix(dataset1[, (ncol_bulk + 1):ncol(dataset1)])

    X <- stats::cor(Expression_bulk, Expression_cell)
    # X <- X / norm(X, "F")
    X <- X / sqrt(sum(X^2))

    # phenotype ranking
    if (method == "survival") {
        ss <- scAB::guanrank(phenotype[, c("time", "status")])
        S <- diag(1 - ss[rownames(phenotype), 3])
    } else {
        S <- diag(1 - phenotype)
    }
    # return
    obj <- list(
        X = Matrix::Matrix(X),
        S = Matrix::Matrix(S),
        L = L,
        D = Dhat,
        A = Ahat,
        phenotype = phenotype,
        method = method
    )
    class(obj) <- "scAB_data"

    obj
}

#' @title Selection of Parameter K for Non-negative Matrix Factorization
#'
#' @description
#' Automatically determines the optimal rank (K) for non-negative matrix
#' factorization using an empirical indicator method. This function evaluates
#' multiple candidate ranks and selects the one that provides the best
#' trade-off between model complexity and reconstruction accuracy.
#'
#' @param Object A scAB_data object containing the data matrix to be factorized.
#' @param K_max The maximum rank value to consider in the search. Must be at
#'              least 2. Defaults to 20.
#' @param repeat_times The number of repeated NMF runs for each candidate rank
#'                     to account for random initialization variability.
#'                     Defaults to 10.
#' @param maxiter The maximum number of iterations for each NMF run.
#'                Defaults to 2000.
#' @param seed Random seed for reproducible results. Defaults to 0.
#' @param verbose Logical indicating whether to print progress messages and
#'                intermediate results. Defaults to FALSE.
#'
#' @return An integer value representing the selected optimal rank K.
#'
#' @note
#' This function is from scAB package,
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

    # dist_all <- norm(X, "F")
    dist_all <- sqrt(sum(X^2))
    # initialize
    eii <- numeric(n_K)
    row_means <- numeric(n_K)

    for (Ki_idx in seq_along(K_all)) {
        Ki <- K_all[Ki_idx]

        for (Kj in seq_len(repeat_times)) {
            res_ij <- NMF.optimized(X = X, K = Ki, maxiter = maxiter)
            diff_matrix <- X - res_ij$W %*% res_ij$H
            dist_K[Ki_idx, Kj] <- sum(diff_matrix^2) # equivalent to `norm(., "F")^2`
        }
        if (verbose) {
            cli::cli_li(sprintf(
                "loss of %d: %.6f",
                Ki,
                mean(dist_K[Ki_idx, ])
            ))
        }
        # Ki == 2
        if (Ki_idx == 1) {
            next
        }
        row_means <- rowMeans(dist_K)
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

    eval(K)
}

#' @title Non-negative Matrix Factorization from scAB
#'
#' @description
#' An R implementation of classical non-negative matrix factorization (NMF).
#'
#' @param X A non-negative numeric matrix (features x samples) to be factorized.
#'          All elements must be non-negative.
#' @param K The rank of factorization, i.e., the number of components/latent
#'          features to extract. Must be a positive integer smaller than
#'          both dimensions of X.
#' @param maxiter The maximum number of iterations for the optimization algorithm.
#'                Defaults to 2000.
#' @param tol The convergence tolerance for the loss function. The algorithm
#'            stops when the absolute change in Euclidean distance between
#'            consecutive iterations is less than this value. Defaults to 1e-5.
#'
#' @return A list containing the factorization results:
#' \itemize{
#'   \item \code{W} - The basis matrix (features x K), representing the
#'                    learned features or components
#'   \item \code{H} - The coefficient matrix (K x samples), representing the
#'                    weights or activations of components for each sample
#'   \item \code{iter} - The number of iterations actually performed
#'   \item \code{loss} - The final value of the objective function (squared
#'                       Euclidean distance)
#' }
#'
#' @seealso
#' This function is from scAB package, and it is not recommended to use it because the computational efficiency of the R language is not very high.
#' For more advanced NMF implementations, see:
#' \code{\link[NMF]{nmf}} from the NMF package, which is written in C++ and is much faster than this function.
#'
#' @keywords internal
#' @family scAB
#'
NMF.optimized <- function(X, K, maxiter = 2000L, tol = 1e-5) {
    eps <- 2.2204e-256

    nr <- nrow(X)
    nc <- ncol(X)

    # initiate
    W <- Matrix::Matrix(stats::runif(nr * K), nrow = nr, ncol = K)
    H <- Matrix::Matrix(stats::runif(K * nc), nrow = K, ncol = nc)

    old_eucl <- Inf

    for (iter in seq_len(maxiter)) {
        WtW <- crossprod(W) # t(W) %*% W
        WtX <- crossprod(W, X) # t(W) %*% X
        # avoid zero with eps
        H <- H * WtX / (WtW %*% H + eps)

        HHt <- tcrossprod(H) # H %*% t(H)
        XHt <- tcrossprod(X, H) #  X %*% t(H)

        W <- W * XHt / (W %*% HHt + eps)

        if (iter != 1) {
            WH <- W %*% H
            eucl_dist <- sum((X - WH)^2)

            if (iter > 1) {
                d_eucl <- abs(eucl_dist - old_eucl)

                if (d_eucl < tol) {
                    break
                }
            }

            old_eucl <- eucl_dist
        }
    }

    final_loss <- sum((X - W %*% H)^2)

    list(
        W = W,
        H = H,
        iter = iter,
        loss = final_loss
    )
}


#' @title Non-negative Matrix Factorization with phenotype and cell-cell similarity regularization.
#' @description
#' Non-negative Matrix Factorization with phenotype and cell-cell similarity regularization,
#' for identifing phenotype-associated cell states at different resolutions.
#'
#' @param Object  a scAB_data object
#' @param K  the rank of matrix factorization
#' @param maxiter the maximum number of iterations
#' @param alpha Coefficient of phenotype regularization
#' @param alpha_2 Coefficient of cell-cell similarity regularization
#'
#' @return a list with the submatrix and loss value
#' @keywords internal
#' @family scAB
#'
#'
scAB.optimized <- function(
    Object,
    K,
    alpha = 0.005,
    alpha_2 = 0.005,
    maxiter = 2000L,
    convergence_threshold = 1e-5
) {
    seed <- ifelse(Object$method == "survival", 7L, 5L)
    if (Object$method != "") {
        set.seed(seed)
    }
    X <- Object$X
    A <- Object$A
    L <- Object$L
    D <- Object$D
    S <- Object$S
    eps <- 2.2204e-256
    nr <- nrow(X)
    nc <- ncol(X)
    W <- Matrix::Matrix(stats::runif(nr * K), nrow = nr, ncol = K)
    H <- Matrix::Matrix(stats::runif(K * nc), nrow = K, ncol = nc)
    SS <- S %*% S

    loss_func <- function(X, W, H, S, L, alpha, alpha_2) {
        # loss <- norm(X - W %*% H, "F")^2 +
        #     alpha * (norm(S %*% W, "F")^2) +
        #     alpha_2 * sum(diag(H %*% L %*% t(H)))
        loss1 <- sum((X - W %*% H)^2)
        loss2 <- alpha * sum((S %*% W)^2)
        loss3 <- alpha_2 * sum(H * (H %*% L)) # diag(H %*% L %*% t(H)) is equivalent to rowSums(H * (H %*% L))

        return(loss1 + loss2 + loss3)
    }
    # tD <- t(D)
    # tA <- t(A)
    tX <- Matrix::t(X)

    for (iter in seq_len(maxiter)) {
        # W_old <- W # no use
        H <- H *
            (crossprod(W, X) + alpha_2 * tcrossprod(H, A)) /
            (crossprod(W, W) %*% H + alpha_2 * tcrossprod(H, D) + eps)
        Pena <- SS %*% W
        W <- W *
            tcrossprod(X, H) /
            (W %*% tcrossprod(H, H) + alpha * Pena + eps)

        # tW <- t(W) # no use

        if (iter != 1) {
            eucl_dist <- loss_func(X, W, H, S, L, alpha, alpha_2)
            d_eucl <- abs(eucl_dist - old_eucl)
            if (d_eucl < convergence_threshold) {
                break
            }
            old_eucl <- eucl_dist
        } else {
            old_eucl <- loss_func(X, W, H, S, L, alpha, alpha_2)
        }
    }
    list(
        W = W,
        H = H,
        iter = iter,
        loss = old_eucl,
        method = Object$method
    )
}

#' @title Subsets identification
#'
#' @param Object a Seurat object
#' @param scAB_Object a scAB_data object
#' @param tred threshold
#'
#' @return a Seurat object
#' @keywords internal
#' @family scAB
#'
#'
findSubset.optimized <- function(Object, scAB_Object, tred = 2L) {
    do.dip <- ifelse(scAB_Object$method == "binary", 1L, 0L)
    H <- as.matrix(scAB_Object$H)
    module <- scAB::findModule(H, tred = tred, do.dip = do.dip)
    scAB_index <- unique(unlist(module))
    # Add scAB column in metadata
    n_cells <- ncol(Object)
    scAB_select <- rep("Other", n_cells)
    scAB_select[scAB_index] <- "Positive"
    Object$scAB <- scAB_select

    for (i in seq_along(module)) {
        M <- rep("Other", n_cells)
        M[as.numeric(module[[i]])] <- "Positive"
        Object <- Seurat::AddMetaData(
            Object,
            metadata = M,
            col.name = paste0("scAB_Subset", i)
        ) %>%
            Seurat::AddMetaData(
                metadata = H[i, ],
                col.name = paste0("Subset", i, "_loading")
            )
    }

    Object
}
