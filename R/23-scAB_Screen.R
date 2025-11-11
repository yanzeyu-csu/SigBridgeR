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
#' @param maxiter Maximum number of iterations for NMF (default=2000).
#' @param tred Z-score threshold in finding subsets (default=2).
#' @param ... Additional arguments. Currently supports:
#'    - `verbose`: Logical indicating whether to print progress messages. Defaults to `TRUE`.
#'    - `seed`: For reproducibility, default is `123L`
#'    - `parallel`: Logical. When `alpha` or `alpha_2` is a numeric vector or NULL, whether to enable parallel mode to search for the optimal `alpha` and `alpha_2`. Default: `FALSE`.
#'    - `workers`: Number of workers to use in parallel mode. Default: `NULL` (use all 4 cores).
#
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
    alpha = c(0.005, NULL),
    alpha_2 = c(5e-05, NULL),
    maxiter = 2000L,
    tred = 2L,
    ...
) {
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

    dots <- rlang::list2(...)
    verbose <- dots$verbose %||% getFuncOption("verbose")
    seed <- dots$seed %||% getFuncOption("seed")
    parallel <- dots$parallel %||% getFuncOption("parallel")
    workers <- dots$workers %||% getFuncOption("workers")

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
        Object = scAB_obj,
        K_max = 20L,
        repeat_times = 10L,
        maxiter = 2000L, # default in scAB
        seed = seed,
        verbose = verbose
    )

    if (verbose) {
        ts_cli$cli_alert_info(
            "Run NMF with phenotype and cell-cell similarity regularization at K = {.val {k}}"
        )
    }

    # Find optimal alpha and alpha_2
    if (
        is.null(alpha) ||
            is.null(alpha_2) ||
            length(alpha) > 1 ||
            length(alpha_2) > 1
    ) {
        para_list <- select_alpha.optimized(
            Object = scAB_obj,
            method = phenotype_class,
            K = k,
            cross_k = 5,
            para_1_list = alpha %||% c(0.01, 0.005, 0.001),
            para_2_list = alpha_2 %||% c(0.01, 0.005, 0.001),
            seed = seed,
            parallel = parallel,
            workers = workers,
            verbose = verbose
        )

        alpha <- para_list$para$alpha_1
        alpha_2 <- para_list$para$alpha_2
    }

    scAB_result <- scAB.optimized(
        Object = scAB_obj,
        K = k,
        alpha = alpha,
        alpha_2 = alpha_2,
        maxiter = maxiter,
        convergence_threshold = 1e-05
    )

    if (verbose) {
        ts_cli$cli_alert_info("Screening cells...")
    }

    sc_data <- findSubset.optimized(
        Object = sc_data,
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
#' @seealso [guanrank2()]
#'
create_scAB.v5 <- function(
    Object,
    bulk_dataset,
    phenotype,
    method = c("survival", "binary"),
    verbose = getFuncOption("verbose"),
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
    X <- X / Matrix::norm(X, "F")

    # phenotype ranking
    if (method == "survival") {
        ss <- guanrank2(phenotype[, c("time", "status")])
        S <- Matrix::diag(1 - ss[rownames(phenotype), 3]) # 3 is the rank column
    } else {
        S <- Matrix::diag(1 - phenotype)
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

#' @title Guanrank for Survival Data
#'
#' @param mat A matrix, data frame, or data table containing survival data.
#'            Must contain at least two columns representing time and status.
#'            The function will use the first two columns as time and status
#'            respectively.
#' @param complete Logical indicating whether to include the survival rate
#'                 column in the output. If \code{FALSE} (default), returns
#'                 only time, status, and rank columns.
#'
#' @return A matrix with the following columns:
#'   \itemize{
#'     \item \code{time}: Original survival times (sorted)
#'     \item \code{status}: Original event indicators (0 = censored, 1 = event)
#'     \item \code{rank}: Computed Guan's rank values (normalized to \[0,1\])
#'     \item \code{survival_rate}: (if \code{complete = TRUE}) Kaplan-Meier survival probabilities from \code{km_curve}
#'   }
#'
#'
#' @examples
#' \dontrun{
#' # Create sample survival data (matrix format)
#' survival_data <- matrix(c(
#'   10, 1,   # Event at time 10
#'   15, 0,   # Censored at time 15
#'   20, 1,   # Event at time 20
#'   25, 0,   # Censored at time 25
#'   30, 1    # Event at time 30
#' ), ncol = 2, byrow = TRUE)
#'
#' # Compute Guan's rank
#' result <- guanrank(survival_data)
#' print(result)
#'
#' # Get complete output with survival rates
#' complete_result <- guanrank(survival_data, complete = TRUE)
#' print(complete_result)
#'
#' # Using data frame input
#' df_data <- data.frame(
#'   time = c(10, 15, 20, 25, 30),
#'   status = c(1, 0, 1, 0, 1),
#'   other_col = c("A", "B", "C", "D", "E")  # This column will be ignored
#' )
#' df_result <- guanrank(df_data)
#' print(df_result)
#' }
#'
#' @family scAB
#' @keywords internal
#'
guanrank2 <- function(mat, complete = FALSE) {
    mat <- as.matrix(mat)
    mat <- mat[, 1:2]
    colnames(mat) = c("time", "status")
    storage.mode(mat) = "numeric"
    mat <- mat[order(mat[, "time"]), ]

    mat_curve <- scAB::km_curve(mat)
    n <- nrow(mat_curve)

    mat_guanrank <- cbind(mat_curve, rank = rep(0, n))
    vect <- mat_guanrank[, "time"]
    vecs <- mat_guanrank[, "status"]
    # vectorized pipeline is not faster than the loop as I expected.
    for (i in seq_len(n)) {
        tA <- mat_guanrank[i, "time"]
        rA <- mat_guanrank[i, "survival_rate"]
        sA <- mat_guanrank[i, "status"]
        if (sA == 1) {
            tBgttA <- mat_guanrank[vect > tA, "survival_rate"]
            tBletA_sBeq0 <- mat_guanrank[
                vect <= tA & vecs == 0,
                "survival_rate"
            ]
            tBeqtA_sBeq1 <- mat_guanrank[
                vect == tA & vecs == 1,
                "survival_rate"
            ]
            mat_guanrank[i, "rank"] = ifelse(
                length(tBgttA) == 0,
                0,
                1 * length(tBgttA)
            ) +
                ifelse(length(tBletA_sBeq0) == 0, 0, sum(rA / tBletA_sBeq0)) +
                ifelse(length(tBeqtA_sBeq1) == 0, 0, 0.5 * length(tBeqtA_sBeq1))
        }
        if (sA == 0) {
            tBgetA_sBeq0 <- mat_guanrank[
                vect >= tA & vecs == 0,
                "survival_rate"
            ]
            tBgetA_sBeq1 <- mat_guanrank[
                vect >= tA & vecs == 1,
                "survival_rate"
            ]
            tBlttA_sBeq0 <- mat_guanrank[vect < tA & vecs == 0, "survival_rate"]
            mat_guanrank[i, "rank"] = ifelse(
                length(tBgetA_sBeq0) == 0,
                0,
                sum(1 - 0.5 * tBgetA_sBeq0 / rA)
            ) +
                ifelse(
                    length(tBgetA_sBeq1) == 0,
                    0,
                    sum(1 - tBgetA_sBeq1 / rA)
                ) +
                ifelse(
                    length(tBlttA_sBeq0) == 0,
                    0,
                    sum(0.5 * rA / tBlttA_sBeq0)
                )
        }
    }
    rank <- mat_guanrank[, "rank"]
    # 0.5 is the correction for self-comparison
    # normalization to [0,1]
    mat_guanrank[, "rank"] <- (rank - 0.5) / max(rank)
    if (!complete) {
        mat_guanrank <- mat_guanrank[, c("time", "status", "rank")]
    }

    mat_guanrank
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
#' @seealso [NMF.optimized()]
#'
select_K.optimized <- function(
    Object,
    K_max = 20L,
    repeat_times = 10L,
    maxiter = 2000L,
    seed = 0L,
    verbose = getFuncOption("verbose")
) {
    X <- Object$X
    set.seed(seed)

    K_all <- 2:K_max
    n_K <- length(K_all)
    dist_K <- matrix(NA_real_, nrow = n_K, ncol = repeat_times)

    dist_all <- Matrix::norm(X, "F")
    # initialize
    eii <- numeric(n_K)
    row_means <- numeric(n_K)

    for (Ki_idx in seq_along(K_all)) {
        Ki <- K_all[Ki_idx]

        for (Kj in seq_len(repeat_times)) {
            res_ij <- NMF.optimized(X = X, K = Ki, maxiter = maxiter)
            diff_matrix <- X - res_ij$W %*% res_ij$H
            dist_K[Ki_idx, Kj] <- Matrix::norm(diff_matrix, "F")^2
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
            eucl_dist <- Matrix::norm(X - W %*% H, "F")^2

            if (iter > 1) {
                d_eucl <- abs(eucl_dist - old_eucl)

                if (d_eucl < tol) {
                    break
                }
            }

            old_eucl <- eucl_dist
        }
    }

    list(
        W = W,
        H = H,
        iter = iter,
        loss = eucl_dist
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
        loss1 <- Matrix::norm(X - W %*% H, "F")^2
        loss2 <- alpha * (Matrix::norm(S %*% W, "F")^2)
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


#' @title Selection of parameter alpha and alpha_2
#'
#' @description
#' Performs cross-validation to select optimal regularization parameters
#' (alpha_1 and alpha_2) for matrix factorization models. This function
#' evaluates a grid of parameter combinations and selects the ones that
#' maximize the cross-validation performance metric.
#'
#' @param Object A structured object containing data and fixed matrices for
#'               the model. Expected to have components:
#'               \itemize{
#'                 \item \code{X}: The main data matrix
#'                 \item \code{A}: Fixed matrix A
#'                 \item \code{L}: Fixed matrix L
#'                 \item \code{D}: Fixed matrix D
#'               }
#' @param K The rank for matrix factorization.
#' @param cross_k Number of folds for k-fold cross-validation. Defaults to 5.
#' @param para_1_list Numeric vector of candidate values for the first
#'                    regularization parameter (alpha_1). Defaults to
#'                    \code{c(0.01, 0.005, 0.001)}.
#' @param para_2_list Numeric vector of candidate values for the second
#'                    regularization parameter (alpha_2). Defaults to
#'                    \code{c(0.01, 0.005, 0.001)}.
#' @param seed Random seed for reproducible cross-validation splitting.
#'             Defaults to 0.
#' @param parallel Logical indicating whether to use parallel processing for
#'                 parameter evaluation. Defaults to \code{FALSE}.
#' @param workers Number of parallel workers to use if \code{parallel = TRUE}.
#'                If \code{NULL}, uses available cores minus one. Defaults to
#'                \code{NULL}.
#' @param verbose Logical indicating whether to print progress messages.
#'                Defaults to \code{TRUE}.
#'
#' @return A list containing:
#'   \itemize{
#'     \item \code{para}: A list with selected parameters:
#'       \itemize{
#'         \item \code{alpha_1}: The optimal first regularization parameter
#'         \item \code{alpha_2}: The optimal second regularization parameter
#'         \item \code{result_cv}: Matrix of cross-validation results for all parameter combinations
#'       }
#'   }
#'
#' @note
#' This function supports both parallel and sequential evaluation to balance
#' computational efficiency and resource usage.
#'
#'
#' @keywords internal
#' @family scAB
#' @family scAB_optimal_param
#'
#'
select_alpha.optimized <- function(
    Object,
    method = c("binary", "survival"),
    K,
    cross_k = 5,
    para_1_list = c(0.01, 0.005, 0.001),
    para_2_list = c(0.01, 0.005, 0.001),
    seed = 0,
    parallel = FALSE,
    workers = getFuncOption("workers"),
    verbose = getFuncOption("verbose")
) {
    if (verbose) {
        ts_cli$cli_alert_info(
            "Selecting optimal {.arg alpha} and {.arg alpha_2}, this would take a while"
        )
    }

    train_phenotype <- PreparePheno(Object)
    train_data <- Object$X

    cvlist <- CVgroup2(k = cross_k, datasize = nrow(train_data), seed = seed)

    train_data_norm <- train_data / Matrix::norm(train_data, "F")
    fixed_matrices <- list(A = Object$A, L = Object$L, D = Object$D)

    param_grid <- expand.grid(
        para_1 = para_1_list,
        para_2 = para_2_list,
        stringsAsFactors = FALSE
    )

    cv_results <- if (parallel) {
        ParallelEvaluate(
            method = method,
            param_grid = param_grid,
            train_data = train_data_norm,
            train_phenotype = train_phenotype,
            cvlist = cvlist,
            fixed_matrices = fixed_matrices,
            K = K,
            cross_k = cross_k,
            workers = workers,
            verbose = verbose
        )
    } else {
        SequentialEvaluate(
            method = method,
            param_grid = param_grid,
            train_data = train_data_norm,
            train_phenotype = train_phenotype,
            cvlist = cvlist,
            fixed_matrices = fixed_matrices,
            K = K,
            cross_k = cross_k,
            verbose = verbose
        )
    }

    result_cv <- matrix(
        cv_results,
        nrow = length(para_1_list),
        ncol = length(para_2_list)
    )

    best_idx <- which(result_cv == max(result_cv), arr.ind = TRUE)[1, ]

    alpha_1 <- para_1_list[best_idx[1]]
    alpha_2 <- para_2_list[best_idx[2]]

    if (verbose) {
        ts_cli$cli_alert_success(
            "Best {.arg alpha} and {.arg alpha_2} are {.val {alpha_1}} and {.val {alpha_2}}"
        )
    }

    list(
        para = list(
            alpha_1 = alpha_1,
            alpha_2 = alpha_2,
            result_cv = result_cv
        )
    )
}

#' @keywords internal
#' @family scAB_optimal_param
PreparePheno <- function(Object) {
    if (Object$method == "survival") {
        return(Object$phenotype)
    }

    data.frame(
        status = as.integer(Object$phenotype),
        time = ifelse(Object$phenotype, 1, 100),
        row.names = rownames(Object$X)
    )
}

#' @keywords internal
#' @family scAB_optimal_param
SequentialEvaluate <- function(
    method = c("binary", "survival"),
    param_grid,
    train_data,
    train_phenotype,
    cvlist,
    fixed_matrices,
    K,
    cross_k,
    verbose = getFuncOption("verbose")
) {
    purrr::map_dbl(
        seq_len(nrow(param_grid)),
        function(i) {
            cv_scores <- vapply(
                seq_len(cross_k),
                function(cv_idx) {
                    EvaluateSingleCV(
                        cv_idx = cv_idx,
                        method = method,
                        train_data = train_data,
                        train_phenotype = train_phenotype,
                        cvlist = cvlist,
                        fixed_matrices = fixed_matrices,
                        K = K,
                        para_1 = param_grid$para_1[i],
                        para_2 = param_grid$para_2[i]
                    )
                },
                numeric(1)
            )
            mean(cv_scores)
        },
        .progress = verbose
    )
}

#' @keywords internal
#' @family scAB_optimal_param
ParallelEvaluate <- function(
    method = c("binary", "survival"),
    param_grid,
    train_data,
    train_phenotype,
    cvlist,
    fixed_matrices,
    K,
    cross_k,
    workers = getFuncOption("workers"),
    verbose = getFuncOption("verbose")
) {
    plan(getFuncOption("parallel.type"), workers = workers)
    on.exit(plan("sequential"))

    if (verbose) {
        ts_cli$cli_alert_info(sprintf(
            "Using parallel processing with %d workers",
            workers
        ))
    }
    EvaluateSingleCV <- SigBridgeR:::EvaluateSingleCV
    ginv2 <- SigBridgeR:::ginv2
    ginv2.default <- SigBridgeR:::ginv2.default
    guanrank2 <- SigBridgeR:::guanrank2
    scAB.optimized <- SigBridgeR:::scAB.optimized
    # Pkg carrier is used to pass arguments and avoid used large objects in the closure
    ScoringAll <- function(i) {
        cv_scores <- vapply(
            seq_len(cross_k),
            function(cv_idx) {
                EvaluateSingleCV(
                    cv_idx = cv_idx,
                    method = method,
                    train_data = train_data,
                    train_phenotype = train_phenotype,
                    cvlist = cvlist,
                    fixed_matrices = fixed_matrices,
                    K = K,
                    para_1 = param_grid$para_1[i],
                    para_2 = param_grid$para_2[i]
                )
            },
            numeric(1)
        )

        mean(cv_scores)
    }

    res <- future_map_dbl(
        seq_len(nrow(param_grid)),
        ScoringAll,
        .progress = verbose,
        .options = furrr_options(
            seed = getFuncOption("seed"),
            packages = c("survival", "SigBridgeR"),
            globals = c(
                "train_data",
                "train_phenotype",
                "cvlist",
                "fixed_matrices",
                "K",
                "cross_k",
                "param_grid",
                "EvaluateSingleCV",
                "guanrank2",
                "ginv2",
                "ginv2.default",
                "scAB.optimized"
            )
        )
    )

    res
}

#' @keywords internal
#' @family scAB_optimal_param
EvaluateSingleCV <- function(
    cv_idx,
    method,
    train_data,
    train_phenotype,
    cvlist,
    fixed_matrices,
    K,
    para_1,
    para_2
) {
    test_idx <- cvlist[[cv_idx]]
    train_subset <- train_data[-test_idx, , drop = FALSE]
    test_subset <- train_data[test_idx, , drop = FALSE]
    train_pheno <- train_phenotype[-test_idx, ]
    test_pheno <- train_phenotype[test_idx, ]

    ss <- guanrank2(train_pheno[, c("time", "status")])
    S <- diag(1 - ss[rownames(train_pheno), 3])

    Object_cv <- structure(
        list(
            X = train_subset,
            S = S,
            phenotype = train_pheno,
            A = fixed_matrices$A,
            L = fixed_matrices$L,
            D = fixed_matrices$D,
            method = method
        ),
        class = "scAB_data"
    )

    s_res <- scAB.optimized(
        Object = Object_cv,
        K = K,
        alpha = para_1,
        alpha_2 = para_2,
        maxiter = 2000
    )

    ginvH <- ginv2(s_res$H)
    new_W <- test_subset %*% ginvH

    clin_km <- data.frame(
        time = train_pheno$time,
        status = train_pheno$status,
        s_res$W
    )

    res.cox <- survival::coxph(survival::Surv(time, status) ~ ., data = clin_km)
    pre_test <- stats::predict(res.cox, data.frame(new_W))

    survival::concordance(
        survival::coxph(
            survival::Surv(test_pheno$time, test_pheno$status) ~ pre_test
        )
    )$concordance
}

#' @title Create subsets of cross-validation
#'
#' @param k  k-fold cross validation
#' @param datasize  the size of samples
#' @param seed random seed
#'
#' @return a list with subsets of cross-validation
#' @keywords internal
#' @family scAB
#' @examples
#' \dontrun{
#' CVgroup(10,10)
#' }
#'
#'
CVgroup2 <- function(k, datasize, seed = 0) {
    set.seed(seed)
    folds <- sample(rep_len(seq_len(k), datasize))
    split(seq_len(datasize), folds)
}
