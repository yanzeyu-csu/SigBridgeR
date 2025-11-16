#' @title Check for Rows with Zero Variance in a Matrix
#'
#' @description
#' This function checks each row of a matrix (including sparse matrices of class `dgCMatrix`)
#' for zero variance. Rows with zero variance or only `NA` values are identified, and an error
#' is thrown listing the names of these rows. This is useful for preprocessing data where
#' constant rows may cause issues in analyses (e.g., PCA, regression).
#'
#' @param mat A numeric matrix or a sparse matrix of class `dgCMatrix` (from the `Matrix` package).
#'   Rows represent features (e.g., genes), and columns represent observations.
#' @param call The environment from which the function was called, used for error reporting.
#'   Defaults to rlang::caller_env(). Most users can ignore this parameter.
#'
#' @return Invisibly returns a numeric vector of row variances. If zero-variance rows are found,
#'   the function throws an error with a message listing the problematic row names.
#'
#' @details
#' For dense matrices, variance is computed using `rowVars()` function that
#' efficiently calculates row variances with proper NA handling. For sparse matrices of
#' class `dgCMatrix`, variance is computed using a mathematical identity that avoids
#' creating large intermediate matrices. Rows with fewer than 2 non-zero observations
#' are treated as zero-variance.
#'
#' @examples
#' \dontrun{
#' # Dense matrix example
#' set.seed(123)
#' mat_dense <- matrix(rnorm(100), nrow = 10)
#' rownames(mat_dense) <- paste0("Gene", 1:10)
#' Check0VarRows(mat_dense) # No error if all rows have variance
#'
#' # Introduce zero variance
#' mat_dense[1, ] <- rep(5, 10) # First row is constant
#' Check0VarRows(mat_dense)     # Throws error listing "Gene1"
#'
#' # Sparse matrix example
#' library(Matrix)
#' mat_sparse <- as(matrix(rpois(100, 0.5), nrow = 10), "dgCMatrix")
#' rownames(mat_sparse) <- paste0("Gene", 1:10)
#' Check0VarRows(mat_sparse)
#' }
#'
#'
#' @keywords internal
#' @family scPP
#'
Check0VarRows <- function(mat, call = rlang::caller_env()) {
    if (!inherits(mat, c("matrix", "Matrix"))) {
        cli::cli_warn("Input data was coerced to a matrix.")
        mat <- Matrix::Matrix(as.matrix(mat))
    }
    if (inherits(mat, "Matrix")) {
        n <- ncol(mat)
        row_sums <- Matrix::rowSums(mat)
        row_sums_sq <- Matrix::rowSums(mat^2)

        # Var = (Σx² - (Σx)²/n) / (n-1)
        vars <- (row_sums_sq - (row_sums^2) / n) / (n - 1)

        zero_counts <- Matrix::rowSums(mat != 0)
        vars[zero_counts <= 1] <- 0
        vars[is.nan(vars)] <- 0
    } else {
        vars <- rowVars(mat, na.rm = TRUE)
    }

    zero_idx <- which(vars == 0 | is.na(vars))
    if (length(zero_idx)) {
        bad_genes <- rownames(mat)[zero_idx]
        cli::cli_abort(
            c(
                "Detected {.val {length(bad_genes)}} gene(s) with zero variance:",
                "i" = "{.val {bad_genes}}"
            ),
            class = "ZeroVarianceError",
            call = call
        )
    }
    invisible(vars)
}


#' @keywords internal
rowVars <- function(x, na.rm = FALSE, ...) {
    if (rlang::is_installed("matrixStats")) {
        return(getExportedValue("matrixStats", "rowVars")(
            x,
            na.rm = na.rm,
            ...
        ))
    }

    if (na.rm) {
        n <- rowSums(!is.na(x))
        rowSums((x - rowMeans(x, na.rm = TRUE))^2, na.rm = TRUE) / (n - 1)
    } else {
        rowSums((x - rowMeans(x))^2) / (ncol(x) - 1)
    }
}

#' @keywords internal
colVars <- function(x, na.rm = FALSE, ...) {
    if (rlang::is_installed("matrixStats")) {
        return(getExportedValue("matrixStats", "colVars")(
            x,
            na.rm = na.rm,
            ...
        ))
    }
    if (na.rm) {
        n <- colSums(!is.na(x))
        colSums((x - colMeans(x, na.rm = TRUE))^2, na.rm = TRUE) / (n - 1)
    } else {
        colSums((x - colMeans(x))^2) / (nrow(x) - 1)
    }
}


#' @keywords internal
rowSds <- function(x, na.rm = FALSE, ...) {
    if (rlang::is_installed("matrixStats")) {
        return(getExportedValue("matrixStats", "rowSds")(x, na.rm = na.rm, ...))
    }
    sqrt(rowVars(x, na.rm = na.rm, ...))
}

#' @keywords internal
colSds <- function(x, na.rm = FALSE, ...) {
    if (rlang::is_installed("matrixStats")) {
        return(getExportedValue("matrixStats", "colSds")(x, na.rm = na.rm))
    }
    sqrt(colVars(x, na.rm = na.rm, ...))
}

#' @keywords internal
colQuantiles <- function(x, probs = seq(0, 1, 0.25), ...) {
    if (rlang::is_installed("matrixStats")) {
        return(getExportedValue("matrixStats", "colQuantiles")(
            x,
            probs = probs,
            ...
        ))
    }
    stats::quantile(x, probs = probs, ...)
}

#' @title Quantile normalizes the rows of a matrix.
#'
#' @description
#' If the package is installed, it will use the `normalize.quantiles()` function from `preprocessCore`,
#' otherwise it will use the R implementation, which is a bit slower.
#'
#' @keywords internal
normalize.quantiles <- function(x, copy = TRUE, keep.names = FALSE, ...) {
    if (rlang::is_installed("preprocessCore")) {
        return(getExportedValue("preprocessCore", "normalize.quantiles")(
            x,
            copy,
            keep.names
        ))
    }

    if (!is.matrix(x)) {
        cli::cli_abort(c(
            "x" = "Matrix expected in normalize.quantiles",
            ">" = "Input is a {.field {class(x)}}"
        ))
    }
    rows <- nrow(x)
    cols <- ncol(x)

    if (copy) {
        mat <- matrix(as.numeric(x), rows, cols)
    } else {
        mat <- x
        if (is.integer(mat)) {
            mat <- matrix(as.double(mat), rows, cols)
        }
    }

    if (keep.names) {
        orig_rownames <- rownames(x)
        orig_colnames <- colnames(x)
    }

    na_positions <- is.na(mat)

    sorted_mat <- apply(mat, 2, sort, na.last = TRUE)

    target_dist <- rowMeans(sorted_mat, na.rm = TRUE)

    rank_mat <- apply(mat, 2, function(col) {
        match(col, sort(col, na.last = NA))
    })

    for (j in seq_len(cols)) {
        valid_idx <- !is.na(rank_mat[, j])
        mat[valid_idx, j] <- target_dist[rank_mat[valid_idx, j]]
    }

    mat[na_positions] <- NA

    if (keep.names) {
        rownames(mat) <- orig_rownames
        colnames(mat) <- orig_colnames
    }

    mat
}


#' @title Generalized Inverse (Moore-Penrose Inverse)
#' @description
#' Compute the Moore-Penrose generalized inverse of a matrix.
#' This is an S3 generic function with methods for base matrices,
#' dense Matrix objects, and sparse Matrix objects.
#'
#' @param X A numeric or complex matrix
#' @param tol Tolerance for determining rank. Default is sqrt(.Machine$double.eps)
#' @param ... Additional arguments passed to methods
#'
#' @return The generalized inverse of X
#'
#' @examples
#' \dontrun{
#' # Base R matrix
#' m <- matrix(c(1, 2, 3, 4, 5, 6), 2, 3)
#' ginv2(m)
#'
#' # Dense Matrix
#' library(Matrix)
#' dm <- Matrix(m, sparse = FALSE)
#' ginv2(dm)
#'
#' # Sparse Matrix
#' sm <- Matrix(m, sparse = TRUE)
#' ginv2(sm)
#' }
#' @keywords internal
#'
ginv2 <- function(X, tol = sqrt(.Machine$double.eps), ...) {
    # if (inherits(X, "sparseMatrix")) {
    #     return(ginv2.sparseMatrix(X, tol = tol, ...))
    # }

    if (inherits(X, "Matrix")) {
        return(Matrix::Matrix(ginv2.default(X, tol = tol, ...)))
    }

    # Base R matrix or coercible to matrix
    if (length(dim(X)) > 2L || !(is.numeric(X) || is.complex(X))) {
        cli::cli_abort(c("x" = "'X' must be a numeric or complex matrix"))
    }
    UseMethod("ginv2")
}


#' @title Generalized Inverse for Base R Matrices
#' @description Default method for base R matrices (from MASS::ginv)
#' @rdname ginv2
#' @keywords internal
ginv2.default <- function(X, tol = sqrt(.Machine$double.eps), ...) {
    if (!is.matrix(X)) {
        X <- as.matrix(X)
    }

    Xsvd <- svd(X)

    if (is.complex(X)) {
        Xsvd$u <- Conj(Xsvd$u)
    }

    Positive <- Xsvd$d > max(tol * Xsvd$d[1L], 0)

    if (all(Positive)) {
        Xsvd$v %*% (1 / Xsvd$d * t(Xsvd$u))
    } else if (!any(Positive)) {
        array(0, dim(X)[c(2L, 1L)])
    } else {
        Xsvd$v[, Positive, drop = FALSE] %*%
            ((1 / Xsvd$d[Positive]) *
                t(Xsvd$u[, Positive, drop = FALSE]))
    }
}

# #' @title Generalized Inverse for Sparse Matrix Objects
# #' @description Method for all sparse Matrix package objects
# #' @param method Method to use: "svd" (default) or "qr"
# #' @param return_sparse Whether to return sparse matrix if appropriate
# #' @rdname ginv2
# #' @keywords internal
# ginv2.sparseMatrix <- function(
#     X,
#     tol = sqrt(.Machine$double.eps),
#     method = c("auto", "svd", "qr"),
#     return_sparse = TRUE,
#     ...
# ) {
#     method <- match.arg(method)
#     if (method == "auto") {
#         density <- Matrix::nnzero(X) / length(X)
#         if (density < 0.3 && ncol(X) <= 1000) {
#             method <- "qr"
#         } else {
#             method <- "svd"
#         }
#     }

#     if (method == "qr") {
#         if (!inherits(X, "CsparseMatrix")) {
#             X <- as(X, "CsparseMatrix")
#         }

#         qr_x <- Matrix::qr(X)
#         R <- Matrix::qr.R(qr_x)

#         d <- abs(Matrix::diag(R))
#         rank <- sum(d > tol * max(d))

#         if (rank < ncol(X)) {
#             cli::cli_warn(
#                 "Matrix is rank deficient, using regularized solution"
#             )
#             XtX <- Matrix::crossprod(X)
#             reg_param <- tol * mean(Matrix::diag(XtX))
#             XtX_reg <- XtX + reg_param * Matrix::Diagonal(ncol(X))
#             result <- Matrix::solve(XtX_reg, Matrix::t(X))
#         } else {
#             result <- Matrix::solve(R, Matrix::t(Matrix::qr.Q(qr_x)))
#         }
#     } else {
#         # SVD
#         X_dense <- as.matrix(X)
#         svd_result <- svd(X_dense)
#         d <- svd_result$d
#         keep <- d > tol * max(d)
#         d_inv <- rep(0, length(d))
#         d_inv[keep] <- 1 / d[keep]

#         result <- svd_result$v[, keep, drop = FALSE] %*%
#             (d_inv[keep] * t(svd_result$u[, keep, drop = FALSE]))
#     }

#     if (return_sparse) {
#         result_max <- max(abs(result))
#         if (result_max > 0) {
#             zero_threshold <- tol * result_max
#             result[abs(result) < zero_threshold] <- 0
#         }
#         return(Matrix::drop0(Matrix::Matrix(result, sparse = TRUE)))
#     }

#     as.matrix(result)
# }
