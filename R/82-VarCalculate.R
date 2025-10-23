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
#' For dense matrices, variance is computed using `matrixStats::rowVars()` function that
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
    if (inherits(mat, "dgCMatrix")) {
        n <- ncol(mat)
        row_sums <- Matrix::rowSums(mat)
        row_sums_sq <- Matrix::rowSums(mat^2)

        # Var = (Σx² - (Σx)²/n) / (n-1)
        vars <- (row_sums_sq - (row_sums^2) / n) / (n - 1)

        zero_counts <- Matrix::rowSums(mat != 0)
        vars[zero_counts <= 1] <- 0
        vars[is.nan(vars)] <- 0
    } else {
        vars <- matrixStats::rowVars(mat, na.rm = TRUE)
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
