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
#' For dense matrices, variance is computed using an optimized `RowVars` function that
#' efficiently calculates row variances with proper NA handling. For sparse matrices of
#' class `dgCMatrix`, variance is computed using a mathematical identity that avoids
#' creating large intermediate matrices. Rows with fewer than 2 non-zero observations
#' are treated as zero-variance.
#'
#' This implementation is memory-efficient and handles large matrices better than
#' the original version.
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
#' @importFrom Matrix rowSums
#' @importFrom rlang caller_env
#' @importFrom cli cli_abort
#'
#' @keywords internal
#' @family scPP
#'
Check0VarRows <- function(mat, call = rlang::caller_env()) {
    if (inherits(mat, "dgCMatrix")) {
        n <- ncol(mat)
        row_sums <- rowSums(mat)
        row_sums_sq <- rowSums(mat^2)

        # Var = (Σx² - (Σx)²/n) / (n-1)
        vars <- (row_sums_sq - (row_sums^2) / n) / (n - 1)

        zero_counts <- rowSums(mat != 0)
        vars[zero_counts <= 1] <- 0
        vars[is.nan(vars)] <- 0
    } else {
        vars <- RowVars(mat, na.rm = TRUE)
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

#' @title Calculate Row Variances
#'
#' @description
#' This function calculates the sample variance for each row of a numeric matrix
#' or data frame. It uses the standard sample variance formula with Bessel's
#' correction (dividing by n-1).
#'
#' @param x A numeric matrix or data frame for which row variances are to be
#'   calculated.
#' @param na.rm A logical value indicating whether missing values (NA) should be
#'   removed before calculating variances. If \code{TRUE} (default), missing
#'   values are excluded and the denominator is adjusted accordingly. If
#'   \code{FALSE}, any row containing missing values will result in \code{NA}.
#'
#' @return A numeric vector of length equal to the number of rows in \code{x},
#'   containing the sample variance for each row.
#'
#' @details
#' The function computes sample variances using the formula:
#' \deqn{\sigma^2 = \frac{\sum_{i=1}^{n}(x_i - \bar{x})^2}{n-1}}
#'
#' When \code{na.rm = TRUE}, the function:
#' \itemize{
#'   \item Excludes missing values from mean and variance calculations
#'   \item Adjusts the degrees of freedom (n-1) based on the number of
#'     non-missing values in each row
#'   \item Returns \code{NaN} for rows with fewer than 2 non-missing values
#' }
#'
#' When \code{na.rm = FALSE}, the function:
#' \itemize{
#'   \item Uses all values including missing ones
#'   \item Returns \code{NA} for any row containing missing values
#'   \item Uses \code{ncol(x) - 1} as the degrees of freedom
#' }
#'
#' @examples
#' \dontrun{
#' # Basic usage with a matrix
#' mat <- matrix(1:12, nrow = 3)
#' rowVars(mat)
#'
#' # With missing values
#' mat[1, 2] <- NA
#' mat[2, 3] <- NA
#' rowVars(mat, na.rm = TRUE)   # Excludes NAs
#' rowVars(mat, na.rm = FALSE)  # Includes NAs (returns NA for affected rows)
#'
#' # With a data frame
#' df <- data.frame(
#'   a = c(1, 4, 7),
#'   b = c(2, 5, 8),
#'   c = c(3, 6, 9)
#' )
#' rowVars(df)
#'
#' # Edge case: single column (variance is 0)
#' single_col <- matrix(1:3, ncol = 1)
#' rowVars(single_col)  # Returns NaN due to division by 0
#' }
#'
#' @seealso
#' \code{\link{var}} for column-wise variance calculation,
#' \code{\link{rowMeans}} for row means,
#' \code{\link{apply}} for applying functions across rows or columns
#'
#' @keywords internal
#'
RowVars <- function(x, na.rm = TRUE) {
    if (na.rm) {
        n <- rowSums(!is.na(x))
        rowSums((x - rowMeans(x, na.rm = TRUE))^2, na.rm = TRUE) / (n - 1)
    } else {
        rowSums((x - rowMeans(x))^2) / (ncol(x) - 1)
    }
}
