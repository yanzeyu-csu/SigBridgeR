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
    quantile(x, probs = probs, ...)
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

    # 核心算法：完全向量化

    # 1. 对每列排序
    sorted_mat <- apply(mat, 2, sort, na.last = TRUE)

    # 2. 计算目标分布（每行的均值）
    target_dist <- rowMeans(sorted_mat, na.rm = TRUE)

    # 3. 创建排名矩阵（完全向量化）
    # 对整个矩阵进行排名操作
    rank_mat <- apply(mat, 2, function(col) {
        match(col, sort(col, na.last = NA))
    })

    # 4. 使用矩阵索引一次性分配所有值
    # 为每个元素分配对应排名的目标分布值
    for (j in 1:cols) {
        valid_idx <- !is.na(rank_mat[, j])
        mat[valid_idx, j] <- target_dist[rank_mat[valid_idx, j]]
    }

    # 5. 恢复NA值
    mat[na_positions] <- NA

    # 恢复行列名
    if (keep.names) {
        rownames(mat) <- orig_rownames
        colnames(mat) <- orig_colnames
    }

    mat
}
