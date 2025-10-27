#' @title Convert a column to rownames
#' @keywords internal
Col2Rownames <- function(.data, var = "rowname") {
    rownames(.data) <- .data[[var]]
    .data[[var]] <- NULL
    .data
}

#' @title Convert rownames to a column
#' @keywords internal
Rownames2Col <- function(.data, var = "rowname") {
    rownames_col <- data.frame(rownames(.data))
    names(rownames_col) <- var
    rownames(.data) <- NULL
    .data <- cbind(rownames_col, .data)
    .data
}

#' @keywords internal
`%||%` <- function(lhs, rhs) {
    if (is.null(lhs)) {
        return(rhs)
    }
    return(lhs)
}

#' @keywords internal
all_identical <- function(..., names = NULL) {
    objs <- list(...)
    n <- length(objs)

    if (is.null(names)) {
        names <- paste0("Obj", seq_len(n))
    } else if (length(names) != n) {
        cli::cli_abort(c(
            "x" = "The length of `names` must be equal to the number of objects."
        ))
    }

    result <- matrix(TRUE, n, n, dimnames = list(names, names))

    if (n > 1) {
        for (i in seq_len(n - 1)) {
            for (j in (i + 1):n) {
                is_identical <- identical(objs[[i]], objs[[j]])
                result[i, j] <- is_identical
                result[j, i] <- is_identical # 对称矩阵
            }
        }
    }

    return(result)
}
