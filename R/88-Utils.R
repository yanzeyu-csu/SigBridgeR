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
