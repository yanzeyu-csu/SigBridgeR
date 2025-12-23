#' @title Convert Ensembles Version IDs & TCGA Version IDs to Genes in Bulk Expression Data
#'
#' @description
#' Preprocess bulk expression data: convert Ensembles version IDs and TCGA version IDs to genes. NA values are replaced with `unknown_k` format (k stands for the position of the NA value in the row).
#'
#' @param data bulk expression data (matrix or data.frame)
#' @param unknown_format A glue pattern containing `{k}` for replace the NA value during conversion.
#'     k must be wrapped in curly braces, stands for the position of the NA value in the row.
#'     Default: `"unknown_{k}"`.
#'
#' @export
#'
SymbolConvert <- function(data, unknown_format = "unknown_{k}") {
  row_names <- gsub("\\..*$", "", rownames(data))
  if (is.null(row_names)) {
    cli::cli_abort(c("x" = "Row names are missing in the data"))
  }
  options(
    IDConverter.datapath = system.file("extdata", package = "IDConverter")
  )
  gene_symbols <- IDConverter::convert_hm_genes(row_names)

  na_count <- sum(is.na(gene_symbols))
  if (na_count > 0) {
    cli::cli_warn(c(
      "Found {.val {na_count}} NA values in gene symbols during conversion."
    ))

    k <- which(is.na(gene_symbols))
    gene_symbols[k] <- glue::glue(unknown_format)
    cli::cli_warn(
      "Replaced {.val {na_count}} NA values with {.code {unknown_format}} format."
    )
  }

  rownames(data) <- gene_symbols
  data
}
