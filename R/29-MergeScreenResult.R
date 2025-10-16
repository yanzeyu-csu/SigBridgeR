#' @title Merge Multiple Screening Analysis Results
#'
#' @description
#' Combines results from multiple single-cell screening analyses (Scissor, scPAS, scPP, or scAB)
#' by merging their metadata and miscellaneous information while preserving the original expression data. Performs an inner join
#' on cell barcodes to ensure only cells present in all inputs are retained.
#'
#' @usage
#' MergeResult(...)
#'
#' @param ... Input objects to merge. Can be:
#'        - Seurat objects
#'        - Lists containing `scRNA_data` (Seurat objects)
#'        - Mixed combinations of the above
#'        - The first one will be used as base object for merging, priority given to first one when duplicate columns are found
#'
#' @return A merged Seurat object containing:
#' \itemize{
#'   \item Expression data from the first input object
#'   \item Combined metadata from all input objects
#'   \item Miscellaneous information from all input objects
#'   \item Only cells present in all input objects (inner join)
#' }
#'
#' @section Processing Details:
#' 1. Input Validation: Checks for valid Seurat objects or lists containing Seurat objects
#' 2. Metadata Extraction: Collects metadata from all objects
#' 3. Cell Intersection: Retains only cells present in all datasets
#' 4. Object Merging: Creates new Seurat object with combined metadata
#' 5. Miscellaneous: Adds miscellaneous information to the merged object
#'
#' @examples
#' \dontrun{
#' # Merge mixed analysis types
#' combined <- MergeResult(scissor_output, scAB_output, scPP_output)
#'
#' # Merge list-containing objects
#' merged_list <- MergeResult(list1, list2, seurat_obj)
#' }
#'
#' @export
#' @importFrom tibble rownames_to_column column_to_rownames
#' @importFrom cli cli_abort cli_alert_warning cli_alert_success
#' @importFrom data.table merge.data.table
#' @importFrom purrr map
#'
MergeResult <- function(...) {
    args <- list(...)

    if (length(args) == 0) {
        cli::cli_abort(
            c(
                "x" = "Input objects must be provided."
            ),
            class = "InputsNotFound"
        )
    }

    seurat_objects <- list()
    seurat_objects <- lapply(args, function(x) {
        if (inherits(x, "Seurat")) {
            return(x)
        } else if (is.list(x) && inherits(x$scRNA_data, "Seurat")) {
            return(x$scRNA_data)
        } else {
            cli::cli_warn(
                "Skipping object of class {.code {class(x)}} - not a Seurat object or a list with Seurat object"
            )
            return(NULL)
        }
    }) %>%
        Filter(Negate(is.null), .)

    if (length(seurat_objects) == 0) {
        cli::cli_abort(c("x" = "No valid Seurat objects found in inputs."))
    }

    # extract metadata
    meta_list <- lapply(seurat_objects, function(x) {
        data.table::as.data.table(x@meta.data, keep.rownames = "cell_id")
    })

    merged_meta <- Reduce(
        function(x, y) {
            duplicate_cols <- setdiff(intersect(names(x), names(y)), "cell_id")

            if (length(duplicate_cols) > 0) {
                y <- y[, !..duplicate_cols]
            }

            data.table::merge.data.table(x, y, by = "cell_id", all = FALSE)
        },
        meta_list
    )

    merged_obj <- seurat_objects[[1]]
    merged_obj@meta.data <- tibble::column_to_rownames(merged_meta, "cell_id")

    # merge misc
    all_keys <- unique(unlist(lapply(seurat_objects, function(obj) {
        if (!is.null(obj@misc)) names(obj@misc) else character(0)
    })))

    misc_list <- stats::setNames(vector("list", length(all_keys)), all_keys)

    for (key in all_keys) {
        values <- lapply(seurat_objects, function(obj) {
            if (!is.null(obj@misc) && key %chin% names(obj@misc)) {
                obj@misc[[key]]
            } else {
                NULL
            }
        })

        values <- values[!vapply(values, is.null, FUN.VALUE = logical(1))]

        misc_list[[key]] <- if (length(values) == 1) values[[1]] else values
    }

    merged_obj <- AddMisc(merged_obj, misc_list, cover = TRUE)

    cli::cli_alert_success(
        "Successfully merged {.val {length(seurat_objects)}} objects."
    )

    return(merged_obj)
}
