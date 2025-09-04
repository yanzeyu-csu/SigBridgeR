# * ------------------- other function ----------------------

#' @title Safely Add Miscellaneous Data to Seurat Object
#'
#' @description
#' Adds arbitrary data to the `@misc` slot of a Seurat object with automatic key
#' conflict resolution. If the key already exists, automatically appends a numeric
#' suffix to ensure unique key naming (e.g., "mykey_1", "mykey_2").
#'
#' @usage
#' AddMisc(
#'   seurat_obj,
#'   ..., # key = value
#'   cover = TRUE # overwrite existing data
#' )
#'
#' @param seurat_obj A Seurat object to modify
#' @param ... key-value pairs to add to the `@misc` slot.
#' @param cover Logical indicating whether to overwrite existing data. If
#'         (default TRUE).
#'
#' @return The modified Seurat object with added `@misc` data. The original object
#'         structure is preserved with no other modifications.
#'
#' @section Key Generation Rules:
#' 1. If `key` doesn't exist: uses as-is
#' 2. If `key` exists: appends the next available number (e.g., "key_1", "key_2")
#' 3. If numbered keys exist (e.g., "key_2"): increments the highest number
#'
#' @examples
#' \dontrun{
#' # Basic usage
#' seurat_obj <- AddMisc(seurat_obj, "QC_stats" = qc_df)
#'
#' # Auto-incrementing example
#' seurat_obj <- AddMisc(seurat_obj, markers = markers1)
#' seurat_obj <- AddMisc(seurat_obj, markers = markers2, cover=FALSE)
#' # Stores as "markers" and "markers_1"
#'
#' }
#'
#'
#' @export
#'
AddMisc <- function(seurat_obj, ..., cover = TRUE) {
    chk::chk_is(seurat_obj, "Seurat")
    chk::chk_flag(cover)
    # Get the key-value pairs from ... arguments
    kv_pairs <- list(...)

    for (key in names(kv_pairs)) {
        value <- kv_pairs[[key]]

        if (key %in% names(seurat_obj@misc) && !cover) {
            pattern <- glue::glue("^{key}(_\\d+)?$")
            existing_keys <- grep(
                pattern,
                names(seurat_obj@misc),
                value = TRUE
            )
            if (length(existing_keys) > 0) {
                nums <- sapply(existing_keys, function(k) {
                    if (k == key) {
                        return(0)
                    }
                    num_str <- sub(glue::glue("^{key}_(\\d+)$"), "\\1", k)
                    if (grepl("^\\d+$", num_str)) as.integer(num_str) else 0
                })

                nums <- max(nums, na.rm = TRUE)
                key <- if (length(nums) > 0) {
                    glue::glue("{key}_{nums + 1}")
                } else {
                    glue::glue("{key}_1")
                }
            } else {
                key <- glue::glue("{key}_1")
            }
        }
        seurat_obj@misc[[key]] <- value
    }

    return(seurat_obj)
}
