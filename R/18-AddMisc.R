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
#' @importFrom data.table `%chin%`
#' @export
#'
AddMisc <- function(seurat_obj, ..., cover = TRUE) {
    chk::chk_is(seurat_obj, "Seurat")
    chk::chk_flag(cover)
    # Get the key-value pairs from ... arguments
    dots <- list(...)
    if (length(dots) == 0L) {
        return(seurat_obj)
    }

    kv_pairs <- if (
        length(dots) == 1L && is.list(dots[[1L]]) && is.null(names(dots))
    ) {
        dots[[1L]] # A list containing key-value pairs, fetch the key-value pairs in the list
    } else {
        dots # fetch the key-value pairs directly
    }

    if (length(kv_pairs) == 0L) {
        return(seurat_obj)
    }

    input_keys <- names(kv_pairs)
    n_input <- length(kv_pairs)

    # Fast path: cover=TRUE, direct assignment
    if (cover) {
        for (i in seq_len(n_input)) {
            seurat_obj@misc[[input_keys[i]]] <- kv_pairs[[i]]
        }
        return(seurat_obj)
    }

    misc_names <- names(seurat_obj@misc)
    n_misc <- length(misc_names)

    # Build suffix lookup using data.table
    max_suffix_lookup <- NULL

    if (n_misc > 0L) {
        underscore_pos <- regexpr("_[0-9]+$", misc_names, perl = TRUE)
        has_suffix <- underscore_pos > 0L

        base_keys <- character(n_misc)
        suffixes <- integer(n_misc)

        base_keys[has_suffix] <- substr(
            misc_names[has_suffix],
            1L,
            underscore_pos[has_suffix] - 1L
        )
        base_keys[!has_suffix] <- misc_names[!has_suffix]

        suffixes[has_suffix] <- as.integer(substr(
            misc_names[has_suffix],
            underscore_pos[has_suffix] + 1L,
            nchar(misc_names[has_suffix])
        ))

        dt <- data.table::data.table(base_key = base_keys, suffix = suffixes)
        data.table::setkey(dt, base_key)
        max_suffix_lookup <- dt[, .(max_suffix = max(suffix)), by = base_key]
        data.table::setkey(max_suffix_lookup, base_key)
    }

    # Track keys added in this batch
    local_added <- character(0L)
    local_max_suffix <- integer(0L)

    # Resolve conflicts and assign values
    for (i in seq_len(n_input)) {
        key <- input_keys[i]

        if (!key %chin% misc_names && !key %chin% local_added) {
            seurat_obj@misc[[key]] <- kv_pairs[[i]]
            local_added <- c(local_added, key)
            local_max_suffix <- c(local_max_suffix, 0L)
            next
        }

        max_num <- if (!is.null(max_suffix_lookup)) {
            result <- max_suffix_lookup[.(key), max_suffix, nomatch = NA]
            if (!is.na(result)) result else NULL
        } else {
            NULL
        }

        # Check for conflicts in locally added keys
        if (length(local_added) > 0L) {
            local_match_idx <- which(startsWith(local_added, key))

            if (length(local_match_idx) > 0L) {
                key_len <- nchar(key)
                local_matches <- local_added[local_match_idx]

                local_suffixes <- vapply(
                    local_matches,
                    function(k) {
                        if (k == key) {
                            return(0L)
                        }
                        if (
                            nchar(k) > key_len + 1L &&
                                substr(k, key_len + 1L, key_len + 1L) == "_"
                        ) {
                            num_str <- substr(k, key_len + 2L, nchar(k))
                            if (grepl("^[0-9]+$", num_str, perl = TRUE)) {
                                return(as.integer(num_str))
                            }
                        }
                        return(-1L)
                    },
                    integer(1L)
                )

                local_max <- max(
                    local_suffixes[local_suffixes >= 0L],
                    na.rm = TRUE
                )
                if (is.finite(local_max)) {
                    max_num <- if (is.null(max_num)) {
                        local_max
                    } else {
                        max(max_num, local_max)
                    }
                }
            }
        }

        # Generate new key with suffix
        new_key <- if (is.null(max_num)) {
            paste0(key, "_1")
        } else {
            paste0(key, "_", max_num + 1L)
        }

        seurat_obj@misc[[new_key]] <- kv_pairs[[i]]

        # Update local tracking
        local_added <- c(local_added, new_key)
        local_max_suffix <- c(
            local_max_suffix,
            if (is.null(max_num)) 1L else max_num + 1L
        )
    }

    seurat_obj
}
