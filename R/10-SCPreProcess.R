# * ---- scRNA-seq preprocessing ----

#' @title Single-Cell RNA-seq Preprocessing Pipeline
#'
#' @description
#' A generic function for standardized preprocessing of single-cell RNA-seq data
#' from multiple sources. Handles data.frame/matrix, AnnData, and Seurat inputs
#' with tumor cell filtering. Implements a complete analysis pipeline
#' from raw data to clustered embeddings.
#'
#' @name SCPreProcess
#' @usage
#' SCPreProcess(sc, ...)
#'
#' @param sc Input data, one of:
#'    - `data.frame/matrix/dgCMatrix`: Raw count matrix (features x cells)
#'    - `AnnDataR6`: Python AnnData object via reticulate
#'    - `Seurat`: Preprocessed Seurat object
#' @param ... Method-specific arguments (see below)
#'
#' @return A Seurat object containing:
#' \itemize{
#'   \item Data filter and quality control
#'   \item Normalized and scaled expression data
#'   \item Variable features
#'   \item PCA/tSNE/UMAP reductions
#'   \item Cluster identities
#'   \item When tumor cells filtered: original dimensions in `@misc$raw_dim`
#'   \item Final dimensions in `@misc$self_dim`
#' }
NULL


#' @rdname SCPreProcess
#' @export
SCPreProcess <- function(sc, ...) {
    UseMethod("SCPreProcess")
}

#' @rdname SCPreProcess
#' @export
SCPreProcess.default <- function(
    sc,
    meta_data = NULL,
    column2only_tumor = NULL,
    project = glue::glue("SC_Screening_Proj"),
    min_cells = 400L,
    min_features = 0L,
    quality_control = TRUE,
    quality_control.pattern = c("^MT-", "^mt-"),
    data_filter = TRUE,
    data_filter.nFeature_RNA_thresh = c(200L, 6000L),
    data_filter.percent.mt = 20L,
    normalization_method = "LogNormalize",
    scale_factor = 10000L,
    scale_features = NULL,
    selection_method = "vst",
    resolution = 0.6,
    dims = seq_len(10L),
    verbose = TRUE,
    ...
) {
    if (!is.null(meta_data)) {
        chk::chk_data(meta_data)
    }
    if (!is.null(column2only_tumor)) {
        chk::chk_character(column2only_tumor)
    }

    sc_seurat <- Seurat::CreateSeuratObject(
        counts = sc,
        project = project,
        meta.data = meta_data,
        min.cells = min_cells,
        min.features = min_features
    )

    if (quality_control) {
        if (length(quality_control.pattern) != 1) {
            quality_control.pattern %<>%
                MatchArg(choices = c("^MT-", "^mt-"), default = "^MT-")
        }

        chk::chk_character(quality_control.pattern)

        sc_seurat[["percent.mt"]] <- Seurat::PercentageFeatureSet(
            sc_seurat,
            pattern = quality_control.pattern
        )
    }
    if (data_filter) {
        chk::chk_length(data_filter.nFeature_RNA_thresh, 2)
        chk::chk_numeric(data_filter.nFeature_RNA_thresh)
        chk::chk_lt(
            data_filter.nFeature_RNA_thresh[1],
            data_filter.nFeature_RNA_thresh[2]
        )
        chk::chk_numeric(data_filter.percent.mt)
        chk::chk_range(data_filter.percent.mt, c(0, 100))

        sc_seurat <- subset(
            x = sc_seurat,
            subset = nFeature_RNA > data_filter.nFeature_RNA_thresh[1] &
                `nFeature_RNA` < data_filter.nFeature_RNA_thresh[2] &
                `percent.mt` < data_filter.percent.mt
        )
    }
    sc_seurat <- ProcessSeuratObject(
        obj = sc_seurat,
        normalization_method = normalization_method,
        scale_factor = scale_factor,
        scale_features = scale_features,
        selection_method = selection_method,
        verbose = verbose
    )

    sc_seurat <- ClusterAndReduce(
        sc_seurat,
        dims = dims,
        resolution = resolution,
        verbose = verbose
    )

    FilterTumorCell(
        obj = sc_seurat,
        column2only_tumor = column2only_tumor,
        verbose = verbose
    )
}

#' @rdname SCPreProcess
#' @export
#'
SCPreProcess.matrix <- function(
    sc,
    ...
) {
    dots <- list(...)
    verbose <- if ("verbose" %in% names(dots)) dots$verbose else TRUE
    # sc is a count matrix
    if (verbose) {
        ts_cli$cli_alert_info("Start from matrix")
    }
    NextMethod(generic = "SCPreProcess", sc = Matrix::Matrix(sc), ...)
}

#' @rdname SCPreProcess
#' @export
#'
SCPreProcess.data.frame <- function(
    sc,
    ...
) {
    dots <- list(...)
    verbose <- if ("verbose" %in% names(dots)) dots$verbose else TRUE
    # sc is a count matrix
    if (verbose) {
        ts_cli$cli_alert_info("Start from data.frame, convert to matrix")
    }
    NextMethod(
        generic = "SCPreProcess",
        sc = Matrix::Matrix(as.matrix(sc)),
        ...
    )
}

#' @rdname SCPreProcess
#' @export
#'
SCPreProcess.dgCMatrix <- function(
    sc,
    ...
) {
    dots <- list(...)
    verbose <- if ("verbose" %in% names(dots)) dots$verbose else TRUE
    # sc is a count matrix
    if (verbose) {
        ts_cli$cli_alert_info("Start from dgCMatrix")
    }
    NextMethod(generic = "SCPreProcess", sc = sc, ...)
}


#' @rdname SCPreProcess
#' @export
SCPreProcess.AnnDataR6 <- function(
    sc,
    meta_data = NULL,
    ...
) {
    dots <- list(...)
    verbose <- if ("verbose" %in% names(dots)) dots$verbose else TRUE
    if (is.null(sc$X)) {
        cli::cli_abort(c("x" = "Input must contain $X matrix"))
    }
    if (verbose) {
        ts_cli$cli_alert_info("Start from anndata object")
    }

    NextMethod(
        generic = "SCPreProcess",
        sc = Matrix::t(sc$X),
        meta_data = sc$obs,
        ...
    )
}

#' @rdname SCPreProcess
#' @export
#'
SCPreProcess.Seurat <- function(
    sc,
    column2only_tumor = NULL,
    ...
) {
    dots <- list(...)
    verbose <- if ("verbose" %in% names(dots)) dots$verbose else TRUE
    if (verbose) {
        ts_cli$cli_alert_info("Start from Seurat object")
    }

    # Validation message can be TRUE or message vector
    valid_msg <- methods::validObject(object = sc, test = TRUE)

    # * Successful validation
    if (is.logical(valid_msg)) {
        return(FilterTumorCell(
            obj = sc,
            column2only_tumor = column2only_tumor,
            verbose = verbose
        ))
    }

    # * Failure to validate the Seurat object
    if (verbose) {
        ts_cli$cli_alert_info(
            "Seurat object validation failed, try repairing it..."
        )
        cli::cli_h3(cli::style_bold("Validation message:"))
        purrr::walk(valid_msg, cli::cli_alert_danger)
        cli::cli_h3(cli::style_bold(
            "Try repairing with {.fn UpdateSeuratObject}:"
        ))
    }

    # Decorator function to wrap the UpdateSeuratObject function
    SafelyUpdateSeuratObject <- purrr::safely(
        SeuratObject::UpdateSeuratObject
    )
    updated_res <- SafelyUpdateSeuratObject(sc)
    # Successful repair
    if (is.null(updated_res$error)) {
        if (verbose) {
            ts_cli$cli_alert_success(cli::col_green(
                "Successfully repaired the Seurat object"
            ))
        }
        return(FilterTumorCell(
            obj = updated_res$result,
            column2only_tumor = column2only_tumor,
            verbose = verbose
        ))
    }
    # Failure to repair
    if (verbose) {
        cli::cli_alert_danger(updated_res$error)
        cli::cli_warn(
            "Seurat object repair failed. It is recommended to rebuild the Seurat object. Filtering is still being performed but may not be reliable."
        )
    }

    return(FilterTumorCell(
        obj = sc,
        column2only_tumor = column2only_tumor,
        verbose = verbose
    ))
}

#' @title Process a Seurat object (internal)
#'
#' @description
#' Normalize, find variable features, scale, and run PCA
#'
#' @param obj Seurat object
#' @param normalization_method Normalization method ("LogNormalize", "CLR", or "RC")
#' @param scale_factor Scaling factor for normalization
#' @param scale_features Features to scale
#' @param selection_method Variable feature selection method ("vst", "mvp", or "disp")
#' @param verbose Print progress messages
#' @return Seurat object
#'
#' @keywords internal
#' @family single_cell_preprocess
#'
ProcessSeuratObject <- function(
    obj,
    normalization_method = "LogNormalize",
    scale_factor = 10000,
    scale_features = NULL,
    selection_method = "vst",
    verbose = TRUE
) {
    Seurat::NormalizeData(
        object = obj,
        normalization.method = normalization_method,
        scale.factor = scale_factor,
        verbose = verbose
    ) %>%
        Seurat::FindVariableFeatures(
            selection.method = selection_method,
            verbose = verbose
        ) %>%
        Seurat::ScaleData(verbose = verbose, features = scale_features) %>%
        Seurat::RunPCA(
            features = Seurat::VariableFeatures(.),
            verbose = verbose
        )
}

#' @title Cluster and reduce dimensions (internal)
#'
#' @description
#' FindNeighbors, FindClusters, RunTSNE, RunUMAP
#'
#' @param obj Seurat object
#' @param dims Dimension to use for clustering and dimension reduction
#' @param resolution Resolution for clustering
#' @param verbose logical, whether to print progress messages
#' @return Seurat object
#'
#' @keywords internal
#' @family single_cell_preprocess
ClusterAndReduce <- function(
    obj,
    dims = seq_len(10),
    resolution = 0.6,
    verbose = TRUE
) {
    n_pcs <- ncol(obj@reductions$pca)
    if (is.null(dims)) {
        dims <- seq_len(FindRobustElbow(obj, verbose = verbose, ndims = 50))
    }
    if (max(dims) > n_pcs) {
        dims <- seq_len(min(max(dims), n_pcs))
        cli::cli_warn(
            "The input dimension is greater than the dimension in PCA. It is now set to the maximum dimension in PCA."
        )
    }

    Seurat::FindNeighbors(object = obj, dims = dims, verbose = verbose) %>%
        Seurat::FindClusters(
            resolution = resolution,
            verbose = verbose
        ) %>%
        Seurat::RunTSNE(dims = dims) %>%
        Seurat::RunUMAP(dims = dims, verbose = verbose)
}

#' @title Filter tumor cells (internal)
#'
#' @description
#' Filter tumor cells from Seurat object.
#'
#' @param obj Seurat object with a column to filter out tumor cells.
#' @param column2only_tumor Name of the column to filter out tumor cells.
#' @param verbose logical, whether to print progress messages
#'
#' @keywords internal
#' @family single_cell_preprocess
#'
FilterTumorCell <- function(
    obj,
    column2only_tumor = NULL,
    verbose = TRUE
) {
    obj %<>% AddMisc(self_dim = dim(.), cover = TRUE)

    if (is.null(column2only_tumor)) {
        return(obj)
    }
    if (!column2only_tumor %chin% colnames(obj[[]])) {
        cli::cli_warn(
            "Column '{.emph column2only_tumor}' not found, skip tumor cell filtering"
        )
        return(obj)
    }
    if (verbose) {
        ts_cli$cli_alert_info(
            "Filtering tumor cells from '{.emph column2only_tumor}'..."
        )
    }

    labels <- obj[[column2only_tumor]][[1]]
    tumor_cells <- grepl(
        "^[Tt]umo.?r|[Cc]ancer[Mm]alignant|[Nn]eoplasm|[Tt]um",
        labels
    )

    tumor_seurat <- obj[, tumor_cells] %>%
        AddMisc(
            raw_dim = dim(obj),
            self_dim = dim(.),
            column2only_tumor = column2only_tumor,
            cover = TRUE
        )

    return(tumor_seurat)
}
