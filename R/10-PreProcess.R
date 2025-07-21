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
#'        - `data.frame/matrix`: Raw count matrix (features x cells)
#'        - `AnnDataR6`: Python AnnData object via reticulate
#'        - `Seurat`: Preprocessed Seurat object
#' @param ... Method-specific arguments (see below)
#'
#' @return A Seurat object containing:
#' \itemize{
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
SCPreProcess.default <- function(sc, ...) {
    stop("Unknown Input type")
}

#' @rdname SCPreProcess
#' @param column2only_tumor Metadata column for tumor cell filtering (regex patterns:
#'        "[Tt]umo.?r", "[Cc]ancer", "[Mm]alignant", "[Nn]eoplasm")
#' @param project Project name for Seurat object
#' @param min_cells Minimum cells per gene to retain
#' @param min_features Minimum features per cell to retain
#' @param normalization_method Normalization method ("LogNormalize", "CLR", or "RC")
#' @param scale_factor Scaling factor for normalization
#' @param selection_method Variable feature selection method ("vst", "mvp", or "disp")
#' @param resolution Cluster resolution (higher for more clusters)
#' @param dims PCA dimensions to use
#' @param verbose Print progress messages
#' @param future_global_maxsize Memory limit for parallelization (bytes)
#' @export
#'
SCPreProcess.data.frame <- function(
    sc,
    column2only_tumor = NULL,
    project = "Scissor_Single_Cell",
    min_cells = 400,
    min_features = 0,
    normalization_method = "LogNormalize",
    scale_factor = 10000,
    selection_method = "vst",
    resolution = 0.6,
    dims = 1:10,
    verbose = TRUE,
    future_global_maxsize = 6 * 1024^3,
    ...
) {
    options(future.globals.maxSize = future_global_maxsize)

    # sc is a count matrix
    cli::cli_alert_info("Start from count matrix")

    sc_seurat <- SeuratObject::CreateSeuratObject(
        counts = sc,
        project = project,
        min.cells = min_cells,
        min.features = min_features
    ) %>%
        ProcessSeuratObject(
            normalization_method = normalization_method,
            scale_factor = scale_factor,
            selection_method = selection_method,
            verbose = verbose
        )

    # Add metadata
    if ("obs" %in% names(sc)) {
        sc_seurat <- sc_seurat %>% Seurat::AddMetaData(sc$obs)
    }

    sc_seurat <- ClusterAndReduce(
        sc_seurat,
        dims = dims,
        resolution = resolution,
        verbose = verbose
    )

    FilterTumorCell(
        sc_seurat,
        column2only_tumor = column2only_tumor,
        verbose = verbose
    )
}

#' @rdname SCPreProcess
#' @export
SCPreProcess.AnnDataR6 <- function(
    sc,
    column2only_tumor = NULL,
    project = "Scissor_Single_Cell",
    min_cells = 400,
    min_features = 0,
    normalization_method = "LogNormalize",
    scale_factor = 10000,
    selection_method = "vst",
    resolution = 0.6,
    dims = 1:10,
    verbose = TRUE,
    future_global_maxsize = 6 * 1024^3,
    ...
) {
    options(future.globals.maxSize = future_global_maxsize)

    if (is.null(sc$X)) {
        stop("Input must contain $X matrix")
    }
    cli::cli_alert_info("Start from anndata object")

    sc_matrix <- if (inherits(sc$X, "sparseMatrix")) {
        if (!requireNamespace("Matrix", quietly = TRUE)) {
            stop("Matrix package required for sparse matrix support")
        }
        Matrix::t(sc$X)
    } else {
        t(sc$X)
    }

    sc_seurat <- SeuratObject::CreateSeuratObject(
        counts = sc_matrix,
        project = project,
        min.cells = min_cells,
        min.features = min_features
    ) %>%
        ProcessSeuratObject(
            normalization_method = normalization_method,
            scale_factor = scale_factor,
            selection_method = selection_method,
            verbose = verbose
        )

    # Add metadata
    if ("obs" %in% names(sc)) {
        sc_seurat <- sc_seurat %>% Seurat::AddMetaData(sc$obs)
    }

    sc_seurat <- ClusterAndReduce(
        sc_seurat,
        dims = dims,
        resolution = resolution,
        verbose = verbose
    )

    FilterTumorCell(
        sc_seurat,
        column2only_tumor = column2only_tumor,
        verbose = verbose
    )
}

#' @rdname SCPreProcess
#' @export
#'
SCPreProcess.Seurat <- function(
    sc,
    column2only_tumor = NULL,
    verbose = TRUE,
    future_global_maxsize = 6 * 1024^3,
    ...
) {
    options(future.globals.maxSize = future_global_maxsize)

    FilterTumorCell(
        obj = sc,
        column2only_tumor = column2only_tumor,
        verbose = verbose
    )
}

#' @title Process a Seurat object (internal)
#'
#' @description
#' Normalize, find variable features, scale, and run PCA
#'
#' @param obj Seurat object
#' @param normalization_method Normalization method ("LogNormalize", "CLR", or "RC")
#' @param scale_factor Scaling factor for normalization
#' @param selection_method Variable feature selection method ("vst", "mvp", or "disp")
#' @param verbose Print progress messages
#' @return Seurat object
#'
#' @keywords internal
#'
ProcessSeuratObject <- function(
    obj,
    normalization_method = "LogNormalize",
    scale_factor = 10000,
    selection_method = "vst",
    verbose = TRUE
) {
    obj %>%
        Seurat::NormalizeData(
            normalization.method = normalization_method,
            scale.factor = scale_factor,
            verbose = verbose
        ) %>%
        Seurat::FindVariableFeatures(
            selection.method = selection_method,
            verbose = verbose
        ) %>%
        Seurat::ScaleData(verbose = verbose) %>%
        Seurat::RunPCA(
            features = SeuratObject::VariableFeatures(.),
            verbose = verbose
        )
}

#' Cluster and reduce dimensions (internal)
#'
#' @description
#' FindNeighbors, FindClusters, RunTSNE, RunUMAP
#'
#'
#' @keywords internal
#'
ClusterAndReduce <- function(
    obj,
    dims = 1:10,
    resolution = 0.6,
    verbose = TRUE
) {
    n_pcs <- ncol(obj@reductions$pca)
    if (is.null(dims) || max(dims) > n_pcs) {
        dims <- 1:min(max(dims), n_pcs)
        cli::cli_alert_warning(crayon::yellow(
            "The input dimension is greater than the dimension in PCA. It is now set to the maximum dimension in PCA."
        ))
    }

    obj %>%
        Seurat::FindNeighbors(dims = dims, verbose = verbose) %>%
        Seurat::FindClusters(
            resolution = resolution,
            verbose = verbose
        ) %>%
        Seurat::RunTSNE(dims = dims) %>%
        Seurat::RunUMAP(dims = dims, verbose = verbose)
}

#' Filter tumor cells (internal)
#'
#' @description
#' Filter tumor cells from Seurat object.
#'
#' @param obj Seurat object with a column to filter out tumor cells.
#' @param name2only_tumor Name of the column to filter out tumor cells.
#' @param verbose Logical. Whether to print messages.
#'
#' @keywords internal
#'
#'
FilterTumorCell <- function(
    obj,
    column2only_tumor = NULL,
    verbose = TRUE
) {
    obj = AddMisc(obj, self_dim = dim(obj), cover = TRUE)

    if (!is.null(column2only_tumor)) {
        ifelse(
            !column2only_tumor %in% colnames(obj@meta.data),
            {
                cli::cli_alert_danger(crayon::red(
                    "Column '{column2only_tumor}' not found, skip tumor cell filtering"
                ))
                return(obj)
            },
            {
                labels <- obj[[column2only_tumor]][[1]]
                tumor_cells <- grepl(
                    "^[Tt]umo.?r|[Cc]ancer[Mm]alignant|[Nn]eoplasm",
                    labels
                )

                tumor_seurat <- obj[, tumor_cells] %>%
                    AddMisc(
                        raw_dim = dim(obj),
                        self_dim = dim(.),
                        column2only_tumor = column2only_tumor,
                        cover = TRUE
                    )

                return(list(tumor_seurat = tumor_seurat, raw_seurat = obj))
            }
        )
    } else {
        return(obj)
    }
}


# * ---- Preprocess bulk expression data ----

#' @title Preprocess bulk expression data
#'
#' @description
#' Preprocess bulk expression data: convert Ensembles version IDs and TCGA version IDs to genes.
#' @param data raw bulk expression data
#'
#' @export
#'
BulkPreProcess = function(data) {
    #   rownames(data) <- substr(rownames(data), 1, 15)
    options(
        IDConverter.datapath = system.file("extdata", package = "IDConverter")
    )
    rownames(data) <- IDConverter::convert_hm_genes(rownames(data))
    return(data)
}


