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
#'        - `data.frame/matrix/dgCMatrix`: Raw count matrix (features x cells)
#'        - `AnnDataR6`: Python AnnData object via reticulate
#'        - `Seurat`: Preprocessed Seurat object
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
SCPreProcess.default <- function(sc, ...) {
    cli::cli_abort(c("x" = "Unsupported input type: {.var {class(sc)}}"))
}

#' @rdname SCPreProcess
#' @param meta_data Optional metadata dataframe (rows = cells, columns = attributes)
#' @param column2only_tumor Metadata column used for filtering tumor cells, matching patterns
#' such as "tumor" (or "tumour"), "cancer", "malignant", or "neoplasm" (case-insensitive).
#' @param project Project name for Seurat object
#' @param min_cells Minimum cells per gene to retain (features present in at least this many cells)
#' @param min_features Minimum features per cell to retain (cells with at least this many features)
#' @param quality_control Logical indicating whether to perform mitochondrial QC
#' @param quality_control.pattern Regex pattern to identify mitochondrial genes (e.g., "^MT-" for human)
#' @param data_filter Logical indicating whether to filter low-quality cells
#' @param data_filter.nFeature_RNA_thresh Numeric vector of length 2 specifying (min, max) features per cell
#' @param data_filter.percent.mt Maximum mitochondrial percentage allowed (0-100)
#' @param normalization_method Normalization method ("LogNormalize", "CLR", or "RC")
#' @param scale_factor Scaling factor for normalization
#' @param scale_features Scale features to unit variance
#' @param selection_method Variable feature selection method ("vst", "mvp", or "disp")
#' @param resolution Cluster resolution (higher for more clusters)
#' @param dims PCA dimensions to use
#' @param verbose Print progress messages
#' @export
#'
SCPreProcess.matrix <- function(
    sc,
    meta_data = NULL,
    column2only_tumor = NULL,
    project = glue::glue("{TimeStamp()}_SC_Screening_Proj"),
    min_cells = 400,
    min_features = 0,
    quality_control = TRUE,
    quality_control.pattern = c("^MT-", "^mt-"),
    data_filter = TRUE,
    data_filter.nFeature_RNA_thresh = c(200, 6000),
    data_filter.percent.mt = 20,
    normalization_method = "LogNormalize",
    scale_factor = 10000,
    scale_features = NULL,
    selection_method = "vst",
    resolution = 0.6,
    dims = 1:10,
    verbose = TRUE,
    ...
) {
    chk::chk_null_or(meta_data, chk::chk_data)
    chk::chk_null_or(column2only_tumor, chk::chk_character)

    # sc is a count matrix
    if (verbose) {
        cli::cli_alert_info("[{TimeStamp()}] Start from count matrix")
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
            cli::cli_abort(c(
                "x" = "{.arg quality_control.pattern} must be specified.",
                "i" = "human: '^MT-', mouse: '^mt-', specify your own or set {.arg quality_control=FALSE}."
            ))
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

        sc_seurat = subset(
            x = sc_seurat,
            subset = nFeature_RNA > data_filter.nFeature_RNA_thresh[1] &
                nFeature_RNA < data_filter.nFeature_RNA_thresh[2] &
                percent.mt < data_filter.percent.mt
        )
    }
    sc_seurat = ProcessSeuratObject(
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
        column2only_tumor = column2only_tumor
    )
}

#' @rdname SCPreProcess
#' @param meta_data Optional metadata dataframe (rows = cells, columns = attributes)
#' @param column2only_tumor Metadata column used for filtering tumor cells, matching patterns
#' such as "tumor" (or "tumour"), "cancer", "malignant", or "neoplasm" (case-insensitive).
#' @param project Project name for Seurat object
#' @param min_cells Minimum cells per gene to retain (features present in at least this many cells)
#' @param min_features Minimum features per cell to retain (cells with at least this many features)
#' @param quality_control Logical indicating whether to perform mitochondrial QC
#' @param quality_control.pattern Regex pattern to identify mitochondrial genes (e.g., "^MT-" for human)
#' @param data_filter Logical indicating whether to filter low-quality cells
#' @param data_filter.nFeature_RNA_thresh Numeric vector of length 2 specifying (min, max) features per cell
#' @param data_filter.percent.mt Maximum mitochondrial percentage allowed (0-100)
#' @param normalization_method Normalization method ("LogNormalize", "CLR", or "RC")
#' @param scale_factor Scaling factor for normalization
#' @param scale_features Scale features to unit variance
#' @param selection_method Variable feature selection method ("vst", "mvp", or "disp")
#' @param resolution Cluster resolution (higher for more clusters)
#' @param dims PCA dimensions to use
#' @param verbose Print progress messages
#' @export
#'
SCPreProcess.data.frame <- function(
    sc,
    meta_data = NULL,
    column2only_tumor = NULL,
    project = glue::glue("{TimeStamp()}_SC_Screening_Proj"),
    min_cells = 400,
    min_features = 0,
    quality_control = TRUE,
    quality_control.pattern = c("^MT-", "^mt-"),
    data_filter = TRUE,
    data_filter.nFeature_RNA_thresh = c(200, 6000),
    data_filter.percent.mt = 20,
    normalization_method = "LogNormalize",
    scale_factor = 10000,
    scale_features = NULL,
    selection_method = "vst",
    resolution = 0.6,
    dims = 1:10,
    verbose = TRUE,
    ...
) {
    SCPreProcess.matrix(
        sc = as.matrix(sc),
        meta_data = meta_data,
        column2only_tumor = column2only_tumor,
        project = project,
        min_cells = min_cells,
        min_features = min_features,
        quality_control = quality_control,
        quality_control.pattern = quality_control.pattern,
        data_filter = data_filter,
        data_filter.nFeature_RNA_thresh = data_filter.nFeature_RNA_thresh,
        data_filter.percent.mt = data_filter.percent.mt,
        normalization_method = normalization_method,
        scale_factor = scale_factor,
        scale_features = scale_features,
        selection_method = selection_method,
        resolution = resolution,
        dims = dims,
        verbose = verbose,
        ...
    )
}

#' @rdname SCPreProcess
#' @export
#'
SCPreProcess.dgCMatrix <- function(
    sc,
    meta_data = NULL,
    column2only_tumor = NULL,
    project = glue::glue("{TimeStamp()}_SC_Screening_Proj"),
    min_cells = 400,
    min_features = 0,
    quality_control = TRUE,
    quality_control.pattern = "^MT-",
    data_filter = TRUE,
    data_filter.nFeature_RNA_thresh = c(200, 6000),
    data_filter.percent.mt = 20,
    normalization_method = "LogNormalize",
    scale_factor = 10000,
    scale_features = NULL,
    selection_method = "vst",
    resolution = 0.6,
    dims = 1:10,
    verbose = TRUE,
    ...
) {
    SCPreProcess.matrix(
        sc = sc,
        meta_data = meta_data,
        column2only_tumor = column2only_tumor,
        project = project,
        min_cells = min_cells,
        min_features = min_features,
        quality_control = quality_control,
        quality_control.pattern = quality_control.pattern,
        data_filter = data_filter,
        data_filter.nFeature_RNA_thresh = data_filter.nFeature_RNA_thresh,
        data_filter.percent.mt = data_filter.percent.mt,
        normalization_method = normalization_method,
        scale_factor = scale_factor,
        scale_features = scale_features,
        selection_method = selection_method,
        resolution = resolution,
        dims = dims,
        verbose = verbose,
        ...
    )
}


#' @rdname SCPreProcess
#' @export
SCPreProcess.AnnDataR6 <- function(
    sc,
    meta_data = NULL,
    column2only_tumor = NULL,
    project = glue::glue("{TimeStamp()}_SC_Screening_Proj"),
    min_cells = 400,
    min_features = 0,
    quality_control = TRUE,
    quality_control.pattern = c("^MT-", "^mt-"),
    data_filter = TRUE,
    data_filter.nFeature_RNA_thresh = c(200, 6000),
    data_filter.percent.mt = 20,
    normalization_method = "LogNormalize",
    scale_factor = 10000,
    scale_features = NULL,
    selection_method = "vst",
    resolution = 0.6,
    dims = 1:10,
    verbose = TRUE,
    ...
) {
    chk::chk_null_or(meta_data, chk::chk_data)
    chk::chk_null_or(column2only_tumor, chk::chk_character)

    if (is.null(sc$X)) {
        cli::cli_abort(c("x" = "Input must contain $X matrix"))
    }
    if (verbose) {
        cli::cli_alert_info("[{TimeStamp()}] Start from anndata object")
    }

    sc_matrix <- t(sc$X)

    sc_seurat <- Seurat::CreateSeuratObject(
        counts = sc_matrix,
        project = project,
        meta.data = meta_data,
        min.cells = min_cells,
        min.features = min_features
    )

    if (quality_control) {
        if (length(quality_control.pattern) != 1) {
            cli::cli_abort(c(
                "x" = "{.arg quality_control.pattern} must be specified",
                "i" = "human: '^MT-', mouse: '^mt-', specify your own or set {.arg quality_control=FALSE}."
            ))
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

        sc_seurat = subset(
            x = sc_seurat,
            subset = nFeature_RNA > data_filter.nFeature_RNA_thresh[1] &
                nFeature_RNA < data_filter.nFeature_RNA_thresh[2] &
                percent.mt < data_filter.percent.mt
        )
    }

    sc_seurat = ProcessSeuratObject(
        obj = sc_seurat,
        normalization_method = normalization_method,
        scale_factor = scale_factor,
        scale_features = scale_features,
        selection_method = selection_method,
        verbose = verbose
    )

    # Add metadata
    if ("obs" %in% names(sc)) {
        sc_seurat <- Seurat::AddMetaData(sc_seurat, sc$obs)
    }

    sc_seurat <- ClusterAndReduce(
        sc_seurat,
        dims = dims,
        resolution = resolution,
        verbose = verbose
    )

    FilterTumorCell(
        obj = sc_seurat,
        column2only_tumor = column2only_tumor
    )
}

#' @rdname SCPreProcess
#' @export
#'
SCPreProcess.Seurat <- function(
    sc,
    column2only_tumor = NULL,
    verbose = TRUE,
    ...
) {
    FilterTumorCell(
        obj = sc,
        column2only_tumor = column2only_tumor
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
#' @param scale_features Features to scale
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
    scale_features = NULL,
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
        Seurat::ScaleData(verbose = verbose, features = scale_features) %>%
        Seurat::RunPCA(
            features = Seurat::VariableFeatures(.),
            verbose = verbose
        )
}

#' Cluster and reduce dimensions (internal)
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
        cli::cli_warn(
            "The input dimension is greater than the dimension in PCA. It is now set to the maximum dimension in PCA."
        )
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
#' @param column2only_tumor Name of the column to filter out tumor cells.
#'
#' @keywords internal
#'
#'
FilterTumorCell <- function(
    obj,
    column2only_tumor = NULL
) {
    obj = AddMisc(obj, self_dim = dim(obj), cover = TRUE)

    if (!is.null(column2only_tumor)) {
        if (!column2only_tumor %in% colnames(obj@meta.data)) {
            cli::cli_alert_danger(crayon::red(
                "Column '{column2only_tumor}' not found, skip tumor cell filtering"
            ))
            return(obj)
        } else {
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
    } else {
        return(obj)
    }
}
