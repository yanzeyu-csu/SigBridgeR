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
#' @param meta_data A data.frame containing metadata for each cell. It will be
#'    added to the Seurat object as `@meta.data`. If `NULL`, it will be
#'    extracted from the input object if possible.
#' @param column2only_tumor A character of column names in `meta_data`, used to
#'    filter the Seurat object to only tumor cells. If `NULL`, no filtering is performed.
#' @param project A character of project name, used to name the Seurat object.
#' @param min_cells Minimum number of cells that must express a feature for it
#'    to be included in the analysis. Defaults to `400`.
#' @param min_features Minimum number of features that must be detected in a
#'    cell for it to be included in the analysis. Defaults to `0`.
#' @param quality_control Logical indicating whether to perform mitochondrial
#'    percentage quality control. Defaults to `TRUE`.
#' @param quality_control.pattern Character pattern to identify mitochondrial
#'    genes, ribosomal protein genes, or other unwanted genes, as well as combinations
#'    of these genes. Customized patterns are supported. Defaults to `"^MT-"`.
#' @param data_filter Logical indicating whether to filter cells based on
#'    quality metrics. Defaults to `TRUE`.
#' @param data_filter.nFeature_RNA_thresh Numeric vector of length 2 specifying
#'    the minimum and maximum number of features per cell. Defaults to `c(200, 6000)`.
#' @param data_filter.percent.mt Maximum mitochondrial percentage allowed.
#'    Defaults to `20`.
#' @param normalization_method Method for normalization: "LogNormalize", "CLR",
#'    or "RC". Defaults to `"LogNormalize"`.
#' @param scale_factor Scaling factor for normalization. Defaults to `10000`.
#' @param scale_features Features to use for scaling. If NULL, uses all variable
#'    features. Defaults to `NULL`.
#' @param selection_method Method for variable feature selection: "vst", "mvp",
#'    or "disp". Defaults to `"vst"`.
#' @param resolution Resolution parameter for clustering. Higher values lead to
#'    more clusters. Defaults to `0.6`.
#' @param dims Dimensions to use for clustering and dimensionality reduction.
#'    If NULL, automatically determined by elbow method. Defaults to `NULL`.
#' @param verbose Logical indicating whether to print progress messages.
#'    Defaults to `TRUE`.
#' @param ... Additional arguments passed to specific methods. Currently unused.
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
#'
#' @examples
#' \dontrun{
#' # Example with matrix input
#' counts_matrix <- matrix(rpois(1000, 5), nrow = 100, ncol = 10)
#' rownames(counts_matrix) <- paste0("Gene", 1:100)
#' colnames(counts_matrix) <- paste0("Cell", 1:10)
#'
#' seurat_obj <- SCPreProcess(
#'   sc = counts_matrix,
#'   project = "TestProject",
#'   min_features = 50,
#'   resolution = 0.8
#' )
#'
#' # Example with tumor cell filtering
#' metadata <- data.frame(
#'   cell_type = c(rep("Tumor", 5), rep("Normal", 5)),
#'   row.names = paste0("Cell", 1:10)
#' )
#'
#' tumor_seurat <- SCPreProcess(
#'   sc = counts_matrix,
#'   meta_data = metadata,
#'   column2only_tumor = "cell_type",
#'   project = "TumorAnalysis"
#' )
#' }
#'
#' @seealso
#' \code{\link[Seurat]{CreateSeuratObject}},
#' \code{\link[Seurat]{NormalizeData}},
#' \code{\link[Seurat]{ScaleData}},
#' \code{\link[Seurat]{FindVariableFeatures}},
#' \code{\link[Seurat]{RunPCA}},
#' \code{\link[Seurat]{RunTSNE}},
#' \code{\link[Seurat]{RunUMAP}},
#' \code{\link[Seurat]{FindNeighbors}},
#' \code{\link[Seurat]{FindClusters}}
#'
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
    project = "SC_Screening_Proj",
    min_cells = 400L,
    min_features = 0L,
    quality_control = TRUE,
    quality_control.pattern = c("^MT-"),
    data_filter = TRUE,
    data_filter.nFeature_RNA_thresh = c(200L, 6000L),
    data_filter.percent.mt = 20L,
    normalization_method = "LogNormalize",
    scale_factor = 10000L,
    scale_features = NULL,
    selection_method = "vst",
    resolution = 0.6,
    dims = NULL,
    verbose = TRUE,
    ...
) {
    if (!is.null(meta_data)) {
        chk::chk_data(meta_data)
    }
    if (!is.null(column2only_tumor)) {
        chk::chk_character(column2only_tumor)
    }

    sc_seurat <- SeuratObject::CreateSeuratObject(
        counts = sc,
        project = project,
        meta.data = meta_data,
        min.cells = min_cells,
        min.features = min_features
    )

    if (quality_control) {
        chk::chk_character(quality_control.pattern)

        sc_seurat <- QCPatternFilter(
            obj = sc_seurat,
            pattern = quality_control.pattern,
            verbose = verbose
        )
    }
    if (data_filter) {
        chk::chk_length(data_filter.nFeature_RNA_thresh, 2)
        chk::chk_numeric(data_filter.nFeature_RNA_thresh)
        chk::chk_lt(
            data_filter.nFeature_RNA_thresh[1],
            data_filter.nFeature_RNA_thresh[2]
        )
        chk::chk_range(data_filter.percent.mt, range = c(0, 100))

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
    meta_data = NULL,
    column2only_tumor = NULL,
    project = "SC_Screening_Proj",
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
    dims = NULL,
    verbose = TRUE,
    ...
) {
    # sc is a count matrix
    if (verbose) {
        ts_cli$cli_alert_info("Start from matrix")
    }
    NextMethod(
        generic = "SCPreProcess",
        sc = Matrix::Matrix(sc),
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
SCPreProcess.data.frame <- function(
    sc,
    meta_data = NULL,
    column2only_tumor = NULL,
    project = "SC_Screening_Proj",
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
    dims = NULL,
    verbose = TRUE,
    ...
) {
    # sc is a count matrix
    if (verbose) {
        ts_cli$cli_alert_info("Start from data.frame, convert to matrix")
    }
    NextMethod(
        generic = "SCPreProcess",
        sc = Matrix::Matrix(as.matrix(sc)),
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
    project = "SC_Screening_Proj",
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
    dims = NULL,
    verbose = TRUE,
    ...
) {
    # sc is a count matrix
    if (verbose) {
        ts_cli$cli_alert_info("Start from dgCMatrix")
    }
    NextMethod(
        generic = "SCPreProcess",
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
    project = "SC_Screening_Proj",
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
    dims = NULL,
    verbose = TRUE,
    ...
) {
    if (is.null(sc$X)) {
        cli::cli_abort(c("x" = "Input must contain $X matrix"))
    }
    if (verbose) {
        ts_cli$cli_alert_info("Start from anndata object")
    }
    meta_data <- if (!is.null(sc$obs) && !is.null(meta_data)) {
        cbind(meta_data, sc$obs)
    } else if (!is.null(sc$obs)) {
        sc$obs
    }

    NextMethod(
        generic = "SCPreProcess",
        sc = Matrix::t(sc$X),
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
SCPreProcess.Seurat <- function(
    sc,
    column2only_tumor = NULL,
    verbose = TRUE,
    ...
) {
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
    if (is.null(updated_res$error) && verbose) {
        ts_cli$cli_alert_success(cli::col_green(
            "Successfully repaired the Seurat object"
        ))
    } else if (verbose) {
        # Failure to repair
        cli::cli_alert_danger(updated_res$error)
        cli::cli_warn(
            "Seurat object repair failed. It is recommended to rebuild the Seurat object. Filtering is still being performed but may not be reliable."
        )
    }

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
    dims = NULL,
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
    obj <- AddMisc(obj, self_dim = dim(.), cover = TRUE)

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

    AddMisc(
        seurat_obj = obj[, tumor_cells],
        raw_dim = dim(obj),
        self_dim = dim(.),
        column2only_tumor = column2only_tumor,
        cover = TRUE
    )
}

#' @title Calculate Percentage of Features Matching Patterns
#'
#' @description
#' This function calculates the percentage of counts coming from features matching
#' specified patterns (e.g., mitochondrial genes, ribosomal genes) and adds them
#' as metadata columns to the Seurat object.
#'
#' @param obj A seurat object.
#' @param pattern Character pattern to identify mitochondrial
#'    genes, ribosomal protein genes, or other unwanted genes, as well as combinations
#'    of these genes. Customized patterns are supported.
#' @param verbose logical, whether to print progress messages
#' @param ... Additional arguments passed to \code{\link[Seurat]{PercentageFeatureSet}}
#'
#' @details
#' The function automatically generates friendly column names based on the patterns:
#' - "mt" for mitochondrial patterns
#' - "rp" for ribosomal patterns
#' - "rrna" for ribosomal RNA patterns
#' - For combined patterns (using |), creates names like "mt_rp"
#' - For other patterns, creates cleaned lowercase names
#'
#' @keywords internal
#' @family single_cell_preprocess
#'
QCPatternFilter <- function(
    obj,
    pattern = c("^MT-", "^mt-", "^RP[SL]", "^MT-|^RP[SL]"),
    verbose = TRUE,
    ...
) {
    Pattern2Colname <- function(pat) {
        pat_lower <- tolower(pat)

        if (grepl("\\|", pat)) {
            # Handle combined patterns (with | separator)
            parts <- strsplit(pat, "|", fixed = TRUE)[[1]]
            names <- purrr::map_chr(parts, function(p) {
                dplyr::case_when(
                    grepl("mt", p_lower) ~ "mt",
                    grepl("rp", p_lower) ~ "rp",
                    grepl("rrna|rna[0-9]", p_lower) ~ "rrna",
                    TRUE ~ tolower(gsub("[^[:alnum:]]", "", p))
                )
            })

            paste(sort(unique(names)), collapse = "_")
        } else {
            # Handle single patterns
            dplyr::case_when(
                grepl("mt", pat_lower) ~ "mt",
                grepl("rp", pat_lower) ~ "rp",
                grepl("rrna|rna[0-9]", pat_lower) ~ "rrna",
                TRUE ~ {
                    clean <- gsub("[^[:alnum:]]", "_", pat)
                    clean <- gsub("_+", "_", clean) # Collapse multiple underscores
                    clean <- gsub("^_+|_+$", "", clean) # Trim leading/trailing underscores
                    tolower(clean)
                }
            )
        }
    }

    for (pat in pattern) {
        col_name <- paste0("percent.", Pattern2Colname(pat))

        if (col_name %chin% colnames(obj[[]])) {
            if (verbose) {
                cli::cli_warn(
                    "Column {.val {col_name}} already exists. Skipping pattern: {.val {pat}}",
                )
            }
            next
        }

        obj[[col_name]] <- Seurat::PercentageFeatureSet(
            obj,
            pattern = pat,
            ...
        )
    }

    obj
}
