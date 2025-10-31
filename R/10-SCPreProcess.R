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
#'    - `R6`: Python AnnData object, obtained via package `anndata` or `anndataR`
#'    - `Seurat`: Preprocessed Seurat object
#' @param meta_data A data.frame containing metadata for each cell. It will be
#'    added to the Seurat object as `@meta.data`. If sc is an anndata object,
#'    `obs` will be automatically used.
#' @param column2only_tumor A character of column names in `meta_data`, used to
#'    filter the Seurat object to only tumor cells. If `NULL`, no filtering is performed.
#' @param project A character of project name, used to name the Seurat object.
#' @param min_cells Minimum number of cells that must express a feature for it
#'    to be included in the analysis. Defaults to `400L`.
#' @param min_features Minimum number of features that must be detected in a
#'    cell for it to be included in the analysis. Defaults to `0L`.
#' @param quality_control Logical, for identification of the proportion of mitochondrial genes,
#'    ribosomal protein genes, or other types of genes (without filtering),
#'    the results will be stored in `meta.data`. Defaults to `TRUE`.
#' @param quality_control.pattern A vector or character containing regex pattern(s) to identify
#'    mitochondrial genes, ribosomal protein genes, or other unwanted genes, as well as combinations
#'    of these genes. Customized patterns are supported. Defaults to `"^MT-"`.
#' @param data_filter Logical indicating whether to filter cells based on
#'    quality metrics. Defaults to `TRUE`.
#' @param data_filter.thresh A list containing filtering thresholds for different quality metrics:
#'    \itemize{
#'      \item `nFeature_RNA_thresh`: Numeric vector of length 2 specifying the minimum and maximum
#'            number of features per cell. Defaults to `c(200L, 6000L)`
#'      \item `percent.mt`: Maximum mitochondrial percentage allowed. Defaults to `20L`
#'      \item `percent.rp`: Maximum ribosomal protein percentage allowed. Not used in default
#'            unless ribosomal protein genes filter pattern (like `^RP[LS]`) is added to
#'            `quality_control.pattern`
#'    }
#' @param data_filter.percent.mt Maximum mitochondrial percentage allowed.
#'    Defaults to `20L`.
#' @param normalization_method Method for normalization: "LogNormalize", "CLR",
#'    or "RC". Defaults to `"LogNormalize"`.
#' @param scale_factor Scaling factor for normalization. Defaults to `10000L`.
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
#'   \item Normalized and scaled expression data
#'   \item Variable features identified by selection method
#'   \item PCA, t-SNE, and UMAP dimensionality reductions
#'   \item Cluster identities at specified resolution
#'   \item Quality control metrics in `@meta.data`
#'   \item When tumor cells filtered: original dimensions in `@misc$raw_dim`
#'   \item Final dimensions in `@misc$self_dim`
#'   \item Quality control column names in `@misc$qc_colnames`
#' }
#'
#' @details
#' \strong{Quality Control Patterns:}
#' The function supports flexible pattern matching for quality control, for example:
#' \itemize{
#'   \item \code{"^MT-"} - Mitochondrial genes (default)
#'   \item \code{"^RP\[LS\]"} - Ribosomal protein genes
#'   \item \code{"^\[rt\]rna"} - rRNA and tRNA genes
#'   \item Custom patterns using regular expressions
#'   \item Combined patterns: \code{"^MT-|^RP\[LS\]"} for both mitochondrial and ribosomal genes
#' }
#'
#' \strong{Flexible Filtering:}
#' The filtering system dynamically adapts to detected quality control patterns:
#' \itemize{
#'   \item Column names are automatically generated from patterns
#'   \item Multiple thresholds can be specified in \code{data_filter.thresh}
#'   \item Use \code{SigBridgeR:::Pattern2Colname} to determine correct column names for custom patterns if still confused
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
    data_filter.thresh = list(
        nFeature_RNA_thresh = c(200L, 6000L),
        # * only used when specifed in `quality_control.pattern`
        percent.mt = 20L, # mitochondrial genes
        percent.rp = 60L # ribosomal protein genes
        # ? When combined pattern is used, like `quality_control.pattern = "^MT-|^RP[LS]"`
        # ? Use `_` to separate different patterns like this:
        # percent.mt_rp = 60L

        # ? When filtering for non-mitochondrial genes and non-ribosomal proteins RNA genes,
        # ? the column names are in lowercase letter form with regular expression symbols removed.
        # `quality_control.pattern = "^[rt]rna"`
        # Correct threshhold setting is `percent.rt_rna = 60L`

        # ? Use `SigBridgeR:::Pattern2Colname` to get the correct colname if still confused.
    ),
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

        sc_seurat <- QCPatternDetect(
            obj = sc_seurat,
            pattern = quality_control.pattern,
            verbose = verbose
        )
    }
    if (data_filter) {
        # first 2 numbers will be used
        chk::chk_lt(
            data_filter.thresh$nFeature_RNA[1],
            data_filter.thresh$nFeature_RNA[2]
        )

        data_filter.thresh <- unlist(data_filter.thresh)
        # filter expr is a string of the form
        # which is used to subset the Seurat object
        nfeat_condition <- expr(
            nFeature_RNA > !!data_filter.thresh[["nFeature_RNA_thresh1"]] &
                nFeature_RNA < !!data_filter.thresh[["nFeature_RNA_thresh2"]]
        )

        # see `QCPatternDetect()` for the column names generation
        qc_colnames <- unlist(sc_seurat@misc$qc_colnames)
        # make sure the exact matching works
        qc_condition <- if (!is.null(qc_colnames)) {
            purrr::map(
                qc_colnames,
                ~ expr(!!sym(.x) < !!data_filter.thresh[[.x]])
            )
        } else {
            NULL
        }
        full_expr <- purrr::reduce(
            .x = c(list(nfeat_condition), qc_condition),
            .f = function(x, y) expr(!!x & !!y)
        )

        sc_seurat <- subset(
            x = sc_seurat,
            subset = full_expr
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
    # sc is a count matrix
    if (verbose) {
        ts_cli$cli_alert_info("Start from data.frame, convert it to matrix")
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

# * ---- anndata and anndataR ----

#' @rdname SCPreProcess
#' @export
SCPreProcess.R6 <- function(
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
    # Both `anndata` and `anndataR` are based on R6
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
    raw_dim <- dim(obj)
    obj <- AddMisc(obj, self_dim = raw_dim, cover = TRUE)

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

    obj = obj[, tumor_cells]

    AddMisc(
        seurat_obj = obj,
        raw_dim = raw_dim,
        self_dim = dim(obj),
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
#' @param pattern A character vector or list containing regex patterns to identify mitochondrial
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
QCPatternDetect <- function(
    obj,
    pattern = c("^MT-", "^mt-", "^RP[SL]", "^MT-|^RP[SL]"),
    verbose = TRUE,
    ...
) {
    # Each pattern will be stored in a list
    colname_mapping <- list()

    for (pat in unlist(pattern)) {
        col_name <- paste0("percent.", Pattern2Colname(pat))
        colname_mapping[[pat]] <- col_name

        if (col_name %chin% colnames(obj[[]])) {
            if (verbose) {
                cli::cli_warn(
                    "Column {.val {col_name}} already exists. Skipping pattern: {.val {pat}}",
                )
            }
            next
        }

        obj <- Seurat::PercentageFeatureSet(
            obj,
            pattern = pat,
            col.name = col_name,
            ...
        )
    }
    # Record these colnames to misc slot for further data filter
    obj@misc$qc_colnames <- colname_mapping

    obj
}

#' @title convert regex patterns to column names (internal)
#'
#' @keywords internal
Pattern2Colname <- function(pat) {
    pat_lower <- tolower(pat)

    if (grepl("\\|", pat_lower)) {
        # Handle combined patterns (with | separator)
        parts <- strsplit(pat_lower, "|", fixed = TRUE)[[1]]
        names <- purrr::map_chr(parts, function(p) {
            dplyr::case_when(
                grepl("mt", p) ~ "mt",
                grepl("rp", p) ~ "rp",
                TRUE ~ tolower(gsub("[^[:alnum:]]", "", p))
            )
        })

        return(paste(sort(unique(names)), collapse = "_"))
    }
    # Handle single patterns
    dplyr::case_when(
        grepl("mt", pat_lower) ~ "mt",
        grepl("rp", pat_lower) ~ "rp",
        TRUE ~ {
            clean <- gsub("[^[:alnum:]]", "_", pat_lower)
            clean <- gsub("_+", "_", clean) # Collapse multiple underscores
            clean <- gsub("^_+|_+$", "", clean) # Trim leading/trailing underscores
            tolower(clean)
        }
    )
}
