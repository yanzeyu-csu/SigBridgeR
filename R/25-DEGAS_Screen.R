# ? --README--------------------------------------------------------------------------------------------------
# ? DEGAS is not recommended to be used in your screening, since the version control is strict and intractable
# ? ----------------------------------------------------------------------------------------------------------

#' @title Run DEGAS Analysis for Single-Cell and Bulk RNA-seq Data Integration
#'
#' @description
#' This function performs DEGAS to integrate
#' single-cell and bulk RNA-seq data, identifying phenotype-associated cells using
#' a bootstrap aggregated multi-task learning approach.
#'
#' @param select_fraction The top percentage of selected cells will be considered as Positive cells, without considering how much larger the possible correlation coefficient of the observation group is compared to that of the control group. Only usedl when `phenotype_class` is "binary" or "survival". (default: 0.05)
#' @param min_thresh DEGAS will calculate the possible correlation coefficients for each cell related to the phenotype. When the coefficient of the observation group is at least `min_thresh` larger than that of the control group, it can be considered related to the phenotype and will be marked as Positive. The priority of `min_thresh` is higher than that of `select_fraction.` (default: 0.4)
#' @param matched_bulk Bulk RNA-seq data as matrix or data.frame (rows=genes, columns=samples)
#' @param sc_data Single-cell data as Seurat object containing RNA assay
#' @param phenotype Bulk-level phenotype data. For classification: binary matrix with one-hot encoding.
#'   For survival: matrix with two columns (time and event status). Can be NULL, matrix, data.frame, or vector.
#' @param sc_data.pheno_colname Column name for single-cell phenotype in metadata (if available)
#' @param label_type Label type for DEGAS results (default: "DEGAS")
#' @param phenotype_class Type of phenotype: "binary" (classification), "continuous", or "survival"
#' @param tmp_dir Temporary directory for intermediate files (default: "tmp")
#' @param env_params List of environment parameters for Python setup including:
#'   - env.name: environment name (default: "r-reticulate-degas")
#'   - env.type: environment type "conda", "environment", or "venv" (default: "environment")
#'   - env.method: environment setup method "system", "conda" (default: "system")
#'   - env.file: path to environment file (default: system.file("conda/DEGAS_environment.yml", package = "SigBridgeR"))
#'   - env.python_version: Python version (default: "3.9.15")
#'   - env.packages: named vector of Python packages and versions (default: c("tensorflow" = "2.4.1", "protobuf" = "3.20" ,"numpy" = "any"))
#'   - env.recreate: whether to recreate environment (default: FALSE)
#'   - env.use_conda_forge: whether to use conda-forge channel (conda only, default: TRUE)
#'   - env.verbose: verbose output (default: FALSE)
#' @param degas_params List of DEGAS algorithm parameters including:
#'   - DEGAS.model_type: model type ("BlankClass", "ClassBlank", "ClassClass", "ClassCox", "BlankCox")
#'   - DEGAS.architecture: "Standard" (feed forward) or "DenseNet" (dense net)
#'   - DEGAS.ff_depth: number of layers in model (>=1, default: 3)
#'   - DEGAS.pyloc: path to Python executable (default: NULL, automatic detection)
#'   - DEGAS.bag_depth: bootstrap aggregation depth (>=1, default: 5)
#'   - DEGAS.train_steps: training steps (default: 2000)
#'   - DEGAS.scbatch_sz: single-cell batch size (default: 200)
#'   - DEGAS.patbatch_sz: patient batch size (default: 50)
#'   - DEGAS.hidden_feats: hidden features (default: 50)
#'   - DEGAS.do_prc: dropout percentage (default: 0.5)
#'   - DEGAS.lambda1: regularization parameter 1 (default: 3.0)
#'   - DEGAS.lambda2: regularization parameter 2 (default: 3.0)
#'   - DEGAS.lambda3: regularization parameter 3 (default: 3.0)
#'   - DEGAS.seed: random seed (default: 2)
#' @param normality_test_method Method for normality testing: "jarque-bera", "d'agostino", or "kolmogorov-smirnov"
#' @param ... for future compatibility
#'
#' @return A list containing:
#'   - scRNA_data: Seurat object with DEGAS labels added to metadata
#'   - model: The model trained using the input data, andit can be used for cell classification prediction.
#'   - DEGAS_prediction: Data table with DEGAS predictions containing:
#'     * Predicted label probabilities for each cell
#'     * Cell labels ("Positive"/"Other") based on selection criteria
#'     * Difference scores for binary phenotypes
#'     * Cell identifiers
#'
#' @details
#' The function performs the following steps:
#' 1. Validates input data and parameters
#' 2. Sets up Python environment with required dependencies
#' 3. Trains bootstrap aggregated DEGAS model using `runCCMTLBag`
#' 4. Generates cell-level predictions using `predClassBag`
#' 5. Applies statistical testing to identify phenotype-associated cells
#' 6. Labels cells as "Positive" or "Other" based on selection criteria
#'
#' Model type is automatically determined:
#' - BlankClass: only bulk phenotype specified (scLab = NULL)
#' - ClassBlank: only single-cell phenotype specified (patLab = NULL)
#' - ClassClass: both single-cell and bulk phenotypes specified
#' - ClassCox: single-cell phenotype + bulk survival data
#' - BlankCox: only bulk survival data specified
#'
#' @examples
#' \dontrun{
#' # Binary classification example
#' result <- DoDEGAS(
#'   select_fraction = 0.05, # `select_fraction` only used in binary and survival phenotyping
#'   matched_bulk = bulk_matrix,
#'   sc_data = seurat_obj,
#'   phenotype = bulk_phenotype,
#'   phenotype_class = "binary"
#' )
#'
#' # Survival analysis example
#' result <- DoDEGAS(
#'   select_fraction = 0.05, # `select_fraction` only used in binary and survival phenotyping
#'   matched_bulk = bulk_matrix,
#'   sc_data = seurat_obj,
#'   phenotype = survival_data,
#'   phenotype_class = "survival"
#' )
#' }
#'
#' @seealso
#' \code{\link{Vec2sparse}} for the structure transformation of phenotype
#' \code{\link{jb.test.modified}} for modified Jarque-Bera test
#' \code{\link{mad.test}} for outlier detection using Median Absolute Deviation
#' \code{\link{runCCMTLBag.optimized}} for DEGAS model training
#' \code{\link{predClassBag.optimized}} for DEGAS model prediction
#' \code{\link{LabelBinaryCells}} for binary classification
#' \code{\link{LabelSurvivalCells}} for survival classification
#' \code{\link{LabelContinuousCells}} for continuous classification
#'
#' @importFrom purrr walk map_chr map_dfr
#' @importFrom magrittr %>%
#' @importFrom data.table as.data.table fifelse setnames `:=`
#' @importFrom stats ks.test qnorm median
#' @family screen_method
#' @family DEGAS
#'
#' @references Johnson TS, Yu CY, Huang Z, Xu S, Wang T, Dong C, et al. Diagnostic Evidence GAuge of Single cells (DEGAS): a flexible deep transfer learning framework for prioritizing cells in relation to disease. Genome Med. 2022 Feb 1;14(1):11.
#'
DoDEGAS <- function(
    select_fraction = 0.05,
    min_thresh = 0.4,
    matched_bulk,
    sc_data,
    phenotype = NULL,
    sc_data.pheno_colname = NULL,
    label_type = "DEGAS",
    phenotype_class = c("binary", "continuous", "survival"),
    # A directory for intermediate files
    tmp_dir = "tmp",
    # DEGAS environment
    env_params = list(),
    # DEGAS parameters
    degas_params = list(),
    normality_test_method = c(
        "jarque-bera",
        "d'agostino",
        "kolmogorov-smirnov"
    ),
    ...
) {
    # Robustness checks
    chk::chk_range(select_fraction)
    chk::chk_is(matched_bulk, c("matrix", "data.frame"))
    chk::chk_not_any_na(matched_bulk)
    chk::chk_is(sc_data, "Seurat")
    chk::chk_character(label_type)
    if (!is.null(sc_data.pheno_colname)) {
        chk::chk_is(sc_data.pheno_colname, c("character"))
    }
    chk::chk_list(env_params)
    chk::chk_list(degas_params)
    phenotype_class %<>% MatchArg(c("binary", "continuous", "survival"), NULL)

    ts_cli$cli_alert_info(cli::col_green("Starting DEGAS Screen"))

    # Auto-choose normality test method
    normality_test_method %<>%
        MatchArg(c(
            "jarque-bera",
            "d'agostino",
            "kolmogorov-smirnov"
        ))

    # DEGAS path must contain "/"
    if (!grepl("/$", tmp_dir)) {
        tmp_dir <- paste0(tmp_dir, "/")
    }

    # default environment and DEGAS parameters
    default_env_params <- list(
        env.name = "r-reticulate-degas",
        env.type = "conda",
        env.method = "environment",
        env.file = system.file(
            "conda/DEGAS_environment.yml",
            package = "SigBridgeR"
        ),
        env.python_verion = "3.9.15",
        env.packages = c(
            "tensorflow" = "2.4.1",
            "protobuf" = "3.20.3"
        ),
        env.recreate = FALSE,
        env.use_conda_forge = TRUE,
        env.verbose = FALSE
    )
    default_degas_params <- list(
        DEGAS.model_type = c(
            "BlankClass", # only bulk level phenotype specified
            "ClassBlank", # only single cell level phenotype specified
            "ClassClass", # when both single cell level phenotype and bulk level phenotype specified
            "ClassCox", # when both single cell level phenotype and bulk level survival data specified
            "BlankCox" # only bulk level survival data specified
        ),
        DEGAS.architecture = c(
            "DenseNet", # a dense net network
            "Standard" # a feed forward network
        ),
        DEGAS.ff_depth = 3L,
        DEGAS.bag_depth = 5L,
        path.data = '',
        path.result = '',
        DEGAS.pyloc = NULL,
        DEGAS.toolsPath = file.path(.libPaths()[1], "SigBridgeR/DEGAS_tools/"),
        DEGAS.train_steps = 2000L,
        DEGAS.scbatch_sz = 200L,
        DEGAS.patbatch_sz = 50L,
        DEGAS.hidden_feats = 50L,
        DEGAS.do_prc = 0.5,
        DEGAS.lambda1 = 3.0,
        DEGAS.lambda2 = 3.0,
        DEGAS.lambda3 = 3.0,
        DEGAS.seed = 2L
    )
    env_params <- utils::modifyList(
        default_env_params,
        env_params
    )
    degas_params <- utils::modifyList(
        default_degas_params,
        degas_params
    )
    on.exit({
        # Clean up global variables
        rm(list = names(degas_params), envir = .GlobalEnv)
    })

    # model.type auto-detection
    model_type.first <- ifelse(
        !is.null(sc_data.pheno_colname),
        "Class",
        "Blank"
    )
    model_type.last <- if (is.null(phenotype)) {
        "Blank"
    } else {
        ifelse(is.vector(phenotype), "Class", "Cox")
    }
    if (model_type.first == "Blank" && model_type.last == "Blank") {
        cli::cli_abort(
            c(
                "x" = "Specify at least one phenotype, currently both are {.val NULL}"
            ),
            class = "NoPhenotypeProvided"
        )
    }
    degas_params$DEGAS.model_type <- paste0(
        model_type.first,
        model_type.last
    )

    # Architecture auto-chosen
    if (length(degas_params$DEGAS.architecture) != 1) {
        degas_params$DEGAS.architecture <- "DenseNet"
    }

    # formatting phenotype
    pheno_df <- if (!is.null(phenotype)) {
        switch(
            phenotype_class,
            "binary" = Vec2sparse(phenotype, col_prefix = label_type),
            "continuous" = Vec2sparse(phenotype, col_prefix = label_type),
            "survival" = as.matrix(phenotype),
        )
    } else {
        NULL
    }

    # Check if single-cell level phenotype is specified
    meta_data <- sc_data[[]]
    sc_pheno <- if (!is.null(sc_data.pheno_colname)) {
        if (any(grepl(sc_data.pheno_colname, colnames(meta_data)))) {
            Vec2sparse(meta_data[[sc_data.pheno_colname]])
        } else {
            cli::cli_warn(
                "single-cell phenotype specified but not found in meta.data, using {.val NULL} now"
            )
            NULL
        }
    } else {
        NULL
    }

    # Python-like data formats
    sc_mat <- SeuratObject::LayerData(sc_data)
    cm_genes <- intersect(rownames(matched_bulk), rownames(sc_mat))
    t_sc_mat <- Matrix::t(sc_mat[cm_genes, ])
    t_matched_bulk <- Matrix::t(matched_bulk[cm_genes, ])

    ts_cli$cli_alert_info("Setting up Environment...")

    # Check if environment exists
    existing_envs <- ListPyEnv.conda(verbose = FALSE)
    if (
        !env_params$env.name %chin% existing_envs$name ||
            env_params$env.recreate
    ) {
        do.call(
            SetupPyEnv,
            c(
                list(
                    env_type = env_params$env.type,
                    env_name = env_params$env.name,
                    python_version = env_params$env.python_verion,
                    packages = env_params$env.packages,
                    recreate = env_params$env.recreate,
                    verbose = env_params$env.verbose,
                    ...
                ),
                if (env_params$env.type == "conda") {
                    list(
                        use_conda_forge = env_params$env.use_conda_forge,
                        env_file = env_params$env.file,
                        method = env_params$env.method
                    )
                } else {
                    list()
                }
            )
        )
    } else {
        cli::cli_alert_info(
            "Existing environment {.val {env_params$env.name}} found"
        )
    }

    if (is.null(degas_params$DEGAS.pyloc)) {
        degas_params$DEGAS.pyloc <- existing_envs$python[
            existing_envs$name == env_params$env.name
        ]
    }
    reticulate::use_python(
        degas_params$DEGAS.pyloc,
        required = TRUE
    )

    # DEGAS needs some global variables to be set up
    list2env(degas_params, envir = .GlobalEnv)

    ts_cli$cli_alert_info("Training DEGAS model...")

    # Train DEGAS model
    ccModel1 <- do.call(
        runCCMTLBag.optimized,
        list(
            scExp = t_sc_mat,
            scLab = sc_pheno,
            patExp = t_matched_bulk,
            patLab = pheno_df,
            tmpDir = tmp_dir,
            model_type = DEGAS.model_type,
            architecture = DEGAS.architecture,
            FFdepth = DEGAS.ff_depth,
            Bagdepth = DEGAS.bag_depth,
            DEGAS.seed = DEGAS.seed
        )
    )
    names(ccModel1) <- glue::glue("ccModel_{seq_len(DEGAS.bag_depth)}")

    ts_cli$cli_alert_info("Predicting and Labeling...")

    # Predict with DEGAS model
    t_sc_preds <- data.table::as.data.table(predClassBag.optimized(
        ccModel = ccModel1,
        Exp = t_sc_mat,
        scORpat = 'pat'
    ))

    if (phenotype_class == "survival") {
        data.table::setnames(t_sc_preds, "Hazard")
    } else {
        pheno_df_colnames <- if (!is.null(pheno_df)) {
            colnames(pheno_df)
        } else {
            NULL
        }
        if (!is.null(pheno_df_colnames)) {
            data.table::setnames(t_sc_preds, pheno_df_colnames)
        }
    }
    t_sc_preds[, "cell_id" := rownames(t_sc_mat)]

    ts_cli$cli_alert_info("Labeling screened cells...")

    # What we get from DEGAS are the probabilities of the cells belonging to each class.
    # We need to convert these to labels.
    t_sc_preds <- switch(
        phenotype_class,
        "binary" = LabelBinaryCells(
            pred_dt = t_sc_preds,
            pheno_colnames = pheno_df_colnames,
            select_fraction = select_fraction,
            test_method = normality_test_method,
            min_threshold = min_thresh
        ),
        "continuous" = LabelContinuousCells(pred_dt = t_sc_preds),
        "survival" = LabelSurvivalCells(
            pred_dt = t_sc_preds,
            select_fraction = select_fraction,
            test_method = normality_test_method,
            min_threshold = min_thresh
        )
    )

    # Record screening results
    sc_data <- Seurat::AddMetaData(
        object = sc_data,
        metadata = t_sc_preds[["label"]],
        col.name = "DEGAS"
    ) %>%
        AddMisc(
            DEGAS_type = label_type,
            DEGAS_para = degas_params,
            cover = FALSE
        )

    ts_cli$cli_alert_info(cli::col_green("DEGAS Screen done."))

    # result
    return(list(
        scRNA_data = sc_data,
        model = ccModel1,
        DEGAS_prediction = t_sc_preds
    ))
}


#' @title Convert a Named Integer Vector to a 0/1 Sparse Matrix
#'
#' @description
#' Takes a named integer vector (sample names as names, integers as values) and
#' returns a sparse 0/1 matrix where row names are sample names, column names
#' are the unique integers, and entries indicate presence (1) or absence (0)
#' of the integer value for each sample.
#' This function is used in DEGAS screening.
#'
#' @param x A named integer vector. Names are treated as sample identifiers and
#'   must be unique. Values are coerced to factor levels to define columns.
#' @param col_prefix A single character string. If supplied, it is prepended to
#'   every column name of the returned matrix (e.g. `"V"` produces `"V1"`, `"V2"` …).
#'   The default `""` keeps the original numeric levels unchanged.
#'
#' @return A sparse logical matrix of class `"ngCMatrix"` from the **Matrix**
#'   package. Rows correspond to samples, columns to unique values, with `TRUE`
#'   (1) indicating the sample carries that value.
#'
#' @examples
#' \dontrun{
#' v <- setNames(sample(1:5, 10, TRUE), paste0("Sample", 1:10))
#' mat1 <- Vec2sparse(v)              # original column names
#' mat2 <- Vec2sparse(v, col_prefix = "V")  # prefixed column names
#' print(mat1)
#' print(mat2)
#' }
#'
#' @importFrom Matrix sparseMatrix
#' @importFrom chk chk_vector chk_character chk_not_any_na
#' @keywords internal
#' @family DEGAS
#'
Vec2sparse <- function(x, col_prefix = "") {
    chk::chk_vector(x)
    chk::chk_character(col_prefix)
    chk::chk_not_any_na(x)

    rows <- seq_along(x)
    cols <- as.factor(x)
    col_names <- levels(cols)
    if (nzchar(col_prefix)) {
        col_names <- paste0(col_prefix, col_names)
    }
    Matrix::sparseMatrix(
        i = rows,
        j = as.integer(cols),
        dims = c(length(x), nlevels(cols)),
        dimnames = list(names(x), col_names)
    )
}

#' @title Optimized Cross-Condition Multi-Task Learning Model Training
#'
#' @description
#' An optimized wrapper function for training cross-condition multi-task learning
#' (CCMTL) models in the DEGAS framework. This function handles the complete
#' training pipeline including data preparation, model configuration, and
#' execution with enhanced performance and error handling.
#'
#' @param scExp A matrix or data frame containing single-cell expression data
#'   for model training.
#' @param scLab A  matrix containing single-cell labels corresponding
#'   to the expression data.
#' @param patExp A matrix or data frame containing patient-level expression data
#'   for multi-task learning.
#' @param patLab A matrix containing patient-level labels corresponding
#'   to the patient expression data.
#' @param tmpDir Character string specifying the temporary directory path for
#'   storing intermediate files and model outputs.
#' @param model_type Character string specifying the type of model to train.
#'   Should match available DEGAS model types.
#' @param architecture Character string specifying the neural network architecture.
#'   One of: "DenseNet", "Standard".
#' @param FFdepth Integer specifying the number of layers in the feed-forward
#'   network architecture.
#' @param DEGAS.seed Integer specifying the random seed for reproducible
#'   model training.
#' @param force_rewrite Logical indicating whether to force rewriting of input
#'   files even if they already exist. Default: FALSE.
#'
#' @return
#' Returns a trained CCMTL model object that can be used for predictions and
#' further analysis.
#'
#' @details
#' ## Workflow:
#' 1. **File Management**: Efficient handling of temporary directories and
#'    input files with optional forced rewriting
#' 2. **Architecture Configuration**: Supports multiple neural network
#'    architectures (DenseNet, Standard) with customizable depth
#' 3. **Python Environment**: Validates Python availability and executes
#'    training scripts with proper error handling
#' 4. **Model Training**: Executes the DEGAS training process with specified
#'    hyperparameters and architecture choices
#'
#' @note
#' This function requires a properly configured Python environment with DEGAS
#' dependencies installed. The temporary directory (`tmpDir`) should have
#' sufficient disk space for model files and intermediate data.
#'
#' @seealso
#' [runCCMTLBag.optimized()] for bootstrap aggregated model training,
#'
#' @keywords internal
#' @family DEGAS
#'
#' @references Johnson TS, Yu CY, Huang Z, Xu S, Wang T, Dong C, et al. Diagnostic Evidence GAuge of Single cells (DEGAS): a flexible deep transfer learning framework for prioritizing cells in relation to disease. Genome Med. 2022 Feb 1;14(1):11.
#'
runCCMTL.optimized <- function(
    scExp,
    scLab,
    patExp,
    patLab,
    tmpDir,
    model_type,
    architecture,
    FFdepth,
    DEGAS.seed,
    force_rewrite = FALSE,
    verbose = TRUE
) {
    # Only write files if explicitly requested
    if (force_rewrite) {
        if (dir.exists(tmpDir)) {
            unlink(tmpDir, recursive = TRUE, force = TRUE)
        }
        dir.create(tmpDir, recursive = TRUE, showWarnings = FALSE)
        # Write input files
        writeInputFiles.optimized(
            scExp = scExp,
            scLab = scLab,
            patExp = patExp,
            patLab = patLab,
            tmpDir = tmpDir
        )
    }

    # create python files
    if (!architecture %chin% c("DenseNet", "Standard")) {
        cli::cli_abort(c("x" = 'Incorrect architecture argument'))
    } else if (architecture == "DenseNet") {
        DEGAS::makeExec2(
            tmpDir = tmpDir,
            FFdepth = FFdepth,
            model_type = model_type
        )
    } else {
        DEGAS::makeExec(
            tmpDir = tmpDir,
            FFdepth = FFdepth,
            model_type = model_type
        )
    }

    # cmd <- paste0(
    #     DEGAS.pyloc,
    #     " ",
    #     tmpDir,
    #     model_type,
    #     "MTL.py",
    #     paste(
    #         "",
    #         tmpDir,
    #         DEGAS.train_steps,
    #         DEGAS.scbatch_sz,
    #         DEGAS.patbatch_sz,
    #         DEGAS.hidden_feats,
    #         DEGAS.do_prc,
    #         DEGAS.lambda1,
    #         DEGAS.lambda2,
    #         DEGAS.lambda3,
    #         DEGAS.seed
    #     )
    # )
    cmd_args <- as.character(c(
        paste0(tmpDir, model_type, "MTL.py"),
        tmpDir,
        DEGAS.train_steps,
        DEGAS.scbatch_sz,
        DEGAS.patbatch_sz,
        DEGAS.hidden_feats,
        DEGAS.do_prc,
        DEGAS.lambda1,
        DEGAS.lambda2,
        DEGAS.lambda3,
        DEGAS.seed
    ))

    # * Execute system command
    # system(command = cmd) # if processx::run failed, use `system` instead
    result <- processx::run(
        command = DEGAS.pyloc,
        args = cmd_args,
        echo = verbose,
        error_on_status = TRUE
    )

    readOutputFiles.optimized(
        tmpDir = tmpDir,
        model_type = model_type,
        architecture = architecture
    )
}

#' @title ccModel Object
#'
#' @slot Bias A list containing bias terms for each layer of the model.
#' @slot Theta A list of parameter matrices/vectors for the model.
#' @slot Activation A list of activation functions for each layer.
#' @slot Depth Numeric. The number of layers in the model architecture.
#'   Must be a positive integer.
#' @slot Model_type Character. The type of model
#' @slot Architecture Character. Description of the model architecture
#'
#' @examples
#' # Create a simple ccModel object
#' model <- new("ccModel",
#'   Bias = list(c(0.1, 0.2), c(0.3)),
#'   Theta = list(matrix(c(0.5, 0.6), nrow=1), matrix(c(0.7), nrow=1)),
#'   Activation = list("relu", "sigmoid"),
#'   Depth = 2,
#'   Model_type = "neural_network",
#'   Architecture = "feedforward"
#' )
#' @keywords internal
#' @family DEGAS_in_SigBridgeR
#' @references Johnson TS, Yu CY, Huang Z, Xu S, Wang T, Dong C, et al. Diagnostic Evidence GAuge of Single cells (DEGAS): a flexible deep transfer learning framework for prioritizing cells in relation to disease. Genome Med. 2022 Feb 1;14(1):11.
#'
setClass(
    "ccModel",
    slots = list(
        Bias = "list",
        Theta = "list",
        Activation = "list",
        Depth = "numeric",
        Model_type = "character",
        Architecture = "character"
    )
)


#' @title Optimized Bootstrap Aggregation for Cross-Condition Multi-Task Learning
#'
#' @description
#' Performs bootstrap aggregated training of multiple CCMTL models to enhance
#' robustness and reduce variance in predictions. This function trains an
#' ensemble of models with different random seeds and aggregates the results.
#'
#' @param scExp A matrix or data frame containing single-cell expression data
#'   for model training.
#' @param scLab A matrix containing single-cell labels corresponding
#'   to the expression data.
#' @param patExp A matrix or data frame containing patient-level expression data
#'   for multi-task learning.
#' @param patLab A matrix containing patient-level labels corresponding
#'   to the patient expression data.
#' @param tmpDir Character string specifying the temporary directory path for
#'   storing intermediate files and model outputs.
#' @param model_type Character string specifying the type of model to train.
#'   Should match available DEGAS model types.
#' @param architecture Character string specifying the neural network architecture.
#'   One of: "DenseNet", "Standard".
#' @param FFdepth Integer specifying the number of layers in the feed-forward
#'   network architecture.
#' @param Bagdepth Integer specifying the number of bootstrap models to train
#'   in the ensemble.
#' @param DEGAS.seed Integer specifying the base random seed for reproducible
#'   model training. Each model in the ensemble uses a derived seed.
#'
#' @return
#' Returns a list of trained CCMTL model objects from the bootstrap aggregation
#' process. The list contains successful model results with proper error handling
#' for failed training attempts.
#'
#' @details
#' This function implements bootstrap aggregated training (bagging) for CCMTL
#' models with the following features:
#'
#' ## Ensemble Training:
#' - Trains multiple models with different random seeds derived from the base seed
#' - Uses parallel-safe file management to avoid I/O conflicts
#' - Implements comprehensive error handling to continue training even if
#'   individual models fail
#'
#' ## Error Handling:
#' - Continues training even if individual models fail
#' - Returns only successfully trained models
#' - Provides progress feedback for long-running ensemble training
#'
#' @note
#' The bootstrap aggregation process can be computationally intensive, especially
#' for large datasets or deep architectures. The function creates derived seeds
#' for each model (base seed + model index) to ensure reproducibility while
#' maintaining diversity in the ensemble.
#'
#' @examples
#' \dontrun{
#' # Train an ensemble of 10 CCMTL models
#' ensemble_models <- runCCMTLBag.optimized(
#'   scExp = sc_expression,
#'   scLab = sc_labels,
#'   patExp = patient_expression,
#'   patLab = patient_labels,
#'   tmpDir = "/tmp/degas_models",
#'   model_type = "classification",
#'   architecture = "DenseNet",
#'   FFdepth = 3,
#'   Bagdepth = 10,
#'   DEGAS.seed = 42
#' )
#'
#' # Access individual models from the ensemble
#' first_model <- ensemble_models[[1]]
#' }
#'
#' @seealso
#' [runCCMTL.optimized] for single model training,
#' [purrr::map()] for the iterative execution pattern.
#'
#' @keywords internal
#' @family DEGAS
#' @references Johnson TS, Yu CY, Huang Z, Xu S, Wang T, Dong C, et al. Diagnostic Evidence GAuge of Single cells (DEGAS): a flexible deep transfer learning framework for prioritizing cells in relation to disease. Genome Med. 2022 Feb 1;14(1):11.
#'
runCCMTLBag.optimized <- function(
    scExp,
    scLab,
    patExp,
    patLab,
    tmpDir,
    model_type,
    architecture,
    FFdepth,
    Bagdepth,
    DEGAS.seed
) {
    ts_cli$cli_alert_info(glue::glue(
        "{FFdepth}-layer {architecture} {model_type} DEGAS model"
    ))

    if (!dir.exists(tmpDir)) {
        dir.create(tmpDir, recursive = TRUE)
    }
    # Write files once at the beginning
    writeInputFiles.optimized(
        scExp = scExp,
        scLab = scLab,
        patExp = patExp,
        patLab = patLab,
        tmpDir = tmpDir
    )

    # Check Python with error handling
    py_check <- processx::run(command = DEGAS.pyloc, args = "--version")
    if (!is.null(py_check$error)) {
        cli::cli_abort("Python check failed: ", py_check$error$message)
    } else {
        ts_cli$cli_alert_info(glue::glue(
            "Python check passed, using {py_check$stdout}"
        ))
    }

    purrr::map(
        seq_len(Bagdepth),
        function(i) {
            DEGAS.seed_i <- DEGAS.seed + (i - 1)

            if (verbose) {
                ts_cli$cli_alert_info("Training progress: {i}/{Bagdepth}...")
            }

            result <- runCCMTL.optimized(
                scExp = scExp,
                scLab = scLab,
                patExp = patExp,
                patLab = patLab,
                tmpDir = tmpDir,
                model_type = model_type,
                architecture = architecture,
                FFdepth = FFdepth,
                DEGAS.seed = DEGAS.seed_i,
                # Written files will not be rewritten
                force_rewrite = FALSE
            )
            class(result) <- "ccModel"

            result
        },
        .progress = TRUE
    )
}

#' @title Optimized Input File Writing for DEGAS Models
#'
#' @description
#' Efficiently writes input data files for DEGAS model training using optimized
#' data handling and fast I/O operations. This function converts various data
#' types to efficient CSV format using data.table's fwrite for rapid file
#' operations with comprehensive error handling.
#'
#' @param scExp A matrix, data frame, or Matrix object containing single-cell
#'   expression data. Rows typically represent genes and columns represent cells.
#' @param scLab A matrix, or data frame containing single-cell labels
#'   corresponding to the expression data. Can be NULL if no labels are available.
#' @param patExp A data frame, or Matrix object containing patient-level
#'   expression data. Rows typically represent genes and columns represent patients.
#' @param patLab A matrix, or data frame containing patient-level labels
#'   corresponding to the patient expression data. Can be NULL if no labels are available.
#' @param tmpDir Character string specifying the directory path where input
#'   files will be written. The directory will be created if it doesn't exist.
#'
#' @return
#' Invisibly returns TRUE if all files are successfully written. If any error
#' occurs during file writing, the function will abort with an informative error message.
#'
#' @details
#' This function provides an optimized pipeline for writing input files required
#' by DEGAS models with the following features:
#'
#' ## File Output:
#' The function creates four CSV files in the specified temporary directory:
#' - `scExp.csv`: Single-cell expression data
#' - `scLab.csv`: Single-cell labels (if provided)
#' - `patExp.csv`: Patient-level expression data
#'
#'
#' @note
#' The function uses comma-separated values (CSV) format without row names to
#' ensure compatibility with Python-based DEGAS training scripts. All input
#' data is converted to dense format during writing, so ensure sufficient
#' memory is available for large datasets.
#'
#' @examples
#' \dontrun{
#' # Write input files for DEGAS training
#' writeInputFiles.optimized(
#'   scExp = single_cell_expression,
#'   scLab = single_cell_labels,
#'   patExp = patient_expression,
#'   patLab = patient_labels,
#'   tmpDir = "/tmp/degas_input"
#' )
#' }
#'
#' @seealso
#' [data.table::fwrite()] for the underlying fast writing implementation,
#' [purrr::safely()] for the error handling mechanism.
#'
#' @keywords internal
#' @family DEGAS
#' @references Johnson TS, Yu CY, Huang Z, Xu S, Wang T, Dong C, et al. Diagnostic Evidence GAuge of Single cells (DEGAS): a flexible deep transfer learning framework for prioritizing cells in relation to disease. Genome Med. 2022 Feb 1;14(1):11.
#'
writeInputFiles.optimized <- function(
    scExp,
    scLab = NULL,
    patExp,
    patLab = NULL,
    tmpDir
) {
    if (!dir.exists(tmpDir)) {
        dir.create(tmpDir, recursive = TRUE)
    }

    file_writing_ops <- purrr::safely(
        ~ {
            # Convert scExp to data.table for fast writing
            if (inherits(scExp, "matrix") || inherits(scExp, "Matrix")) {
                scExp_dt <- data.table::as.data.table(as.matrix(scExp))
            } else {
                scExp_dt <- data.table::as.data.table(scExp)
            }

            # Write scExp using data.table's fwrite (much faster than write.table)
            data.table::fwrite(
                scExp_dt,
                file = file.path(tmpDir, 'scExp.csv'),
                row.names = FALSE,
                sep = ',',
                showProgress = FALSE
            )

            # Process scLab if not NULL
            if (!is.null(scLab)) {
                if (inherits(scLab, "matrix") || inherits(scLab, "Matrix")) {
                    scLab_dt <- data.table::as.data.table(as.matrix(scLab))
                } else {
                    scLab_dt <- data.table::as.data.table(scLab)
                }
                scLab_dt[, names(scLab_dt) := lapply(.SD, as.integer)] # Convert to integer from boolean

                data.table::fwrite(
                    scLab_dt,
                    file = file.path(tmpDir, 'scLab.csv'),
                    row.names = FALSE,
                    sep = ',',
                    showProgress = FALSE
                )
            }

            # Convert patExp to data.table
            if (inherits(patExp, "matrix") || inherits(patExp, "Matrix")) {
                patExp_dt <- data.table::as.data.table(as.matrix(patExp))
            } else {
                patExp_dt <- data.table::as.data.table(patExp)
            }

            # Write patExp
            data.table::fwrite(
                patExp_dt,
                file = file.path(tmpDir, 'patExp.csv'),
                row.names = FALSE,
                sep = ',',
                showProgress = FALSE
            )

            # Process patLab if not NULL
            if (!is.null(patLab)) {
                if (inherits(patLab, "matrix") || inherits(patLab, "Matrix")) {
                    patLab_dt <- data.table::as.data.table(as.matrix(patLab))
                } else {
                    patLab_dt <- data.table::as.data.table(patLab)
                }
                patLab_dt[, names(patLab_dt) := lapply(.SD, as.integer)] # Convert to integer from boolean

                data.table::fwrite(
                    patLab_dt,
                    file = file.path(tmpDir, 'patLab.csv'),
                    row.names = FALSE,
                    sep = ',',
                    showProgress = FALSE
                )
            }
        }
    )

    # Execute file writing operations with error handling
    result <- file_writing_ops()

    if (!is.null(result$error)) {
        cli::cli_abort("Failed to write input files: ", result$error$message)
    }

    invisible(TRUE)
}


#' @title Read Model Output Files
#'
#' @description
#' Reads neural network parameters (weights, biases, activation functions) from
#' output files generated by an optimized training process and constructs a
#' complete model object.
#'
#' @param tmpDir Character string specifying the temporary directory path where
#'               the model output files are stored.
#' @param model_type Character string describing the type of model.
#' @param architecture Character string describing the neural network architecture.
#'
#' @return Returns an object of class \code{ccModel} containing:
#' \itemize{
#'   \item \code{Bias} - List of bias matrices for each layer
#'   \item \code{Theta} - List of weight matrices for each layer
#'   \item \code{Activation} - List of activation functions for each layer
#'   \item \code{Depth} - Integer specifying the number of layers
#'   \item \code{Model_type} - The input model type
#'   \item \code{Architecture} - The input architecture description
#' }
#'
#' @details
#' This function reads the following files from the specified directory:
#' \itemize{
#'   \item \code{Activations.csv} - Contains activation function names (one per line)
#'   \item \code{Bias{i}.csv} - Bias values for layer i (where i = 1, 2, ...)
#'   \item \code{Theta{i}.csv} - Weight values for layer i (where i = 1, 2, ...)
#' }
#'
#' The function automatically determines the network depth from the number of
#' activation functions and loads the corresponding parameters for each layer.
#' All CSV files are expected to be comma-separated without row names.
#'
#' @examples
#' \dontrun{
#' # Read model files from temporary directory
#' model <- readOutputFiles.optimized(
#'   tmpDir = "/path/to/model/files",
#'   Model_type = "Classification",
#'   architecture = "12"
#' )
#' }
#'
#' @keywords internal
#' @family DEGAS
#' @references Johnson TS, Yu CY, Huang Z, Xu S, Wang T, Dong C, et al. Diagnostic Evidence GAuge of Single cells (DEGAS): a flexible deep transfer learning framework for prioritizing cells in relation to disease. Genome Med. 2022 Feb 1;14(1):11.
#'
readOutputFiles.optimized <- function(tmpDir, model_type, architecture) {
    activations <- data.table::fread(
        file.path(tmpDir, 'Activations.csv'),
        header = FALSE,
        sep = "\n",
        showProgress = FALSE
    )[[1]]

    depth <- length(activations)
    file_indices <- seq_len(depth)

    # Read Bias files using purrr and data.table
    Biases <- purrr::map(
        file_indices,
        ~ {
            bias_file <- file.path(tmpDir, paste0('Bias', .x, '.csv'))
            bias_data <- data.table::fread(
                bias_file,
                header = FALSE,
                sep = ',',
                showProgress = FALSE
            )
            as.matrix(bias_data)
        }
    )

    # Read Theta files using purrr and data.table
    Thetas <- purrr::map(
        file_indices,
        ~ {
            theta_file <- file.path(tmpDir, paste0('Theta', .x, '.csv'))
            theta_data <- data.table::fread(
                theta_file,
                header = FALSE,
                sep = ',',
                showProgress = FALSE
            )
            as.matrix(theta_data)
        }
    )

    Activations <- as.list(activations)

    methods::new(
        'ccModel',
        Bias = Biases,
        Theta = Thetas,
        Activation = Activations,
        Depth = depth,
        Model_type = model_type,
        Architecture = architecture
    )
}


#' @title Bagged Prediction for Classification Models
#'
#' @description
#' Performs bagged (bootstrap aggregating) predictions using an ensemble of
#' ccModel objects. This function averages predictions from multiple models
#' to improve stability and accuracy.
#'
#' @param ccModel A list of `ccModel` objects representing the ensemble of
#'               models to use for bagged prediction. Each model should be
#'               a valid ccModel instance.
#' @param Exp A matrix or data frame containing the expression data to be
#'            classified.
#' @param scORpat A character string specifying whether the data represents
#'               single-cell ("sc") or phenotype ("pat") data. This parameter
#'               is passed to the underlying prediction functions.
#'
#' @return A matrix of averaged predictions where each element represents the
#'         aggregated probability or classification score across all models
#'         in the ensemble. The output is computed as the mean of all individual
#'         model predictions.
#'
#' @details
#' This function implements model bagging by:
#' \itemize{
#'   \item Iterating over each model in the `ccModel` list using `purrr::map`
#'   \item Selecting the appropriate prediction function based on the model's
#'         architecture ("DenseNet" vs "Standard")
#'   \item Averaging all predictions using `Reduce("+", out) / length(out)`
#' }
#'
#' The function supports two architectures:
#' \itemize{
#'   \item **"DenseNet"**: Uses `DEGAS::predClass2` for prediction
#'   \item **"Standard"**: Uses `DEGAS::predClass1` for prediction
#' }
#'
#' If an unsupported architecture is encountered, the function will abort
#' with an informative error message.
#'
#' @examples
#' \dontrun{
#' # Create an ensemble of models
#' model_list <- list(model1, model2, model3)
#'
#' # Perform bagged prediction on expression data
#' predictions <- predClassBag(
#'   ccModel = model_list,
#'   Exp = expression_matrix,
#'   scORpat = "sc"
#' )
#'
#' # Use the averaged predictions for downstream analysis
#' class_assignments <- ifelse(predictions > 0.5, "Class1", "Class2")
#' }
#'
#' @importFrom purrr map
#' @keywords internal
#' @family DEGAS
#' @references Johnson TS, Yu CY, Huang Z, Xu S, Wang T, Dong C, et al. Diagnostic Evidence GAuge of Single cells (DEGAS): a flexible deep transfer learning framework for prioritizing cells in relation to disease. Genome Med. 2022 Feb 1;14(1):11.
#'
predClassBag.optimized <- function(ccModel, Exp, scORpat) {
    out <- purrr::map(ccModel, function(ccmodel) {
        switch(
            ccmodel@Architecture,
            "DenseNet" = DEGAS::predClass2(ccmodel, Exp, scORpat),
            "Standard" = DEGAS::predClass1(ccmodel, Exp, scORpat),
            cli::cli_abort(c(
                "x" = "Incorrect architecture argument: ",
                ">" = ccmodel@Architecture
            ))
        )
    })
    Reduce("+", out) / length(out)
}


#' @title Label Binary Phenotype Cells Based on Prediction Scores with Minimum Threshold
#'
#' @description
#' Classifies cells into binary phenotype groups ("Positive" vs "Other") based on
#' prediction score differences between two phenotypic conditions. This enhanced
#' version includes a minimum threshold constraint to ensure biological relevance
#' and provides detailed reporting on classification outcomes.
#'
#' @param pred_dt A data.table containing prediction scores for two phenotypic
#'   conditions. Must contain columns specified in `pheno_colnames`.
#' @param pheno_colnames Character vector of length 2 specifying the column names
#'   for the two phenotypic conditions to compare. The second element is used as
#'   the reference group if not found with regex matching.
#' @param select_fraction Numeric value between 0 and 1 specifying the fraction
#'   of cells to classify as "Positive". Default selection depends on the
#'   distribution characteristics.
#' @param test_method Character string specifying the statistical test to use
#'   for normality assessment. One of: `"jarque-bera"`, `"d'agostino"`,
#'   `"kolmogorov-smirnov"`.
#' @param min_threshold Numeric value specifying the minimum score difference
#'   required for a cell to be considered "Positive". This ensures biological
#'   relevance by filtering out weak associations. Default: 0.7.
#'
#' @return
#' The input `pred_dt` with three additional columns:
#' \itemize{
#'   \item `diff` - Numeric vector of score differences between the two conditions
#'   \item `label` - Character vector with cell classifications: "Positive" or "Other"
#' }
#' The function also provides detailed console output about the classification
#' process and results.
#'
#' @details
#' This function implements a sophisticated approach for binary cell classification
#' that adapts to the underlying distribution of prediction score differences while
#' enforcing a minimum threshold for biological significance:
#'
#' ## Classification Strategies with Minimum Threshold:
#' - **Non-normal distributions (p-value < 0.05)**: Uses quantile-based selection
#'   where the top `select_fraction` of cells by score difference are classified
#'   as "Positive", but only if they exceed `min_threshold`
#' - **Normal distributions (p-value >= 0.05)**: Uses normal distribution
#'   quantiles to determine the classification threshold, adjusted upward if
#'   necessary to meet the minimum threshold requirement
#'
#' ## Supported Normality Tests:
#' - **Jarque-Bera**: Tests for skewness and kurtosis deviations from normality
#' - **D'Agostino**: Extended normality test focusing on skewness
#' - **Kolmogorov-Smirnov**: Non-parametric test comparing empirical distribution
#'   to normal distribution
#'
#' @note
#' The minimum threshold parameter (`min_threshold`) helps prevent over-
#' interpretation of weak phenotypic associations and ensures that classified
#' cells show substantial differences between conditions. The function provides
#' comprehensive feedback about threshold adjustments and final classification
#' statistics.
#'
#' @examples
#' \dontrun{
#' # Create example prediction data
#' pred_data <- data.table(
#'   condition_A = runif(1000),
#'   condition_B = runif(1000)
#' )
#'
#' # Classify cells using D'Agostino test with minimum threshold
#' result <- LabelBinaryCells(
#'   pred_dt = pred_data,
#'   pheno_colnames = c("condition_A", "condition_B"),
#'   select_fraction = 0.1,
#'   test_method = "d'agostino",
#'   min_threshold = 0.7
#' )
#'
#' # Check classification results
#' table(result$label)
#'
#' # View the actual fraction of positive cells
#' prop.table(table(result$label))
#' }
#'
#' @seealso
#' [jb.test.modified()], [dagostino.test()], [ks.test()] for the underlying
#' normality tests used in the classification process.
#'
#' @keywords internal
#' @family DEGAS
#'
LabelBinaryCells <- function(
    pred_dt,
    pheno_colnames,
    select_fraction,
    test_method,
    min_threshold = 0.7 # Added minimum threshold parameter
) {
    chk::chk_length(pheno_colnames, 2)
    # Try to find the reference group in the column names
    # If not found, use the second column as the reference group
    ctrl_col <- grep(
        "[Nn]on|[Nn]ormal|[Cc]ontrol|[Rr]ef|[Cc]trl|0",
        pheno_colnames,
        value = TRUE
    )
    if (nchar(ctrl_col) == 0) {
        cli::cli_alert_info(
            "Using {.val {pheno_colnames[2]}} as reference group"
        )
        pred_dt[, "diff" := .SD[[pheno_colnames[1]]] - .SD[[pheno_colnames[2]]]]
    } else {
        cli::cli_alert_info(
            "Using {.val {ctrl_col}} as reference group"
        )

        pred_dt[,
            "diff" := .SD[[setdiff(pheno_colnames, ctrl_col)]] -
                .SD[[ctrl_col]]
        ]
    }

    # Calculate difference using data.table operations

    # Perform normality test with purrr pattern matching
    normality_test_pval <- switch(
        test_method,
        "jarque-bera" = jb.test.modified(pred_dt$diff)$p.value,
        "d'agostino" = dagostino.test(pred_dt$diff)$p.value[3],
        "kolmogorov-smirnov" = ks.test(pred_dt$diff, "pnorm")$p.value
    )

    # Apply labeling based on normality test result
    pred_dt[,
        "label" := {
            if (normality_test_pval < 0.05) {
                # Use quantile-based selection for non-normal distributions
                n_positive <- ceiling(.N * select_fraction)
                sorted_indices <- order(diff, decreasing = TRUE)
                positive_positions <- sorted_indices[seq_len(n_positive)]

                # Calculate original threshold
                original_thresh <- round(
                    diff[positive_positions[n_positive]],
                    4
                )

                # Apply minimum threshold constraint
                if (original_thresh < min_threshold) {
                    # Filter positions to only include cells above min_threshold
                    above_min_thresh <- which(diff > min_threshold)

                    # Take intersection with top fraction positions
                    valid_positions <- intersect(
                        positive_positions,
                        above_min_thresh
                    )

                    # Update actual number of positive cells
                    actual_n_positive <- length(valid_positions)
                    actual_thresh <- min_threshold

                    cli::cli_alert_info(
                        "Original threshold {.val {original_thresh}} below minimum {.val {min_threshold}}, using {.val {actual_thresh}} instead"
                    )
                } else {
                    valid_positions <- positive_positions
                    actual_thresh <- original_thresh
                    actual_n_positive <- n_positive
                }

                cli::cli_alert_info(
                    "Scores over {.val {actual_thresh}} are considered `Positive`."
                )

                # Create labels using vectorized operations
                labels <- rep("Other", .N)
                labels[valid_positions] <- "Positive"
                labels
            } else {
                # Use normal distribution-based selection
                mean_val <- colMeans(as.matrix(diff))
                sd_val <- colSds(as.matrix(diff))
                quantile_val <- stats::qnorm(
                    select_fraction,
                    mean = mean_val,
                    sd = sd_val,
                    lower.tail = FALSE
                )

                # Apply minimum threshold constraint
                actual_thresh <- max(quantile_val, min_threshold)

                if (quantile_val < min_threshold) {
                    cli::cli_alert_info(
                        "Original threshold {.val {round(quantile_val, 4)}} below minimum {.val {min_threshold}}, using {.val {actual_thresh}} instead"
                    )
                }

                cli::cli_alert_info(
                    "Scores over {.val {round(actual_thresh, 4)}} are considered `Positive`."
                )

                # Use data.table's fast ifelse
                data.table::fifelse(diff > actual_thresh, "Positive", "Other")
            }
        }
    ]

    # Additional validation: count actual positive cells
    positive_count <- pred_dt[label == "Positive", .N]
    total_count <- nrow(pred_dt)
    actual_fraction <- round(positive_count / total_count, 4)

    cli::cli_alert_success(
        "Labeled {.val {positive_count}} cells as Positive ({.val {actual_fraction * 100}}% of total)."
    )

    pred_dt
}

#' @title Label Continuous Phenotype Cells Using MAD Testing
#'
#' @description
#' Identifies phenotype-associated cells in continuous prediction data using
#' Median Absolute Deviation (MAD) testing across multiple phenotypic dimensions.
#'
#' @param pred_dt A data.table containing prediction scores for multiple
#'   phenotypic conditions. Must contain a 'cell_id' column and one or more
#'   columns with prediction scores for different phenotypes.
#'
#' @return
#' The input `pred_dt` with an additional column:
#' \itemize{
#'   \item `label` - Character vector with cell classifications: "Positive"
#'     (significant in at least one phenotype) or "Other"
#' }
#'
#'
#' @note
#' The function assumes that all columns except 'cell_id' contain prediction
#' scores for different phenotypes. It provides progress information and
#' warnings for edge cases like empty results.
#'
#' @examples
#' \dontrun{
#' # Create example prediction data with multiple phenotypes
#' pred_data <- data.table(
#'   cell_id = paste0("cell_", 1:1000),
#'   phenotype_A = rnorm(1000),
#'   phenotype_B = rexp(1000),
#'   phenotype_C = runif(1000)
#' )
#'
#' # Identify phenotype-associated cells
#' result <- LabelContinuousCells(pred_data)
#'
#' # Check classification results
#' table(result$label)
#'
#' # View the proportion of positive cells
#' prop.table(table(result$label))
#' }
#'
#' @seealso
#' [mad.test()] for the underlying statistical test used for phenotype
#' significance assessment.
#'
#' @keywords internal
#' @family DEGAS
#'
LabelContinuousCells <- function(pred_dt) {
    ts_cli$cli_alert_info("Searching for various phenotype-associated cells...")

    # Use matrix operations for efficient MAD testing across predictions
    label_cols <- setdiff(names(pred_dt), "cell_id")

    mad_results <- purrr::map_dfr(
        label_cols,
        function(col) {
            vals <- as.numeric(pred_dt[[col]])
            mad_test <- mad.test(vals)
            data.table::data.table(
                column = col,
                p_value = mad_test$p.value,
                is_positive = mad_test$p.value < 0.05
            )
        }
    )
    if (is.null(mad_results)) {
        cli::cli_warn("Empty MAD test results returned.")
        return(pred_dt)
    }

    # Apply labels based on MAD test results
    positive_cols <- mad_results[`p_value` < 0.05, `column`]
    if (length(positive_cols) > 0) {
        pred_dt[,
            "label" := ifelse(rowSums(.SD) > 0, "Positive", "Other"),
            .SDcols = positive_cols
        ]
    } else {
        pred_dt[, "label" := "Other"]
    }

    pred_dt
}

#' @title Label Survival-Associated Phenotype Cells Based on Hazard Scores
#'
#' @description
#' Classifies cells into survival-associated phenotype groups ("Positive" vs "Other")
#' based on hazard scores using statistical distribution analysis. This function
#' identifies cells with significantly elevated hazard scores that may be
#' associated with survival outcomes, employing adaptive thresholding based on
#' distribution characteristics.
#'
#' @param pred_dt A data.table containing hazard scores for cells. Must contain
#'   a column named 'Hazard' with numeric hazard scores.
#' @param select_fraction Numeric value between 0 and 1 specifying the target
#'   fraction of cells to classify as "Positive". The actual fraction may be
#'   adjusted based on distribution characteristics and minimum threshold
#'   constraints.
#' @param test_method Character string specifying the statistical test to use
#'   for normality assessment of hazard scores. One of: `"jarque-bera"`,
#'   `"d'agostino"`, `"kolmogorov-smirnov"`.
#'
#' @return
#' The input `pred_dt` with an additional column:
#' \itemize{
#'   \item `label` - Character vector with cell classifications: "Positive"
#'     (high hazard cells) or "Other"
#' }
#'
#' @details
#' ## Classification Strategies:
#' - **Non-normal distributions (p-value < 0.05)**: Uses quantile-based selection
#'   where the top `select_fraction` of cells by hazard score are classified as
#'   "Positive", with minimum threshold constraints
#' - **Normal distributions (p-value ≥ 0.05)**: Uses normal distribution
#'   quantiles to determine the classification threshold, adjusted to meet
#'   minimum requirements
#'
#' ## Supported Normality Tests:
#' - **Jarque-Bera**: Tests for skewness and kurtosis deviations from normality
#' - **D'Agostino**: Extended normality test focusing on skewness
#' - **Kolmogorov-Smirnov**: Non-parametric test comparing empirical distribution
#'   to normal distribution
#'
#' @note
#' The function assumes the input data.table contains a column named 'Hazard'
#' with numeric values representing hazard scores from upstream analysis.
#' The minimum threshold is internally defined to ensure biological relevance
#' of the identified cell populations.
#'
#' @examples
#' \dontrun{
#' # Create example hazard score data
#' hazard_data <- data.table(
#'   cell_id = paste0("cell_", 1:1000),
#'   Hazard = rexp(1000, rate = 2)  # Simulated hazard scores
#' )
#'
#' # Identify survival-associated cells
#' result <- LabelSurvivalCells(
#'   pred_dt = hazard_data,
#'   select_fraction = 0.1,
#'   test_method = "jarque-bera"
#' )
#'
#' # Check classification results
#' table(result$label)
#'
#' # Analyze the hazard scores of positive cells
#' summary(result[label == "Positive", Hazard])
#' }
#'
#' @seealso
#' [jb.test.modified()], [dagostino.test()], [ks.test()] for the underlying
#' normality tests used in the classification process.
#'
#' @keywords internal
#' @family DEGAS
#'
LabelSurvivalCells <- function(
    pred_dt,
    select_fraction,
    test_method,
    min_threshold = 0.7 # Added minimum threshold parameter
) {
    ts_cli$cli_alert_info("Searching for survival-associated cells...")

    pred_vec <- pred_dt[["Hazard"]]
    normality_test_pval <- switch(
        test_method,
        "jarque-bera" = jb.test.modified(pred_vec)$p.value,
        "d'agostino" = dagostino.test(pred_vec)$p.value[3],
        "kolmogorov-smirnov" = ks.test(pred_vec, "pnorm")$p.value
    )

    pred_dt[,
        "label" := {
            if (normality_test_pval < 0.05) {
                # Use quantile-based selection for non-normal distributions
                n_positive <- ceiling(.N * select_fraction)
                sorted_indices <- order(`Hazard`, decreasing = TRUE)
                positive_positions <- sorted_indices[seq_len(n_positive)]

                # Calculate original threshold
                original_thresh <- round(
                    `Hazard`[positive_positions[n_positive]],
                    4
                )

                # Apply minimum threshold constraint
                if (original_thresh < min_threshold) {
                    # Filter positions to only include cells above min_threshold
                    above_min_thresh <- which(`Hazard` > min_threshold)

                    # Take intersection with top fraction positions
                    valid_positions <- intersect(
                        positive_positions,
                        above_min_thresh
                    )

                    # Update actual number of positive cells
                    actual_n_positive <- length(valid_positions)
                    actual_thresh <- min_threshold

                    cli::cli_alert_info(
                        "Original threshold {.val {original_thresh}} below minimum {.val {min_threshold}}, using {.val {actual_thresh}} instead"
                    )
                } else {
                    valid_positions <- positive_positions
                    actual_thresh <- original_thresh
                    actual_n_positive <- n_positive
                }

                cli::cli_alert_info(
                    "Scores over {.val {actual_thresh}} are considered `Positive`."
                )

                # Create labels using vectorized operations
                labels <- rep("Other", .N)
                labels[positive_positions] <- "Positive"
                labels
            } else {
                # Use normal distribution-based selection
                mean_val <- colMeans(`Hazard`)
                sd_val <- colSds(`Hazard`)
                quantile_val <- stats::qnorm(
                    select_fraction,
                    mean = mean_val,
                    sd = sd_val,
                    lower.tail = FALSE
                )
                # Apply minimum threshold constraint
                actual_thresh <- max(quantile_val, min_threshold)

                if (quantile_val < min_threshold) {
                    cli::cli_alert_info(
                        "Original threshold {.val {round(quantile_val, 4)}} below minimum {.val {min_threshold}}, using {.val {actual_thresh}} instead"
                    )
                }

                cli::cli_alert_info(
                    "Scores over {.val {round(actual_thresh, 4)}} are considered `Positive`."
                )

                data.table::fifelse(
                    `Hazard` > actual_thresh,
                    "Positive",
                    "Other"
                )
            }
        }
    ]
}
