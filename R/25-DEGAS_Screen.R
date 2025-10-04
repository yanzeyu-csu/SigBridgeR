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
#' @param select_fraction Fraction of cells to select as positive when `phenotype_class` is "binary" or "survival" (default: 0.05)
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
#' \code{\link{Vec2sparse}} for the structure of phenotype
#' \code{\link{jb.test.modified}} for modified Jarque-Bera test
#' \code{\link{mad.test}} for outlier detection using Median Absolute Deviation
#'
#' @importFrom purrr walk map_chr map_dfr
#' @importFrom dplyr %>%
#' @importFrom data.table as.data.table fifelse setnames `:=`
#' @importFrom matrixStats colMeans2 colSds
#' @importFrom stats ks.test qnorm median

#' @keywords internal
#' @family screen_method
#'
DoDEGAS <- function(
    select_fraction = 0.05,
    matched_bulk,
    sc_data,
    phenotype = NULL,
    sc_data.pheno_colname = NULL,
    label_type = "DEGAS",
    phenotype_class = c("binary", "continuous", "survival"),
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
    result <- local({
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
        phenotype_class <- match.arg(phenotype_class)
        chk::chk_subset(
            normality_test_method,
            c("kolmogorov-smirnov", "d'agostino", "jarque-bera")
        )
        # DEGAS path must contain "/"
        if (!grepl("/$", tmp_dir)) {
            tmp_dir <- paste0(tmp_dir, "/")
        }

        cli::cli_alert_info(c(
            "[{TimeStamp()}] ",
            crayon::green("Starting DEGAS Screen")
        ))

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
                "numpy" = "any"
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
            DEGAS.ff_depth = 3,
            DEGAS.bag_depth = 5,
            path.data = '',
            path.result = '',
            DEGAS.pyloc = NULL,
            DEGAS.toolsPath = paste0(.libPaths()[1], "/DEGAS/tools/"),
            DEGAS.train_steps = 2000,
            DEGAS.scbatch_sz = 200,
            DEGAS.patbatch_sz = 50,
            DEGAS.hidden_feats = 50,
            DEGAS.do_prc = 0.5,
            DEGAS.lambda1 = 3.0,
            DEGAS.lambda2 = 3.0,
            DEGAS.lambda3 = 3.0,
            DEGAS.seed = 2
        )
        env_params <- utils::modifyList(default_env_params, env_params)
        degas_params <- utils::modifyList(default_degas_params, degas_params)

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
            cli::cli_abort(c(
                "x" = "Specify at least on phenotype"
            ))
        }
        degas_params$DEGAS.model_type = paste0(
            model_type.first,
            model_type.last
        )

        if (length(degas_params$DEGAS.architecture) != 1) {
            degas_params$DEGAS.architecture <- "DenseNet"
        }

        pheno_df <- if (!is.null(phenotype)) {
            switch(
                phenotype_class,
                "binary" = Vec2sparse(phenotype),
                "continuous" = Vec2sparse(phenotype),
                "survival" = as.matrix(phenotype),
            )
        } else {
            NULL
        }

        meta_data = sc_data@meta.data
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

        sc_mat <- if (utils::packageVersion("Seurat") >= "5.0.0") {
            sc_data@assays$RNA$data
        } else {
            sc_data@assays$RNA@data
        }
        cm_genes <- intersect(rownames(matched_bulk), rownames(sc_mat))
        t_sc_mat <- t(sc_mat[cm_genes, ])
        t_matched_bulk <- t(matched_bulk[cm_genes, ])

        cli::cli_alert_info(c("[{TimeStamp()}] ", "Setting up Environment..."))

        existing_envs <- ListPyEnv()

        if (
            !env_params$env.name %in% existing_envs$name ||
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

        DEGAS::set_seed_term(2)

        list2env(degas_params, envir = environment())

        cli::cli_alert_info(c("[{TimeStamp()}] ", "Training DEGAS model..."))

        ccModel1 <- DEGAS::runCCMTLBag(
            scExp = t_sc_mat,
            scLab = sc_pheno,
            patExp = t_matched_bulk,
            patLab = pheno_df,
            tmpDir = tmp_dir,
            model_type = DEGAS.model_type,
            architecture = DEGAS.architecture,
            FFdepth = DEGAS.ff_depth,
            Bagdepth = DEGAS.bag_depth
        )

        cli::cli_alert_info(c(
            "[{TimeStamp()}] ",
            "Predicting with DEGAS model..."
        ))

        t_sc_preds <- data.table::as.data.table(DEGAS::predClassBag(
            ccModel1,
            t_sc_mat,
            'pat'
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

        if (length(normality_test_method) != 1) {
            normality_test_method = "jarque-bera"
        }

        cli::cli_alert_info(c("[{TimeStamp()}] ", "Labeling screened cells..."))

        t_sc_preds <- switch(
            phenotype_class,
            "binary" = LabelBinaryCells(
                pred_dt = t_sc_preds,
                pheno_colnames = pheno_df_colnames,
                select_fraction = select_fraction,
                test_method = normality_test_method
            ),
            "continuous" = LabelContinuousCells(pred_dt = t_sc_preds),
            "survival" = LabelSurvivalCells(
                pred_dt = t_sc_preds,
                select_fraction = select_fraction,
                test_method = normality_test_method
            )
        )

        sc_data <- Seurat::AddMetaData(
            object = sc_data,
            metadata = t_sc_preds[["label"]],
            col.name = "DEGAS"
        ) %>%
            AddMisc(DEGAS_type = label_type)

        cli::cli_alert_info(c(
            "[{TimeStamp()}] ",
            crayon::green("DEGAS Screen done.")
        ))

        list(
            scRNA_data = sc_data,
            DEGAS_prediction = t_sc_preds
        )
    })

    return(result)
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
#' @import Matrix
#' @keywords internal
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
    sparseMatrix(
        i = rows,
        j = as.integer(cols),
        dims = c(length(x), nlevels(cols)),
        dimnames = list(names(x), col_names)
    )
}

#' @title Modified Jarque-Bera Test for Normailty
#'
#' @description
#' Performs a modified version of the Jarque-Bera test for normality that
#' provides enhanced sensitivity to departures from normal distribution.
#' This implementation includes four different test statistics based on
#' whether population mean and standard deviation are known or estimated.
#'
#' @param x A numeric vector of data values to test for normality.
#' @param mean An optional numeric value specifying the known population mean.
#'   If `NA` (default), the sample mean is used.
#' @param sd An optional numeric value specifying the known population standard
#'   deviation. If `NA` (default), the sample standard deviation is used.
#'
#' @return
#' An object of class `"htest"` containing the following components:
#' \itemize{
#'   \item `statistic` - The test statistic (chi-squared distributed)
#'   \item `parameter` - Degrees of freedom (always 2)
#'   \item `p.value` - The p-value for the test
#'   \item `method` - Description of the test method
#'   \item `data.name` - Name of the data vector
#' }
#'
#' @details
#' The modified Jarque-Bera test evaluates normality by testing whether
#' the sample data has the skewness and kurtosis matching a normal distribution.
#' The test statistic is based on the sample size, skewness, and kurtosis,
#' and follows a chi-squared distribution with 2 degrees of freedom.
#'
#' ## Test Variants:
#' The function implements four different test statistics depending on
#' whether the population mean and standard deviation are specified:
#'
#' - **Both unknown**: Uses sample estimates for mean and variance
#' - **Mean known, SD unknown**: Uses known mean but estimates variance
#' - **Mean unknown, SD known**: Uses sample mean with known variance
#' - **Both known**: Uses known population parameters
#'
#' Each variant uses different coefficients in the test statistic calculation
#' to account for the uncertainty in parameter estimation.
#'
#' ## Mathematical Formulation:
#' The test statistic is calculated as:
#' \deqn{JB = n \times \left( \frac{S^2}{6} + \frac{(K - 3)^2}{24} \right)}
#' where \eqn{S} is skewness, \eqn{K} is kurtosis, and \eqn{n} is sample size,
#' with modified coefficients for different parameter knowledge scenarios.
#'
#' @section Reference:
#' KhrushchevSergey/modified_jarque_bera_test Internet. cited 2025 Sep 28.
#' Available from: \url{https://github.com/KhrushchevSergey/Modified-Jarque-Bera-test/blob/main/modified_jarque_bera_test.R}
#'
#' @examples
#' \dontrun{
#' # Test with sample data (both parameters unknown)
#' x <- rnorm(100)
#' test_result <- jb.test.modified(x)
#' print(test_result)
#'
#' # Test with known population parameters
#' x <- rnorm(100, mean = 5, sd = 2)
#' test_result <- jb.test.modified(x, mean = 5, sd = 2)
#' print(test_result)
#'
#' # Test with only known mean
#' test_result <- jb.test.modified(x, mean = 5)
#' print(test_result)
#' }
#'
#' @seealso
#' [stats::shapiro.test()] for the Shapiro-Wilk test
#'
#' @keywords internal
jb.test.modified <- function(x, mean = NA, sd = NA) {
    if ((NCOL(x) > 1) || is.data.frame(x)) {
        cli::cli_abort(c("x" = "x is not a vector or univariate time series"))
    }
    if (anyNA(x)) {
        cli::cli_abort(c("x" = "x contains {.val NA}"))
    }
    DNAME <- deparse(substitute(x))
    n <- length(x)

    if (is.na(mean) & is.na(sd)) {
        m1 <- sum(x) / n
        m2 <- sum((x - m1)^2) / n
        m3 <- sum((x - m1)^3) / n
        m4 <- sum((x - m1)^4) / n
        b1 <- (m3 / m2^(3 / 2))^2
        b2 <- (m4 / m2^2)
        STATISTIC <- n * (b1 / 6 + (b2 - 3)^2 / 24)
    }

    if (!is.na(mean) & is.na(sd)) {
        m1 <- mean
        m2 <- sum((x - m1)^2) / n
        m3 <- sum((x - m1)^3) / n
        m4 <- sum((x - m1)^4) / n
        b1 <- (m3 / m2^(3 / 2))^2
        b2 <- (m4 / m2^2)
        STATISTIC <- n * (b1 / 15 + (b2 - 3)^2 / 24)
    }

    if (is.na(mean) & !is.na(sd)) {
        m1 <- mean(x)
        m2 <- sd^2
        m3 <- sum((x - m1)^3) / n
        m4 <- sum((x - m1)^4) / n
        b1 <- (m3 / m2^(3 / 2))^2
        b2 <- (m4 / m2^2)
        STATISTIC <- n * (b1 / 6 + (b2 - 3)^2 / 96)
    }

    if (!is.na(mean) & !is.na(sd)) {
        m1 <- mean
        m2 <- sd^2
        m3 <- sum((x - m1)^3) / n
        m4 <- sum((x - m1)^4) / n
        b1 <- (m3 / m2^(3 / 2))^2
        b2 <- (m4 / m2^2)
        STATISTIC <- n * (b1 / 15 + (b2 - 3)^2 / 96)
    }

    PVAL <- 1 - stats::pchisq(STATISTIC, df = 2)
    PARAMETER <- 2
    METHOD <- "Modified Jarque Bera Test"
    names(STATISTIC) <- "X-squared"
    names(PARAMETER) <- "df"
    structure(
        list(
            statistic = STATISTIC,
            parameter = PARAMETER,
            p.value = PVAL,
            method = METHOD,
            data.name = DNAME
        ),
        class = "htest"
    )
}


#' @title D'Agostino Test of Normality
#'
#' @description
#' Performs the D'Agostino test for normality, which is particularly suitable
#' for large sample sizes where other tests like Shapiro-Wilk may be too sensitive.
#'
#' @param x a numeric vector of data values. Missing values are allowed and will be removed.
#'
#' @return A list with class `"htest"` containing the following components:
#' \item{statistic}{a named vector containing three test statistics:
#'   \describe{
#'     \item{Skewness}{the Z statistic for skewness test}
#'     \item{Kurtosis}{the Z statistic for kurtosis test}
#'     \item{Omnibus}{the chi-squared statistic combining both tests}
#'   }
#' }
#' \item{parameter}{the degrees of freedom for the omnibus test}
#' \item{p.value}{a named vector containing three p-values corresponding to the statistics}
#' \item{method}{the character string "D'Agostino Normality Test"}
#' \item{data.name}{a character string giving the name(s) of the data}
#' \item{alternative}{a character string describing the alternative hypothesis}
#' \item{estimates}{a named vector containing sample statistics:
#'   \describe{
#'     \item{n}{sample size}
#'     \item{Skewness}{sample skewness coefficient}
#'     \item{Kurtosis}{sample kurtosis coefficient (excess kurtosis)}
#'     \item{Mean}{sample mean}
#'     \item{SD}{sample standard deviation}
#'   }
#' }
#'
#' @details
#' The D'Agostino test is a powerful omnibus test for normality that combines
#' separate tests for skewness and kurtosis. It is based on the transformation
#' of the sample skewness and kurtosis to approximate normal distributions,
#' which are then combined into a chi-squared statistic.
#'
#' This test is particularly recommended for:
#' \itemize{
#'   \item Large sample sizes (n > 50)
#'   \item Situations where Shapiro-Wilk test is overly sensitive
#'   \item Testing both tail behavior and symmetry simultaneously
#' }
#'
#' The test has good power against a wide range of alternative distributions
#' and is less sensitive to sample size fluctuations than many other normality tests.
#'
#' @examples
#' \dontrun{
#' # Test normally distributed data
#' set.seed(123)
#' normal_data <- rnorm(100)
#' result <- dagostino.test(normal_data)
#' print(result)
#'
#' # Test non-normal data (exponential distribution)
#' non_normal_data <- rexp(100)
#' result2 <- dagostino.test(non_normal_data)
#' print(result2)
#'
#' # Access specific components
#' result$p.value["Omnibus"]  # Overall test p-value
#' result$estimates["Skewness"]  # Sample skewness
#' result$statistic["Kurtosis"]  # Kurtosis test statistic
#'
#' # Large sample performance
#' large_sample <- rnorm(5000)
#' dagostino.test(large_sample)
#' }
#'
#'
#' @keywords htest
#' @keywords internal
#'
dagostino.test <- function(x) {
    if (!is.numeric(x)) {
        stop("'x' must be a numeric vector")
    }
    if (length(x) < 8) {
        stop("'x' must have at least 8 elements")
    }

    x <- x[!is.na(x)]
    n <- as.double(length(x))

    x_mean <- mean(x)
    x_centered <- x - x_mean
    m2 <- sum(x_centered^2) / n
    m3 <- sum(x_centered^3) / n
    m4 <- sum(x_centered^4) / n

    g1 <- m3 / (m2^(3 / 2))

    n_dbl <- as.double(n)
    n2 <- n_dbl * n_dbl
    n3 <- n2 * n_dbl

    Y <- g1 * sqrt((n_dbl + 1) * (n_dbl + 3) / (6 * (n_dbl - 2)))
    beta2 <- 3 *
        (n2 + 27 * n_dbl - 70) *
        (n_dbl + 1) *
        (n_dbl + 3) /
        ((n_dbl - 2) * (n_dbl + 5) * (n_dbl + 7) * (n_dbl + 9))

    if (beta2 <= 1) {
        W_sq <- 1.0
    } else {
        W_sq <- sqrt(2 * beta2 - 2)
    }

    if (W_sq <= 1) {
        Z_g1 <- 0
    } else {
        delta <- 1 / sqrt(log(W_sq))
        alpha <- sqrt(2 / (W_sq - 1))

        ratio <- Y / alpha
        Z_g1 <- delta * log(ratio + sqrt(ratio^2 + 1))
    }

    g2 <- m4 / (m2^2) - 3

    E_g2 <- -6 / (n_dbl + 1)
    Var_g2 <- 24 *
        n_dbl *
        (n_dbl - 2) *
        (n_dbl - 3) /
        ((n_dbl + 1)^2 * (n_dbl + 3) * (n_dbl + 5))

    if (Var_g2 <= 0) {
        standardized_g2 <- 0
    } else {
        standardized_g2 <- (g2 - E_g2) / sqrt(Var_g2)
    }

    if (n_dbl <= 3) {
        beta2_kurt <- 1.0
    } else {
        beta2_kurt <- 6 *
            (n2 - 5 * n_dbl + 2) /
            ((n_dbl + 7) * (n_dbl + 9)) *
            sqrt(
                6 *
                    (n_dbl + 3) *
                    (n_dbl + 5) /
                    (n_dbl * (n_dbl - 2) * (n_dbl - 3))
            )
    }

    if (beta2_kurt <= 0 || is.infinite(beta2_kurt)) {
        Z_g2 <- 0
    } else {
        A <- 6 +
            (8 / beta2_kurt) * (2 / beta2_kurt + sqrt(1 + 4 / (beta2_kurt^2)))

        if (A <= 4) {
            Z_g2 <- 0
        } else {
            term1 <- 1 - 2 / (9 * A)
            term2 <- (1 - 2 / A) / (1 + standardized_g2 * sqrt(2 / (A - 4)))
            if (term2 <= 0) {
                Z_g2 <- 0
            } else {
                Z_g2 <- (term1 - term2^(1 / 3)) / sqrt(2 / (9 * A))
            }
        }
    }

    K_sq <- Z_g1^2 + Z_g2^2

    p_value_skew <- 2 * stats::pnorm(-abs(Z_g1))
    p_value_kurt <- 2 * stats::pnorm(-abs(Z_g2))
    p_value_combined <- stats::pchisq(K_sq, df = 2, lower.tail = FALSE)

    p_value_skew <- max(0, min(1, p_value_skew))
    p_value_kurt <- max(0, min(1, p_value_kurt))
    p_value_combined <- max(0, min(1, p_value_combined))

    DNAME <- deparse(substitute(x))
    METHOD <- "D'Agostino Normality Test"

    result <- structure(
        list(
            statistic = c(
                Skewness = Z_g1,
                Kurtosis = Z_g2,
                Omnibus = K_sq
            ),
            parameter = c(df = 2),
            p.value = c(
                Skewness = p_value_skew,
                Kurtosis = p_value_kurt,
                Omnibus = p_value_combined
            ),
            method = METHOD,
            data.name = DNAME,
            alternative = "data are not normally distributed",
            estimates = c(
                n = n,
                Skewness = g1,
                Kurtosis = g2,
                Mean = x_mean,
                SD = sqrt(m2 * n / (n - 1))
            )
        ),
        class = "htest"
    )
    return(result)
}


#' @title outlier detection using Median Absolute Deviation (MAD)
#'
#' @description
#' This function identifies outliers in a numeric vector using Median Absolute
#' Deviation (MAD) method, with specialized handling for very small samples.
#' Returns results in htest format for consistency with statistical tests.
#'
#' @param x A numeric vector. Missing values (`NA`) are allowed and will be
#'          excluded from calculations.
#' @param na.rm Logical. Should missing values be removed? Default is `TRUE`.
#'
#' @return An object of class "htest" containing:
#'   \item{statistic}{The test statistic (MAD-based z-score of most extreme outlier)}
#'   \item{parameter}{Named vector with: n, median, mad, threshold}
#'   \item{p.value}{Approximate p-value based on outlier probability}
#'   \item{method}{Description of the method used}
#'   \item{data.name}{Name of the data vector}
#'   \item{alternative}{Character string describing the alternative hypothesis}
#'   \item{outlier.indices}{Indices of detected outliers (as additional element)}
#'
#' @details
#' For samples with n > 5, uses standard MAD method with threshold of 3.
#' For very small samples (n <= 5), uses MAD with adjusted thresholds:
#' - n = 3-4: threshold = 2.5
#' - n = 5: threshold = 2.8
#'
#' The test statistic is calculated as the maximum absolute deviation from median
#' divided by adjusted MAD.
#'
#' @examples
#' \dontrun{
#' # Basic usage
#' x <- c(1, 2, 3, 100, 4, 5, 3)
#' mad.test(x)
#'
#' # Small sample
#' x3 <- c(1, 2, 10)
#' mad.test(x3)
#' }
#'
#' @import data.table
#' @importFrom chk chk_length chk_numeric chk_not_any_na
#'
#' @keywords internal
#'
mad.test <- function(x, na.rm = TRUE) {
    # Store data name for htest output
    dna <- deparse(substitute(x))

    # Input validation
    chk::chk_numeric(x)

    # Handle missing values
    if (na.rm) {
        complete_cases <- !is.na(x)
        x_clean <- x[complete_cases]
        original_indices <- which(complete_cases)
    } else {
        chk::chk_not_any_na(x)
        x_clean <- x
        original_indices <- seq_along(x)
    }

    n <- length(x_clean)

    chk::chk_length(n)

    # Create data.table for efficient computation
    dt <- data.table::data.table(
        index = original_indices,
        value = x_clean
    )

    # Calculate median and MAD
    median_val <- dt[, median(`value`)]
    abs_dev <- abs(dt$`value` - median_val)
    mad_val <- median(abs_dev)
    mad_adjusted <- mad_val * 1.4826

    # Determine threshold based on sample size
    if (n <= 4) {
        threshold <- 2.5
        method_used <- "MAD-based outlier detection (small sample n<=4, threshold=2.5)"
    } else if (n == 5) {
        threshold <- 2.8
        method_used <- "MAD-based outlier detection (small sample n=5, threshold=2.8)"
    } else {
        threshold <- 3.0
        method_used <- "MAD-based outlier detection (threshold=3.0)"
    }

    # Calculate MAD-based z-scores
    dt[, `z_score` := abs(`value` - `median_val`) / `mad_adjusted`] # id mad_adjusted == 0, z_score will be NA

    # Identify outliers
    dt[, `is_outlier` := `z_score` > `threshold`]

    # Get outlier information
    outliers_dt <- dt[`is_outlier` == TRUE]
    n_outliers <- nrow(outliers_dt)

    # Calculate test statistic (z-score of most extreme outlier)
    if (n_outliers > 0) {
        statistic_val <- max(dt$z_score)
        # Get indices of outliers in original data (before NA removal)
        outlier_indices <- dt[`is_outlier` == TRUE, `index`]
    } else {
        statistic_val <- max(dt$z_score)
        outlier_indices <- integer(0)
    }

    # Calculate approximate p-value
    # Based on probability of observing such extreme value in normal distribution
    if (mad_adjusted > 0) {
        pval <- 2 * (1 - stats::pnorm(statistic_val))
        # Adjust for multiple testing (Bonferroni correction)
        pval <- min(pval * n, 1.0)
    } else {
        pval <- NA_real_
    }

    # Create htest structure
    result <- list(
        statistic = c("MAD z-score" = statistic_val),
        parameter = c(
            "n" = n,
            "median" = median_val,
            "mad" = mad_adjusted,
            "threshold" = threshold
        ),
        p.value = pval,
        method = method_used,
        data.name = dna,
        alternative = "at least one value is an outlier",
        outlier.indices = outlier_indices
    )

    class(result) <- "htest"
    return(result)
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
        "[Nn]on|[Nn]ormal|[Cc]ontrol|[Rr]ef|[Cc]trl",
        pheno_colnames
    )
    if (is.null(ctrl_col)) {
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
                mean_val <- matrixStats::colMeans2(as.matrix(diff))
                sd_val <- matrixStats::colSds(as.matrix(diff))
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

    return(pred_dt)
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
#'
LabelContinuousCells <- function(pred_dt) {
    cli::cli_alert_info(c(
        "[{TimeStamp()}] ",
        "Searching for various phenotype-associated cells..."
    ))

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

    return(pred_dt)
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
#'
LabelSurvivalCells <- function(
    pred_dt,
    select_fraction,
    test_method,
    min_threshold = 0.7 # Added minimum threshold parameter
) {
    cli::cli_alert_info(c(
        "[{TimeStamp()}] ",
        "Searching for phenotype-associated cells..."
    ))

    pred_vec = pred_dt[["Hazard"]]
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
                mean_val <- matrixStats::colMeans2(`Hazard`)
                sd_val <- matrixStats::colSds(`Hazard`)
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

    return(pred_dt)
}
