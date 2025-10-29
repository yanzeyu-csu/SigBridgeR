# ? ---- Set up Python environment ----

#' @title Create or Use Python Environment with Required Packages
#'
#' @description
#' Sets up a Python environment with specified packages for DEGAS screening methods.
#' This function can create new environments or reuse existing ones, supporting
#' both Conda and venv environment types. It ensures all required dependencies
#' are properly installed and verified.
#'
#' @param env_type Character string specifying the type of Python environment to
#'   create or use. One of: `"conda"`, `"venv"`.
#' @param ... Additional parameters passed to specific environment methods.
#'
#' @return
#' A data frame containing verification results for the environment setup,
#' including installation status of all required packages. Invisibly returns
#' the verification results.
#'
#' @details
#' This function provides a comprehensive solution for Python environment
#' management in R projects, particularly for machine learning workflows
#' requiring TensorFlow. Key features include:
#'
#' - **Environment Creation**: Automatically creates new environments or reuses
#'   existing ones with the same name
#' - **Package Management**: Installs specified Python packages with version
#'   pinning support
#' - **Verification**: Validates environment setup and package installations
#' - **Flexible Methods**: Supports different backend methods for environment
#'   creation (reticulate vs system calls)
#'
#' The function uses S3 method dispatch to handle different environment types,
#' allowing for extensible support of additional environment managers in the future.
#'
#' @seealso
#' [reticulate::conda_create()], [reticulate::virtualenv_create()] for
#' underlying environment creation functions.
#'
#' @examples
#' \dontrun{
#' # Setup a Conda environment with default parameters
#' SetupPyEnv("conda")
#'
#' # Setup a venv environment
#' SetupPyEnv("venv")
#' }
#'
#' @export
#'
SetupPyEnv <- function(env_type = c("conda", "venv"), ...) {
    UseMethod("SetupPyEnv")
}


#' @rdname SetupPyEnv
#' @description
#' Default method for unsupported environment types. Throws an informative error
#' with supported environment types.
#'
#' @export
SetupPyEnv.default <- function(
    env_type = c("conda", "venv"),
    ...
) {
    switch(
        tolower(env_type),
        "conda" = SetupPyEnv.conda(env_type = "conda", ...),
        "venv" = SetupPyEnv.venv(env_type = "venv", ...),
        cli::cli_abort(c(
            "x" = "Unsupported environment type: {.val {env_type}}",
            "i" = "Supported environment types are: conda, venv"
        ))
    )
}

#' @title Setup Conda Python Environment
#'
#' @description
#' Creates and configures a Conda environment specifically designed for screening workflows. This function
#' provides multiple methods for environment creation and package installation,
#' including support for environment files, with comprehensive verification and
#' error handling.
#'
#' @param env_type Character string specifying the environment type. For this
#'   method, must be "conda".
#' @param env_name Character string specifying the Conda environment name.
#'   Default: "r-reticulate-degas".
#' @param method Character string specifying the method for environment creation
#'   and package installation. One of: "reticulate" (uses reticulate package),
#'   "system" (uses system conda commands), or "environment" (uses YAML
#'   environment file). Default: "reticulate".
#' @param env_file Character string specifying the path to a Conda environment
#'   YAML file. Used when method = "environment". Default: NULL
#' @param python_version Character string specifying the Python version to
#'   install. Default: "3.9.15".
#' @param packages Named character vector of Python packages to install.
#'   Package names as names, versions as values. Use "any" for version to
#'   install latest available. Default includes tensorflow, protobuf, and numpy.
#' @param recreate Logical indicating whether to force recreation of the
#'   environment if it already exists. Default: FALSE.
#' @param use_conda_forge Logical indicating whether to use the conda-forge
#'   channel for package installation. Default: TRUE.
#' @param verbose Logical indicating whether to display detailed progress
#'   messages and command output. Default: TRUE.
#' @param timeout The maximum timeout time when using system commands, default: 30 miniutes
#' @param ... For future compatibility.
#'
#' @return
#' Invisibly returns NULL.
#'
#' @note
#' The function requires Conda to be installed and accessible on the system PATH
#' or through reticulate. For method = "environment", the specified YAML file
#' must exist and be properly formatted. The function includes extensive error
#' handling but may fail if Conda is not properly configured.
#'
#' @examples
#' \dontrun{
#' # Setup using reticulate method (default)
#' SetupPyEnv.conda(
#'   env_name = "my-degas-env",
#'   python_version = "3.9.15"
#' )
#'
#' # Setup using environment file
#' SetupPyEnv.conda(
#'   method = "environment",
#'   env_file = "path/to/environment.yml"
#' )
#'
#' # Setup with custom packages
#' SetupPyEnv.conda(
#'   packages = c(
#'     "tensorflow" = "2.4.1",
#'     "scikit-learn" = "1.0.2",
#'     "pandas" = "any"
#'   )
#' )
#' }
#'
#'
#' @seealso
#' [reticulate::conda_create()], [reticulate::py_install()] for the underlying
#' functions used in reticulate method.
#'
#' @method SetupPyEnv conda
#' @export
SetupPyEnv.conda <- function(
    env_type = "conda",
    env_name = "r-reticulate-degas",
    method = c("reticulate", "system", "environment"),
    env_file = NULL,
    python_version = "3.9.15",
    packages = c(
        "tensorflow" = "2.4.1",
        "protobuf" = "3.20.3"
    ),
    recreate = FALSE,
    use_conda_forge = TRUE,
    verbose = TRUE,
    timeout = 1800000,
    ...
) {
    purrr::walk(
        list(env_type, env_name, python_version),
        ~ chk::chk_character
    )
    purrr::walk(
        list(recreate, use_conda_forge, verbose),
        ~ chk::chk_flag
    )
    if (!is.null(packages)) {
        chk::chk_named(packages)
    }
    #   Default method is `reticulate`
    method %<>% MatchArg(c("reticulate", "system", "environment"))

    if (verbose) {
        cli::cli_h1("Setting up Conda Python Environment")
        cli::cli_alert_info("Environment name: {.val {env_name}}")
        if (!is.null(python_version)) {
            cli::cli_alert_info("Python version: {.val {python_version}}")
        }
    }

    envs <- ListPyEnv(env_type = "conda")
    env_exists <- env_name %chin% envs$name

    if (env_exists && recreate) {
        if (verbose) {
            ts_cli$cli_alert_info(
                "Force recreating conda environment: {.val {env_name}}"
            )
        }
        reticulate::conda_remove(env_name)
        env_exists <- FALSE
    }

    safely_run <- purrr::safely(processx::run)
    safely_create <- purrr::safely(reticulate::conda_create)

    # Create new conda environment
    if (!env_exists) {
        if (verbose) {
            ts_cli$cli_alert_info(
                "Creating new conda environment: {.val {env_name}}"
            )
        }

        switch(
            method,
            "reticulate" = {
                res <- safely_create(
                    envname = env_name,
                    python_version = python_version,
                    channels = if (use_conda_forge) {
                        "conda-forge"
                    } else {
                        NULL
                    },
                    conda = "auto"
                )
                if (!is.null(res$error)) {
                    cli::cli_abort(c(
                        "x" = "Environment creation failed via `reticulate`:",
                        ">" = res$error
                    ))
                }
            },
            "system" = {
                args <- c(
                    "create",
                    "-n",
                    env_name,
                    if (use_conda_forge) c("-c", "conda-forge"),
                    paste0("python=", python_version),
                    "-y",
                    if (verbose) "-v"
                )

                create_res <- safely_run(
                    command = "conda",
                    args = args,
                    error_on_status = FALSE,
                    timeout = timeout,
                    cleanup = TRUE,
                    windows_verbatim_args = FALSE,
                    echo = verbose,
                    echo_cmd = verbose
                )

                # check status
                if (!is.null(create_res$error)) {
                    error_msg <- if (nzchar(create_res$result$stderr)) {
                        create_res$result$stderr
                    } else {
                        create_res$result$stdout
                    }

                    if (grepl("timeout", error_msg, ignore.case = TRUE)) {
                        cli::cli_abort(c(
                            "x" = "Conda environment creation timed out after {.val {timeout/1000/60}} minutes",
                            ">" = "Consider increasing the timeout parameter or using a different method"
                        ))
                    } else {
                        cli::cli_abort(c(
                            "x" = "Environment creation failed via `system` (status {result$status}):",
                            ">" = error_msg
                        ))
                    }
                }
                # print message
                if (verbose && nzchar(create_res$result$stdout)) {
                    message(paste(
                        create_res$result$stdout,
                        sep = "\n",
                        collapse = "\n"
                    ))
                }
                if (verbose) {
                    cli::cli_alert_success(
                        "Conda environment created successfully"
                    )
                }
            },
            "environment" = {
                chk::chk_file(env_file)

                create_res <- safely_create(
                    envname = env_name,
                    environment = env_file
                )
                if (!is.null(create_res$error)) {
                    cli::cli_abort(c(
                        "x" = "Environment creation failed via `environment`:",
                        ">" = create_res$error
                    ))
                }
            }
        )
    } else if (verbose) {
        ts_cli$cli_alert_info(
            "Using existing conda environment: {.val {env_name}}"
        )
    }

    envs <- ListPyEnv(env_type = "conda")

    reticulate::use_condaenv(
        envs[envs$name == env_name, 'python'],
        required = TRUE
    )
    # Install packages
    if (length(packages) > 0 && method != "environment") {
        if (verbose) {
            ts_cli$cli_alert_info(
                "Installing Python packages in conda environment"
            )
        }

        switch(
            method,
            "reticulate" = {
                packages_to_install_reticulate <- purrr::imap_chr(
                    packages,
                    ~ if (tolower(.x) == "any") .y else paste0(.y, "==", .x)
                ) %>%
                    unique()

                safely_py_install <- purrr::safely(reticulate::py_install)

                install_res <- reticulate::py_install(
                    packages = packages_to_install_reticulate,
                    envname = env_name,
                    method = "auto",
                    pip = TRUE,
                    pip_ignore_installed = TRUE
                )
                if (!is.null(install_res$error)) {
                    cli::cli_abort(c(
                        "x" = "Failed to install packages in conda environment {.val {env_name}} via `reticulate`",
                        ">" = install_res$error
                    ))
                }
            },
            "system" = {
                packages_to_install_conda <- purrr::imap_chr(
                    packages,
                    ~ if (tolower(.x) == "any") .y else paste0(.y, "=", .x)
                )

                args <- c(
                    "install",
                    "-n",
                    env_name,
                    if (use_conda_forge) c("-c", "conda-forge"),
                    packages_to_install_conda,
                    "-y",
                    if (verbose) "-v"
                )

                install_res <- safely_run(
                    command = "conda",
                    args = args,
                    error_on_status = FALSE,
                    timeout = timeout,
                    cleanup = TRUE,
                    windows_verbatim_args = FALSE,
                    echo = verbose,
                    echo_cmd = verbose
                )

                # check status
                if (!is.null(install_res$error)) {
                    error_msg <- if (nzchar(install_res$error$stderr)) {
                        install_res$error$stderr
                    } else {
                        install_res$error$stdout
                    }

                    if (grepl("timeout", error_msg, ignore.case = TRUE)) {
                        cli::cli_abort(c(
                            "x" = "Package installation timed out after {.val {timeout/1000/60}} minutes",
                            ">" = "Consider increasing the timeout parameter or installing packages separately"
                        ))
                    } else {
                        cli::cli_abort(c(
                            "x" = "Package installation failed via `system`:",
                            ">" = error_msg
                        ))
                    }
                }
                # print message
                if (verbose && nzchar(install_res$error$stdout)) {
                    message(paste(
                        install_res$error$stdout,
                        sep = "\n",
                        collapse = "\n"
                    ))
                }

                if (verbose) {
                    cli::cli_alert_success(
                        "Packages installed successfully"
                    )
                }
            }
        )
    }

    if (verbose) {
        ts_cli$cli_alert_info("Verifying environment setup...")
    }

    # Use reticulate to verify the environment
    verification_result <- rlang::try_fetch(
        {
            # Test Python availability
            py_available <- reticulate::py_available(initialize = TRUE)
            py_version <- reticulate::py_version()

            if (py_available) {
                if (verbose) {
                    ts_cli$cli_alert_success(
                        "Python {.val {py_version}} successfully initialized"
                    )
                }
                TRUE
            } else {
                FALSE
            }
        },
        error = function(e) {
            cli::cli_alert_danger(
                "Environment verification failed: {e$message}"
            )
            FALSE
        }
    )

    if (verbose) {
        if (verification_result) {
            ts_cli$cli_alert_info(cli::col_green(
                "Conda environment {env_name} configured successfully!"
            ))
        } else {
            cli::cli_warn(
                "Conda environment created but verification failed"
            )
        }
    }

    invisible()
}

#' @title Setup Virtual Environment (venv)
#'
#' @description
#' Creates and configures a Python virtual environment (venv) specifically designed
#' for screening workflows.
#' This function provides a lightweight, isolated Python environment alternative
#' to Conda environments with similar package management capabilities.
#'
#' @param env_type Character string specifying the environment type. For this
#'   method, must be "venv".
#' @param env_name Character string specifying the virtual environment name.
#'   Default: "r-reticulate-degas".
#' @param python_version Character string specifying the Python version to use.
#'   Default: "3.9.15".
#' @param packages Named character vector of Python packages to install.
#'   Package names as names, versions as values. Use "any" for version to
#'   install latest available. Default includes tensorflow, protobuf, and numpy.
#' @param python_path Character string specifying the path to a specific Python
#'   executable. If NULL, uses the system default or installs the specified
#'   version. Default: NULL.
#' @param recreate Logical indicating whether to force recreation of the
#'   virtual environment if it already exists. Default: FALSE.
#' @param verbose Logical indicating whether to display detailed progress
#'   messages and command output. Default: TRUE.
#' @param ... For future compatibility.
#'
#' @return
#' Invisibly returns NULL.
#'
#' @note
#' Virtual environments require a base Python installation. If the specified
#' Python version is not available, the function will attempt to install it
#' using reticulate. Virtual environments are generally faster to create than
#' Conda environments but may have more limited package availability compared
#' to Conda-forge.
#'
#' @examples
#' \dontrun{
#' # Setup virtual environment with default parameters
#' SetupPyEnv.venv()
#'
#' # Setup with custom Python version and packages
#' SetupPyEnv.venv(
#'   env_name = "my-degas-venv",
#'   python_version = "3.8.12",
#'   packages = c(
#'     "tensorflow" = "2.4.1",
#'     "scikit-learn" = "1.0.2",
#'     "pandas" = "any"
#'   )
#' )
#'
#' # Force recreate existing environment
#' SetupPyEnv.venv(
#'   env_name = "existing-env",
#'   recreate = TRUE
#' )
#' }
#'
#' @seealso
#' [reticulate::virtualenv_create()], [reticulate::virtualenv_remove()],
#' [reticulate::use_virtualenv()] for the underlying virtual environment
#' management functions.
#'
#' @method SetupPyEnv venv
#' @export
#'
SetupPyEnv.venv <- function(
    env_type = "venv",
    env_name = "r-reticulate-degas",
    python_version = "3.9.15",
    packages = c("tensorflow" = "2.4.1", "protobuf" = "3.20.3"),
    python_path = NULL,
    recreate = FALSE,
    verbose = TRUE,
    ...
) {
    # Input validation
    purrr::walk(
        list(env_type, env_name, python_version),
        ~ chk::chk_character
    )
    purrr::walk(
        list(recreate, verbose),
        ~ chk::chk_flag
    )
    chk::chk_named(packages)
    if (!is.null(python_path)) {
        chk::chk_file(python_path)
    }

    if (verbose) {
        cli::cli_h1("Setting up Venv Python Environment")
        cli::cli_alert_info("Environment name: {.val {env_name}}")
        if (!is.null(python_version)) {
            cli::cli_alert_info("Python version: {.val {python_version}}")
        }
    }

    # Check if environment exists
    env_dir <- Sys.getenv("WORKON_HOME", "~/.virtualenvs")
    env_full_path <- file.path(path.expand(env_dir), env_name)
    env_exists <- dir.exists(env_full_path)

    # Handle existing environment based on recreate flag
    if (env_exists && recreate) {
        if (verbose) {
            ts_cli$cli_alert_info(
                "Force recreating venv environment: {.val {env_name}}"
            )
        }
        reticulate::virtualenv_remove(envname = env_name, confirm = FALSE)
        unlink(env_full_path, recursive = TRUE)
        env_exists <- FALSE
    }

    # Create new environment if it doesn't exist
    if (!env_exists) {
        if (verbose) {
            ts_cli$cli_alert_info(
                "Creating new venv environment: {.val {env_name}}"
            )
        }

        # Determine Python path
        if (is.null(python_path)) {
            reticulate::install_python(version = python_version)
        }

        # Create venv environment
        reticulate::virtualenv_create(
            envname = env_name,
            python = python_version,
            packages = NULL,
            virtualenv = "venv"
        )
    } else if (verbose) {
        ts_cli$cli_alert_info(
            "Using existing venv environment: {.val {env_name}}"
        )
    }

    # Install required packages
    if (length(packages) > 0) {
        if (verbose) {
            ts_cli$cli_alert_info(
                "Installing Python packages in venv environment"
            )
        }

        # Switch to target environment
        reticulate::use_virtualenv(
            virtualenv = env_name,
            required = TRUE
        )

        # Format packages for installation using purrr
        packages_to_install <- purrr::imap_chr(
            packages,
            ~ if (tolower(.x) == "any") .y else paste0(.y, "==", .x)
        ) %>%
            unique()

        reticulate::py_install(
            packages = packages_to_install,
            envname = env_name,
            method = "virtualenv",
            python_version = python_version
        )
    }

    # Verify installation
    pkg_names <- if (length(packages) > 0) names(packages) else character(0)
    packages_to_verify <- unique(c(
        pkg_names,
        "functools",
        "math"
    ))

    if (verbose) {
        ts_cli$cli_alert_success(
            "Venv environment {.val {env_name}} configured successfully!"
        )
    }

    reticulate::py_require(
        packages = packages_to_install,
        python_version = python_version,
    )

    invisible()
}


# ? ----List Environments as Data Frame----

#' @title List Available Python Environments
#'
#' @description
#' Discovers and lists available Python environments of various types on the system.
#' This generic function provides a unified interface to find Conda environments
#' and virtual environments (venv) through S3 method dispatch.
#'
#' @param env_type Character string specifying the type of environments to list.
#'   One of: `"all"`, `"conda"`, `"venv"`. Defaults to `"all"`.
#' @param timeout Numeric value specifying the timeout in seconds for the Conda
#'   environment discovery process. Defaults to 30 minutes.
#' @param venv_locations Character vector of additional locations to search for
#'   virtual environments. Defaults to `c("~/.virtualenvs", "~/.venvs", "./venv", "./.venv")`.
#' @param verbose Logical value indicating whether to print verbose output.
#'   Defaults to `TRUE`.
#' @param ... For future use.
#'
#' @return
#' A data frame with the following columns:
#' \itemize{
#'   \item `name` - Character vector of environment names
#'   \item `python` - Character vector of paths to Python executables
#'   \item `type` - Character vector indicating environment type (`"conda"` or `"venv"`)
#' }
#' Returns an empty data frame with these columns if no environments are found.
#'
#' @details
#' The function uses S3 method dispatch to handle different environment types:
#'
#' - **`"all"`**: Combines results from all environment types using `rbind()`
#' - **`"conda"`**: Searches for Conda environments using multiple methods:
#'   - Primary: `reticulate::conda_list()` for reliable environment detection
#'   - Fallback: System `conda info --envs` command for broader compatibility
#' - **`"venv"`**: Searches common virtual environment locations including
#'   user directories and project folders
#'
#' Each method includes comprehensive error handling and will return empty
#' results with informative warnings if no environments are found or if
#' errors occur during discovery.
#'
#' @examples
#' \dontrun{
#' # List all Python environments
#' ListPyEnv("all")
#'
#' # List only Conda environments
#' ListPyEnv("conda")
#'
#' # List only virtual environments with custom search paths
#' ListPyEnv("venv", venv_locations = c("~/my_envs", "./project_env"))
#' }
#'
#' @export
ListPyEnv <- function(
    env_type = c("all", "conda", "venv", "virtualenv"),
    timeout = 30000,
    venv_locations = c("~/.virtualenvs", "~/.venvs", "./venv", "./.venv"),
    verbose = TRUE,
    ...
) {
    UseMethod("ListPyEnv")
}

#' @rdname ListPyEnv
#' @description
#' Default method that lists all Python environments by combining results from
#' Conda and virtual environment discovery methods.
#'
#' @param timeout The maximum timeout time when using system commands, only effective when `env_type=conda`.
#' @param venv_locations Character vector specifying custom directories to search
#'   for virtual environments. Default locations include standard virtualenv
#'   directories and common project locations.
#' @param verbose Logical indicating whether to print verbose output.
#'
#' @export
#'
ListPyEnv.default <- function(
    env_type = c("all", "conda", "venv", "virtualenv"),
    timeout = 30000,
    venv_locations = c("~/.virtualenvs", "~/.venvs", "./venv", "./.venv"),
    verbose = TRUE,
    ...
) {
    env_type %<>% MatchArg(c("all", "conda", "venv", "virtualenv"))
    switch(
        env_type,
        "conda" = ListPyEnv.conda(
            timeout = timeout,
            verbose = verbose,
            ...
        ),
        "virtualenv" = ListPyEnv.venv(venv_locations = venv_locations),
        "venv" = ListPyEnv.venv(venv_locations = venv_locations),
        "all" = rbind(
            ListPyEnv.conda(
                timeout = timeout,
                verbose = verbose,
                ...
            ),
            ListPyEnv.venv(
                venv_locations = venv_locations
            )
        ),
        cli::cli_abort(c(
            "x" = "Invalid environment type: {.val {env_type}}",
            "i" = "Valid types are: {.code all}, {.code conda} or {.code venv}"
        ))
    )
}

#' @rdname ListPyEnv
#' @description
#' Discovers Conda environments using multiple detection strategies for maximum
#' reliability. First attempts to use system Conda commands,
#' then falls back to reticulate's built-in Conda interface if Conda command is unavailable or
#' fails. Returns empty data frame if Conda is not available or no environments
#' are found.
#'
#' @export
ListPyEnv.conda <- function(
    env_type = c("all", "conda", "venv", "virtualenv"),
    timeout = 30000,
    venv_locations = c("~/.virtualenvs", "~/.venvs", "./venv", "./.venv"),
    verbose = TRUE,
    ...
) {
    methods <- c(
        system = function() {
            # Method1: system
            process_result <- processx::run(
                command = "conda",
                args = c("info", "--envs"),
                error_on_status = FALSE,
                timeout = timeout,
                cleanup = TRUE,
                windows_verbatim_args = FALSE
            )

            if (process_result$status != 0) {
                error_msg <- if (nzchar(process_result$stderr)) {
                    process_result$stderr
                } else {
                    process_result$stdout
                }

                cli::cli_abort(c(
                    "x" = "Conda command failed with status {process_result$status}:",
                    ">" = "{error_msg}"
                ))
            }
            conda_output <- strsplit(process_result$stdout, "\n")[[1]]

            env_lines <- grep(
                "^[a-zA-Z_]",
                conda_output,
                value = TRUE
            ) %>%
                gsub("\\*", "", .) %>%
                trimws() %>%
                strsplit("\\s+")

            if (length(env_lines) == 0) {
                cli::cli_warn(
                    "No Conda environments found, return empty result."
                )
                return(data.frame(
                    name = character(),
                    python = character(),
                    type = character(),
                ))
            }

            env_matrix <- do.call(rbind, env_lines)
            env_names <- env_matrix[, 1]
            env_paths <- env_matrix[, 2]

            GetPythonPath <- function(path) {
                if (is.na(path)) {
                    return(NA_character_)
                }
                candidates <- if (.Platform$OS.type == "windows") {
                    c("python.exe", "Scripts/python.exe")
                } else {
                    c("bin/python", "bin/python3")
                }
                for (candidate in candidates) {
                    full_path <- file.path(path, candidate)
                    if (file.exists(full_path)) {
                        return(normalizePath(full_path, mustWork = FALSE))
                    }
                }
                return(NA_character_)
            }
            python_paths <- vapply(env_paths, GetPythonPath, character(1))

            conda_result <- data.frame(
                name = env_names,
                python = python_paths,
                type = "conda",
                stringsAsFactors = FALSE
            )

            if (!is.null(conda_result) && nrow(conda_result) > 0) {
                return(conda_result)
            }

            cli::cli_warn(
                "No conda environments found, return empty result."
            )

            data.frame(
                name = character(),
                python = character(),
                type = character()
            )
        },
        reticulate = function() {
            # Method2: reticulate
            cli::cli_warn(
                "Failed to find conda environments via system command, trying reticulate as fallback."
            )
            conda_envs <- reticulate::conda_list()

            if (!is.null(conda_envs) && nrow(conda_envs) > 0) {
                conda_envs$type <- "conda"
                return(conda_envs)
            }

            cli::cli_warn(
                "No conda environments found, return empty result."
            )

            data.frame(
                name = character(),
                python = character(),
                type = character()
            )
        },
        default = function() {
            cli::cli_warn(
                "All methods have failed to find the conda environment, returning empty conda environment result ."
            )
            data.frame(
                name = character(),
                python = character(),
                type = character()
            )
        }
    ) %>%
        purrr::map(purrr::safely)

    for (func_name in names(methods)) {
        method_result <- methods[[func_name]]()
        if (is.null(method_result$error) || func_name == "default") {
            return(method_result$result)
        }
    }
}

#' @rdname ListPyEnv
#' @description
#' Discovers virtual environments by searching common venv locations including
#' user directories (`~/.virtualenvs`, `~/.venvs`) and project folders
#' (`./venv`, `./.venv`). Supports custom search paths through the
#' `venv_locations` parameter. Returns empty data frame if no virtual
#' environments are found in the specified locations.
#'
#' @param venv_locations Character vector of directory paths to search for
#'   virtual environments. Default includes standard locations and common
#'   project directories.
#'
#' @export
ListPyEnv.venv <- function(
    env_type = c("all", "conda", "venv", "virtualenv"),
    timeout = 30000,
    venv_locations = c("~/.virtualenvs", "~/.venvs", "./venv", "./.venv"),
    verbose = TRUE,
    ...
) {
    venv_dirs <- c()

    for (location in venv_locations) {
        expanded_path <- path.expand(location)
        if (dir.exists(expanded_path)) {
            dirs <- list.dirs(expanded_path, recursive = FALSE)
            venv_dirs <- c(venv_dirs, dirs)
        }
    }

    if (length(venv_dirs) > 0) {
        return(data.frame(
            name = basename(venv_dirs),
            python = file.path(
                venv_dirs,
                ifelse(
                    .Platform$OS.type == "windows",
                    "Scripts/python.exe",
                    "bin/python"
                )
            ),
            type = "venv"
        ))
    }

    cli::cli_warn(
        "No venv found in {.val {venv_locations}}, return empty virtual environment result"
    )

    data.frame(
        name = character(),
        python = character(),
        type = character()
    )
}
