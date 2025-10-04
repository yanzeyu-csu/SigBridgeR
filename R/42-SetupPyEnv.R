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
        TRUE ~
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
        "numpy" = "any"
    ),
    recreate = FALSE,
    use_conda_forge = TRUE,
    verbose = TRUE,
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
    chk::chk_named(packages)
    #   Default method is `reticulate`
    method = tolower(method)
    if (length(method) > 1) {
        method <- "reticulate"
    }

    if (verbose) {
        cli::cli_h1("Setting up Conda Python Environment")
        cli::cli_alert_info("Environment name: {.val {env_name}}")
        cli::cli_alert_info("Python version: {.val {python_version}}")
    }

    envs <- ListPyEnv(env_type = "conda")
    env_exists <- env_name %in% envs$name

    if (env_exists && recreate) {
        if (verbose) {
            cli::cli_alert_info(c(
                "[{TimeStamp()}] ",
                "Force recreating conda environment: {.val {env_name}}"
            ))
        }
        reticulate::conda_remove(env_name)
        env_exists <- FALSE
    }

    # Create new conda environment
    if (!env_exists) {
        if (verbose) {
            cli::cli_alert_info(c(
                "[{TimeStamp()}] ",
                "Creating new conda environment: {.val {env_name}}"
            ))
        }

        env_created <- switch(
            method,
            "reticulate" = {
                tryCatch(
                    {
                        reticulate::conda_create(
                            envname = env_name,
                            python_version = python_version,
                            channels = if (use_conda_forge) {
                                "conda-forge"
                            } else {
                                NULL
                            },
                            conda = "auto"
                        )
                        TRUE
                    },
                    error = function(e) {
                        cli::cli_alert_danger(
                            "Reticulate creation failed: {e$message}"
                        )
                        FALSE
                    }
                )
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

                result <- tryCatch(
                    {
                        system2(
                            "conda",
                            args = args,
                            stdout = verbose,
                            stderr = verbose
                        )
                        TRUE
                    },
                    error = function(e) {
                        cli::cli_alert_danger(
                            "System creation failed: {e$message}"
                        )
                        FALSE
                    }
                )
                result
            },
            "environment" = {
                chk::chk_file(env_file)

                tryCatch(
                    {
                        system2(
                            "conda",
                            args = c(
                                "env",
                                "create",
                                "-f",
                                env_file,
                                "-y",
                                if (verbose) "-v"
                            ),
                            stdout = verbose,
                            stderr = verbose
                        )

                        # reticulate::conda_create(
                        #     envname = env_name,
                        #     environment = env_file,
                        #     channel = if (use_conda_forge) {
                        #         "conda-forge"
                        #     } else {
                        #         NULL
                        #     },
                        #     conda = "auto"
                        # )
                        TRUE
                    },
                    error = function(e) {
                        cli::cli_alert_danger(
                            "Environment creation from file failed: {e$message}"
                        )
                        FALSE
                    }
                )
            }
        )

        if (!env_created) {
            cli::cli_abort(c(
                "x" = "Failed to create conda environment: {.val {env_name}}"
            ))
        }
    } else if (verbose) {
        cli::cli_alert_info(c(
            "[{TimeStamp()}] ",
            "Using existing conda environment: {.val {env_name}}"
        ))
    }

    envs <- ListPyEnv(env_type = "conda")

    reticulate::use_condaenv(
        envs[envs$name == env_name, 'python'],
        required = TRUE
    )
    # Install packages
    if (length(packages) > 0) {
        if (verbose) {
            cli::cli_alert_info(c(
                "[{TimeStamp()}] ",
                "Installing Python packages in conda environment"
            ))
        }

        switch(
            method,
            "reticulate" = {
                packages_to_install_reticulate <- purrr::imap_chr(
                    packages,
                    ~ if (tolower(.x) == "any") .y else paste0(.y, "==", .x)
                ) %>%
                    unique()

                reticulate::py_install(
                    packages = packages_to_install_reticulate,
                    envname = env_name,
                    method = "auto",
                    pip = TRUE,
                    pip_ignore_installed = TRUE
                )
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

                system2(
                    'conda',
                    args = args,
                    stdout = verbose,
                    stderr = verbose
                )
            },
            "environment" = {
                chk::chk_file(env_file)

                args <- c(
                    "env",
                    "update",
                    if (use_conda_forge) c("-c", "conda-forge"),
                    "-n",
                    env_name,
                    "-f",
                    env_file,
                    "-y",
                    if (verbose) "-v"
                )

                system2(
                    "conda",
                    args = args,
                    stdout = verbose,
                    stderr = verbose
                )
            }
        )
    }

    if (verbose) {
        cli::cli_alert_info(c(
            "[{TimeStamp()}] ",
            "Verifying environment setup..."
        ))
    }

    # Use reticulate to verify the environment
    verification_result <- tryCatch(
        {
            # Test Python availability
            py_available <- reticulate::py_available(initialize = TRUE)
            py_version <- reticulate::py_version()

            if (py_available) {
                if (verbose) {
                    cli::cli_alert_success(c(
                        "[{TimeStamp()}] ",
                        "Python {.val {py_version}} successfully initialized"
                    ))
                }
                TRUE
            } else {
                FALSE
            }
        },
        error = function(e) {
            cli::cli_alert_warning(
                "Environment verification failed: {e$message}"
            )
            FALSE
        }
    )

    if (verbose) {
        if (verification_result) {
            cli::cli_alert_success(c(
                "[{TimeStamp()}] ",
                crayon::green(
                    "Conda environment {env_name} configured successfully!"
                )
            ))
        } else {
            cli::cli_warn(c(
                "[{TimeStamp()}] ",
                "Conda environment created but verification failed"
            ))
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
SetupPyEnv.venv <- function(
    env_type = "venv",
    env_name = "r-reticulate-degas",
    python_version = "3.9.15",
    packages = c("tensorflow" = "2.4.1", "numpy" = "any"),
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
        cli::cli_alert_info("Python version: {.val {python_version}}")
    }

    # Check if environment exists
    env_dir <- Sys.getenv("WORKON_HOME", "~/.virtualenvs")
    env_full_path <- file.path(path.expand(env_dir), env_name)
    env_exists <- dir.exists(env_full_path)

    # Handle existing environment based on recreate flag
    if (env_exists && recreate) {
        if (verbose) {
            cli::cli_alert_info(c(
                "[{TimeStamp()}] ",
                "Force recreating venv environment: {.val {env_name}}"
            ))
        }
        reticulate::virtualenv_remove(envname = env_name, confirm = FALSE)
        unlink(env_full_path, recursive = TRUE)
        env_exists <- FALSE
    }

    # Create new environment if it doesn't exist
    if (!env_exists) {
        if (verbose) {
            cli::cli_alert_info(c(
                "[{TimeStamp()}] ",
                "Creating new venv environment: {.val {env_name}}"
            ))
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
        cli::cli_alert_info(c(
            "[{TimeStamp()}] ",
            "Using existing venv environment: {.val {env_name}}"
        ))
    }

    # Install required packages
    if (length(packages) > 0) {
        if (verbose) {
            cli::cli_alert_info(c(
                "[{TimeStamp()}] ",
                "Installing Python packages in venv environment"
            ))
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
            python_version = python_version,
            pip = TRUE,
            pip_ignore_installed = TRUE
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
        cli::cli_alert_success(c(
            "[{TimeStamp()}] ",
            "Venv environment {.val {env_name}} configured successfully!"
        ))
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
#' @param ... Additional arguments passed to specific methods. Currently supports
#'   `venv_locations` for custom virtual environment search paths.
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
    env_type = c("all", "conda", "venv"),
    ...
) {
    env_type <- match.arg(env_type)
    UseMethod("ListPyEnv")
}

#' @rdname ListPyEnv
#' @description
#' Default method that lists all Python environments by combining results from
#' Conda and virtual environment discovery methods.
#'
#' @param venv_locations Character vector specifying custom directories to search
#'   for virtual environments. Default locations include standard virtualenv
#'   directories and common project locations.
#'
#' @export
ListPyEnv.default <- function(
    env_type = "all",
    venv_locations = c("~/.virtualenvs", "~/.venvs", "./venv", "./.venv"),
    ...
) {
    rbind(
        ListPyEnv.conda(),
        ListPyEnv.venv(
            venv_locations = venv_locations
        )
    )
}

#' @rdname ListPyEnv
#' @description
#' Discovers Conda environments using multiple detection strategies for maximum
#' reliability. First attempts to use reticulate's built-in Conda interface,
#' then falls back to system Conda commands if reticulate is unavailable or
#' fails. Returns empty data frame if Conda is not available or no environments
#' are found.
#'
#' @export
ListPyEnv.conda <- function(env_type = "conda", ...) {
    tryCatch(
        {
            # Method1: reticulate
            conda_envs <- reticulate::conda_list()

            if (!is.null(conda_envs)) {
                conda_envs$type <- "conda"
                return(conda_envs)
            }

            # Method2: system2
            conda_output <- system2(
                "conda",
                c("info", "--envs"),
                stdout = TRUE,
                stderr = FALSE
            )

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

            if (.Platform$OS.type == "windows") {
                python_paths <- file.path(env_paths, "python.exe")
            } else {
                python_paths <- file.path(env_paths, "bin", "python")
            }

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
                "Failed to find conda environments, return empty result."
            )

            data.frame(
                name = character(),
                python = character(),
                type = character()
            )
        },
        error = function(e) {
            cli::cli_warn(c(
                "Error while discovering Conda environments, return empty result:",
                e$message
            ))
            return(data.frame(
                name = character(),
                python = character(),
                type = character()
            ))
        }
    )
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
    env_type = "venv",
    venv_locations = c("~/.virtualenvs", "~/.venvs", "./venv", "./.venv"),
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
    } else {
        cli::cli_warn(
            "No venv found in {.val {venv_locations}}, return empty result"
        )
        return(data.frame(
            name = character(),
            python = character(),
            type = character()
        ))
    }
}
