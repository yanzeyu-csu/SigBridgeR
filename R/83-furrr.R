#' @title Future plan for furrr
#' @keywords internal
#' @inheritParams future::plan
plan <- function(
    strategy = NULL,
    ...,
    substitute = TRUE,
    .skip = FALSE,
    .call = TRUE,
    .cleanup = NA,
    .init = TRUE
) {
    if (!rlang::is_installed("future")) {
        cli::cli_abort(c(
            "x" = "Package {.pkg future} is required. Please install it to use this function."
        ))
    }
    getExportedValue("future", "plan")(
        strategy = strategy,
        ...,
        substitute = substitute,
        .skip = .skip,
        .call = .call,
        .cleanup = .cleanup,
        .init = .init
    )
}

#' @title Available cores for furrr
#' @keywords internal
availableCores <- function() {
    if (!rlang::is_installed("future")) {
        cli::cli_abort(c(
            "x" = "Package {.pkg future} is required. Please install it to use this function."
        ))
    }
    getExportedValue("future", "availableCores")()
}


#' @title Future map for furrr
#' @keywords internal
#' @inheritParams furrr::future_map_dbl
future_map_dbl <- function(
    .x,
    .f,
    ...,
    .options = furrr_options(),
    .env_globals = parent.frame(),
    .progress = FALSE
) {
    if (!rlang::is_installed("furrr")) {
        cli::cli_abort(c(
            "x" = "Package {.pkg furrr} is required. Please install it to use this function."
        ))
    }
    getExportedValue("furrr", "future_map_dbl")(
        .x = .x,
        .f = .f,
        ...,
        .options = .options,
        .env_globals = .env_globals,
        .progress = .progress
    )
}

#' @title Future map for furrr
#' @keywords internal
#' @inheritParams furrr::future_map
future_map <- function(
    .x,
    .f,
    ...,
    .options = furrr_options(),
    .env_globals = parent.frame(),
    .progress = FALSE
) {
    if (!rlang::is_installed("furrr")) {
        cli::cli_abort(c(
            "x" = "Package {.pkg furrr} is required. Please install it to use this function."
        ))
    }
    getExportedValue("furrr", "future_map")(
        .x = .x,
        .f = .f,
        ...,
        .options = .options,
        .env_globals = .env_globals,
        .progress = .progress
    )
}


#' @title Future options for furrr
#' @keywords internal
#' @inheritParams furrr::furrr_options
furrr_options <- function(
    ...,
    stdout = TRUE,
    conditions = "condition",
    globals = TRUE,
    packages = NULL,
    seed = FALSE,
    scheduling = 1,
    chunk_size = NULL,
    prefix = NULL
) {
    if (!rlang::is_installed("furrr")) {
        cli::cli_abort(c(
            "x" = "Package {.pkg furrr} is required. Please install it to use this function."
        ))
    }
    getExportedValue("furrr", "furrr_options")(
        ...,
        stdout = stdout,
        conditions = conditions,
        globals = globals,
        packages = packages,
        seed = seed,
        scheduling = scheduling,
        chunk_size = chunk_size,
        prefix = prefix
    )
}
