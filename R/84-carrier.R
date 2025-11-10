#' @title Function from package carrier
#' @keywords internal
#' @inheritParams carrier::crate
crate <- function(
    .fn,
    ...,
    .parent_env = baseenv(),
    .error_arg = ".fn",
    .error_call = environment()
) {
    if (!rlang::is_installed("carrier")) {
        cli::cli_abort(c(
            "x" = "Package {.pkg carrier} is required. Please install it to use this function."
        ))
    }
    getExportedValue("carrier", "crate")(
        .fn = .fn,
        ...,
        .parent_env = .parent_env,
        .error_arg = .error_arg,
        .error_call = .error_call
    )
}
