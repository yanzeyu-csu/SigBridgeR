#' @title Generate Timestamp String
#'
#' @description
#' Creates a formatted character string representing the current system time.
#' The format is "YYYY/mm/DD HH:MM:SS" (year/month/day hour:minute:second).
#'
#' @return Character string with current time in "YYYY/mm/DD HH:MM:SS" format.
#'         If system time is unavailable, returns a fixed timestamp.
#'
#' @examples
#' \dontrun{
#' # Current time as formatted string
#' TimeStamp()
#' # Returns something like: "2025/06/15 16:04:00"
#' }
#'
#' @keywords internal
#' @family TimeStamp
#'
TimeStamp <- function() {
    safely_time <- purrr::safely(Sys.time)
    time_res <- safely_time()
    time <- if (!is.null(time_res$error)) {
        as.POSIXct("1970/01/01 00:00:00", tz = "UTC")
    } else {
        time_res$result
    }
    format(time, "%Y/%m/%d %H:%M:%S")
}

#' @title A Decorator for Adding Timestamp to CLI Functions
#'
#' @description
#' A higher-order function that wraps CLI functions to automatically prepend
#' timestamps to their output messages. This creates a modified version of
#' any CLI function that includes timestamp information in its output.
#'
#' @param cli_func A CLI function from the \code{cli} package (e.g., \code{cli_alert_info},
#'                 \code{cli_warn}) that will be wrapped with timestamp functionality.
#'
#' @return Returns a modified version of the input function that automatically
#'         adds a timestamp in the format \code{"[{timestamp}]"} to the beginning
#'         of all character messages passed to it.
#'
#' @details
#' This function uses \code{\link{force}} to ensure the CLI function is evaluated
#' at creation time. The timestamp is generated using a \code{TimeStamp()} function
#' which should be available in the execution environment and is inserted using
#' cli's glue-like syntax \code{"{ }"}.
#'
#' @examples
#' \dontrun{
#' # Create a timestamp-enabled version of cli_alert_info
#' timestamp_alert <- AddTimeStamp2cli(cli::cli_alert_info)
#' timestamp_alert("This message will have a timestamp")
#' }
#'
#' @seealso \code{\link{CreateTimeStampCliEnv}} for creating a complete environment
#'          of timestamped CLI functions.
#'
#' @keywords internal
#' @family TimeStamp
#'
AddTimeStamp2cli <- function(cli_func) {
    force(cli_func)

    function(...) {
        messages <- list(...)

        if (length(messages) > 0) {
            if (is.character(messages[[1]])) {
                messages[[1]] <- paste0("[{TimeStamp()}] ", messages[[1]])
            }
        }

        do.call(cli_func, messages)
    }
}

#' @title Create Environment with Timestamped CLI Functions
#'
#' @description
#' Generates an environment containing wrapped versions of common CLI functions
#' that automatically include timestamps in their output. This provides a
#' convenient way to use multiple CLI functions with consistent timestamping.
#'
#' @return Returns an environment containing timestamp-wrapped versions of:
#' \itemize{
#'   \item \code{cli_warn}
#'   \item \code{cli_alert_info}
#'   \item \code{cli_alert_success}
#'   \item \code{cli_alert_warning}
#'   \item \code{cli_alert_danger}
#' }
#' Each function in the environment will automatically prepend timestamps to
#' its output messages.
#'
#' @details
#' This function creates a new environment and populates it with timestamped
#' versions of commonly used CLI functions from the \code{cli} package. Only
#' functions that exist in the loaded \code{cli} package are added to the
#' environment. The function uses \code{\link[purrr]{walk}} for side-effect
#' iteration over the function names.
#'
#' @examples
#' \dontrun{
#' # Create environment with timestamped CLI functions
#' cli_env <- CreateTimeStampCliEnv()
#'
#' # Use timestamped functions
#' cli_env$cli_alert_info("System started")
#' cli_env$cli_alert_success("Operation completed")
#' }
#'
#' @seealso \code{\link{AddTimeStamp2cli}} for the wrapper function used internally.
#' @importFrom purrr walk
#'
#' @keywords internal
#' @family TimeStamp
#'
CreateTimeStampCliEnv <- function(
    cli_functions = c(
        "cli_alert_info",
        "cli_alert_success",
        "cli_alert_warning",
        "cli_alert_danger"
    )
) {
    cli_env <- new.env()

    purrr::walk(cli_functions, function(func_name) {
        if (exists(func_name, envir = asNamespace("cli"))) {
            orig_func <- get(func_name, envir = asNamespace("cli"))

            new_func <- eval(substitute(
                function(..., .envir = parent.frame()) {
                    args <- list(...)
                    if (length(args) > 0 && is.character(args[[1]])) {
                        args[[1]] <- paste0(
                            "[{.dim ",
                            TimeStamp(),
                            "}] ",
                            args[[1]]
                        )
                    }
                    args$.envir <- .envir
                    do.call(ORIG_FUNC, args)
                },
                list(ORIG_FUNC = orig_func)
            ))

            assign(func_name, new_func, envir = cli_env)
        }
    })

    invisible(cli_env)
}
