#' @title Match Argument with Enhanced Feedback and Efficiency
#'
#' @description
#' Matches an argument against a set of valid choices, with support for partial matching,
#' case-insensitive matching, and multiple selections. Provides rich feedback for user-friendly error and success messages.
#'
#' @param arg The argument to match, typically a character string or vector.
#' @param choices A character vector of valid choices to match against.
#' @param multiple Logical; if `TRUE`, allows multiple matches. Defaults to `FALSE`.
#' @param case_sensitive Logical; if `TRUE`, matching is case-sensitive. Defaults to `FALSE`.
#' @param arg_name Character; name of the argument for error messages. Defaults to the
#'   name of `arg` as derived from `deparse(substitute(arg))`.
#'
#' @return A character vector of matched choices. If `multiple = FALSE`, returns a single
#'   character string. If no match is found, an error is thrown with detailed feedback.
#'
#' @details
#' This function is similar to `match.arg()` but provides enhanced features:
#' - Supports multiple selections when `multiple = TRUE`.
#' - Uses a hash table for efficient matching.
#' - Provides rich, user-friendly feedback.
#' - Supports case-insensitive matching (controlled by `case_sensitive`).
#'
#' @examples
#' # Single selection
#' MatchArg("gau", c("gaussian", "poisson", "binomial"))
#'
#' # Multiple selections
#' MatchArg(c("gau", "poi"), c("gaussian", "poisson", "binomial"), multiple = TRUE)
#'
#' # Case-sensitive matching
#' MatchArg("Gau", c("gaussian", "poisson", "binomial"), case_sensitive = TRUE)
#'
#' # Invalid input (will throw an error)
#' try(MatchArg("xyz", c("gaussian", "poisson", "binomial")))
#'
#' @importFrom cli cli_abort cli_alert_info cli_alert_success cli_alert_danger
#' @importFrom glue glue_collapse
#'
#' @keywords internal
#'
#' @export
#'
MatchArg <- function(
    arg,
    choices,
    multiple = FALSE,
    case_sensitive = FALSE,
    arg_name = deparse(substitute(arg))
) {
    arg_name <- as.character(arg_name)

    if (missing(choices)) {
        call_env <- parent.frame()
        call_fn <- sys.function(-1)
        call_formals <- formals(call_fn)

        if (as.character(substitute(arg)) %in% names(call_formals)) {
            choices <- eval(
                call_formals[[as.character(substitute(arg))]],
                call_env
            )
        } else {
            cli::cli_abort(c(
                "x" = "Cannot determine choices for argument: {.arg {arg_name}}",
                "i" = "Please provide the choices explicitly"
            ))
        }
    }

    if (length(choices) == 0) {
        cli::cli_abort(c(
            "x" = "Choices cannot be empty for argument: {.arg {arg_name}}"
        ))
    }
    if (any(nchar(choices) == 0) && is.character(choices)) {
        cli::cli_abort(c(
            "x" = "Choices cannot contain empty strings for argument: {.arg {arg_name}}"
        ))
    }

    # Handle NULL or empty input
    if (is.null(arg) || length(arg) == 0) {
        if (!multiple) {
            cli::cli_alert_info(
                "No value provided for {.arg {arg_name}}, defaulting to: {.val {choices[1]}}"
            )
            return(choices[1])
        } else {
            cli::cli_abort(c(
                "x" = "Multiple selections require non-empty input for: {.arg {arg_name}}"
            ))
        }
    }

    choices <- as.character(choices)
    arg <- as.character(arg)

    if (!case_sensitive) {
        choices <- tolower(choices)
        arg <- tolower(arg)
    }

    choice_set <- unique(choices)
    choice_map <- setNames(seq_along(choice_set), choice_set)

    # Match logic
    if (multiple) {
        matches <- vapply(
            arg,
            function(x) {
                if (nchar(x) == 0) {
                    return(NA_integer_)
                }
                exact <- choice_map[x]
                if (!is.na(exact)) {
                    return(exact)
                }
                partial <- which(startsWith(names(choice_map), x))
                if (length(partial) == 1) {
                    return(partial)
                }
                NA_integer_
            },
            integer(1)
        )
    } else {
        if (length(arg) > 1 && !is.null(arg)) {
            cli::cli_abort(c(
                "x" = "Only one value allowed for {.arg {arg_name}}.",
                "i" = "Valid choices are: {.val {glue::glue_collapse(choice_set, ', ', last = ' or ')}}"
            ))
        }
        exact <- choice_map[arg[1]]
        if (!is.na(exact)) {
            matches <- exact
        } else {
            partial <- which(startsWith(names(choice_map), arg[1]))
            matches <- if (length(partial) == 1) partial else NA_integer_
        }
    }
    # Handle non-matches
    if (anyNA(matches)) {
        cli::cli_abort(
            c(
                "x" = "Matching failed for argument: {.arg {arg_name}} = {.val {arg}}",
                "i" = "Valid choices are: {.val {glue::glue_collapse(choice_set, ', ', last = ' or ')}}"
            ),
            class = "MatchArgOptionsNotFound"
        )
    }

    return(choices[matches])
}
