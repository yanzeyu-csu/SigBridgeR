#' @title Configuration Functions for SigBridgeR Package
#' @description
#' These functions provide a centralized configuration system for the SigBridgeR
#' package, allowing users to set and retrieve package-specific options with
#' automatic naming conventions.
#'
#' @section Available Opions:
#' - `verbose`: A logical value indicating whether to print verbose messages, defaults to `TRUE`
#' - `parallel`: A logical value indicating whether to use parallel processing, defaults to `FALSE``
#' - `parallel.type`: A character string specifying the type of parallel processing to use, defaults to `"multisession"`
#' - `workers`: An integer specifying the number of workers to use for parallel processing, defaults to `4L`
#' - `timeout`: An integer specifying the timeout in seconds for parallel processing, defaults to `180L``
#' - `seed`: An integer specifying the random seed for reproducible results, defaults to `123L``
#'
#'
#' @name SigBridgeR_Function_Setting
NULL

#' @rdname SigBridgeR_Function_Setting
#' @description
#' setFunOption sets one or more configuration options for the SigBridgeR package.
#' Options are automatically prefixed with "SigBridgeR." if not already present,
#' ensuring proper namespace isolation.
#'
#' @param ... Named arguments representing option-value pairs. Options can be
#' provided with or without the "SigBridgeR." prefix. If the prefix is missing,
#' it will be automatically added.
#'
#' @return Invisibly returns NULL. The function is called for its side effects
#' of setting package options.
#'
#' @examples
#' # Set options with automatic prefixing
#' setFuncOption(verbose = TRUE)
#'
#' @export
setFuncOption <- function(...) {
    opts <- rlang::list2(...)
    if (length(opts) > 0) {
        opt_names <- names(opts)
        needs_prefix <- !startsWith(opt_names, "SigBridgeR.")
        opt_names[needs_prefix] <- paste0(
            "SigBridgeR.",
            opt_names[needs_prefix]
        )
        names(opts) <- opt_names

        options(opts)
    }

    invisible()
}

#' @rdname SigBridgeR_Function_Setting
#' @description
#' getFuncOption retrieves configuration options for the SigBridgeR package.
#' The function automatically handles the "SigBridgeR." prefix, allowing users
#' to reference options with or without the explicit prefix.
#'
#' @param option Character string specifying the option name to retrieve.
#' The "SigBridgeR." prefix is optional and will be added if missing.
#' @param default The default value to return if the option is not set.
#' Defaults to NULL.
#'
#' @return The value of the specified option if it exists, otherwise the
#' provided default value.
#'
#' @examples
#' # Retrieve options with automatic prefixing
#' getFuncOption("verbose")
#'
#' @export
getFuncOption <- function(option, default = NULL) {
    if (!startsWith(option, "SigBridgeR.")) {
        option <- paste0("SigBridgeR.", option)
    }
    getOption(option, default = default)
}
