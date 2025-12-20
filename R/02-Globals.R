# ? Global variables

#' @importFrom data.table `:=` `%chin%`
#' @importFrom dplyr %>%
NULL

# These symbols are used inside the package
utils::globalVariables(c(
    # environment variables
    "ts_cli",
    # pkg options
    "verbose",
    "seed"
))
