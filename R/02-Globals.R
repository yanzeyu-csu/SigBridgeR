# ? Global variables

#' @importFrom data.table `:=` `%chin%`
#' @importFrom dplyr %>%
#' @importFrom SigBridgeRUtils %||%
NULL

# Silence R CMD check NOTES
# These symbols are used inside dplyr/tidyeval or data.table calls
utils::globalVariables(c(
    # environment variables
    "ts_cli",
    # pkg options
    "verbose",
    "seed",
    "parallel",
    "workers"
))
