# ? Global variables

#' @importFrom data.table `:=` `%chin%`
#' @importFrom dplyr %>%
#' @import SigBridgeRUtils
#' @importFrom DEGAS SetupPyEnv ListPyEnv
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
    "workers",
    # dplyr/tidyeval and data.table variables
    ".",
    "..duplicate_cols",
    "PC",
    "PC1",
    "PC2",
    "Variance",
    "condition",
    "batch",
    "Feature",
    "n",
    "nFeature_RNA",
    "Total",
    "Fraction",
    "sets",
    "count",
    # DEGAS-related variables
    "DEGAS.model_type",
    "DEGAS.architecture",
    "DEGAS.ff_depth",
    "DEGAS.bag_depth",
    "DEGAS.seed"
))
