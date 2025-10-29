# ? Global variables

#' @importFrom data.table `:=` `%chin%`
#' @importFrom magrittr %>% %<>%
#' @importFrom Matrix Matrix
NULL

# Silence R CMD check NOTES
# These symbols are used inside dplyr/tidyeval or data.table calls
utils::globalVariables(c(
    "!!!",
    # environment variables
    "ts_cli",
    # dplyr pipe & tidy-select
    ".",
    "..duplicate_cols",
    # data.table internals
    "sym",
    # column names and suffixes
    "scPP_AUCdown",
    "scPP_AUCup",
    "scPP",
    "base_key",
    "batch",
    "cell_id",
    "cell_label",
    "column",
    "condition",
    "DEGAS.architecture",
    "DEGAS.bag_depth",
    "DEGAS.do_prc",
    "DEGAS.ff_depth",
    "DEGAS.hidden_feats",
    "DEGAS.lambda1",
    "DEGAS.lambda2",
    "DEGAS.lambda3",
    "DEGAS.model_type",
    "DEGAS.patbatch_sz",
    "DEGAS.pyloc",
    "DEGAS.scbatch_sz",
    "DEGAS.seed",
    "DEGAS.train_steps",
    "Feature",
    "Fraction",
    "Hazard",
    "index",
    "is_outlier",
    "label",
    "MatchArg",
    "max_suffix",
    "model_results",
    "n",
    "nFeature_RNA",
    "p_value",
    "percent.mt",
    "PC",
    "PC1",
    "PC2",
    "scAB_select",
    "sets",
    "suffix",
    "Total",
    "value",
    "Variance",
    "z_score"
))
