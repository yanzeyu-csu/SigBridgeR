#' @title Perform scPP screening analysis
#'
#' @description
#' This function performs scPP screening on single-cell data using matched bulk data and phenotype information.
#' It supports binary, continuous, and survival phenotype types.
#'
#' @usage
#' DoscPP(
#'   matched_bulk,
#'   sc_data,
#'   phenotype,
#'   label_type = "scPP",
#'   phenotype_class = c("Binary", "Continuous", "Survival"),
#'   ref_group = 1,
#'   Log2FC_cutoff = 0.585,
#'   estimate_cutoff = 0.2,
#'   probs = 0.2
#' )
#'
#' @param matched_bulk Bulk expression data (genes × samples) where:
#'        - Column names must match `phenotype` row names
#' @param sc_data Seurat object containing preprocessed single-cell data:
#'        - Normalized counts in `RNA` assay
#' @param phenotype Data frame or tibble or named vector with:
#'        - Rownames matching `matched_bulk` columns
#'        - For survival: must contain time and status columns
#' @param label_type Character specifying phenotype label type (e.g., "SBS1"), stored in `scRNA_data@misc`
#' @param phenotype_class Analysis type (case-sensitive):
#'        - `"Binary"`: Case-control studies (e.g., tumor/normal)
#'        - `"Continuous"`: Quantitative traits (e.g., drug response)
#'        - `"Survival"`: Time-to-event data (requires time/status columns)
#' @param ref_group Reference group for binary comparisons (default: 1)
#' @param Log2FC_cutoff Minimum log2 fold-change for binary markers (default: 0.585)
#' @param estimate_cutoff Effect size threshold for continuous traits (default: 0.2)
#' @param probs Quantile cutoff for cell classification (default: 0.2)
#'
#' @return A list containing:
#' \describe{
#'   \item{scRNA_data}{Seurat object with added metadata:
#'     \describe{
#'       \item{ScPP}{"Positive"/"Negative"/"Neutral" classification}
#'       \item{label_type}{Outcome label used}
#'     }
#'   }
#' }
#'
#' @section Algorithm Steps:
#' 1. Data Validation: Checks sample alignment between bulk and phenotype data
#' 2. Marker Selection: Identifies phenotype-associated genes from bulk data
#' 3. Single-cell Screening: Projects bulk markers onto single-cell data
#' 4. Cell Classification: Categorizes cells based on phenotype association
#'
#' @section Reference:
#' Wang X Lab (2023). "scPP: Single-cell Phenotype Prediction".
#' GitHub: \url{https://github.com/WangX-Lab/ScPP}
#'
#' @examples
#' \dontrun{
#' # Binary phenotype analysis
#' res <- DoscPP(
#'   matched_bulk = bulk_data,
#'   sc_data = seurat_obj,
#'   phenotype = ms_data,
#'   label_type = "SBS1",
#'   phenotype_class = "Binary"
#' )
#'
#' # Survival analysis
#' surv_res <- DoscPP(
#'   sc_data = seurat_obj,
#'   matched_bulk = bulk_data,
#'   phenotype = surv_df,
#'   label_type = "OS_status",
#'   phenotype_class = "Survival"
#' )
#' }
#'
#' @importFrom ScPP marker_Binary marker_Continuous marker_Survival ScPP
#' @importFrom Seurat AddMetaData
#' @importFrom dplyr %>% rename
#' @importFrom glue glue
#' @importFrom crayon green
#' @importFrom cli cli_abort cli_alert_info cli_alert_success
#' @importFrom tibble rownames_to_column
#' @importFrom data.table as.data.table fifelse
#' @importFrom stats quantile
#'
#' @family screen method
#'
#' @keywords SigBridgeR_internal
#' @export
#'
DoscPP = function(
    matched_bulk,
    sc_data,
    phenotype,
    label_type = "scPP",
    phenotype_class = c("Binary", "Continuous", "Survival"),
    ref_group = 1,
    Log2FC_cutoff = 0.585,
    estimate_cutoff = 0.2,
    probs = 0.2
) {
    # robust
    if (!all(rownames(phenotype) == colnames(matched_bulk))) {
        cli::cli_abort(c(
            "x" = "Please check the rownames of {.var phenotype} and colnames of {.var bulk_dataset}, they should be the same."
        ))
    }

    cli::cli_alert_info(c(
        "[{TimeStamp()}]",
        crayon::green(" Start scPP screening.")
    ))

    cli::cli_alert_info(c(
        "[{TimeStamp()}]",
        " Finding markers..."
    ))

    matched_bulk = as.data.frame(matched_bulk)
    # decide which type of phenotype data is used
    if (is.vector(phenotype)) {
        phenotype = as.data.frame(phenotype) %>%
            tibble::rownames_to_column("Sample") %>%
            dplyr::rename("Feature" = 2) %>%
            dplyr::mutate(Feature = as.numeric(Feature))
    }
    if (tolower(phenotype_class) == "binary") {
        gene_list = ScPP::marker_Binary(
            bulk_data = matched_bulk,
            features = phenotype,
            ref_group = ref_group,
            Log2FC_cutoff = Log2FC_cutoff
        )
    } else if (tolower(phenotype_class) == "continuous") {
        gene_list = ScPP::marker_Continuous(
            bulk_data = matched_bulk,
            features = phenotype$Feature,
            estimate_cutoff = estimate_cutoff
        )
    } else if (tolower(phenotype_class) == "survival") {
        gene_list = ScPP::marker_Survival(
            bulk_data = matched_bulk,
            survival_data = phenotype
        )
    } else {
        cli::cli_abort(c(
            "x" = "Unknown phenotype type, please check the `phenotype_class` and `phenotype`"
        ))
    }

    l = lapply(gene_list, length)
    pos_null = FALSE
    neg_null = FALSE
    if ("gene_pos" %in% names(l)) {
        # Cannot combine the conditions due to the feature of `gene_list`
        if (l[["gene_pos"]] == 0) {
            cli::cli_alert_info(c(
                "[{TimeStamp()}]",
                " No significant positive genes found"
            ))
            pos_null = TRUE
        }
    }
    if ("gene_neg" %in% names(l)) {
        if (l[["gene_neg"]] == 0) {
            cli::cli_alert_info(c(
                "[{TimeStamp()}]",
                " No significant negative genes found"
            ))
            neg_null = TRUE
        }
    }
    if (pos_null & neg_null) {
        cli::cli_alert_info(c(
            "[{TimeStamp()}]",
            crayon::green(" scPP screening done.")
        ))
        return(list(scRNA_data = NULL))
    }

    cli::cli_alert_info(c(
        "[{TimeStamp()}]",
        " Screening..."
    ))

    # *Start screen
    tryCatch(
        scPP_result <- ScPP::ScPP(sc_data, gene_list, probs = probs),
        error = function(e) {
            cli::cli_alert_danger(c("[{TimeStamp()}]", e$message))
            cli::cli_alert_danger(c("[{TimeStamp()}]", conditionMessage(e)))

            cli::cli_alert_info(c(
                "[{TimeStamp()}]",
                crayon::yellow(" scPP screening exit.")
            ))

            return(list(scRNA_data = NULL))
        }
    )

    sc_data@meta.data[, "scPP"] <- data.table::as.data.table(
        scPP_result$metadata
    )[,
        .(
            scPP = data.table::fifelse(
                ScPP == "Phenotype+",
                "Positive",
                data.table::fifelse(
                    ScPP == "Phenotype-",
                    "Negative",
                    "Neutral"
                )
            )
        )
    ]$scPP

    sc_data <- sc_data %>%
        AddMisc(scPP_type = label_type, cover = FALSE)

    cli::cli_alert_success(c(
        "[{TimeStamp()}]",
        crayon::green(" scPP screening done.")
    ))

    return(list(scRNA_data = sc_data))
}
