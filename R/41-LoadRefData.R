#' @title Download and load the reference data
#'
#' @description
#' This function checks if the data already exists in a local cache.
#' If not, it downloads the file from the remote repository.
#'
#' @param data_type The type of data to download. Must be one of "survival", "binary", or "continuous".
#' @param path Optional path to save the downloaded file.
#' @param cache Logical. If TRUE (default), saves the data for future sessions.
#' @return The requested datasets
#' @export
#'
#' @examples
#' \dontrun{
#' ref_data <- LoadRefData()
#' }
#'
LoadRefData <- function(
    data_type = c("survival", "binary", "continuous"),
    path = tempdir(),
    cache = TRUE
) {
    chk::chk_subset(data_type, c("survival", "binary", "continuous"))
    chk::chk_length(data_type, 1)

    local_file <- file.path(path, glue::glue("{data_type}_ref_data.rds"))

    if (!file.exists(local_file)) {
        cli::cli_alert_info("Downloading reference data...")
        tryCatch(
            {
                data_url <- switch(
                    data_type,
                    survival = "https://raw.githubusercontent.com/WangLabCSU/SigBridgeR/refs/heads/main/vignettes/example_data/survival_example_data.rds",
                    binary = "https://raw.githubusercontent.com/WangLabCSU/SigBridgeR/refs/heads/main/vignettes/example_data/binary_example_data.rds",
                    continuous = "https://raw.githubusercontent.com/WangLabCSU/SigBridgeR/refs/heads/main/vignettes/example_data/continuous_example_data.rds"
                )

                utils::download.file(
                    url = data_url,
                    destfile = local_file,
                    mode = "wb"
                )
            },
            error = function(e) {
                if (file.exists(local_file)) {
                    unlink(local_file)
                }
                cli::cli_abort(c(
                    "x" = crayon::red("Download failed. Error: "),
                    "i" = e$message
                ))
            }
        )
    } else {
        cli::cli_alert_info("Found cached data.")
    }

    data <- readRDS(local_file)

    cli::cli_alert_success(crayon::green("Data loaded successfully."))

    if (!cache) {
        on.exit(unlink(local_file))
    }

    return(data)
}
