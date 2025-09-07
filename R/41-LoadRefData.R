#' @title Download and load the reference data
#'
#' @description
#' This function checks if the data already exists in a local cache.
#' If not, it downloads the file from the remote repository.
#'
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
LoadRefData <- function(path = tempdir(), cache = TRUE) {
    data_url <- "https://zenodo.org/record/XXXXXX/files/ref_data.rds"

    local_file <- file.path(path, "ref_data.rds")

    if (!file.exists(local_file)) {
        cli::cli_alert_info("Downloading reference data...")
        tryCatch(
            {
                utils::download.file(
                    url = data_url,
                    destfile = local_file,
                    mode = "wb"
                )
            },
            error = function(e) {
                cli::cli_abort(c(
                    "x" = crayon::red("Download failed. Error: "),
                    e$message
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
