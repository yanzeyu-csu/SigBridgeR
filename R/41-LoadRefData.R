#' @title Download & Load Reference Data
#'
#' @description
#' This function checks if the data already exists in a local cache.
#' If not, it downloads the file from the remote repository using multiple sources with fallback.
#'
#' @param data_type The type of data to download. Must be one of "survival", "binary", or "continuous".
#' @param path Optional path to save the downloaded file.
#' @param cache Logical. If TRUE (default), saves the data for future sessions.
#' @param timeout Integer. Connection timeout in seconds (default: 60).
#' @return The requested datasets, stored in a list.
#' @export
#'
#' @examples
#' \dontrun{
#' ref_data <- LoadRefData(data_type = c("survival"))
#' }
#'
LoadRefData <- function(
    data_type = c("survival", "binary", "continuous"),
    path = tempdir(),
    cache = TRUE,
    timeout = 60
) {
    data_type <- match.arg(data_type)
    chk::chk_dir(path)
    chk::chk_flag(cache)
    chk::chk_whole_number(timeout)

    local_file <- file.path(path, glue::glue("{data_type}_ref_data.rds"))

    if (!file.exists(local_file)) {
        cli::cli_alert_info("Downloading reference data...")

        # Define multiple sources with priority order
        data_urls <- list(
            survival = list(
                "CDN (jsDelivr)" = "https://cdn.jsdelivr.net/gh/WangLabCSU/SigBridgeR@main/vignettes/example_data/survival_example_data.rds",
                "GitHub Raw" = "https://raw.githubusercontent.com/WangLabCSU/SigBridgeR/refs/heads/main/vignettes/example_data/survival_example_data.rds",
                "GitHub API" = "https://raw.githubusercontent.com/WangLabCSU/SigBridgeR/main/vignettes/example_data/survival_example_data.rds"
            ),
            binary = list(
                "CDN (jsDelivr)" = "https://cdn.jsdelivr.net/gh/WangLabCSU/SigBridgeR@main/vignettes/example_data/binary_example_data.rds",
                "GitHub Raw" = "https://raw.githubusercontent.com/WangLabCSU/SigBridgeR/refs/heads/main/vignettes/example_data/binary_example_data.rds",
                "GitHub API" = "https://raw.githubusercontent.com/WangLabCSU/SigBridgeR/main/vignettes/example_data/binary_example_data.rds"
            ),
            continuous = list(
                "CDN (jsDelivr)" = "https://cdn.jsdelivr.net/gh/WangLabCSU/SigBridgeR@main/vignettes/example_data/continuous_example_data.rds",
                "GitHub Raw" = "https://raw.githubusercontent.com/WangLabCSU/SigBridgeR/refs/heads/main/vignettes/example_data/continuous_example_data.rds",
                "GitHub API" = "https://raw.githubusercontent.com/WangLabCSU/SigBridgeR/main/vignettes/example_data/continuous_example_data.rds"
            )
        )

        available_urls <- data_urls[[data_type]]
        success <- FALSE

        # Try each source in priority order
        for (i in seq_along(available_urls)) {
            source_name <- names(available_urls)[i]
            data_url <- available_urls[[i]]

            cli::cli_alert_info(glue::glue("Trying {source_name}..."))

            tryCatch(
                {
                    # Set timeout for the download
                    old_timeout <- getOption("timeout")
                    options(timeout = timeout)

                    utils::download.file(
                        url = data_url,
                        destfile = local_file,
                        mode = "wb",
                        quiet = FALSE # Show progress
                    )

                    cli::cli_alert_success(glue::glue(
                        "Successfully downloaded from {source_name}"
                    ))
                    success <- TRUE
                    break # Exit loop on success
                },
                error = function(e) {
                    if (file.exists(local_file)) {
                        unlink(local_file)
                    }

                    cli::cli_warn(glue::glue(
                        "Failed from {source_name}: {e$message}"
                    ))

                    # If this was the last source, show final error
                    if (i == length(available_urls)) {
                        options(timeout = old_timeout)
                        cli::cli_abort(c(
                            "x" = cli::col_red("All download attempts failed."),
                            "i" = "Please check your internet connection or try again later.",
                            "i" = "Error from last attempt: {e$message}"
                        ))
                    }
                }
            )

            if (success) break
        }
    } else {
        cli::cli_alert_info("Found cached data.")
    }

    data <- tryCatch(
        readRDS(local_file),
        error = function(e) {
            if (file.exists(local_file)) {
                unlink(local_file) # Clean up corrupted file
            }
            cli::cli_abort(c(
                "x" = cli::col_red("Downloaded file appears to be corrupted."),
                "i" = "Please try again. Error: {e$message}"
            ))
        }
    )
    cli::cli_alert_success(cli::col_green("Data loaded successfully."))

    if (!cache) {
        on.exit(
            {
                options(timeout = old_timeout)

                if (file.exists(local_file)) {
                    unlink(local_file)
                }
            },
            add = TRUE
        )
    }

    return(data)
}
