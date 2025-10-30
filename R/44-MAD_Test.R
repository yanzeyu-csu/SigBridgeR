#' @title outlier detection using Median Absolute Deviation (MAD)
#'
#' @description
#' This function identifies outliers in a numeric vector using Median Absolute
#' Deviation (MAD) method, with specialized handling for very small samples.
#' Returns results in htest format for consistency with statistical tests.
#'
#' @param x A numeric vector. Missing values (`NA`) are allowed and will be
#'          excluded from calculations.
#' @param na.rm Logical. Should missing values be removed? Default is `TRUE`.
#'
#' @return An object of class "htest" containing:
#'   \item{statistic}{The test statistic (MAD-based z-score of most extreme outlier)}
#'   \item{parameter}{Named vector with: n, median, mad, threshold}
#'   \item{p.value}{Approximate p-value based on outlier probability}
#'   \item{method}{Description of the method used}
#'   \item{data.name}{Name of the data vector}
#'   \item{alternative}{Character string describing the alternative hypothesis}
#'   \item{outlier.indices}{Indices of detected outliers (as additional element)}
#'
#' @details
#' For samples with n > 5, uses standard MAD method with threshold of 3.
#' For very small samples (n <= 5), uses MAD with adjusted thresholds:
#' - n = 3-4: threshold = 2.5
#' - n = 5: threshold = 2.8
#'
#' The test statistic is calculated as the maximum absolute deviation from median
#' divided by adjusted MAD.
#'
#' @examples
#' \dontrun{
#' # Basic usage
#' x <- c(1, 2, 3, 100, 4, 5, 3)
#' mad.test(x)
#'
#' # Small sample
#' x3 <- c(1, 2, 10)
#' mad.test(x3)
#' }
#'
#' @import data.table
#' @importFrom chk chk_length chk_numeric chk_not_any_na
#'
#' @keywords internal
#'
mad.test <- function(x, na.rm = TRUE) {
    # Store data name for htest output
    data_name <- deparse(substitute(x))

    # Input validation
    chk::chk_numeric(x)

    # Handle missing values
    if (na.rm) {
        complete_cases <- !is.na(x)
        x_clean <- x[complete_cases]
        original_indices <- which(complete_cases)
    } else {
        chk::chk_not_any_na(x)
        x_clean <- x
        original_indices <- seq_along(x)
    }

    n <- length(x_clean)

    chk::chk_length(n)

    # Create data.table for efficient computation
    dt <- data.table::data.table(
        index = original_indices,
        value = x_clean
    )

    # Calculate median and MAD
    median_val <- dt[, stats::median(`value`)]
    abs_dev <- abs(dt$`value` - median_val)
    mad_val <- stats::median(abs_dev)
    mad_adjusted <- mad_val * 1.4826

    # Determine threshold based on sample size
    if (n <= 4) {
        threshold <- 2.5
        method_used <- "MAD-based outlier detection (small sample n<=4, threshold=2.5)"
    } else if (n == 5) {
        threshold <- 2.8
        method_used <- "MAD-based outlier detection (small sample n=5, threshold=2.8)"
    } else {
        threshold <- 3.0
        method_used <- "MAD-based outlier detection (threshold=3.0)"
    }

    # Calculate MAD-based z-scores
    dt[, `z_score` := abs(`value` - `median_val`) / `mad_adjusted`] # id mad_adjusted == 0, z_score will be NA

    # Identify outliers
    dt[, `is_outlier` := `z_score` > `threshold`]

    # Get outlier information
    outliers_dt <- dt[`is_outlier` == TRUE]
    n_outliers <- nrow(outliers_dt)

    # Calculate test statistic (z-score of most extreme outlier)
    if (n_outliers > 0) {
        statistic_val <- max(dt$z_score)
        # Get indices of outliers in original data (before NA removal)
        outlier_indices <- dt[`is_outlier` == TRUE, `index`]
    } else {
        statistic_val <- max(dt$z_score)
        outlier_indices <- integer(0)
    }

    # Calculate approximate p-value
    # Based on probability of observing such extreme value in normal distribution
    if (mad_adjusted > 0) {
        pval <- 2 * (1 - stats::pnorm(statistic_val))
        # Adjust for multiple testing (Bonferroni correction)
        pval <- min(pval * n, 1.0)
    } else {
        pval <- NA_real_
    }

    # Create htest structure
    result <- list(
        statistic = c("MAD z-score" = statistic_val),
        parameter = c(
            "n" = n,
            "median" = median_val,
            "mad" = mad_adjusted,
            "threshold" = threshold
        ),
        p.value = pval,
        method = method_used,
        data.name = data_name,
        alternative = "at least one value is an outlier",
        outlier.indices = outlier_indices
    )

    class(result) <- "htest"
    result
}
