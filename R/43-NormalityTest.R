#' @title Modified Jarque-Bera Test for Normailty
#'
#' @description
#' Performs a modified version of the Jarque-Bera test for normality that
#' provides enhanced sensitivity to departures from normal distribution.
#' This implementation includes four different test statistics based on
#' whether population mean and standard deviation are known or estimated.
#'
#' @param x A numeric vector of data values to test for normality.
#' @param mean An optional numeric value specifying the known population mean.
#'   If `NA` (default), the sample mean is used.
#' @param sd An optional numeric value specifying the known population standard
#'   deviation. If `NA` (default), the sample standard deviation is used.
#'
#' @return
#' An object of class `"htest"` containing the following components:
#' \itemize{
#'   \item `statistic` - The test statistic (chi-squared distributed)
#'   \item `parameter` - Degrees of freedom (always 2)
#'   \item `p.value` - The p-value for the test
#'   \item `method` - Description of the test method
#'   \item `data.name` - Name of the data vector
#' }
#'
#' @details
#' The modified Jarque-Bera test evaluates normality by testing whether
#' the sample data has the skewness and kurtosis matching a normal distribution.
#' The test statistic is based on the sample size, skewness, and kurtosis,
#' and follows a chi-squared distribution with 2 degrees of freedom.
#'
#' ## Test Variants:
#' The function implements four different test statistics depending on
#' whether the population mean and standard deviation are specified:
#'
#' - **Both unknown**: Uses sample estimates for mean and variance
#' - **Mean known, SD unknown**: Uses known mean but estimates variance
#' - **Mean unknown, SD known**: Uses sample mean with known variance
#' - **Both known**: Uses known population parameters
#'
#' Each variant uses different coefficients in the test statistic calculation
#' to account for the uncertainty in parameter estimation.
#'
#' ## Mathematical Formulation:
#' The test statistic is calculated as:
#' \deqn{JB = n \times \left( \frac{S^2}{6} + \frac{(K - 3)^2}{24} \right)}
#' where \eqn{S} is skewness, \eqn{K} is kurtosis, and \eqn{n} is sample size,
#' with modified coefficients for different parameter knowledge scenarios.
#'
#' @section Reference:
#' KhrushchevSergey/modified_jarque_bera_test Internet. cited 2025 Sep 28.
#' Available from: \url{https://github.com/KhrushchevSergey/Modified-Jarque-Bera-test/blob/main/modified_jarque_bera_test.R}
#'
#' @examples
#' \dontrun{
#' # Test with sample data (both parameters unknown)
#' x <- rnorm(100)
#' test_result <- jb.test.modified(x)
#' print(test_result)
#'
#' # Test with known population parameters
#' x <- rnorm(100, mean = 5, sd = 2)
#' test_result <- jb.test.modified(x, mean = 5, sd = 2)
#' print(test_result)
#'
#' # Test with only known mean
#' test_result <- jb.test.modified(x, mean = 5)
#' print(test_result)
#' }
#'
#' @seealso
#' [stats::shapiro.test()] for the Shapiro-Wilk test
#'
#' @keywords internal
jb.test.modified <- function(x, mean = NA, sd = NA) {
    if ((NCOL(x) > 1) || is.data.frame(x)) {
        cli::cli_abort(c("x" = "x is not a vector or univariate time series"))
    }
    if (anyNA(x)) {
        cli::cli_abort(c("x" = "x contains {.val NA}"))
    }
    DNAME <- deparse(substitute(x))
    n <- length(x)

    if (is.na(mean) & is.na(sd)) {
        m1 <- sum(x) / n
        m2 <- sum((x - m1)^2) / n
        m3 <- sum((x - m1)^3) / n
        m4 <- sum((x - m1)^4) / n
        b1 <- (m3 / m2^(3 / 2))^2
        b2 <- (m4 / m2^2)
        STATISTIC <- n * (b1 / 6 + (b2 - 3)^2 / 24)
    }

    if (!is.na(mean) & is.na(sd)) {
        m1 <- mean
        m2 <- sum((x - m1)^2) / n
        m3 <- sum((x - m1)^3) / n
        m4 <- sum((x - m1)^4) / n
        b1 <- (m3 / m2^(3 / 2))^2
        b2 <- (m4 / m2^2)
        STATISTIC <- n * (b1 / 15 + (b2 - 3)^2 / 24)
    }

    if (is.na(mean) & !is.na(sd)) {
        m1 <- mean(x)
        m2 <- sd^2
        m3 <- sum((x - m1)^3) / n
        m4 <- sum((x - m1)^4) / n
        b1 <- (m3 / m2^(3 / 2))^2
        b2 <- (m4 / m2^2)
        STATISTIC <- n * (b1 / 6 + (b2 - 3)^2 / 96)
    }

    if (!is.na(mean) & !is.na(sd)) {
        m1 <- mean
        m2 <- sd^2
        m3 <- sum((x - m1)^3) / n
        m4 <- sum((x - m1)^4) / n
        b1 <- (m3 / m2^(3 / 2))^2
        b2 <- (m4 / m2^2)
        STATISTIC <- n * (b1 / 15 + (b2 - 3)^2 / 96)
    }

    PVAL <- 1 - stats::pchisq(STATISTIC, df = 2)
    PARAMETER <- 2
    METHOD <- "Modified Jarque Bera Test"
    names(STATISTIC) <- "X-squared"
    names(PARAMETER) <- "df"
    structure(
        list(
            statistic = STATISTIC,
            parameter = PARAMETER,
            p.value = PVAL,
            method = METHOD,
            data.name = DNAME
        ),
        class = "htest"
    )
}


#' @title D'Agostino Test of Normality
#'
#' @description
#' Performs the D'Agostino test for normality, which is particularly suitable
#' for large sample sizes where other tests like Shapiro-Wilk may be too sensitive.
#'
#' @param x a numeric vector of data values. Missing values are allowed and will be removed.
#'
#' @return A list with class `"htest"` containing the following components:
#' \item{statistic}{a named vector containing three test statistics:
#'   \describe{
#'     \item{Skewness}{the Z statistic for skewness test}
#'     \item{Kurtosis}{the Z statistic for kurtosis test}
#'     \item{Omnibus}{the chi-squared statistic combining both tests}
#'   }
#' }
#' \item{parameter}{the degrees of freedom for the omnibus test}
#' \item{p.value}{a named vector containing three p-values corresponding to the statistics}
#' \item{method}{the character string "D'Agostino Normality Test"}
#' \item{data.name}{a character string giving the name(s) of the data}
#' \item{alternative}{a character string describing the alternative hypothesis}
#' \item{estimates}{a named vector containing sample statistics:
#'   \describe{
#'     \item{n}{sample size}
#'     \item{Skewness}{sample skewness coefficient}
#'     \item{Kurtosis}{sample kurtosis coefficient (excess kurtosis)}
#'     \item{Mean}{sample mean}
#'     \item{SD}{sample standard deviation}
#'   }
#' }
#'
#' @details
#' The D'Agostino test is a powerful omnibus test for normality that combines
#' separate tests for skewness and kurtosis. It is based on the transformation
#' of the sample skewness and kurtosis to approximate normal distributions,
#' which are then combined into a chi-squared statistic.
#'
#' This test is particularly recommended for:
#' \itemize{
#'   \item Large sample sizes (n > 50)
#'   \item Situations where Shapiro-Wilk test is overly sensitive
#'   \item Testing both tail behavior and symmetry simultaneously
#' }
#'
#' The test has good power against a wide range of alternative distributions
#' and is less sensitive to sample size fluctuations than many other normality tests.
#'
#' @examples
#' \dontrun{
#' # Test normally distributed data
#' set.seed(123)
#' normal_data <- rnorm(100)
#' result <- dagostino.test(normal_data)
#' print(result)
#'
#' # Test non-normal data (exponential distribution)
#' non_normal_data <- rexp(100)
#' result2 <- dagostino.test(non_normal_data)
#' print(result2)
#'
#' # Access specific components
#' result$p.value["Omnibus"]  # Overall test p-value
#' result$estimates["Skewness"]  # Sample skewness
#' result$statistic["Kurtosis"]  # Kurtosis test statistic
#'
#' # Large sample performance
#' large_sample <- rnorm(5000)
#' dagostino.test(large_sample)
#' }
#'
#'
#' @keywords htest
#' @keywords internal
#'
dagostino.test <- function(x) {
    if (!is.numeric(x)) {
        stop("'x' must be a numeric vector")
    }
    if (length(x) < 8) {
        stop("'x' must have at least 8 elements")
    }

    x <- x[!is.na(x)]
    n <- as.double(length(x))

    x_mean <- mean(x)
    x_centered <- x - x_mean
    m2 <- sum(x_centered^2) / n
    m3 <- sum(x_centered^3) / n
    m4 <- sum(x_centered^4) / n

    g1 <- m3 / (m2^(3 / 2))

    n_dbl <- as.double(n)
    n2 <- n_dbl * n_dbl
    n3 <- n2 * n_dbl

    Y <- g1 * sqrt((n_dbl + 1) * (n_dbl + 3) / (6 * (n_dbl - 2)))
    beta2 <- 3 *
        (n2 + 27 * n_dbl - 70) *
        (n_dbl + 1) *
        (n_dbl + 3) /
        ((n_dbl - 2) * (n_dbl + 5) * (n_dbl + 7) * (n_dbl + 9))

    if (beta2 <= 1) {
        W_sq <- 1.0
    } else {
        W_sq <- sqrt(2 * beta2 - 2)
    }

    if (W_sq <= 1) {
        Z_g1 <- 0
    } else {
        delta <- 1 / sqrt(log(W_sq))
        alpha <- sqrt(2 / (W_sq - 1))

        ratio <- Y / alpha
        Z_g1 <- delta * log(ratio + sqrt(ratio^2 + 1))
    }

    g2 <- m4 / (m2^2) - 3

    E_g2 <- -6 / (n_dbl + 1)
    Var_g2 <- 24 *
        n_dbl *
        (n_dbl - 2) *
        (n_dbl - 3) /
        ((n_dbl + 1)^2 * (n_dbl + 3) * (n_dbl + 5))

    if (Var_g2 <= 0) {
        standardized_g2 <- 0
    } else {
        standardized_g2 <- (g2 - E_g2) / sqrt(Var_g2)
    }

    if (n_dbl <= 3) {
        beta2_kurt <- 1.0
    } else {
        beta2_kurt <- 6 *
            (n2 - 5 * n_dbl + 2) /
            ((n_dbl + 7) * (n_dbl + 9)) *
            sqrt(
                6 *
                    (n_dbl + 3) *
                    (n_dbl + 5) /
                    (n_dbl * (n_dbl - 2) * (n_dbl - 3))
            )
    }

    if (beta2_kurt <= 0 || is.infinite(beta2_kurt)) {
        Z_g2 <- 0
    } else {
        A <- 6 +
            (8 / beta2_kurt) * (2 / beta2_kurt + sqrt(1 + 4 / (beta2_kurt^2)))

        if (A <= 4) {
            Z_g2 <- 0
        } else {
            term1 <- 1 - 2 / (9 * A)
            term2 <- (1 - 2 / A) / (1 + standardized_g2 * sqrt(2 / (A - 4)))
            if (term2 <= 0) {
                Z_g2 <- 0
            } else {
                Z_g2 <- (term1 - term2^(1 / 3)) / sqrt(2 / (9 * A))
            }
        }
    }

    K_sq <- Z_g1^2 + Z_g2^2

    p_value_skew <- 2 * stats::pnorm(-abs(Z_g1))
    p_value_kurt <- 2 * stats::pnorm(-abs(Z_g2))
    p_value_combined <- stats::pchisq(K_sq, df = 2, lower.tail = FALSE)

    p_value_skew <- max(0, min(1, p_value_skew))
    p_value_kurt <- max(0, min(1, p_value_kurt))
    p_value_combined <- max(0, min(1, p_value_combined))

    DNAME <- deparse(substitute(x))
    METHOD <- "D'Agostino Normality Test"

    result <- structure(
        list(
            statistic = c(
                Skewness = Z_g1,
                Kurtosis = Z_g2,
                Omnibus = K_sq
            ),
            parameter = c(df = 2),
            p.value = c(
                Skewness = p_value_skew,
                Kurtosis = p_value_kurt,
                Omnibus = p_value_combined
            ),
            method = METHOD,
            data.name = DNAME,
            alternative = "data are not normally distributed",
            estimates = c(
                n = n,
                Skewness = g1,
                Kurtosis = g2,
                Mean = x_mean,
                SD = sqrt(m2 * n / (n - 1))
            )
        ),
        class = "htest"
    )
    return(result)
}
