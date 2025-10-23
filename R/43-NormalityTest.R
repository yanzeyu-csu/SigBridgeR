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
#' @note
#' This test cannot be used for small sample sizes (n < 500) due to the assumption
#' of normality. For smaller samples, consider using the Shapiro-Wilk test.
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
#' @keywords htest
#' @family normality_test
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
        m1 <- matrixStats::sum2(x) / n
        m2 <- matrixStats::sum2((x - m1)^2) / n
        m3 <- matrixStats::sum2((x - m1)^3) / n
        m4 <- matrixStats::sum2((x - m1)^4) / n
        b1 <- (m3 / m2^(3 / 2))^2
        b2 <- (m4 / m2^2)
        STATISTIC <- n * (b1 / 6 + (b2 - 3)^2 / 24)
    }

    if (!is.na(mean) & is.na(sd)) {
        m1 <- mean
        m2 <- matrixStats::sum2((x - m1)^2) / n
        m3 <- matrixStats::sum2((x - m1)^3) / n
        m4 <- matrixStats::sum2((x - m1)^4) / n
        b1 <- (m3 / m2^(3 / 2))^2
        b2 <- (m4 / m2^2)
        STATISTIC <- n * (b1 / 15 + (b2 - 3)^2 / 24)
    }

    if (is.na(mean) & !is.na(sd)) {
        m1 <- matrixStats::mean2(x)
        m2 <- sd^2
        m3 <- matrixStats::sum2((x - m1)^3) / n
        m4 <- matrixStats::sum2((x - m1)^4) / n
        b1 <- (m3 / m2^(3 / 2))^2
        b2 <- (m4 / m2^2)
        STATISTIC <- n * (b1 / 6 + (b2 - 3)^2 / 96)
    }

    if (!is.na(mean) & !is.na(sd)) {
        m1 <- mean
        m2 <- sd^2
        m3 <- matrixStats::sum2((x - m1)^3) / n
        m4 <- matrixStats::sum2((x - m1)^4) / n
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
#' @keywords htest
#' @family normality_test
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

    x_mean <- matrixStats::mean2(x)
    x_centered <- x - x_mean
    m2 <- matrixStats::sum2(x_centered^2) / n
    m3 <- matrixStats::sum2(x_centered^3) / n
    m4 <- matrixStats::sum2(x_centered^4) / n

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

#' @title Anderson-Darling Normality Test
#'
#' @description
#' Performs the Anderson-Darling test for normality to assess whether a given
#' sample comes from a normal distribution. This implementation follows the
#' methodology from the `nortest` R package.
#'
#' @param x A numeric vector of data values. Missing values will be
#'          automatically removed.
#'
#' @return A list with class `"htest"` containing the following components:
#' \itemize{
#'   \item `statistic` - The Anderson-Darling test statistic (A)
#'   \item `p.value` - The p-value for the test
#'   \item `method` - The name of the method ("Anderson-Darling normality test")
#'   \item `data.name` - The name of the data used in the test
#' }
#'
#' @details
#' The Anderson-Darling test is a statistical test of whether a given sample
#' of data is drawn from a normal distribution. The test is a modification
#' of the Kolmogorov-Smirnov (K-S) test and gives more weight to the tails
#' than the K-S test.
#'
#' This implementation requires a minimum sample size of 8 observations.
#' The test statistic is calculated as:
#' \deqn{A^2 = -n - \frac{1}{n}\sum_{i=1}^n(2i-1)[\ln(p_i) + \ln(1-p_{n+1-i})]}
#' where \eqn{p_i = \Phi((x_i - \mu)/\sigma)} and \eqn{\Phi} is the cumulative
#' distribution function of the standard normal distribution.
#'
#' The p-value is calculated using different approximation formulas depending
#' on the value of the adjusted test statistic:
#' \itemize{
#'   \item For AA < 0.2: \eqn{p = 1 - \exp(-13.436 + 101.14 \times AA - 223.73 \times AA^2)}
#'   \item For 0.2 ≤ AA < 0.34: \eqn{p = 1 - \exp(-8.318 + 42.796 \times AA - 59.938 \times AA^2)}
#'   \item For 0.34 ≤ AA < 0.6: \eqn{p = \exp(0.9177 - 4.279 \times AA - 1.38 \times AA^2)}
#'   \item For 0.6 ≤ AA < 10: \eqn{p = \exp(1.2937 - 5.709 \times AA + 0.0186 \times AA^2)}
#'   \item For AA ≥ 10: \eqn{p = 3.7 \times 10^{-24}}
#' }
#' where \eqn{AA = A^2 \times (1 + 0.75/n + 2.25/n^2)}.
#'
#' @references
#' This implementation is based on the `nortest` package:
#'
#' Gross, J. and Ligges, U. (2015). nortest: Tests for Normality. R package
#' version 1.0-4. https://CRAN.R-project.org/package=nortest
#'
#' The original test is described in:
#'
#' Anderson, T.W. and Darling, D.A. (1954). A Test of Goodness of Fit.
#' Journal of the American Statistical Association, 49(268), 765-769.
#'
#' Stephens, M.A. (1974). EDF Statistics for Goodness of Fit and Some
#' Comparisons. Journal of the American Statistical Association, 69(347), 730-737.
#'
#' @examples
#' # Test a sample from normal distribution
#' set.seed(123)
#' normal_data <- rnorm(100)
#' ad.test(normal_data)
#'
#' # Test a sample from non-normal distribution
#' exponential_data <- rexp(50)
#' ad.test(exponential_data)
#'
#' # Test with real data
#' if (require("datasets")) {
#'   ad.test(iris$Sepal.Length)
#' }
#'
#' @seealso
#' \code{\link[nortest]{ad.test}} for the original implementation in the nortest package.
#' \code{\link{shapiro.test}} for the Shapiro-Wilk normality test.
#' \code{\link{ks.test}} for the Kolmogorov-Smirnov test.
#'
#' @keywords htest
#' @family normality_test
#' @keywords internal
ad.test <- function(x) {
    DNAME <- deparse(substitute(x))
    x <- sort(x[complete.cases(x)])
    n <- length(x)
    if (n < 8) {
        stop("sample size must be greater than 7")
    }
    logp1 <- pnorm((x - mean(x)) / sd(x), log.p = TRUE)
    logp2 <- pnorm(-(x - mean(x)) / sd(x), log.p = TRUE)
    h <- (2 * seq(1:n) - 1) * (logp1 + rev(logp2))
    A <- -n - mean(h)
    AA <- (1 + 0.75 / n + 2.25 / n^2) * A
    if (AA < 0.2) {
        pval <- 1 - exp(-13.436 + 101.14 * AA - 223.73 * AA^2)
    } else if (AA < 0.34) {
        pval <- 1 - exp(-8.318 + 42.796 * AA - 59.938 * AA^2)
    } else if (AA < 0.6) {
        pval <- exp(0.9177 - 4.279 * AA - 1.38 * AA^2)
    } else if (AA < 10) {
        pval <- exp(1.2937 - 5.709 * AA + 0.0186 * AA^2)
    } else {
        pval <- 3.7e-24
    }
    RVAL <- list(
        statistic = c(A = A),
        p.value = pval,
        method = "Anderson-Darling normality test",
        data.name = DNAME
    )
    class(RVAL) <- "htest"
    return(RVAL)
}

#' @title Cramer-von Mises Normality Test
#'
#' @description
#' Performs the Cramer-von Mises test for normality to assess whether a given
#' sample comes from a normal distribution. This implementation follows the
#' methodology from the `nortest` R package.
#'
#' @param x A numeric vector of data values. Missing values will be
#'          automatically removed.
#'
#' @return A list with class `"htest"` containing the following components:
#' \itemize{
#'   \item `statistic` - The Cramer-von Mises test statistic (W)
#'   \item `p.value` - The p-value for the test
#'   \item `method` - The name of the method ("Cramer-von Mises normality test")
#'   \item `data.name` - The name of the data used in the test
#' }
#'
#' @details
#' The Cramer-von Mises test is a statistical test of whether a given sample
#' of data is drawn from a normal distribution. It is an empirical distribution
#' function (EDF) test that compares the empirical distribution function of the
#' sample with the cumulative distribution function of the normal distribution.
#'
#' This implementation requires a minimum sample size of 8 observations.
#' The test statistic is calculated as:
#' \deqn{W = \frac{1}{12n} + \sum_{i=1}^n \left(p_i - \frac{2i-1}{2n}\right)^2}
#' where \eqn{p_i = \Phi((x_i - \mu)/\sigma)} and \eqn{\Phi} is the cumulative
#' distribution function of the standard normal distribution.
#'
#' The p-value is calculated using different approximation formulas depending
#' on the value of the adjusted test statistic WW = W × (1 + 0.5/n):
#' \itemize{
#'   \item For WW < 0.0275: \eqn{p = 1 - \exp(-13.953 + 775.5 \times WW - 12542.61 \times WW^2)}
#'   \item For 0.0275 ≤ WW < 0.051: \eqn{p = 1 - \exp(-5.903 + 179.546 \times WW - 1515.29 \times WW^2)}
#'   \item For 0.051 ≤ WW < 0.092: \eqn{p = \exp(0.886 - 31.62 \times WW + 10.897 \times WW^2)}
#'   \item For 0.092 ≤ WW < 1.1: \eqn{p = \exp(1.111 - 34.242 \times WW + 12.832 \times WW^2)}
#'   \item For WW ≥ 1.1: \eqn{p = 7.37 \times 10^{-10}} with a warning
#' }
#'
#' @references
#' This implementation is based on the `nortest` package:
#'
#' Gross, J. and Ligges, U. (2015). nortest: Tests for Normality. R package
#' version 1.0-4. https://CRAN.R-project.org/package=nortest
#'
#' The original test is described in:
#'
#' Cramer, H. (1928). On the Composition of Elementary Errors.
#' Scandinavian Actuarial Journal, 1928(1), 13-74.
#'
#' Von Mises, R. (1931). Wahrscheinlichkeitsrechnung und ihre Anwendung
#' in der Statistik und theoretischen Physik. Leipzig: Deuticke.
#'
#' Stephens, M.A. (1970). Use of the Kolmogorov-Smirnov, Cramer-von Mises
#' and Related Statistics Without Extensive Tables. Journal of the Royal
#' Statistical Society. Series B (Methodological), 32(1), 115-122.
#'
#' @examples
#' # Test a sample from normal distribution
#' set.seed(123)
#' normal_data <- rnorm(100)
#' cvm.test(normal_data)
#'
#' # Test a sample from non-normal distribution
#' exponential_data <- rexp(50)
#' cvm.test(exponential_data)
#'
#' # Test with real data
#' if (require("datasets")) {
#'   cvm.test(iris$Sepal.Length)
#' }
#'
#' @seealso
#' \code{\link[nortest]{cvm.test}} for the original implementation in the nortest package.
#' \code{\link{ad.test}} for the Anderson-Darling normality test.
#' \code{\link{shapiro.test}} for the Shapiro-Wilk normality test.
#' \code{\link{ks.test}} for the Kolmogorov-Smirnov test.
#'
#' @keywords htest
#' @family normality_test
#' @keywords internal
cvm.test <- function(x) {
    DNAME <- deparse(substitute(x))
    x <- sort(x[complete.cases(x)])
    n <- length(x)
    if (n < 8) {
        stop("sample size must be greater than 7")
    }
    p <- pnorm((x - mean(x)) / sd(x))
    W <- (1 / (12 * n) + sum((p - (2 * seq(1:n) - 1) / (2 * n))^2))
    WW <- (1 + 0.5 / n) * W
    if (WW < 0.0275) {
        pval <- 1 - exp(-13.953 + 775.5 * WW - 12542.61 * WW^2)
    } else if (WW < 0.051) {
        pval <- 1 - exp(-5.903 + 179.546 * WW - 1515.29 * WW^2)
    } else if (WW < 0.092) {
        pval <- exp(0.886 - 31.62 * WW + 10.897 * WW^2)
    } else if (WW < 1.1) {
        pval <- exp(1.111 - 34.242 * WW + 12.832 * WW^2)
    } else {
        warning(
            "p-value is smaller than 7.37e-10, cannot be computed more accurately"
        )
        pval <- 7.37e-10
    }
    RVAL <- list(
        statistic = c(W = W),
        p.value = pval,
        method = "Cramer-von Mises normality test",
        data.name = DNAME
    )
    class(RVAL) <- "htest"
    return(RVAL)
}

#' @title Pearson Chi-Square Normality Test
#'
#' @description
#' Performs the Pearson chi-square test for normality to assess whether a given
#' sample comes from a normal distribution. This test divides the data into
#' classes and compares observed frequencies with expected frequencies under
#' the normal distribution.
#'
#' @param x A numeric vector of data values. Missing values will be
#'          automatically removed.
#' @param n.classes Number of classes to use for the chi-square test.
#'                  Defaults to \code{ceiling(2 * (n^(2 / 5)))}, where n is the
#'                  sample size. This formula provides a data-driven approach
#'                  to determining the number of classes.
#' @param adjust Logical indicating whether to adjust the degrees of freedom.
#'               If \code{TRUE} (default), 2 degrees of freedom are subtracted
#'               to account for estimated parameters (mean and standard deviation).
#'
#' @return A list with class \code{"htest"} containing the following components:
#' \itemize{
#'   \item \code{statistic} - The Pearson chi-square test statistic (P)
#'   \item \code{p.value} - The p-value for the test
#'   \item \code{method} - The name of the method ("Pearson chi-square normality test")
#'   \item \code{data.name} - The name of the data used in the test
#'   \item \code{n.classes} - The number of classes used in the test
#'   \item \code{df} - The degrees of freedom used for the chi-square distribution
#' }
#'
#' @details
#' The Pearson chi-square normality test is a classical goodness-of-fit test
#' that compares the observed frequency distribution with the expected frequency
#' distribution under the normal distribution assumption.
#'
#' The test procedure involves:
#' \enumerate{
#'   \item Estimating the mean and standard deviation from the sample
#'   \item Dividing the range of the normal distribution into \code{n.classes} equal-probability intervals
#'   \item Counting the number of observations falling into each interval
#'   \item Calculating the chi-square statistic:
#'         \deqn{P = \sum \frac{(O_i - E_i)^2}{E_i}}
#'         where \eqn{O_i} is the observed frequency and \eqn{E_i} is the expected frequency
#'   \item Comparing the test statistic to a chi-square distribution with
#'         \code{n.classes - 1 - dfd} degrees of freedom, where \code{dfd} is 2 if
#'         \code{adjust = TRUE} (accounting for estimated parameters) or 0 otherwise
#' }
#'
#' The default number of classes follows the recommendation:
#' \code{ceiling(2 * (n^(2 / 5)))}, which adapts to the sample size.
#'
#' @references
#' This implementation is based on the `nortest` package:
#'
#' Pearson, K. (1900). On the Criterion that a Given System of Deviations from
#' the Probable in the Case of a Correlated System of Variables is such that
#' it can be Reasonably Supposed to have Arisen from Random Sampling.
#' Philosophical Magazine, 50(302), 157-175.
#'
#' Moore, D.S. (1986). Tests of the chi-squared type. In: D'Agostino, R.B.
#' and Stephens, M.A. (eds.), Goodness-of-Fit Techniques, Marcel Dekker, New York.
#'
#' @examples
#' # Test a sample from normal distribution
#' set.seed(123)
#' normal_data <- rnorm(100)
#' pearson.test(normal_data)
#'
#' # Test with custom number of classes
#' pearson.test(normal_data, n.classes = 10)
#'
#' # Test without degrees of freedom adjustment
#' pearson.test(normal_data, adjust = FALSE)
#'
#' # Test a sample from non-normal distribution
#' exponential_data <- rexp(50)
#' pearson.test(exponential_data)
#'
#' @seealso
#' \code{\link{chisq.test}} for the general chi-square test of independence.
#' \code{\link{ad.test}} for the Anderson-Darling normality test.
#' \code{\link{cvm.test}} for the Cramer-von Mises normality test.
#' \code{\link{shapiro.test}} for the Shapiro-Wilk normality test.
#'
#' @keywords htest
#' @family normality_test
#' @keywords internal
#'
pearson.test <- function(
    x,
    n.classes = ceiling(2 * (n^(2 / 5))),
    adjust = TRUE
) {
    DNAME <- deparse(substitute(x))
    x <- x[complete.cases(x)]
    n <- length(x)
    if (adjust) {
        dfd <- 2
    } else {
        dfd <- 0
    }
    num <- floor(1 + n.classes * pnorm(x, mean(x), sd(x)))
    count <- tabulate(num, n.classes)
    prob <- rep(1 / n.classes, n.classes)
    xpec <- n * prob
    h <- ((count - xpec)^2) / xpec
    P <- sum(h)
    pvalue <- pchisq(P, n.classes - dfd - 1, lower.tail = FALSE)
    RVAL <- list(
        statistic = c(P = P),
        p.value = pvalue,
        method = "Pearson chi-square normality test",
        data.name = DNAME,
        n.classes = n.classes,
        df = n.classes -
            1 -
            dfd
    )
    class(RVAL) <- "htest"
    return(RVAL)
}

#' @title Shapiro-Francia Normality Test
#'
#' @description
#' Performs the Shapiro-Francia test for normality, which is similar to the
#' Shapiro-Wilk test but designed to handle larger sample sizes (up to 5000
#' observations).
#'
#' @param x A numeric vector of data values. Missing values will be
#'          automatically removed.
#'
#' @return A list with class \code{"htest"} containing the following components:
#' \itemize{
#'   \item \code{statistic} - The Shapiro-Francia test statistic (W)
#'   \item \code{p.value} - The p-value for the test
#'   \item \code{method} - The name of the method ("Shapiro-Francia normality test")
#'   \item \code{data.name} - The name of the data used in the test
#' }
#'
#' @details
#' The Shapiro-Francia test is a powerful test for normality that is particularly
#' useful for larger sample sizes where the Shapiro-Wilk test may not be applicable.
#' The test is based on the correlation between the ordered sample values and the
#' corresponding normal quantiles.
#'
#' The test procedure involves:
#' \enumerate{
#'   \item Sorting the data and computing normal quantiles using \code{\link{qnorm}}
#'         with \code{\link{ppoints}} (Blom's method with a = 3/8)
#'   \item Calculating the squared correlation coefficient between the ordered
#'         data and the normal quantiles: \eqn{W = [cor(x, y)]^2}
#'   \item Transforming the statistic using: \eqn{z = (\log(1 - W) - \mu) / \sigma}
#'   \item Computing the p-value from the standard normal distribution
#' }
#'
#' The parameters for the transformation are:
#' \deqn{\mu = -1.2725 + 1.0521 \times (\log(u) - u)}
#' \deqn{\sigma = 1.0308 - 0.26758 \times (\log(u) + 2/u)}
#' where \eqn{u = \log(n)} and \eqn{n} is the sample size.
#'
#' @note
#' The Shapiro-Francia test has the following limitations:
#' \itemize{
#'   \item Sample size must be between 5 and 5000
#'   \item The test may be less powerful than the Shapiro-Wilk test for small samples
#'   \item For sample sizes > 5000, consider using other tests like
#'         Anderson-Darling or Cramer-von Mises
#' }
#'
#' @references
#' This implementation is based on the `nortest` package:
#'
#' Shapiro, S.S. and Francia, R.S. (1972). An Approximate Analysis of Variance
#' Test for Normality. Journal of the American Statistical Association, 67(337), 215-216.
#'
#' Royston, P. (1993). A Pocket-Calculator Algorithm for the Shapiro-Francia Test
#' for Non-Normality: An Application to Medicine. Statistics in Medicine, 12(2), 181-184.
#'
#' Thode, H.C. (2002). Testing for Normality. Marcel Dekker, New York.
#'
#' @examples
#' # Test a sample from normal distribution
#' set.seed(123)
#' normal_data <- rnorm(100)
#' sf.test(normal_data)
#'
#' # Test a sample from non-normal distribution
#' exponential_data <- rexp(50)
#' sf.test(exponential_data)
#'
#' # Test with real data
#' if (require("datasets")) {
#'   sf.test(iris$Sepal.Length)
#' }
#'
#' @seealso
#' \code{\link{shapiro.test}} for the Shapiro-Wilk normality test (limited to 5000 observations).
#' \code{\link{ad.test}} for the Anderson-Darling normality test.
#' \code{\link{cvm.test}} for the Cramer-von Mises normality test.
#' \code{\link{pearson.test}} for the Pearson chi-square normality test.
#' \code{\link{ppoints}} for the generation of probability points.
#' \code{\link{qnorm}} for the normal quantile function.
#'
#' @keywords htest
#' @family normality_test
#' @keywords internal
#'
sf.test <- function(x) {
    DNAME <- deparse(substitute(x))
    x <- sort(x[complete.cases(x)])
    n <- length(x)
    if ((n < 5 || n > 5000)) {
        stop("sample size must be between 5 and 5000")
    }
    y <- qnorm(ppoints(n, a = 3 / 8))
    W <- cor(x, y)^2
    u <- log(n)
    v <- log(u)
    mu <- -1.2725 + 1.0521 * (v - u)
    sig <- 1.0308 - 0.26758 * (v + 2 / u)
    z <- (log(1 - W) - mu) / sig
    pval <- pnorm(z, lower.tail = FALSE)
    RVAL <- list(
        statistic = c(W = W),
        p.value = pval,
        method = "Shapiro-Francia normality test",
        data.name = DNAME
    )
    class(RVAL) <- "htest"
    return(RVAL)
}
