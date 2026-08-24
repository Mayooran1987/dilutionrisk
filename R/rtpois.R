#' Truncated Poisson Random Number Generation
#'
#' Generates random deviates from a truncated Poisson distribution.
#' The function uses rejection sampling to generate observations from a Poisson
#' distribution conditioned to lie within the interval [a, b].
#'
#' @param n Number of observations. Must be a positive integer.
#' @param lambda Expected count parameter. Must be non-negative.
#' @param a Lower truncation point. Must be non-negative.
#' @param b Upper truncation point. Must be greater than \code{a}.
#'
#' @return A vector of random numbers from the truncated Poisson distribution.
#'   Returns \code{NA} if rejection sampling fails to converge.
#'
#' @details
#' The truncated Poisson distribution has probability mass function:
#' \deqn{P(X = x) = \frac{\lambda^x e^{-\lambda}}{x! \sum_{k=a}^{b} \frac{\lambda^k e^{-\lambda}}{k!}}}
#' for \eqn{a \le x \le b}.
#'
#' The function uses rejection sampling with a maximum of 1,000,000 iterations
#' to ensure convergence. For lambda = 0, all observations are 0 if within bounds.
#'
#' @references
#' \itemize{
#'   \item Johnson, N. L., Kemp, A. W., & Kotz, S. (2005). Univariate Discrete
#'     Distributions (3rd ed.). Wiley.
#' }
#'
#' @examples
#' # Generate 10 observations from truncated Poisson with lambda = 5
#' # bounded between 0 and 300
#' set.seed(123)
#' rtpois(n = 10, lambda = 5, a = 0, b = 300)
#'
#' # Generate 5 observations with lower bound 10 and upper bound 50
#' rtpois(n = 5, lambda = 20, a = 10, b = 50)
#'
#' # Lambda = 0 case (all zeros)
#' rtpois(n = 5, lambda = 0, a = 0, b = 300)
#'
#' # Generate for different lambda values
#' lambdas <- c(1, 5, 10, 20)
#' set.seed(456)
#' for (lam in lambdas) {
#'   cat("Lambda =", lam, ":",
#'       rtpois(n = 3, lambda = lam, a = 0, b = 100), "\n")
#' }
#'
#' # Compare with regular Poisson (when bounds are wide)
#' set.seed(789)
#' trunc_sample <- rtpois(n = 1000, lambda = 10, a = 0, b = 300)
#' regular_sample <- rpois(n = 1000, lambda = 10)
#'
#' # Both should be similar when bounds are not restrictive
#' cat("Truncated mean:", mean(trunc_sample), "\n")
#' cat("Regular mean:", mean(regular_sample), "\n")
#'
#' # When bounds are restrictive (0-10 only)
#' set.seed(321)
#' trunc_restricted <- rtpois(n = 1000, lambda = 20, a = 0, b = 10)
#' cat("Restricted truncated mean:", mean(trunc_restricted), "\n")
#'
#' # Note: Rejection sampling may be slow for very restrictive bounds
#' # or very large lambda values
#'
#' @export
rtpois <- function(n, lambda, a, b) {
  if (!is.numeric(n) || length(n) == 0) {
    stop("'n' must be a positive integer", call. = FALSE)
  }

  if (length(n) > 1) {
    n <- length(n)
  } else {
    if (n < 1 || n != round(n)) {
      stop("'n' must be a positive integer", call. = FALSE)
    }
  }

  if (!is.numeric(lambda) || length(lambda) != 1 || lambda < 0) {
    stop("'lambda' must be a non-negative numeric scalar", call. = FALSE)
  }

  if (!is.numeric(a) || length(a) != 1 || a < 0) {
    stop("'a' must be a non-negative numeric scalar", call. = FALSE)
  }
  if (!is.numeric(b) || length(b) != 1 || b <= a) {
    stop("'b' must be a numeric scalar greater than 'a'", call. = FALSE)
  }

  if (lambda == 0) {
    if (a <= 0 && b >= 0) {
      return(rep(0, n))
    } else {
      warning("lambda = 0 but 0 is not within bounds [", a, ", ", b, "]", call. = FALSE)
      return(rep(NA_real_, n))
    }
  }

  result <- numeric(n)
  for (i in seq_len(n)) {
    x <- stats::rpois(1, lambda)
    iter <- 0
    while ((x < a || x > b) && iter < 1000000) {
      iter <- iter + 1
      x <- stats::rpois(1, lambda)
    }
    if (iter >= 1000000) {
      warning("Rejection sampling reached maximum iterations for lambda = ", lambda,
              ", bounds = [", a, ", ", b, "]", call. = FALSE)
      result[i] <- NA
    } else {
      result[i] <- x
    }
  }
  return(result)
}

#' Truncated Poisson-Lognormal Random Number Generation
#'
#' Generates random deviates from a truncated Poisson-lognormal distribution.
#' This is a compound distribution where the Poisson mean parameter follows
#' a lognormal distribution, conditioned on the Poisson outcome being within
#' the interval [a, b].
#'
#' @param n Number of observations. Must be a positive integer.
#' @param mu Mean of the lognormal distribution on the log scale.
#' @param sd Standard deviation of the lognormal distribution on the log scale.
#'   Must be positive.
#' @param a Lower truncation point (lower bound of cell count).
#'   Must be non-negative.
#' @param b Upper truncation point (upper bound of cell count).
#'   Must be greater than \code{a}.
#'
#' @return A vector of random numbers from the truncated Poisson-lognormal
#'   distribution. Returns \code{NA} if rejection sampling fails to converge.
#'
#' @details
#' The Poisson-lognormal distribution is a mixture distribution where:
#' \itemize{
#'   \item \eqn{\Lambda \sim \text{Lognormal}(\mu, \sigma^2)}
#'   \item \eqn{X | \Lambda = \lambda \sim \text{Poisson}(\lambda)}
#' }
#'
#' The truncation limits the support to \eqn{a \le X \le b}.
#'
#' The probability mass function is given by:
#' \deqn{P(X = x) = \frac{1}{x! \sqrt{2\pi}\sigma}
#'       \int_0^\infty \lambda^{x-1} e^{-\lambda}
#'       \exp\left(-\frac{(\log(\lambda)-\mu)^2}{2\sigma^2}\right) d\lambda}
#'
#' This distribution is particularly useful for modeling microbial counts
#' in heterogeneous (non-homogeneous) batches where the contamination level
#' varies across the batch.
#'
#' @references
#' \itemize{
#'   \item Gonzales-Barron, U.A., Pilão Cadavez, V.A., Butler, F., 2013.
#'     Statistical approaches for the design of sampling plans for
#'     microbiological monitoring of foods. In: Mathematical and Statistical
#'     Methods in Food Science and Technology. Wiley, pp. 363–384.
#'   \item Bulmer, M. G. (1974). On fitting the Poisson lognormal distribution
#'     to species-abundance data. Biometrics, 30(1), 101-110.
#' }
#'
#' @examples
#' # Generate 10 observations from truncated Poisson-lognormal
#' # with mu = 2, sd = 0.5, bounded between 0 and 300
#' set.seed(123)
#' rtrunpoilog(n = 10, mu = 2, sd = 0.5, a = 0, b = 300)
#'
#' # Generate with different parameters
#' set.seed(456)
#' rtrunpoilog(n = 5, mu = 1, sd = 0.3, a = 0, b = 100)
#'
#' # Compare with regular Poisson-lognormal (when bounds are wide)
#' set.seed(789)
#' n_sim <- 1000
#' trunc_sample <- rtrunpoilog(n = n_sim, mu = 2, sd = 0.5, a = 0, b = 300)
#' regular_sample <- rpois(n = n_sim, lambda = rlnorm(n_sim, 2, 0.5))
#'
#' cat("Truncated mean:", mean(trunc_sample), "\n")
#' cat("Regular mean:", mean(regular_sample), "\n")
#' cat("Truncated variance:", var(trunc_sample), "\n")
#' cat("Regular variance:", var(regular_sample), "\n")
#'
#' # When bounds are restrictive (0-10 only)
#' set.seed(321)
#' trunc_restricted <- rtrunpoilog(n = n_sim, mu = 3, sd = 0.8, a = 0, b = 10)
#' cat("Restricted truncated mean:", mean(trunc_restricted), "\n")
#'
#' # Generate samples with different sd values
#' mu_val <- 2
#' sd_values <- c(0.2, 0.5, 1.0)
#' set.seed(654)
#' for (sd_val in sd_values) {
#'   sample <- rtrunpoilog(n = 10, mu = mu_val, sd = sd_val, a = 0, b = 300)
#'   cat("mu =", mu_val, "sd =", sd_val, ":",
#'       summary(sample), "\n")
#' }
#'
#' # Note: Rejection sampling may be slow for very restrictive bounds
#' # or when sd is very large
#'
#' @export
rtrunpoilog <- function(n, mu, sd, a, b) {
  if (!is.numeric(n) || length(n) == 0) {
    stop("'n' must be a positive integer", call. = FALSE)
  }

  if (length(n) > 1) {
    n <- length(n)
  } else {
    if (n < 1 || n != round(n)) {
      stop("'n' must be a positive integer", call. = FALSE)
    }
  }

  if (!is.numeric(mu) || length(mu) != 1) {
    stop("'mu' must be a numeric scalar", call. = FALSE)
  }
  if (!is.numeric(sd) || length(sd) != 1 || sd <= 0) {
    stop("'sd' must be a positive numeric scalar", call. = FALSE)
  }
  if (!is.numeric(a) || length(a) != 1 || a < 0) {
    stop("'a' must be a non-negative numeric scalar", call. = FALSE)
  }
  if (!is.numeric(b) || length(b) != 1 || b <= a) {
    stop("'b' must be a numeric scalar greater than 'a'", call. = FALSE)
  }

  result <- numeric(n)
  for (i in seq_len(n)) {
    lambda <- stats::rlnorm(1, meanlog = mu, sdlog = sd)
    x <- stats::rpois(1, lambda)
    iter <- 0
    while ((x < a || x > b) && iter < 1000000) {
      iter <- iter + 1
      lambda <- stats::rlnorm(1, meanlog = mu, sdlog = sd)
      x <- stats::rpois(1, lambda)
    }
    if (iter >= 1000000) {
      warning("Rejection sampling reached maximum iterations for mu = ", mu,
              ", sd = ", sd, ", bounds = [", a, ", ", b, "]", call. = FALSE)
      result[i] <- NA
    } else {
      result[i] <- x
    }
  }
  return(result)
}

#' Compare Truncated and Regular Distributions
#'
#' Helper function to visualize the effect of truncation on Poisson
#' and Poisson-lognormal distributions.
#'
#' @param n Number of observations for simulation.
#' @param lambda Poisson parameter for comparison.
#' @param mu Mean parameter for Poisson-lognormal.
#' @param sd Standard deviation parameter for Poisson-lognormal.
#' @param a Lower bound for truncation.
#' @param b Upper bound for truncation.
#' @param distribution Type of distribution: "poisson" or "poilog".
#' @param ... Additional arguments passed to plot.
#'
#' @return Invisibly returns a list of generated samples.
#' @examples
#' # Compare truncated vs regular Poisson
#' \dontrun{
#' compare_truncation(n = 1000, lambda = 10,
#'                    a = 0, b = 20,
#'                    distribution = "poisson")
#'
#' # Compare truncated vs regular Poisson-lognormal
#' compare_truncation(n = 1000, mu = 2, sd = 0.5,
#'                    a = 0, b = 50,
#'                    distribution = "poilog")
#' }
#' @keywords internal
compare_truncation <- function(n, lambda = NULL, mu = NULL, sd = NULL,
                               a = 0, b = 300,
                               distribution = c("poisson", "poilog"),
                               ...) {
  distribution <- match.arg(distribution)

  if (distribution == "poisson") {
    if (is.null(lambda)) stop("lambda must be provided for Poisson")
    regular <- stats::rpois(n, lambda)
    truncated <- rtpois(n, lambda, a, b)
    title <- "Truncated vs Regular Poisson"
    param <- paste("lambda =", lambda)
  } else {
    if (is.null(mu) || is.null(sd)) {
      stop("mu and sd must be provided for Poisson-lognormal")
    }
    regular <- stats::rpois(n, stats::rlnorm(n, mu, sd))
    truncated <- rtrunpoilog(n, mu, sd, a, b)
    title <- "Truncated vs Regular Poisson-Lognormal"
    param <- paste("mu =", mu, ", sd =", sd)
  }

  par(mfrow = c(1, 2))
  hist(regular, main = paste("Regular", param),
       xlab = "Count", col = "lightblue", ...)
  hist(truncated, main = paste("Truncated", param),
       xlab = "Count", col = "lightgreen", ...)

  cat("Regular mean:", mean(regular), "var:", var(regular), "\n")
  cat("Truncated mean:", mean(truncated), "var:", var(truncated), "\n")

  invisible(list(regular = regular, truncated = truncated))
}
