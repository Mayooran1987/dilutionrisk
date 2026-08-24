#' Probability of Acceptance for Heterogeneous Batches
#'
#' Calculates the probability of acceptance (PA) in the original sample when
#' diluted samples are collected from a heterogeneous (non-homogeneous) batch.
#'
#' @param c Acceptance number (maximum number of defective units allowed).
#' @param mu Mean microbial count on the log scale. Must satisfy
#'   \code{log(a) - log(f*u) <= mu <= log(b) - log(f*u)}.
#' @param sd Standard deviation of the normal distribution on the log scale.
#'   Must be positive.
#' @param a Lower bound of cell count domain. Must be non-negative.
#' @param b Upper bound of cell count domain. Must be greater than \code{a}.
#' @param f Final dilution factor. Must be between 0 and 1.
#' @param u Amount placed on the plate. Must be positive.
#' @param USL Upper specification limit for microbial count.
#' @param n Number of samples inspected. Must be positive integer.
#' @param type Type of calculation: "theory" (default) or "simulation".
#' @param n_sim Number of simulations. Required when \code{type = "simulation"}.
#'
#' @return Numeric value representing the probability of acceptance.
#'
#' @details
#' The probability of acceptance is calculated using the binomial distribution:
#' \deqn{P_a = P(X \le c) = \sum_{k=0}^{c} \binom{n}{k} p_d^k (1-p_d)^{n-k}}
#' where \eqn{p_d} is the probability of detection and \eqn{n} is the number
#' of samples inspected.
#'
#' @examples
#' # Basic usage
#' prob_acceptance_heterogeneous(c = 2, mu = 7, sd = 0.2, a = 0, b = 300,
#'                               f = 0.01, u = 0.1, USL = 1000, n = 5)
#'
#' # Multiple dilution schemes
#' prob_acceptance_heterogeneous_multiple(c = 2, mu = 7, sd = 0.2,
#'                                        a = 0, b = 300,
#'                                        f = c(0.01, 0.1), u = c(0.1, 0.1),
#'                                        USL = 1000, n = 5)
#'
#' @name prob_acceptance_heterogeneous
#' @aliases prob_acceptance_heterogeneous prob_acceptance_heterogeneous_multiple
#' @rdname prob_acceptance_heterogeneous
#' @export
prob_acceptance_heterogeneous <- function(c, mu, sd, a, b, f, u, USL, n,
                                          type = "theory", n_sim = NA) {
  # Input validation
  if (!is.numeric(c) || length(c) != 1 || c < 0 || c != round(c)) {
    stop("'c' must be a non-negative integer", call. = FALSE)
  }
  if (!is.numeric(n) || length(n) != 1 || n < 1 || n != round(n)) {
    stop("'n' must be a positive integer", call. = FALSE)
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
  if (!is.numeric(f) || length(f) != 1 || f <= 0 || f >= 1) {
    stop("'f' must be a numeric scalar between 0 and 1", call. = FALSE)
  }
  if (!is.numeric(u) || length(u) != 1 || u <= 0) {
    stop("'u' must be a positive numeric scalar", call. = FALSE)
  }
  if (!is.numeric(USL) || length(USL) != 1 || USL <= 0) {
    stop("'USL' must be a positive numeric scalar", call. = FALSE)
  }

  # Calculate probability of detection
  pd <- prob_detection_heterogeneous(mu, sd, a, b, f, u, USL, type, n_sim)

  # Calculate probability of acceptance using binomial distribution
  pa <- stats::pbinom(c, n, pd)

  return(pa)
}

#' @rdname prob_acceptance_heterogeneous
#' @export
prob_acceptance_heterogeneous_multiple <- function(c, mu, sd, a, b, f, u, USL, n,
                                                   type = "theory", n_sim = NA) {
  # Input validation
  if (length(f) != length(u)) {
    stop("'f' and 'u' must have equal length", call. = FALSE)
  }
  if (length(f) == 0) {
    stop("'f' must have at least one element", call. = FALSE)
  }

  results <- vapply(seq_along(f), function(i) {
    prob_acceptance_heterogeneous(c, mu, sd, a, b, f[i], u[i], USL, n,
                                  type, n_sim)
  }, numeric(1))

  return(matrix(results, nrow = 1))
}
