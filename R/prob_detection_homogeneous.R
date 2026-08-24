#' Probability of Detection for Homogeneous Batches
#'
#' Calculates the probability of detection (PD) in the original sample when
#' diluted samples are collected from a homogeneous batch.
#'
#' @param lambda Expected microbial count (\eqn{\lambda}). Must be positive.
#' @param a Lower bound of cell count domain. Must be non-negative.
#' @param b Upper bound of cell count domain. Must be greater than \code{a}.
#' @param f Final dilution factor. Must be between 0 and 1.
#' @param u Amount placed on the plate. Must be positive.
#' @param USL Upper specification limit for microbial count.
#' @param type Type of calculation: "theory" (default) or "simulation".
#' @param n_sim Number of simulations. Required when \code{type = "simulation"}.
#'
#' @return Numeric value representing the probability of detection.
#'
#' @details
#' The probability of detection is defined as the probability that the estimated
#' microbial count exceeds the upper specification limit (USL). For homogeneous
#' batches, the count follows a truncated Poisson distribution.
#'
#' The theoretical probability is calculated as:
#' \deqn{P_d = 1 - \sum_{i=1}^{USL \times f \times u} \frac{(\lambda f u)^i}{i!(e^{\lambda f u} - 1)}}
#'
#' @references
#' \itemize{
#'   \item Schothorst, M. van, Zwietering, M.H., Ross, T., Buchanan, R.L.,
#'     Cole, M.B., 2009. Relating microbiological criteria to food safety
#'     objectives and performance objectives. Food Control 20, 967–979.
#' }
#'
#' @examples
#' # Theoretical calculation
#' prob_detection_homogeneous(lambda = 2000, a = 0, b = 300,
#'                            f = 0.01, u = 0.1, USL = 1000)
#'
#' # Simulation-based calculation
#' prob_detection_homogeneous(lambda = 2000, a = 0, b = 300,
#'                            f = 0.01, u = 0.1, USL = 1000,
#'                            type = "simulation", n_sim = 50000)
#'
#' @name prob_detection_homogeneous
#' @aliases prob_detection_homogeneous prob_detection_homogeneous_multiple
#' @rdname prob_detection_homogeneous
#' @export
prob_detection_homogeneous <- function(lambda, a, b, f, u, USL,
                                       type = "theory", n_sim = NA) {
  # Input validation
  if (!is.numeric(lambda) || length(lambda) != 1 || lambda <= 0) {
    stop("'lambda' must be a positive numeric scalar", call. = FALSE)
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

  # Check type argument
  type <- match.arg(type, c("theory", "simulation"))

  if (type == "theory") {
    # Theoretical calculation
    USL1 <- USL * f * u

    # Handle case where USL1 < 1
    if (USL1 < 1) {
      warning("USL * f * u is less than 1, probability of detection may be 0",
              call. = FALSE)
      return(0)
    }

    # Calculate using stable method (log-sum)
    log_lambda_d <- log(lambda * f * u)

    # Sum probabilities using log-sum-exp for numerical stability
    log_probs <- vapply(seq_len(USL1), function(i) {
      log_lambda_d * i - lfactorial(i) - log(exp(lambda * f * u) - 1)
    }, numeric(1))

    # Use log-sum-exp for numerical stability
    max_log <- max(log_probs)
    sum_exp <- sum(exp(log_probs - max_log))
    sum_prob <- max_log + log(sum_exp)

    Pd <- 1 - exp(sum_prob)

    # Ensure result is within [0, 1]
    return(max(0, min(1, Pd)))

  } else if (type == "simulation") {
    # Simulation-based calculation
    if (is.na(n_sim) || !is.numeric(n_sim) || n_sim < 100) {
      stop("When type = 'simulation', 'n_sim' must be a numeric value >= 100",
           call. = FALSE)
    }

    # Check if lambda * f * u is within bounds
    lambda_d <- lambda * f * u
    if (lambda_d < a || lambda_d > b) {
      warning(sprintf(
        "lambda * f * u = %.2f is outside the bounds [%.0f, %.0f].
        Consider adjusting parameters.",
        lambda_d, a, b
      ), call. = FALSE)
    }

    # Simulation
    sim1 <- matrix(NA, nrow = n_sim, ncol = 1)

    for (j in seq_len(n_sim)) {
      sim1[j, ] <- rtpois(1, lambda * f * u, a, b) * (1 / (f * u))
    }

    pd <- length(which(sim1[, 1] > USL)) / n_sim
    return(pd)
  }
}

#' @rdname prob_detection_homogeneous
#' @export
prob_detection_homogeneous_multiple <- function(lambda, a, b, f, u, USL,
                                                type = "theory", n_sim = NA) {
  # Input validation
  if (length(f) != length(u)) {
    stop("'f' and 'u' must have equal length", call. = FALSE)
  }
  if (length(f) == 0) {
    stop("'f' must have at least one element", call. = FALSE)
  }

  results <- vapply(seq_along(f), function(i) {
    prob_detection_homogeneous(lambda, a, b, f[i], u[i], USL, type, n_sim)
  }, numeric(1))

  return(matrix(results, nrow = 1))
}
