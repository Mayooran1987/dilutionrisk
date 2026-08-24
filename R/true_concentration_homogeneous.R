#' True Concentration Estimation for Homogeneous Batches
#'
#' Estimates the true microbial concentration in the original sample when
#' diluted samples are collected from a homogeneous batch.
#'
#' @param lambda Expected microbial count (\eqn{\lambda}). Must be positive.
#' @param lambda_low Lower bound of expected microbial count for x-axis.
#' @param lambda_high Upper bound of expected microbial count for x-axis.
#' @param a Lower bound of cell count domain. Must be non-negative.
#' @param b Upper bound of cell count domain. Must be greater than \code{a}.
#' @param f Final dilution factor. Must be between 0 and 1.
#' @param u Amount placed on the plate. Must be positive.
#' @param USL Upper specification limit for microbial count.
#' @param n_sim Number of simulations. Larger values provide more precise
#'   estimates. Recommended minimum: 10,000.
#'
#' @return Numeric value representing the estimated true concentration.
#'
#' @details
#' The true concentration is estimated as:
#' \deqn{C = \frac{X}{f \times u}}
#' where X is the count of microorganisms on a plate following a truncated
#' Poisson distribution, f is the final dilution factor, and u is the amount
#' placed on the plate.
#'
#' @examples
#' # Basic usage
#' true_concentration_homogeneous(lambda = 2000, a = 0, b = 300,
#'                                f = 0.01, u = 0.1, USL = 1000, n_sim = 5000)
#'
#' # Multiple dilution schemes
#' true_concentration_homogeneous_multiple(lambda = 2000, a = 0, b = 300,
#'                                         f = c(0.01, 0.1), u = c(0.1, 0.1),
#'                                         USL = 1000, n_sim = 5000)
#'
#' # Plot concentration curves
#' true_concentration_curves_homogeneous(lambda_low = 0, lambda_high = 5000,
#'                                       a = 0, b = 300, f = c(0.01, 0.1),
#'                                       u = c(0.1, 0.1), USL = 1000,
#'                                       n_sim = 5000)
#'
#'
#' @name true_concentration_homogeneous
#' @aliases true_concentration_homogeneous true_concentration_homogeneous_multiple
#' @rdname true_concentration_homogeneous
#' @export
true_concentration_homogeneous <- function(lambda, a, b, f, u, USL, n_sim) {
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
  if (!is.numeric(n_sim) || length(n_sim) != 1 || n_sim < 100) {
    warning("'n_sim' should be at least 100 for reliable results", call. = FALSE)
  }

  # Simulation
  sim1 <- matrix(NA, nrow = n_sim, ncol = 1)

  for (j in seq_len(n_sim)) {
    sim1[j, ] <- rtpois(1, lambda * f * u, a, b) * (1 / (f * u))
  }

  mean_concentr <- mean(sim1)
  return(mean_concentr)
}

#' @rdname true_concentration_homogeneous
#' @export
true_concentration_homogeneous_multiple <- function(lambda, a, b, f, u, USL, n_sim) {
  # Input validation
  if (length(f) != length(u)) {
    stop("'f' and 'u' must have equal length", call. = FALSE)
  }
  if (length(f) == 0) {
    stop("'f' must have at least one element", call. = FALSE)
  }

  results <- vapply(seq_along(f), function(i) {
    true_concentration_homogeneous(lambda, a, b, f[i], u[i], USL, n_sim)
  }, numeric(1))

  return(matrix(results, nrow = 1))
}

#' @rdname true_concentration_homogeneous
#' @export
true_concentration_curves_homogeneous <- function(lambda_low, lambda_high,
                                                  a, b, f, u, USL, n_sim) {
  # Input validation
  if (!is.numeric(lambda_low) || length(lambda_low) != 1 || lambda_low < 0) {
    stop("'lambda_low' must be a non-negative numeric scalar", call. = FALSE)
  }
  if (!is.numeric(lambda_high) || length(lambda_high) != 1 ||
      lambda_high <= lambda_low) {
    stop("'lambda_high' must be greater than 'lambda_low'", call. = FALSE)
  }
  if (length(f) != length(u)) {
    stop("'f' and 'u' must have equal length", call. = FALSE)
  }

  # Generate lambda sequence
  lambda_seq <- seq(lambda_low, lambda_high, by = 0.1)

  # Calculate concentration for each lambda and scheme
  c_matrix <- matrix(NA, nrow = length(lambda_seq), ncol = length(f))

  # Show progress for long calculations
  if (length(lambda_seq) > 100) {
    pb <- utils::txtProgressBar(min = 0, max = length(lambda_seq), style = 3)
  }

  for (i in seq_along(lambda_seq)) {
    c_matrix[i, ] <- true_concentration_homogeneous_multiple(
      lambda_seq[i], a, b, f, u, USL, n_sim
    )

    if (exists("pb") && is(pb, "txtProgressBar")) {
      utils::setTxtProgressBar(pb, i)
    }
  }

  if (exists("pb") && is(pb, "txtProgressBar")) {
    close(pb)
  }

  # Prepare data for plotting
  f_spr <- function(f, u) {
    sprintf("Scheme (f=%.3f, u=%.1f)", f, u)
  }

  Prob <- data.frame(lambda = lambda_seq, c_matrix)
  colnames(Prob) <- c("lambda", f_spr(f, u))

  melten.Prob <- reshape2::melt(
    Prob,
    id = "lambda",
    variable.name = "Dilution_scheme",
    value.name = "C"
  )

  # Create plot
  plot_sam <- ggplot2::ggplot(melten.Prob) +
    ggplot2::geom_line(ggplot2::aes(
      x = lambda,
      y = C,
      group = Dilution_scheme,
      colour = Dilution_scheme
    ), size = 1.2) +
    ggplot2::theme_classic() +
    ggplot2::xlab(expression("Expected microbial count (" ~ lambda * ~ ")")) +
    ggplot2::ylab(expression("Estimated true concentration (" ~ widehat(C) ~ ")")) +
    ggthemes::scale_colour_colorblind() +
    ggplot2::geom_vline(xintercept = USL, linetype = "dashed",
                        color = "gray40", size = 0.8) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(hjust = 0.5, size = 14, face = "bold"),
      legend.position = c(0.85, 0.25),
      legend.title = ggplot2::element_text(size = 10, face = "bold"),
      legend.background = ggplot2::element_rect(fill = "white",
                                                color = "gray80",
                                                linetype = "solid"),
      axis.line.x.top = ggplot2::element_line(color = "red"),
      axis.ticks.x.top = ggplot2::element_line(color = "red"),
      axis.text.x.top = ggplot2::element_text(color = "red", size = 8),
      axis.title.x.top = ggplot2::element_text(color = "red", size = 9)
    )

  return(plot_sam)
}
