#' True Concentration Estimation for Heterogeneous Batches
#'
#' Estimates the true microbial concentration in the original sample when
#' diluted samples are collected from a heterogeneous (non-homogeneous) batch.
#'
#' @param mu Mean microbial count on the log scale. Must satisfy
#'   \code{log(a) - log(f*u) <= mu <= log(b) - log(f*u)}.
#' @param mu_low Lower bound of mean microbial count for x-axis (log scale).
#' @param mu_high Upper bound of mean microbial count for x-axis (log scale).
#' @param sd Standard deviation of the normal distribution on the log scale.
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
#' Poisson-lognormal distribution, f is the final dilution factor, and u is
#' the amount placed on the plate.
#'
#' @references
#' \itemize{
#'   \item Gonzales-Barron, U.A., Pilão Cadavez, V.A., Butler, F., 2013.
#'     Statistical approaches for the design of sampling plans for
#'     microbiological monitoring of foods. In: Mathematical and Statistical
#'     Methods in Food Science and Technology. Wiley, pp. 363–384.
#' }
#'
#' @examples
#' \dontrun{
#' # Basic usage
#' true_concentration_heterogeneous(mu = 2, sd = 0.2, a = 0, b = 300,
#'                                  f = 0.01, u = 0.1, USL = 1000, n_sim = 50000)
#'
#' # Multiple dilution schemes
#' true_concentration_heterogeneous_multiple(mu = 2, sd = 0.2, a = 0, b = 300,
#'                                           f = c(0.01, 0.1), u = c(0.1, 0.1),
#'                                           USL = 1000, n_sim = 50000)
#'
#' # Plot concentration curves
#' true_concentration_curves_heterogeneous(mu_low = 0, mu_high = 10, sd = 0.2,
#'                                         a = 0, b = 300, f = c(0.01, 0.1),
#'                                         u = c(0.1, 0.1), USL = 1000,
#'                                         n_sim = 50000)
#' }
#'
#' @name true_concentration_heterogeneous
#' @aliases true_concentration_heterogeneous true_concentration_heterogeneous_multiple
#' @rdname true_concentration_heterogeneous
#' @export
true_concentration_heterogeneous <- function(mu, sd, a, b, f, u, USL, n_sim) {
  # Input validation
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
  if (!is.numeric(n_sim) || length(n_sim) != 1 || n_sim < 100) {
    warning("'n_sim' should be at least 100 for reliable results", call. = FALSE)
  }

  # Check bounds
  fu <- f * u
  log_fu <- log(fu, exp(1))
  log_a <- log(a, exp(1))
  log_b <- log(b, exp(1))

  if (mu > (log_b - log_fu) || mu < (log_a - log_fu)) {
    stop("The truncated Poisson lognormal (TPLN) random variable must be bounded by a and b,\n",
         "which means that mu must be between log(a)-log(fu) and log(b)-log(fu).\n",
         sprintf("Current mu = %.4f, bounds: [%.4f, %.4f]",
                 mu, log_a - log_fu, log_b - log_fu),
         call. = FALSE)
  }

  # Simulation
  sim1 <- matrix(NA, nrow = n_sim, ncol = 1)

  for (j in seq_len(n_sim)) {
    sim1[j, ] <- rtrunpoilog(1, (mu + log(fu, exp(1))), sd, a, b) *
      (1 / (f * u))
  }

  mean_concentr <- mean(sim1)
  return(mean_concentr)
}

#' @rdname true_concentration_heterogeneous
#' @export
true_concentration_heterogeneous_multiple <- function(mu, sd, a, b, f, u, USL, n_sim) {
  # Input validation
  if (length(f) != length(u)) {
    stop("'f' and 'u' must have equal length", call. = FALSE)
  }
  if (length(f) == 0) {
    stop("'f' must have at least one element", call. = FALSE)
  }

  results <- vapply(seq_along(f), function(i) {
    true_concentration_heterogeneous(mu, sd, a, b, f[i], u[i], USL, n_sim)
  }, numeric(1))

  return(matrix(results, nrow = 1))
}

#' @rdname true_concentration_heterogeneous
#' @export
true_concentration_curves_heterogeneous <- function(mu_low, mu_high, sd,
                                                    a, b, f, u, USL, n_sim) {
  # Input validation
  if (!is.numeric(mu_low) || length(mu_low) != 1) {
    stop("'mu_low' must be a numeric scalar", call. = FALSE)
  }
  if (!is.numeric(mu_high) || length(mu_high) != 1 || mu_high <= mu_low) {
    stop("'mu_high' must be greater than 'mu_low'", call. = FALSE)
  }
  if (!is.numeric(sd) || length(sd) != 1 || sd <= 0) {
    stop("'sd' must be a positive numeric scalar", call. = FALSE)
  }
  if (length(f) != length(u)) {
    stop("'f' and 'u' must have equal length", call. = FALSE)
  }

  # Generate mu sequence
  mu_seq <- seq(mu_low, mu_high, by = 0.1)

  # Calculate concentration for each mu and scheme
  c_matrix <- matrix(NA, nrow = length(mu_seq), ncol = length(f))

  # Show progress for long calculations
  if (length(mu_seq) > 100) {
    pb <- utils::txtProgressBar(min = 0, max = length(mu_seq), style = 3)
  }

  for (i in seq_along(mu_seq)) {
    c_matrix[i, ] <- true_concentration_heterogeneous_multiple(
      mu_seq[i], sd, a, b, f, u, USL, n_sim
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

  Prob <- data.frame(mu = mu_seq, c_matrix)
  colnames(Prob) <- c("mu", f_spr(f, u))

  melten.Prob <- reshape2::melt(
    Prob,
    id = "mu",
    variable.name = "Dilution_scheme",
    value.name = "C"
  )

  # Create plot
  log_USL <- log(USL, exp(1))

  plot_sam <- ggplot2::ggplot(melten.Prob) +
    ggplot2::geom_line(ggplot2::aes(
      x = mu,
      y = C,
      group = Dilution_scheme,
      colour = Dilution_scheme
    ), size = 1.2) +
    ggplot2::theme_classic() +
    ggplot2::xlab(expression("Log mean microbial count (" ~ mu * ~ ")")) +
    ggplot2::ylab(expression("Estimated true concentration (" ~ widehat(C) ~ ")")) +
    ggthemes::scale_colour_colorblind() +
    ggplot2::geom_vline(xintercept = log_USL, linetype = "dashed",
                        color = "gray40", size = 0.8) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(hjust = 0.5, size = 14, face = "bold"),
      legend.position = c(0.25, 0.85),
      legend.title = ggplot2::element_text(size = 10, face = "bold"),
      legend.background = ggplot2::element_rect(fill = "white",
                                                color = "gray80",
                                                linetype = "solid"),
      axis.line.x.top = ggplot2::element_line(color = "red"),
      axis.ticks.x.top = ggplot2::element_line(color = "red"),
      axis.text.x.top = ggplot2::element_text(color = "red", size = 8),
      axis.title.x.top = ggplot2::element_text(color = "red", size = 9)
    ) +
    ggplot2::scale_x_continuous(
      sec.axis = ggplot2::sec_axis(
        ~.,
        name = "Mean microbial count",
        breaks = seq(floor(mu_low), ceiling(mu_high), by = 1),
        labels = sprintf("%0.2f", exp(seq(floor(mu_low), ceiling(mu_high), by = 1)))
      )
    )

  return(plot_sam)
}
