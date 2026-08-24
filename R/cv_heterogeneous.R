#' Coefficient of Variation for Heterogeneous Batches
#'
#' Estimates the coefficient of variation (CV) in the original sample when
#' diluted samples are collected from a heterogeneous (non-homogeneous) batch.
#' The CV measures the relative variability of the estimated microbial concentration.
#'
#' @param mu Mean microbial count on the log scale. Must satisfy
#'   \code{log(a) - log(f*u) <= mu <= log(b) - log(f*u)}.
#' @param sd Standard deviation of the normal distribution on the log scale.
#'   Must be positive.
#' @param a Lower bound of cell count domain. Must be non-negative.
#' @param b Upper bound of cell count domain. Must be greater than \code{a}.
#' @param f Final dilution factor. Must be between 0 and 1.
#' @param u Amount placed on the plate. Must be positive.
#' @param USL Upper specification limit for microbial count.
#' @param n_sim Number of simulations. Larger values provide more precise
#'   estimates. Recommended minimum: 10,000.
#'
#' @return Numeric value representing the coefficient of variation for the
#'   heterogeneous batch.
#'
#' @details
#' The coefficient of variation is defined as:
#' \deqn{CV = \frac{\sigma}{\mu}}
#' where \eqn{\sigma} is the standard deviation and \eqn{\mu} is the mean
#' of the estimated microbial concentration.
#'
#' For heterogeneous batches, the microbial count follows a truncated
#' Poisson-lognormal distribution, accounting for variability in contamination
#' levels across the batch.
#'
#' @references
#' \itemize{
#'   \item Gonzales-Barron, U.A., Pilão Cadavez, V.A., Butler, F., 2013.
#'     Statistical approaches for the design of sampling plans for
#'     microbiological monitoring of foods. In: Mathematical and Statistical
#'     Methods in Food Science and Technology. Wiley, pp. 363–384.
#'   \item Schothorst, M. van, Zwietering, M.H., Ross, T., Buchanan, R.L.,
#'     Cole, M.B., 2009. Relating microbiological criteria to food safety
#'     objectives and performance objectives. Food Control 20, 967–979.
#' }
#'
#' @examples
#' \dontrun{
#' # Basic usage - estimate CV for a single dilution scheme
#' cv_heterogeneous(mu = 2, sd = 0.2, a = 0, b = 300,
#'                  f = 0.01, u = 0.1, USL = 1000, n_sim = 10000)
#'
#' # Compare multiple dilution schemes
#' cv_heterogeneous_multiple(mu = 2, sd = 0.2, a = 0, b = 300,
#'                           f = c(0.01, 0.1), u = c(0.1, 0.1),
#'                           USL = 1000, n_sim = 10000)
#'
#' # Generate CV curves across a range of mean microbial counts
#' cv_curves_heterogeneous(mu_low = -5, mu_high = 10, sd = 0.2,
#'                         a = 0, b = 300, f = c(0.01, 0.1),
#'                         u = c(0.1, 0.1), USL = 1000, n_sim = 10000)
#' }
#'
#' # Quick example with small simulation size for testing
#' cv_heterogeneous(mu = 2, sd = 0.2, a = 0, b = 300,
#'                  f = 0.01, u = 0.1, USL = 1000, n_sim = 500)
#'
#' @export
cv_heterogeneous <- function(mu, sd, a, b, f, u, USL, n_sim) {
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
  n_rep <- max(1, round(1 / (f * u)))

  sim1 <- matrix(NA, nrow = n_sim, ncol = n_rep)

  for (j in seq_len(n_sim)) {
    sim1[j, ] <- rtrunpoilog(n_rep, (mu + log_fu), sd, a, b) * (1 / (f * u))
  }

  sim2 <- apply(sim1, 2, mean, na.rm = TRUE)
  cv <- sqrt(var(sim2, na.rm = TRUE)) / mean(sim2, na.rm = TRUE)

  return(cv)
}

#' @rdname cv_heterogeneous
#' @export
cv_heterogeneous_multiple <- function(mu, sd, a, b, f, u, USL, n_sim) {
  if (length(f) != length(u)) {
    stop("'f' and 'u' must have equal length", call. = FALSE)
  }
  if (length(f) == 0) {
    stop("'f' must have at least one element", call. = FALSE)
  }

  results <- vapply(seq_along(f), function(i) {
    cv_heterogeneous(mu, sd, a, b, f[i], u[i], USL, n_sim)
  }, numeric(1))

  return(matrix(results, nrow = 1))
}

#' @rdname cv_heterogeneous
#' @export
cv_curves_heterogeneous <- function(mu_low, mu_high, sd, a, b, f, u, USL, n_sim) {
  if (!is.numeric(mu_low) || length(mu_low) != 1) {
    stop("'mu_low' must be a numeric scalar", call. = FALSE)
  }
  if (!is.numeric(mu_high) || length(mu_high) != 1 || mu_high <= mu_low) {
    stop("'mu_high' must be greater than 'mu_low'", call. = FALSE)
  }
  if (length(f) != length(u)) {
    stop("'f' and 'u' must have equal length", call. = FALSE)
  }

  # Add validation for f and u
  if (any(f <= 0 | f >= 1)) {
    stop("All 'f' values must be between 0 and 1", call. = FALSE)
  }
  if (any(u <= 0)) {
    stop("All 'u' values must be positive", call. = FALSE)
  }

  mu_seq <- seq(mu_low, mu_high, by = 0.1)
  cv_matrix <- matrix(NA, nrow = length(mu_seq), ncol = length(f))

  # Show progress for long calculations
  if (length(mu_seq) > 100) {
    pb <- utils::txtProgressBar(min = 0, max = length(mu_seq), style = 3)
  }

  for (i in seq_along(mu_seq)) {
    cv_matrix[i, ] <- cv_heterogeneous_multiple(
      mu_seq[i], sd, a, b, f, u, USL, n_sim
    )

    if (exists("pb") && is(pb, "txtProgressBar")) {
      utils::setTxtProgressBar(pb, i)
    }
  }

  if (exists("pb") && is(pb, "txtProgressBar")) {
    close(pb)
  }

  f_spr <- function(f, u) {
    sprintf("Scheme (f=%.3f, u=%.1f)", f, u)
  }

  Prob <- data.frame(mu = mu_seq, cv_matrix)
  colnames(Prob) <- c("mu", f_spr(f, u))

  melten.Prob <- reshape2::melt(
    Prob,
    id = "mu",
    variable.name = "Dilution_scheme",
    value.name = "cv"
  )

  # Add USL vertical line
  log_USL <- log(USL, exp(1))

  plot_sam <- ggplot2::ggplot(melten.Prob) +
    ggplot2::geom_line(ggplot2::aes(
      x = mu,
      y = cv,
      group = Dilution_scheme,
      colour = Dilution_scheme
    ), size = 1.2) +
    ggplot2::theme_classic() +
    ggplot2::xlab(expression("Log mean microbial count (" ~ mu * ~ ")")) +
    ggplot2::ylab(expression("Coefficient of variation (CV)")) +
    ggthemes::scale_colour_colorblind() +
    ggplot2::geom_vline(xintercept = log_USL, linetype = "dashed",
                        color = "gray40", size = 0.8) +
    ggplot2::annotate("text", x = log_USL, y = 0.05,
                      label = sprintf("log(USL) = %0.4f", log_USL),
                      size = 3.5, hjust = -0.1) +
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

#' Helper function to understand CV behavior
#'
#' Compares CV for different dilution schemes and parameters.
#' Useful for sensitivity analysis.
#'
#' @param mu_seq Vector of mu values to evaluate.
#' @param sd_values Vector of sd values to evaluate.
#' @param a Lower bound.
#' @param b Upper bound.
#' @param f Dilution factor.
#' @param u Plate amount.
#' @param USL Upper specification limit.
#' @param n_sim Number of simulations.
#' @param ... Additional arguments passed to cv_heterogeneous.
#'
#' @return A data frame with CV values.
#' @examples
#' \dontrun{
#' # Compare CV across different mu and sd values
#' cv_sensitivity <- compare_cv_heterogeneous(
#'   mu_seq = seq(0, 5, by = 0.5),
#'   sd_values = c(0.1, 0.3, 0.5),
#'   a = 0, b = 300,
#'   f = 0.01, u = 0.1,
#'   USL = 1000, n_sim = 10000
#' )
#' print(cv_sensitivity)
#' }
#' @export
compare_cv_heterogeneous <- function(mu_seq, sd_values, a, b, f, u, USL, n_sim = 10000, ...) {
  result <- expand.grid(mu = mu_seq, sd = sd_values)
  result$cv <- NA

  for (i in seq_len(nrow(result))) {
    result$cv[i] <- cv_heterogeneous(
      mu = result$mu[i],
      sd = result$sd[i],
      a = a, b = b,
      f = f, u = u,
      USL = USL,
      n_sim = n_sim,
      ...
    )
  }

  return(result)
}

#' Visualize CV Comparison
#'
#' Creates a heatmap or contour plot of CV values across mu and sd parameters.
#'
#' @param cv_data Data frame from compare_cv_heterogeneous.
#' @param type Type of plot: "heatmap" or "contour".
#'
#' @return A ggplot object.
#' @examples
#' \dontrun{
#' cv_data <- compare_cv_heterogeneous(
#'   mu_seq = seq(0, 5, by = 0.25),
#'   sd_values = seq(0.1, 0.8, by = 0.05),
#'   a = 0, b = 300,
#'   f = 0.01, u = 0.1,
#'   USL = 1000, n_sim = 5000
#' )
#' plot_cv_comparison(cv_data, type = "heatmap")
#' }
#' @export
plot_cv_comparison <- function(cv_data, type = "heatmap") {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("ggplot2 package is required for this function", call. = FALSE)
  }

  if (type == "heatmap") {
    p <- ggplot2::ggplot(cv_data, ggplot2::aes(x = mu, y = sd, fill = cv)) +
      ggplot2::geom_tile() +
      ggplot2::scale_fill_gradient2(
        low = "blue", mid = "white", high = "red",
        midpoint = mean(cv_data$cv, na.rm = TRUE),
        name = "CV"
      ) +
      ggplot2::theme_classic() +
      ggplot2::xlab(expression("Log mean microbial count (" ~ mu * ~ ")")) +
      ggplot2::ylab(expression("Standard deviation (" ~ sigma * ~ ")")) +
      ggplot2::theme(
        plot.title = ggplot2::element_text(hjust = 0.5, size = 14, face = "bold"),
        legend.position = "right"
      )
  } else {
    # Contour plot
    p <- ggplot2::ggplot(cv_data, ggplot2::aes(x = mu, y = sd, z = cv)) +
      ggplot2::geom_contour_filled() +
      ggplot2::theme_classic() +
      ggplot2::xlab(expression("Log mean microbial count (" ~ mu * ~ ")")) +
      ggplot2::ylab(expression("Standard deviation (" ~ sigma * ~ ")")) +
      ggplot2::theme(
        plot.title = ggplot2::element_text(hjust = 0.5, size = 14, face = "bold"),
        legend.position = "right"
      )
  }

  return(p)
}
