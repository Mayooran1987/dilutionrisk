#' Operating Characteristic (OC) Curves for Heterogeneous Batches
#'
#' Generates operating characteristic curves comparing different dilution
#' schemes when samples are collected from a heterogeneous (non-homogeneous) batch.
#'
#' @param c Acceptance number (maximum number of defective units allowed).
#' @param mu_low Lower bound of mean microbial count for x-axis (log scale).
#' @param mu_high Upper bound of mean microbial count for x-axis (log scale).
#' @param sd Standard deviation of the normal distribution on the log scale.
#' @param a Lower bound of cell count domain. Must be non-negative.
#' @param b Upper bound of cell count domain. Must be greater than \code{a}.
#' @param f Vector of final dilution factors.
#' @param u Vector of amounts placed on the plate.
#' @param USL Upper specification limit for microbial count.
#' @param n Number of samples inspected. Must be positive integer.
#' @param type Type of calculation: "theory" (default) or "simulation".
#' @param n_sim Number of simulations. Required when \code{type = "simulation"}.
#'
#' @return A ggplot object showing OC curves for different dilution schemes.
#'
#' @examples
#' \dontrun{
#' OC_curves_heterogeneous(c = 2, mu_low = 4, mu_high = 9, sd = 0.2,
#'                         a = 0, b = 300, f = c(0.01, 0.1), u = c(0.1, 0.1),
#'                         USL = 1000, n = 5)
#' }
#'
#' @export
OC_curves_heterogeneous <- function(c, mu_low, mu_high, sd, a, b, f, u,
                                    USL, n, type = "theory", n_sim = NA) {
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
  if (length(f) == 0) {
    stop("'f' must have at least one element", call. = FALSE)
  }

  # Generate mu sequence
  mu_seq <- seq(mu_low, mu_high, by = 0.1)

  # Calculate probability of acceptance for each mu and scheme
  pa_matrix <- matrix(NA, nrow = length(mu_seq), ncol = length(f))

  # Show progress for long calculations
  if (length(mu_seq) > 100) {
    pb <- utils::txtProgressBar(min = 0, max = length(mu_seq), style = 3)
  }

  for (i in seq_along(mu_seq)) {
    pa_matrix[i, ] <- prob_acceptance_heterogeneous_multiple(
      c, mu_seq[i], sd, a, b, f, u, USL, n, type, n_sim
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
    sprintf("Scheme (f=%.4f, u=%.1f)", f, u)
  }

  Prob <- data.frame(mu = mu_seq, pa_matrix)
  colnames(Prob) <- c("mu", f_spr(f, u))

  melten.Prob <- reshape2::melt(
    Prob,
    id = "mu",
    variable.name = "Dilution_scheme",
    value.name = "P_a"
  )

  # Create plot
  log_USL <- log(USL, exp(1))

  plot_sam <- ggplot2::ggplot(melten.Prob) +
    ggplot2::geom_line(ggplot2::aes(
      x = mu,
      y = P_a,
      group = Dilution_scheme,
      colour = Dilution_scheme
    ), size = 1.2) +
    ggplot2::theme_classic() +
    ggplot2::xlab(expression("Log mean microbial count (" ~ mu * ~ ")")) +
    ggplot2::ylab(expression("Probability of acceptance (" ~ P[a] ~ ")")) +
    ggthemes::scale_colour_colorblind() +
    ggplot2::geom_vline(xintercept = log_USL, linetype = "dashed",
                        color = "gray40", size = 0.8) +
    ggplot2::annotate("text", x = log_USL, y = 0.05,
                      label = sprintf("log(USL) = %0.4f", log_USL),
                      size = 3.5, hjust = -0.1) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(hjust = 0.5, size = 14, face = "bold"),
      legend.position = c(0.85, 0.75),
      legend.title = ggplot2::element_text(size = 10, face = "bold"),
      legend.background = ggplot2::element_rect(fill = "white",
                                                color = "gray80",
                                                linetype = "solid"),
      legend.key.size = ggplot2::unit(0.8, "cm"),
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
