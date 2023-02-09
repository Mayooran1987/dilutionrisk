##' \code{\link{pd_curves_heterogeneous}} provides the probability of detection curves when samples collected from a heterogeneous batch.
##' @title Comparison based on probability of detection curves for different dilution schemes when the diluted samples collected from a heterogeneous batch.
##' @param mu_low the lower value of the mean microbial count (\eqn{\mu}) for use in the graphical display's x-axis (on the log scale).
##' @param mu_high the upper value of the mean microbial count (\eqn{\mu}) for use in the graphical display's x-axis (on the log scale).
##' @param sd the standard deviation of the normal distribution (on the log scale).
##' @param a lower domain of the number of cell counts.
##' @param b upper domain of the number of cell counts.
##' @param f final dilution factor.
##' @param u amount put on the plate.
##' @param USL upper specification limit.
##' @param type what type of the results you would like to consider such as "theory" or "simulation" (default "theory").
##' @param n_sim number of simulations (large simulations provide more precise estimations).
##' @details \code{\link{pd_curves_heterogeneous}} provides probability of detection curves for different dilution schemes when the diluted samples collected from a heterogeneous batch (this section will be updated later on).
##' @return Probability of detection curves when samples collected from a heterogeneous batch.
##' @examples
##' mu_low <- 0
##' mu_high <- 10
##' sd <- 0.2
##' a <- 0
##' b <- 300
##' f <- c(0.01,0.1)
##' u <- c(0.1,0.1)
##' USL <- 1000
##' pd_curves_heterogeneous(mu_low, mu_high, sd, a, b, f, u, USL)
##' @usage  pd_curves_heterogeneous(mu_low, mu_high, sd, a, b, f, u, USL, type, n_sim)
##' @export
pd_curves_heterogeneous <- function(mu_low, mu_high, sd = 0.8, a, b, f, u, USL, type = "theory", n_sim = NA) {
  p_d <- NULL
  Dilution_scheme <- NULL
  f_spr <- function(f, u) {
    sprintf("Scheme (f=%.4f, u=%.1f)", f, u)
  }
  mu <- seq(mu_low, mu_high, 0.1)
  # lambda <- 10^(mu + (sd^2/2) * log(10, exp(1)))
  Pd <- matrix(NA, nrow = length(mu), ncol = length(f))
  for (i in 1:length(mu)) {
    Pd[i, ] <- cbind(prob_detection_heterogeneous_multiple(mu[i], sd, a, b, f, u, USL, type, n_sim))
  }
  Prob <- data.frame(mu, Pd)
  colnames(Prob) <- c("mu", f_spr(f, u))
  melten.Prob <- reshape2::melt(Prob, id = "mu", variable.name = "Dilution_scheme", value.name = "p_d")
  plot_sam <- ggplot2::ggplot(melten.Prob) +
    ggplot2::geom_line(ggplot2::aes(x = mu, y = p_d, group = Dilution_scheme, colour = Dilution_scheme)) +
    # ggplot2::ggtitle("OC curve based on Lognormal distribution") +
    ggplot2::theme_classic() +
    ggplot2::xlab(expression(" log mean microbial count  (" ~ mu * ~")")) +
    ggplot2::ylab(expression("Probability of detection" ~ (P[d]))) +
    ggthemes::scale_colour_colorblind() +
    ggplot2::geom_vline(xintercept = log(USL, exp(1)), linetype = "dashed") +
    ggplot2::annotate("text",
      x = log(USL, exp(1)),
      y = 0, label = sprintf("log(USL) = %0.4f", log(USL, exp(1))), size = 3
    ) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(hjust = 0.5), legend.position = c(0.25, 0.85), axis.line.x.top = ggplot2::element_line(color = "red"),
      axis.ticks.x.top = ggplot2::element_line(color = "red"), axis.text.x.top = ggplot2::element_text(color = "red"), axis.title.x.top = ggplot2::element_text(color = "red")
    ) +
    ggplot2::scale_x_continuous(sec.axis = ggplot2::sec_axis(~.,
      name = "mean microbial count", breaks = seq(min(mu), max(mu), 1),
      labels = c(sprintf("%0.2f", exp(seq(min(mu), max(mu), 1))))
    ))
  # ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5), legend.position = c(0.25, 0.85), axis.line.x.top = ggplot2::element_line(color = "red"),
  #                axis.ticks.x.top = ggplot2::element_line(color = "red"), axis.text.x.top = ggplot2::element_text(color = "red"), axis.title.x.top = ggplot2::element_text(color = "red")) +
  # ggplot2::scale_x_continuous(sec.axis = ggplot2::sec_axis(~., name = "expected cell counts", breaks = seq(min(mu), max(mu),1),
  #                                                          labels = c(sprintf("%.4f", 10^(seq(min(mu),max(mu),1) + (sd^2/2) * log(10, exp(1)))))))
  # plot_sam
  return(plot_sam)
}
