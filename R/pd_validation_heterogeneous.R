##' \code{\link{pd_validation_heterogeneous}} provides the probability of detection curves for validate the results when samples collected from a heterogeneous batch.
##' @title Comparison based on probability of detection curves for different dilution schemes when diluted samples collected from a heterogeneous batch.
##' @param mu_low the lower value of the mean microbial count (\eqn{\mu}) for use in the graphical display's x-axis (on the log scale).
##' @param mu_high the upper value of the mean microbial count (\eqn{\mu}) for use in the graphical display's x-axis (on the log scale).
##' @param sd the standard deviation of the normal distribution (on the log scale).
##' @param a lower domain of the number of cell counts.
##' @param b upper domain of the number of cell counts.
##' @param f final dilution factor.
##' @param u amount put on the plate.
##' @param USL upper specification limit.
##' @param n_sim number of simulations (large simulations provide more precise estimations).
##' @details \code{\link{pd_curves_heterogeneous}} provides probability of detection curves for different dilution schemes when samples collected from a heterogeneous batch (this section will be updated later on).
##' @return Probability of detection curves when diluted samples collected from a heterogeneous batch.
##' @examples
##' mu_low <- 1
##' mu_high <- 12
##' sd <- 0.2
##' a <- 0
##' b <- 300
##' f <- 0.01
##' u <- 0.1
##' USL <- 1000
##' n_sim <- 50000
##' pd_validation_heterogeneous(mu_low,  mu_high, sd, a, b, f, u, USL, n_sim)
##' @usage  pd_validation_heterogeneous(mu_low,  mu_high, sd, a, b, f, u, USL, n_sim)
##' @export
pd_validation_heterogeneous <- function(mu_low, mu_high, sd, a, b, f, u, USL, n_sim) {
  p_d <- NULL
  Methods <- NULL
  mu <- seq(mu_low, mu_high, 0.1)
  # pd_theory_heterogeneous <- function(mu, sd, f, u, USL){
  #   USL1 <- USL*f*u
  #   pd <- matrix(NA, nrow = USL1, ncol = 1)
  #   for (i in 1:USL1) {
  #     mu_d <- mu + log(f*u, exp(1))
  #     # pd[i,1] <- integrate(function(t) exp(t*i - 0.5 * ((t - mu_d)/sd)^2)/(exp(exp(t))-1), lower = 0, upper = 500)$value/(factorial(i)*sqrt(2 * pi) * sd)
  #     pd[i,1] <- (integrate(function(t) t^(i - 1)*exp( -0.5 * ((log(t,exp(1)) - mu_d)/sd)^2)/(exp(t) - 1), lower = 0, upper = Inf )$value)/(exp(lgamma(i + 1))*sqrt(2 * pi) * sd)
  #
  #   }
  #   result <- 1 - sum(pd)
  #   return(result)
  # }
  # pd_theory(lambda = 1, f, u, USL)
  Pd <- matrix(NA, nrow = length(mu), ncol = 2)
  for (i in 1:length(mu)) {
    Pd[i, 1] <- prob_detection_heterogeneous(mu[i], sd, a, b, f, u, USL, type = "simulation", n_sim)
    Pd[i, 2] <- prob_detection_heterogeneous(mu[i], sd, a, b, f, u, USL, type = "theory")
  }
  Prob <- data.frame(mu, Pd)
  # colnames(Prob ) <- c("lambda", f_spr(f,u))
  colnames(Prob) <- c("mu", "simulation", "theory")
  melten.Prob <- reshape2::melt(Prob, id = "mu", variable.name = "Methods", value.name = "p_d")
  plot_sam <- ggplot2::ggplot(melten.Prob) +
    ggplot2::geom_line(ggplot2::aes(x = mu, y = p_d, group = Methods, colour = Methods)) +
    # ggplot2::ggtitle("OC curve based on Lognormal distribution") +
    ggplot2::theme_classic() +
    ggplot2::xlab(expression(" log mean microbial count  (" ~ mu * ~")")) +
    ggplot2::ylab(expression("Probability of detection" ~ (P[d]))) +
    ggthemes::scale_colour_colorblind() +
    ggplot2::theme(
      plot.title = ggplot2::element_text(hjust = 0.5), legend.position = c(0.25, 0.85), axis.line.x.top = ggplot2::element_line(color = "red"),
      axis.ticks.x.top = ggplot2::element_line(color = "red"), axis.text.x.top = ggplot2::element_text(color = "red"), axis.title.x.top = ggplot2::element_text(color = "red")
    ) +
    ggplot2::scale_x_continuous(sec.axis = ggplot2::sec_axis(~.,
      name = "mean microbial count", breaks = seq(min(mu), max(mu), 1),
      labels = c(sprintf("%0.2f", exp(seq(min(mu), max(mu), 1))))
    ))
  return(plot_sam)
}
