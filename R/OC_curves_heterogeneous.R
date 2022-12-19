##' \code{\link{OC_curves_heterogeneous}} provides the operating characteristic(OC) curves when samples collected from a heterogeneous batch.
##' @title Comparison based on OC curves for different dilution schemes when the diluted samples collected from a heterogeneous batch.
##' @param c acceptance number
##' @param mu_low the lower value of the mean microbial count(\eqn{\mu}) for use in the graphical display's x-axis.
##' @param mu_high the upper value of the mean microbial count(\eqn{\mu}) for use in the graphical display's x-axis.
##' @param sd the standard deviation of the normal distribution (on the log scale).
##' @param a lower domain of the number of cell counts.
##' @param b upper domain of the number of cell counts.
##' @param f final dilution factor.
##' @param u amount put on the plate.
##' @param USL upper specification limit.
##' @param n number of samples which are used for inspection.
##' @param type what type of the results you would like to consider such as "theory" or "simulation" (default "theory").
##' @param n_sim number of simulations (large simulations provide more precise estimations).
##' @details \code{\link{OC_curves_heterogeneous}} provides OC curves for different dilution schemes when the diluted samples collected from a heterogeneous batch (this section will be updated later on).
##' @return OC curves when samples collected from a heterogeneous batch.
##' @examples
##' c <- 2
##' mu_low <- 4
##' mu_high <- 9
##' sd <- 0.2
##' a <- 0
##' b <- 300
##' f <- c(0.01,0.1)
##' u <- c(0.1,0.1)
##' USL <- 1000
##' n <- 5
##' OC_curves_heterogeneous(c, mu_low, mu_high, sd, a, b, f, u, USL, n)
##' @usage  OC_curves_heterogeneous(c, mu_low, mu_high, sd, a, b, f, u, USL, n, type, n_sim)
##' @export
OC_curves_heterogeneous <- function(c, mu_low, mu_high, sd, a, b, f, u, USL, n, type = "theory", n_sim = NA){
  P_a <- NULL
  Dilution_scheme <- NULL
  mu <- seq(mu_low, mu_high, 0.001)
  f_spr_1 <- function(f, u) {
    sprintf("Scheme (f=%.4f, u=%.1f)", f, u)
  }

  pa <- matrix(NA, nrow = length(mu), ncol = length(f))
  for (i in 1:length(mu)) {
    pa[i,] <-  cbind(prob_acceptance_heterogeneous_multiple(c, mu[i], sd, a, b, f, u, USL, n, type, n_sim))
  }
  Prob <- data.frame(mu, pa)
  colnames(Prob ) <- c("mu", f_spr_1(f,u))
  melten.Prob <- reshape2::melt(Prob, id = "mu", variable.name = "Dilution_scheme", value.name = "P_a")
  plot_sam <- ggplot2::ggplot(melten.Prob) + ggplot2::geom_line(ggplot2::aes(x = mu, y = P_a, group = Dilution_scheme, colour = Dilution_scheme)) +
    # ggplot2::ggtitle("OC curve based on Poisson Lognormal distribution") +
    ggplot2::theme_classic() + ggplot2::xlab(expression("log mean microbial count  (" ~ mu*~")")) + ggplot2::ylab(expression("Probability of acceptance"~(P[a]))) + ggthemes::scale_colour_colorblind() +
    ggplot2::geom_vline(xintercept = log(USL,exp(1)), linetype = "dashed") +
    ggplot2::annotate("text", x = log(USL,exp(1)),
                      y = 0, label = sprintf("log(USL) = %0.4f", log(USL,exp(1))), size = 3) +
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5), legend.position = c(0.85, 0.75), axis.line.x.top = ggplot2::element_line(color = "red"),
                   axis.ticks.x.top = ggplot2::element_line(color = "red"), axis.text.x.top = ggplot2::element_text(color = "red"), axis.title.x.top = ggplot2::element_text(color = "red")) +
    ggplot2::scale_x_continuous(sec.axis = ggplot2::sec_axis(~., name = "mean microbial count ", breaks = seq(min(mu),max(mu),1),
                                                             labels = c(sprintf("%0.2f", exp(seq(min(mu),max(mu),1))))))
  # ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5), legend.position = c(0.25, 0.85), axis.line.x.top = ggplot2::element_line(color = "red"),
  #                axis.ticks.x.top = ggplot2::element_line(color = "red"), axis.text.x.top = ggplot2::element_text(color = "red"), axis.title.x.top = ggplot2::element_text(color = "red")) +
  #   ggplot2::scale_x_continuous(sec.axis = ggplot2::sec_axis(~., name = "expected cell counts", breaks = seq(min(mu), max(mu),1),
  #                                                            labels = c(sprintf("%.4f", 10^(seq(min(mu),max(mu),1) + (sd^2/2) * log(10, exp(1)))))))
  # plot_sam
  return(plot_sam)
}
