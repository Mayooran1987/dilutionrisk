##' \code{\link{OC_curves_heterogeneous}} provides the operating characteristic(OC) curves when samples collected from a heterogeneous batch.
##' @title Comparison based on OC curves for different dilution schemes when the diluted samples collected from a heterogeneous batch.
##' @param c acceptance number
##' @param meanlog_low the lower value of the mean concentration (\eqn{\mu}) for use in the graphical display's x-axis.
##' @param meanlog_high the upper value of the mean concentration (\eqn{\mu}) for use in the graphical display's x-axis.
##' @param sdlog the standard deviation of the normal distribution (on the log scale).
##' @param a lower domain of the number of cell counts.
##' @param b upper domain of the number of cell counts.
##' @param f final dilution factor.
##' @param u amount put on the plate.
##' @param USL upper specification limit.
##' @param n number of samples which are used for inspection.
##' @param n_sim number of simulations (large simulations provide more precise estimations).
##' @details \code{\link{OC_curves_heterogeneous}} provides OC curves for different dilution schemes when the diluted samples collected from a heterogeneous batch (this section will be updated later on).
##' @return OC curves when samples collected from a heterogeneous batch.
##' @examples
##' c <- 2
##' meanlog_low <- 4
##' meanlog_high <- 9
##' sdlog <- 0.2
##' a <- 0
##' b <- 300
##' f <- c(0.01,0.1)
##' u <- c(0.1,0.1)
##' USL <- 1000
##' n <- 5
##' n_sim <- 50000
##' OC_curves_heterogeneous(c, meanlog_low, meanlog_high, sdlog, a, b, f, USL, u, n, n_sim)
##' @usage  OC_curves_heterogeneous(c, meanlog_low, meanlog_high, sdlog, a, b, f, u, USL, n, n_sim)
##' @export
OC_curves_heterogeneous <- function(c, meanlog_low, meanlog_high, sdlog, a, b, f, u, USL, n, n_sim){
  P_a <- NULL
  Dilution_scheme <- NULL
  meanlog <- seq(meanlog_low, meanlog_high, 0.1)
  f_spr <- function(f, u) {
    sprintf("Scheme (f=%.3f, u=%.1f)", f, u)
  }
  pa <- matrix(NA, nrow = length(meanlog), ncol = length(f))
  for (i in 1:length(meanlog)) {
    pa[i,] <-  cbind(prob_acceptance_heterogeneous_multiple(c, meanlog[i], sdlog, a, b, f, u, USL, n, n_sim))
  }
  Prob <- data.frame(meanlog, pa)
  colnames(Prob ) <- c("meanlog", f_spr(f,u))
  melten.Prob <- reshape2::melt(Prob, id = "meanlog", variable.name = "Dilution_scheme", value.name = "P_a")
  plot_sam <- ggplot2::ggplot(melten.Prob) + ggplot2::geom_line(ggplot2::aes(x = meanlog, y = P_a, group = Dilution_scheme, colour = Dilution_scheme)) +
    # ggplot2::ggtitle("OC curve based on Poisson Lognormal distribution") +
    ggplot2::theme_classic() + ggplot2::xlab(expression("log mean concentrations  (" ~ mu*~")")) + ggplot2::ylab(expression(P[a])) + ggthemes::scale_colour_colorblind() +
    ggplot2::geom_vline(xintercept = log(USL,exp(1)), linetype = "dashed") +
    ggplot2::annotate("text", x = log(USL,exp(1)),
                      y = 0, label = sprintf("log(USL) = %0.4f", log(USL,exp(1))), size = 3) +
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5), legend.position = c(0.85, 0.75), axis.line.x.top = ggplot2::element_line(color = "red"),
                   axis.ticks.x.top = ggplot2::element_line(color = "red"), axis.text.x.top = ggplot2::element_text(color = "red"), axis.title.x.top = ggplot2::element_text(color = "red")) +
    ggplot2::scale_x_continuous(sec.axis = ggplot2::sec_axis(~., name = "mean concentrations", breaks = seq(min(meanlog),max(meanlog),1),
                                                             labels = c(sprintf("%f", exp(seq(min(meanlog),max(meanlog),1))))))
  # ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5), legend.position = c(0.25, 0.85), axis.line.x.top = ggplot2::element_line(color = "red"),
  #                axis.ticks.x.top = ggplot2::element_line(color = "red"), axis.text.x.top = ggplot2::element_text(color = "red"), axis.title.x.top = ggplot2::element_text(color = "red")) +
  #   ggplot2::scale_x_continuous(sec.axis = ggplot2::sec_axis(~., name = "expected cell counts", breaks = seq(min(meanlog), max(meanlog),1),
  #                                                            labels = c(sprintf("%.4f", 10^(seq(min(meanlog),max(meanlog),1) + (sdlog^2/2) * log(10, exp(1)))))))
  # plot_sam
  return(plot_sam)
}
