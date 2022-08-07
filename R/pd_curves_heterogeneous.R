##' \code{\link{pd_curves_heterogeneous}} provides the probability of detection curves when samples collected from a heterogeneous batch.
##' @title comparison based on probability of detection curves for different dilution schemes when the diluted samples collected from a heterogeneous batch.
##' @param meanlog_low the lower value of the mean concentration (\eqn{\mu}) for use in the graphical display's x-axis (on the log scale).
##' @param meanlog_high the upper value of the mean concentration (\eqn{\mu}) for use in the graphical display's x-axis (on the log scale).
##' @param sdlog the standard deviation of the normal distribution (on the log scale).
##' @param a lower domain of the number of cell counts.
##' @param b upper domain of the number of cell counts.
##' @param FDF final dilution factor.
##' @param u amount put on the plate.
##' @param USL upper specification limit.
##' @param n_sim number of simulations (large simulations provide more precise estimations).
##' @details \code{\link{pd_curves_heterogeneous}} provides probability of detection curves for different dilution schemes when the diluted samples collected from a heterogeneous batch (this section will be updated later on).
##' @return Probability of detection curves when samples collected from a heterogeneous batch.
##' @examples
##' meanlog_low <- 0
##' meanlog_high <- 10
##' sdlog <- 0.2
##' a <- 0
##' b <- 300
##' FDF <- c(0.01,0.1)
##' u <- c(0.1,0.1)
##' USL <- 1000
##' n_sim <- 50000
##' pd_curves_heterogeneous(meanlog_low, meanlog_high, sdlog, a, b, FDF, u, USL, n_sim)
##' @usage  pd_curves_heterogeneous(meanlog_low, meanlog_high, sdlog, a, b, FDF, u, USL, n_sim)
pd_curves_heterogeneous <- function(meanlog_low, meanlog_high, sdlog = 0.8, a, b, FDF, u, USL, n_sim){
  p_d <- NULL
  Dilution_scheme <- NULL
  f_spr <- function(FDF, u) {
    sprintf("Scheme (FDF=%.3f, u=%.1f)", FDF, u)
  }
  meanlog <- seq(meanlog_low, meanlog_high, 0.1)
  # lambda <- 10^(mu + (sd^2/2) * log(10, exp(1)))
  Pd <- matrix(NA, nrow = length(meanlog), ncol = length(FDF))
  for (i in 1:length(meanlog)) {
    Pd[i,] <-  cbind(prob_detection_heterogeneous_multiple(meanlog[i], sdlog, a, b, FDF, u, USL, n_sim))
  }
  Prob <- data.frame(meanlog, Pd)
  colnames(Prob ) <- c("meanlog", f_spr(FDF,u))
  melten.Prob <- reshape2::melt(Prob, id = "meanlog", variable.name = "Dilution_scheme", value.name = "p_d")
  plot_sam <- ggplot2::ggplot(melten.Prob) + ggplot2::geom_line(ggplot2::aes(x = meanlog, y = p_d, group = Dilution_scheme, colour = Dilution_scheme)) +
    # ggplot2::ggtitle("OC curve based on Lognormal distribution") +
    ggplot2::theme_classic() + ggplot2::xlab(expression(" log mean concentrations  (" ~ mu*~")")) + ggplot2::ylab(expression(p[d])) + ggthemes::scale_colour_colorblind() +
    ggplot2::geom_vline(xintercept = log(USL,exp(1)), linetype = "dashed") +
    ggplot2::annotate("text", x = log(USL,exp(1)),
                      y = 0, label = sprintf("log(USL) = %0.4f", log(USL,exp(1))), size = 3) +
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5), legend.position = c(0.25, 0.85), axis.line.x.top = ggplot2::element_line(color = "red"),
                   axis.ticks.x.top = ggplot2::element_line(color = "red"), axis.text.x.top = ggplot2::element_text(color = "red"), axis.title.x.top = ggplot2::element_text(color = "red")) +
    ggplot2::scale_x_continuous(sec.axis = ggplot2::sec_axis(~., name = "expected cell counts", breaks = seq(min(meanlog), max(meanlog),1),
                                                             labels = c(sprintf("%.4f", 10^(seq(min(meanlog),max(meanlog),1) + (sdlog^2/2) * log(10, exp(1)))))))
  # plot_sam
  return(plot_sam)
}


