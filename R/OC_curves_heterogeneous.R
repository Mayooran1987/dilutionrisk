##' \code{\link{OC_curves_heterogeneous}} provides the operating characteristic(OC) curves when diluted sample has heterogeneous contaminants.
##' @title Comparison based on OC curves for different dilution schemes when the diluted sample has heterogeneous contaminants.
##' @param c acceptance number
##' @param meanlog_low the lower value of the mean concentration (\eqn{\mu}) for use in the graphical display's x-axis.
##' @param meanlog_high the upper value of the mean concentration (\eqn{\mu}) for use in the graphical display's x-axis.
##' @param sdlog the standard deviation of the normal distribution (on the log scale).
##' @param a lower domain of the number of cell counts.
##' @param b upper domain of the number of cell counts.
##' @param FDF final dilution factor.
##' @param USL upper specification limit.
##' @param n number of samples which are used for inspection.
##' @param n_sim number of simulations (large simulations provide more precise estimation).
##' @details \code{\link{OC_curves_heterogeneous}} provides OC curves for different dilution schemes when the diluted sample has heterogeneous contaminants (this section will be updated later on).
##' @return OC curves when the diluted sample has heterogeneous contaminants.
##' @examples
##' c <- 0
##' meanlog_low <- -6
##' meanlog_high <- 3
##' sdlog <- 0.8
##' a <- 0
##' b <- 300
##' FDF <- 0.001
##' USL <- c(1000, 2000)
##' n <- 5
##' n_sim <- 1000
##' OC_curves_heterogeneous(c, meanlog_low, meanlog_high, sdlog, a, b, FDF, USL, n, n_sim)
##' @usage  OC_curves_heterogeneous(c, meanlog_low, meanlog_high, sdlog, a, b, FDF, USL, n, n_sim)
OC_curves_heterogeneous <- function(c, meanlog_low, meanlog_high, sdlog = 0.8, a, b, FDF, USL, n, n_sim){
  P_a <- NULL
  Dilution_scheme <- NULL
  meanlog <- seq(meanlog_low, meanlog_high, 0.1)
  f_spr <- function(USL) {
    sprintf("Scheme (USL=%.0f)", USL)
  }
  pa <- matrix(NA, nrow = length(meanlog), ncol = 2)
  for (i in 1:length(meanlog)) {
    # for (j in 1:2) {
    pa[i,1] <-  prob_acceptance_heterogeneous(c, meanlog[i], sdlog, a, b, FDF, USL[1], n, n_sim)
    pa[i,2] <-  prob_acceptance_heterogeneous(c, meanlog[i], sdlog, a, b, FDF, USL[2], n, n_sim)
    # }
  }
  Prob <- data.frame(meanlog, pa)
  colnames(Prob ) <- c("meanlog", f_spr(USL))
  melten.Prob <- reshape2::melt(Prob, id = "meanlog", variable.name = "Dilution_scheme", value.name = "P_a")
  plot_sam <- ggplot2::ggplot(melten.Prob) + ggplot2::geom_line(ggplot2::aes(x = meanlog, y = P_a, group = Dilution_scheme, colour = Dilution_scheme)) +
    # ggplot2::ggtitle("OC curve based on Lognormal distribution") +
    ggplot2::theme_classic() + ggplot2::xlab(expression("log mean concentrations  (" ~ mu*~")")) + ggplot2::ylab(expression(P[a])) + ggthemes::scale_colour_colorblind() +
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5), legend.position = c(0.85, 0.75), axis.line.x.top = ggplot2::element_line(color = "red"),
                   axis.ticks.x.top = ggplot2::element_line(color = "red"), axis.text.x.top = ggplot2::element_text(color = "red"), axis.title.x.top = ggplot2::element_text(color = "red")) +
    ggplot2::scale_x_continuous(sec.axis = ggplot2::sec_axis(~., name = "expected cell counts on plate", breaks = seq(min(meanlog),max(meanlog),1),
                                                       labels = c(sprintf("%f", 10^(seq(min(meanlog),max(meanlog),1) + (sdlog^2/2) * log(10, exp(1)))))))
  # plot_sam
  return(plot_sam)
}
