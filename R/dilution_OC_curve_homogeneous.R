##' \code{\link{dilution_OC_curve_homogeneous}} provides the operating characteristic(OC) curves when diluted sample has homogeneous contaminants.
##' @title Comparison based on OC curves for different dilution schemes when the diluted sample has homogeneous contaminants.
##' @param c acceptance number
##' @param lambda_low the lower value of the expected cell count (\eqn{\lambda}) for use in the graphical display's x-axis.
##' @param lambda_high the upper value of the expected cell count (\eqn{\lambda}) for use in the graphical display's x-axis.
##' @param a lower domain of the number of cell counts.
##' @param b upper domain of the number of cell counts.
##' @param pf plating factor (pf = 1/final dilution factor).
##' @param USL upper specification limit.
##' @param n number of samples which are used for inspection.
##' @param n_sim number of simulations (large simulations provide more precise estimation).
##' @details \code{\link{dilution_OC_curve_homogeneous}} provides OC curves for different dilution schemes when the diluted sample has homogeneous contaminants (this section will be updated later on).
##' @return OC curves when the diluted sample has homogeneous contaminants.
##' @examples
##' c <- 0
##' lambda_low <- 0
##' lambda_high <- 5
##' a <- 0
##' b <- 300
##' pf <- 1000
##' USL <- c(1000, 2000)
##' n <- 5
##' n_sim <- 10000
##' dilution_OC_curve_homogeneous(c, lambda_low, lambda_high, a, b, pf, USL, n, n_sim)
##' @usage  dilution_OC_curve_homogeneous(c, lambda_low, lambda_high, a, b, pf, USL, n, n_sim)
dilution_OC_curve_homogeneous <- function(c, lambda_low, lambda_high, a, b, pf, USL, n, n_sim){
  P_a <- NULL
  Dilution_scheme <- NULL
  lambda <- seq(lambda_low, lambda_high, 0.1)
  f_spr <- function(pf) {
      sprintf("Scheme (USL=%.0f)", USL)
  }
  pa <- matrix(NA, nrow = length(lambda), ncol = 2)
  for (i in 1:length(lambda)) {
    # for (j in 1:2) {
    pa[i,1] <-  prob_acceptance_homogeneous(c, lambda[i], a, b, pf, USL[1], n, n_sim)
    pa[i,2] <-  prob_acceptance_homogeneous(c, lambda[i], a, b, pf, USL[2], n, n_sim)
    # }
  }
  Prob <- data.frame(lambda, pa)
  colnames(Prob ) <- c("lambda", f_spr(USL))
  melten.Prob <- reshape2::melt(Prob, id = "lambda", variable.name = "Dilution_scheme", value.name = "P_a")
  plot_sam <- ggplot2::ggplot(melten.Prob) + ggplot2::geom_line(ggplot2::aes(x = lambda, y = P_a, group = Dilution_scheme, colour = Dilution_scheme)) +
    # ggplot2::ggtitle("OC curve based on Lognormal distribution") +
    ggplot2::theme_classic() + ggplot2::xlab(expression("expected cell counts  (" ~ lambda*~")")) + ggplot2::ylab(expression(P[a])) + ggthemes::scale_colour_colorblind() +
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5), legend.position = c(0.85, 0.75), axis.line.x.top = ggplot2::element_line(color = "red"),
                   axis.ticks.x.top = ggplot2::element_line(color = "red"), axis.text.x.top = ggplot2::element_text(color = "red"), axis.title.x.top = ggplot2::element_text(color = "red"))
  # +
  #   ggplot2::scale_x_continuous(sec.axis = ggplot2::sec_axis(~., name = "expected cell counts (cfu/g)", breaks = seq(min(mu),max(mu),1),
  #                                                            labels = c(sprintf("%f", 10^(seq(min(mu),max(mu),1) + (sd^2/2) * log(10, exp(1)))))))
  # plot_sam
  return(plot_sam)
}

