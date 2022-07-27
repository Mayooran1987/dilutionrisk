##' \code{\link{pd_curves_homogeneous}} provides the probability of detection curves when the diluted sample has homogeneous contaminants.
##' @title comparison based on probability of detection curves for different dilution schemes when the diluted sample has homogeneous contaminants.
##' @param lambda_low the lower value of the expected cell count (\eqn{\lambda}) for use in the graphical display's x-axis.
##' @param lambda_high the upper value of the expected cell count (\eqn{\lambda}) for use in the graphical display's x-axis.
##' @param a lower domain of the number of cell counts.
##' @param b upper domain of the number of cell counts.
##' @param FDF final dilution factor.
##' @param USL upper specification limit.
##' @param n_sim number of simulations (large simulations provide more precise estimation).
##' @details \code{\link{pd_curves_homogeneous}} provides probability of detection curves for different dilution schemes when the diluted sample has homogeneous contaminants (this section will be updated later on).
##' @return Probability of detection curves when the diluted sample has homogeneous contaminants.
##' @examples
##' lambda_low <- 0
##' lambda_high <- 10
##' a <- 0
##' b <- 300
##' FDF <- 0.001
##' USL <- c(1000, 2000)
##' n_sim <- 1000
##' pd_curves_homogeneous(lambda_low, lambda_high, a, b, FDF, USL, n_sim)
##' @usage  pd_curves_homogeneous(lambda_low, lambda_high, a, b, FDF, USL, n_sim)
pd_curves_homogeneous <- function(lambda_low, lambda_high, a, b, FDF, USL, n_sim){
  p_d <- NULL
  Sampling_scheme <- NULL
  f_spr <- function(USL) {
    sprintf("Scheme (USL=%.0f)", USL)
  }
  lambda <- seq(lambda_low, lambda_high, 0.1)
  # lambda <- 10^(mu + (sd^2/2) * log(10, exp(1)))
  Pd <- matrix(NA, nrow = length(lambda), ncol = 2)
  for (i in 1:length(lambda)) {
    Pd[i,1] <-  prob_detection_homogeneous(lambda[i], a, b, FDF, USL[1], n_sim)
    Pd[i,2] <-  prob_detection_homogeneous(lambda[i], a, b, FDF, USL[2], n_sim)
  }
  Prob <- data.frame(lambda, Pd)
  colnames(Prob ) <- c("lambda", f_spr(USL))
  melten.Prob <- reshape2::melt(Prob, id = "lambda", variable.name = "Sampling_scheme", value.name = "p_d")
  plot_sam <- ggplot2::ggplot(melten.Prob) + ggplot2::geom_line(ggplot2::aes(x = lambda, y = p_d, group = Sampling_scheme, colour = Sampling_scheme)) +
    # ggplot2::ggtitle("OC curve based on Lognormal distribution") +
    ggplot2::theme_classic() + ggplot2::xlab(expression("expected cell counts  (" ~ lambda*~")")) + ggplot2::ylab(expression(p[d])) + ggthemes::scale_colour_colorblind() +
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5), legend.position = c(0.85, 0.75), axis.line.x.top = ggplot2::element_line(color = "red"),
                   axis.ticks.x.top = ggplot2::element_line(color = "red"), axis.text.x.top = ggplot2::element_text(color = "red"), axis.title.x.top = ggplot2::element_text(color = "red"))
  # +
  #   ggplot2::scale_x_continuous(sec.axis = ggplot2::sec_axis(~., name = "expected cell counts (cfu/g)", breaks = seq(min(mu),max(mu),1),
  #                                                            labels = c(sprintf("%f", 10^(seq(min(mu),max(mu),1) + (sd^2/2) * log(10, exp(1)))))))
  return(plot_sam)
}

