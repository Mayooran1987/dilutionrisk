##' \code{\link{pd_curves_homogeneous}} provides the probability of detection curves when samples collected from a homogeneous batch.
##' @title Comparison based on probability of detection curves for different dilution schemes when diluted samples collected from a homogeneous batch.
##' @param lambda_low the lower value of the expected microbial count (\eqn{\lambda}) for use in the graphical display's x-axis.
##' @param lambda_high the upper value of the expected microbial count (\eqn{\lambda}) for use in the graphical display's x-axis.
##' @param a lower domain of the number of microbial count.
##' @param b upper domain of the number of microbial count.
##' @param f final dilution factor.
##' @param u amount put on the plate.
##' @param USL upper specification limit.
##' @param n_sim number of simulations (large simulations provide more precise estimations).
##' @details \code{\link{pd_curves_homogeneous}} provides probability of detection curves for different dilution schemes when samples collected from a homogeneous batch (this section will be updated later on).
##' @return Probability of detection curves when diluted samples collected from a homogeneous batch.
##' @examples
##' lambda_low <- 0
##' lambda_high <- 5000
##' a <- 0
##' b <- 300
##' f <- c(0.01,0.1)
##' u <- c(0.1,0.1)
##' USL <- 1000
##' n_sim <- 50000
##' pd_curves_homogeneous(lambda_low, lambda_high, a, b, f, u, USL, n_sim)
##' @usage  pd_curves_homogeneous(lambda_low, lambda_high, a, b, f, u, USL, n_sim)
##' @export
pd_curves_homogeneous <- function(lambda_low, lambda_high, a, b, f, u, USL, n_sim){
  p_d <- NULL
  Dilution_scheme <- NULL
  f_spr <- function(f, u) {
    sprintf("Scheme (f=%.3f, u=%.1f)", f, u)
  }
  lambda <- seq(lambda_low, lambda_high, 0.1)
  # lambda <- 10^(mu + (sd^2/2) * log(10, exp(1)))
  Pd <- matrix(NA, nrow = length(lambda), ncol = length(f))
  for (i in 1:length(lambda)) {
    Pd[i,] <-  cbind(prob_detection_homogeneous_multiple(lambda[i], a, b, f, u, USL, n_sim))
  }
  # Pd <- matrix(NA, nrow = length(lambda), ncol = 2)
  # for (i in 1:length(lambda)) {
  #   Pd[i,1] <-  prob_detection_homogeneous(lambda[i], a, b, f[1], u[1], USL, n_sim)
  #   Pd[i,2] <-  prob_detection_homogeneous(lambda[i], a, b, f[2], u[2], USL, n_sim)
  # }
  Prob <- data.frame(lambda, Pd)
  colnames(Prob ) <- c("lambda", f_spr(f,u))
  melten.Prob <- reshape2::melt(Prob, id = "lambda", variable.name = "Dilution_scheme", value.name = "p_d")
  plot_sam <- ggplot2::ggplot(melten.Prob) + ggplot2::geom_line(ggplot2::aes(x = lambda, y = p_d, group = Dilution_scheme, colour = Dilution_scheme)) +
    ggplot2::theme_classic() + ggplot2::xlab(expression("expected microbial count  (" ~ lambda*~")")) + ggplot2::ylab(expression("Probability of detection"~(P[d]))) + ggthemes::scale_colour_colorblind() +
    ggplot2::geom_vline(xintercept = USL, linetype = "dashed") +
    ggplot2::annotate("text", x = USL,
                      y = 0, label = sprintf("USL = %0.0f", USL), size = 3) +
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5), legend.position = c(0.85, 0.25), axis.line.x.top = ggplot2::element_line(color = "red"),
                   axis.ticks.x.top = ggplot2::element_line(color = "red"), axis.text.x.top = ggplot2::element_text(color = "red"), axis.title.x.top = ggplot2::element_text(color = "red"))
  # +
  #   ggplot2::scale_x_continuous(sec.axis = ggplot2::sec_axis(~., name = "expected microbial count (cfu/g)", breaks = seq(min(mu),max(mu),1),
  #                                                            labels = c(sprintf("%f", 10^(seq(min(mu),max(mu),1) + (sd^2/2) * log(10, exp(1)))))))
  return(plot_sam)
}

