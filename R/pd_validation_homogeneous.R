##' \code{\link{pd_validation_homogeneous}} provides the probability of detection curves for validate the results when samples collected from a homogeneous batch.
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
##' f <- 0.01
##' u <- 0.1
##' USL <- 1000
##' n_sim <- 50000
##' pd_validation_homogeneous(lambda_low, lambda_high, a, b, f, u, USL, n_sim = 500)
##' pd_validation_homogeneous(lambda_low, lambda_high, a, b, f, u, USL, n_sim = 50000)
##' @usage  pd_validation_homogeneous(lambda_low, lambda_high, a, b, f, u, USL, n_sim)
##' @export
pd_validation_homogeneous <- function(lambda_low, lambda_high, a, b, f, u, USL, n_sim){
  p_d <- NULL
  Methods <- NULL
  lambda <- seq(lambda_low, lambda_high, 1)
  # pd_theory <- function(lambda, f, u, USL){
  #   USL1 <- USL*f*u
  #   pd <- matrix(NA, nrow = USL1, ncol = 1)
  #   for (i in 1:USL1) {
  #     pd[i,1] <-  (u*f*lambda)^i/(factorial(i)*(exp(u*f*lambda) - 1))
  #   }
  #   pd <- 1 - sum(pd)
  #   return(pd)
  # }
  # pd_theory(lambda = 1, f, u, USL)
  Pd <- matrix(NA, nrow = length(lambda), ncol = 2)
  for (i in 1:length(lambda)) {
    Pd[i,1] <-  prob_detection_homogeneous(lambda[i], a, b, f, u, USL, type = "simulation", n_sim)
    Pd[i,2] <-  prob_detection_homogeneous(lambda[i], a, b, f, u, USL, type = "theory")
  }
  Prob <- data.frame(lambda, Pd)
  # colnames(Prob ) <- c("lambda", f_spr(f,u))
  colnames(Prob ) <- c("lambda", "simulation", "theory")

  melten.Prob <- reshape2::melt(Prob, id = "lambda", variable.name = "Methods", value.name = "p_d")
  plot_sam <- ggplot2::ggplot(melten.Prob) + ggplot2::geom_line(ggplot2::aes(x = lambda, y = p_d, group = Methods, colour = Methods)) +
    ggplot2::theme_classic() + ggplot2::xlab(expression("expected microbial count (" ~ lambda*~")")) + ggplot2::ylab(expression("Probability of detection"~(P[d]))) + ggthemes::scale_colour_colorblind() +
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5), legend.position = c(0.85, 0.25), axis.line.x.top = ggplot2::element_line(color = "red"),
                   axis.ticks.x.top = ggplot2::element_line(color = "red"), axis.text.x.top = ggplot2::element_text(color = "red"), axis.title.x.top = ggplot2::element_text(color = "red"))
  return(plot_sam)
}
