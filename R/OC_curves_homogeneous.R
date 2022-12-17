##' \code{\link{OC_curves_homogeneous}} provides the operating characteristic(OC) curves when diluted sample has homogeneous contaminants.
##' @title Comparison based on OC curves for different dilution schemes when diluted samples collected from a homogeneous batch.
##' @param c acceptance number
##' @param lambda_low the lower value of the expected microbial count(\eqn{\lambda}) for use in the graphical display's x-axis.
##' @param lambda_high the upper value of the expected microbial count(\eqn{\lambda}) for use in the graphical display's x-axis.
##' @param a lower domain of the number of cell counts.
##' @param b upper domain of the number of cell counts.
##' @param f final dilution factor.
##' @param u amount put on the plate.
##' @param USL upper specification limit.
##' @param n number of samples which are used for inspection.
##' @param type what type of the results you would like to consider such as "theory" or "simulation" (default "theory").
##' @param n_sim number of simulations (large simulations provide more precise estimations).
##' @details \code{\link{OC_curves_homogeneous}} provides OC curves for different dilution schemes when samples collected from a homogeneous batch (this section will be updated later on).
##' @return OC curves when diluted samples collected from a homogeneous batch.
##' @examples
##' c <- 2
##' lambda_low <- 1
##' lambda_high <- 5000
##' a <- 0
##' b <- 300
##' f <- c(0.01,0.1)
##' u <- c(0.1,0.1)
##' USL <- 1000
##' n <- 5
##' OC_curves_homogeneous(c, lambda_low, lambda_high, a, b, f, u, USL, n)
##' @usage  OC_curves_homogeneous(c, lambda_low, lambda_high, a, b, f, u, USL, n, type, n_sim)
##' @export
OC_curves_homogeneous <- function(c, lambda_low, lambda_high, a, b, f, u, USL, n, type = "theory", n_sim = NA){
  P_a <- NULL
  Dilution_scheme <- NULL
  lambda <- seq(lambda_low, lambda_high, 0.1)
  f_spr <- function(f, u) {
    sprintf("Scheme (f=%.3f, u=%.1f)", f, u)
  }
  pa <- matrix(NA, nrow = length(lambda), ncol = length(f))
  for (i in 1:length(lambda)) {
    pa[i,] <-  cbind(prob_acceptance_homogeneous_multiple(c, lambda[i], a, b, f, u, USL, n, type, n_sim))
  }
  Prob <- data.frame(lambda, pa)
  colnames(Prob ) <- c("lambda", f_spr(f,u))
  melten.Prob <- reshape2::melt(Prob, id = "lambda", variable.name = "Dilution_scheme", value.name = "P_a")
  plot_sam <- ggplot2::ggplot(melten.Prob) + ggplot2::geom_line(ggplot2::aes(x = lambda, y = P_a, group = Dilution_scheme, colour = Dilution_scheme)) +
    # ggplot2::ggtitle("OC curve based on Lognormal distribution") +
    ggplot2::theme_classic() + ggplot2::xlab(expression("expected microbial count  (" ~ lambda*~")")) + ggplot2::ylab(expression("Probability of acceptance"~(P[a]))) + ggthemes::scale_colour_colorblind() +
    ggplot2::geom_vline(xintercept = USL, linetype = "dashed") +
    ggplot2::annotate("text", x = USL,
                      y = 0, label = sprintf("USL = %0.0f", USL), size = 3) +
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5), legend.position = c(0.85, 0.75), axis.line.x.top = ggplot2::element_line(color = "red"),
                   axis.ticks.x.top = ggplot2::element_line(color = "red"), axis.text.x.top = ggplot2::element_text(color = "red"), axis.title.x.top = ggplot2::element_text(color = "red"))
  # +
  #   ggplot2::scale_x_continuous(sec.axis = ggplot2::sec_axis(~., name = "expected cell counts (cfu/g)", breaks = seq(min(mu),max(mu),1),
  #                                                            labels = c(sprintf("%f", 10^(seq(min(mu),max(mu),1) + (sd^2/2) * log(10, exp(1)))))))
  # plot_sam
  return(plot_sam)
}
