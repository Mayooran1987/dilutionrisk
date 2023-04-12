##' \code{\link{compare_plans_dilution_1_pd_binom_pois}} provides graphical displays of the probability of the detection curves for dilution schemes in the first dilution stage based on the binomial distribution, while the count of microorganisms is modelled by Poisson distribution.
##' @title Comparison based on binomial distribution-based probability of detection curves for different dilution schemes in the first dilution stage.
##' @param S amount of sample (in grams) used for diluted solution preparation
##' @param lambda_lower an lower value of the expected number of microorganisms count per gram.
##' @param lambda_upper an upper value of the expected number of microorganisms count per gram.
##' @param V0 dilution volume in the first dilution stage testing.
##' @param V1 the volume of the diluted solution used for testing (which is equal to the plated amount in this study)
##' @details \code{\link{compare_plans_dilution_1_pd_binom_pois}} provides graphical displays of the probability of the detection curves for dilution schemes in the first dilution stage based on the binomial distribution, while the count of microorganisms is modelled by Poisson distribution.(This section will be updated later on.)
##' @return Comparison based on probability of detection curves for different dilution schemes in the first dilution stage.
##' @examples
##' S <- c(25,125,250)
##' V0 <- 100
##' V1 <- 1
##' lambda_lower <- 0
##' lambda_upper <- 30
##' compare_plans_dilution_1_pd_binom_pois (S, V0, V1, lambda_lower, lambda_upper)
##' @usage  compare_plans_dilution_1_pd_binom_pois (S, V0, V1, lambda_lower, lambda_upper)
##' @export
compare_plans_dilution_1_pd_binom_pois <- function(S, V0, V1, lambda_lower, lambda_upper) {
  dilution_scheme <- NULL # Initalizing
  p_d <- NULL
  C <- NULL
  lambda <- seq(lambda_lower, lambda_upper, by = 0.1)
  f_spr <- function(S) {
    sprintf("sample weight (%.0f gram)", S)
  }
  pd <- matrix(NA, nrow = length(lambda), ncol = length(S))
  for (i in 1:length(lambda)) {
    pd[i, ] <- prob_detect_dilution_1_binom_pois(S, lambda[i], V0, V1)
  }
  Prob <- data.frame(lambda, pd)
  colnames(Prob) <- c("lambda", f_spr(S))
  melten.Prob <- reshape2::melt(Prob, id = "lambda", variable.name = "dilution_scheme", value.name = "p_d")
  plot_sam <- ggplot2::ggplot(melten.Prob) +
    ggplot2::geom_line(ggplot2::aes(x = lambda, y = p_d, group = dilution_scheme, colour = dilution_scheme)) +
    ggplot2::ylab(expression(P[D])) +
    ggplot2::xlab(expression(lambda ~ ("cell count per gram"))) +
    ggplot2::theme_classic() +
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5, size = 20), legend.position = c(0.75, 0.25), legend.text = ggplot2::element_text(size = 12)) +
    ggthemes::scale_colour_colorblind()
  return(plot_sam)
}
