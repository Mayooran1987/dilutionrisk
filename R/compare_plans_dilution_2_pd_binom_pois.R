##' \code{\link{compare_plans_dilution_2_pd_binom_pois}} provides graphical displays of the probability of the detection curves for dilution schemes in the second dilution stage based on the binomial distribution, while the count of microorganisms is modelled by Poisson distribution.
##' @title Comparison based on binomial distribution-based probability of detection curves for different dilution schemes in the second dilution stage.
##' @param S amount of sample (in grams) used for diluted solution preparation
##' @param V0 dilution volume in the first dilution stage testing.
##' @param V1 the volume of the diluted solution used for testing (which is equal to the plated amount in this study)
##' @param V2 dilution volume in the first dilution stage testing.
##' @param V3 the volume of the diluted solution used for testing (which is equal to the plated amount in this study)
##' @param n_sim number of simulations (large simulations provide a more precise estimation).
##' @param lambda_lower an lower value of the expected number of microorganisms count per gram.
##' @param lambda_upper an upper value of the expected number of microorganisms count per gram.
##' @details \code{\link{compare_plans_dilution_2_pd_binom_pois}} provides graphical displays of the probability of the detection curves for dilution schemes in the second dilution stage based on the binomial distribution, while the count of microorganisms is modelled by Poisson distribution. (This section will be updated later on.)
##' @return Comparison based on probability of detection curves for different dilution schemes in the second dilution stage.
##' @examples
##' S <- c(25,125,250)
##' V0 <- 100
##' V1 <- 1
##' V2 <- 10
##' V3 <- 1
##' n_sim <- 20000
##' lambda_lower <- 1
##' lambda_upper <- 100
##' compare_plans_dilution_2_pd_binom_pois(S, V0, V1, V2, V3, n_sim, lambda_lower, lambda_upper)
##' @usage  compare_plans_dilution_2_pd_binom_pois(S, V0, V1, V2, V3, n_sim, lambda_lower, lambda_upper)
##' @export
compare_plans_dilution_2_pd_binom_pois <- function(S, V0, V1, V2, V3, n_sim, lambda_lower, lambda_upper) {
  Sampling_scheme <- NULL # Initalizing
  p_d <- NULL
  C <- NULL
  lambda <- seq(lambda_lower, lambda_upper, by = 1)
  f_spr <- function(S) {
    sprintf("sample weight (%.0f gram)", S)
  }
  prob_detect_dilution_2_multi_binom_pois <- function(S,lambda, V0, V1, V2, V3, n_sim) {
    pd <- matrix(NA, nrow = 1, ncol = length(lambda))
    for (j in 1:length(lambda)) {
      pd[1, j] <- prob_detect_dilution_2_binom_pois(S,lambda[j], V0, V1, V2, V3, n_sim)
    }
    result <- as.numeric(pd)
    return(result)
  }
  Pd <- matrix(NA, nrow = length(lambda), ncol = length(S))
  for (j in 1:length(S)) {
    Pd[, j] <- prob_detect_dilution_2_multi_binom_pois(S[j], lambda,V0, V1, V2, V3, n_sim)
  }
  Prob <- data.frame(lambda, Pd)
  colnames(Prob) <- c("lambda", f_spr(S))
  melten.Prob <- reshape2::melt(Prob, id = "lambda", variable.name = "Sampling_scheme", value.name = "p_d")
  plot_sam <- ggplot2::ggplot(melten.Prob) +
    ggplot2::geom_line(ggplot2::aes(x = lambda, y = p_d, group = Sampling_scheme, colour = Sampling_scheme)) +
    ggplot2::ylab(expression(P[D])) +
    ggplot2::xlab(expression(lambda ("cell count per gram"))) +
    ggplot2::theme_classic() +
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5, size = 20), legend.position = c(0.75, 0.25), legend.text = ggplot2::element_text(size=12)) +
    ggthemes::scale_colour_colorblind()
  return(plot_sam)
}
