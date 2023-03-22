##' \code{\link{compare_plans_dilution_2_pd_betabinom_pois}} provides graphical displays of the probability of the detection curves for dilution schemes in the second dilution stage based on the beta binomial distribution, while the count of microorganisms is modelled by Poisson distribution.
##' @title Comparison based on beta binomial distribution-based probability of detection curves for different dilution schemes in the second dilution stage.
##' @param S amount of sample (in grams) used for diluted solution preparation
##' @param V0 dilution volume in the first dilution stage testing.
##' @param V1 the volume of the diluted solution used for testing (which is equal to the plated amount in this study)
##' @param V2 dilution volume in the first dilution stage testing.
##' @param V3 the volume of the diluted solution used for testing (which is equal to the plated amount in this study)
##' @param n_sim number of simulations (large simulations provide a more precise estimation).
##' @param lambda_lower an lower value of the expected number of microorganisms count per gram.
##' @param lambda_upper an upper value of the expected number of microorganisms count per gram
##' @param alpha non-negative parameters of the beta distribution.
##' @param beta non-negative parameters of the beta distribution.
##' @details \code{\link{compare_plans_dilution_2_pd_betabinom_pois}} provides graphical displays of the probability of the detection curves for dilution schemes in the second dilution stage based on the beta binomial distribution, while the count of microorganisms is modelled by Poisson distribution.(This section will be updated later on.)
##' @return Comparison based on probability of detection curves for different dilution schemes in the second dilution stage.
##' @examples
##' S <- c(1,10,20)
##' V0 <- 100
##' V1 <- 1
##' V2 <- 10
##' V3 <- 1
##' n_sim <- 20000
##' lambda_lower <- 0
##' lambda_upper <- 1000
##' alpha <- 1
##' beta <- 5
##' compare_plans_dilution_2_pd_betabinom_pois(S, V0, V1, V2, V3, n_sim, lambda_lower, lambda_upper, alpha, beta)
##' @usage  compare_plans_dilution_2_pd_betabinom_pois(S, V0, V1, V2, V3, n_sim, lambda_lower, lambda_upper, alpha, beta)
##' @export
compare_plans_dilution_2_pd_betabinom_pois <- function(S, V0, V1, V2, V3, n_sim, lambda_lower, lambda_upper, alpha, beta) {
  message("\033[1;31m","This function takes a few hours to produce the output! Thanks for your patience.")
  Sampling_scheme <- NULL # Initalizing
  p_d <- NULL
  C <- NULL
  lambda <- seq(lambda_lower, lambda_upper, by = 0.01)
  f_spr <- function(S) {
    sprintf("sample weight (%.0f gram)", S)
  }
  prob_detect_dilution_2_multi_betabinom_pois <- function(S,lambda, V0, V1, V2, V3, n_sim, alpha, beta) {
    pd <- matrix(NA, nrow = 1, ncol = length(lambda))
    for (j in 1:length(lambda)) {
      pd[1, j] <- prob_detect_dilution_2_betabinom_pois(S,lambda[j], V0, V1, V2, V3, n_sim, alpha, beta)
    }
    result <- as.numeric(pd)
    return(result)
  }
  Pd <- matrix(NA, nrow = length(lambda), ncol = length(S))
  for (j in 1:length(S)) {
    Pd[, j] <- prob_detect_dilution_2_multi_betabinom_pois(S[j], lambda,V0, V1, V2, V3, n_sim, alpha, beta)
  }
  Prob <- data.frame(lambda, Pd)
  colnames(Prob) <- c("lambda", f_spr(S))
  melten.Prob <- reshape2::melt(Prob, id = "lambda", variable.name = "Sampling_scheme", value.name = "p_d")
  plot_sam <- ggplot2::ggplot(melten.Prob) +
    ggplot2::geom_line(ggplot2::aes(x = lambda, y = p_d, group = Sampling_scheme, colour = Sampling_scheme)) +
    # ggplot2::stat_smooth(data = melten.Prob,size = 0.5, method = 'gam', formula = y ~ s(x, bs = "cs") , ggplot2::aes(x = C, y = p_d, group = Sampling_scheme, colour = Sampling_scheme),se = FALSE, na.rm = TRUE)+
    ggplot2::ylab(expression(P[D])) +
    ggplot2::xlab(expression(lambda ("cell count per gram"))) +
    ggplot2::theme_classic() +
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5, size = 20), legend.position = c(0.75, 0.25), legend.text = ggplot2::element_text(size=12)) +
    ggthemes::scale_colour_colorblind()
  return(plot_sam)

}
