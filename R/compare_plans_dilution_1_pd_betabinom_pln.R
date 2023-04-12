##' \code{\link{compare_plans_dilution_1_pd_betabinom_pln}} provides graphical displays of the probability of the detection curves for dilution schemes in the first dilution stage based on the binomial distribution, while the count of microorganisms is modelled by Poisson lognormal distribution.
##' @title Comparison based on beta binomial distribution-based  probability of detection curves for different dilution schemes in the first dilution stage.
##' @param S amount of sample (in grams) used for diluted solution preparation
##' @param sd the standard deviation of the normal distribution (on the log scale).
##' @param V0 dilution volume in the first dilution stage testing.
##' @param V1 the volume of the diluted solution used for testing (which is equal to the plated amount in this study)
##' @param mu_lower the lower value of the mean microbial count(\eqn{\mu}) for use in the graphical display's x-axis.
##' @param mu_upper the upper value of the mean microbial count(\eqn{\mu}) for use in the graphical display's x-axis.
##' @param alpha non-negative parameters of the beta distribution.
##' @param beta non-negative parameters of the beta distribution.
##' @details \code{\link{compare_plans_dilution_1_pd_betabinom_pln}} provides graphical displays of the probability of the detection curves for dilution schemes in the first dilution stage based on the binomial distribution, while the count of microorganisms is modelled by Poisson lognormal distribution. (This section will be updated later on.)
##' @return Comparison based on probability of detection curves for different dilution schemes in the first dilution stage.
##' @examples
##' S <- c(25,125,250)
##' sd <- 0.8
##' V0 <- 100
##' V1 <- 1
##' mu_lower <- -8
##' mu_upper <- 8
##' alpha <- 1
##' beta <- 5
##' compare_plans_dilution_1_pd_betabinom_pln(S, sd, V0, V1, mu_lower, mu_upper, alpha, beta)
##' @usage  compare_plans_dilution_1_pd_betabinom_pln(S, sd, V0, V1, mu_lower, mu_upper, alpha, beta)
##' @export
compare_plans_dilution_1_pd_betabinom_pln <- function(S, sd, V0, V1, mu_lower, mu_upper, alpha, beta) {
  dilution_scheme <- NULL # Initalizing
  p_d <- NULL
  mu <- seq(mu_lower, mu_upper, by = 0.1)
  f_spr <- function(S) {
    sprintf("sample weight (%.0f gram)", S)
  }
  pd <- matrix(NA, nrow = length(mu), ncol = length(S))
  for (i in 1:length(mu)) {
    pd[i, ] <- prob_detect_dilution_1_betabinom_pln(S, mu[i], sd, V0, V1, alpha, beta)
  }
  Prob <- data.frame(mu, pd)
  colnames(Prob) <- c("mu", f_spr(S))
  melten.Prob <- reshape2::melt(Prob, id = "mu", variable.name = "dilution_scheme", value.name = "p_d")
  plot_sam <- ggplot2::ggplot(melten.Prob) +
    # ggplot2::geom_line(ggplot2::aes(x = mu, y = p_d, group = dilution_scheme, colour = dilution_scheme)) +
    ggplot2::stat_smooth(geom = "smooth", method = "gam", mapping = ggplot2::aes(x = mu, y = p_d, group = dilution_scheme, colour = dilution_scheme), se = FALSE, n = 1000) +
    # ggplot2::ylim(0,1)+
    ggplot2::ylab(expression(P[D])) +
    ggplot2::xlab(expression(mu ~ ("log mean concentration"))) +
    ggplot2::theme_classic() +
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5, size = 20), legend.position = c(0.75, 0.25), legend.text = ggplot2::element_text(size = 12)) +
    ggthemes::scale_colour_colorblind()
  return(plot_sam)
}
