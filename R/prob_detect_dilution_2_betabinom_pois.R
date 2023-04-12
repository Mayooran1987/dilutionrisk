##' \code{\link{prob_detect_dilution_2_betabinom_pois}} provides a probability of detection for a dilution scheme in the second dilution stage based on the beta binomial distribution, while the count of microorganisms is modelled by Poisson distribution.
##' @title Beta binomial distribution-based probability of detection for a dilution scheme in the second dilution stage testing.
##' @param S amount of sample (in grams) used for diluted solution preparation
##' @param lambda the expected number of microorganisms count per gram
##' @param V0 dilution volume in the first dilution stage testing.
##' @param V1 the volume of the diluted solution used for testing (which is equal to the plated amount in this study)
##' @param V2 dilution volume in the first dilution stage testing.
##' @param V3 the volume of the diluted solution used for testing (which is equal to the plated amount in this study)
##' @param n_sim number of simulations (large simulations provide a more precise estimation).
##' @param alpha non-negative parameters of the beta distribution.
##' @param beta non-negative parameters of the beta distribution.
##' @details \code{\link{prob_detect_dilution_2_betabinom_pois}} provides a probability of detection for a dilution scheme in the second dilution stage based on the beta binomial distribution, while the count of microorganisms is modelled by Poisson distribution.(This section will be updated later on.)
##' @return Probability of detection for a dilution scheme in the second dilution stage testing.
##' @examples
##' S <- 20
##' lambda  <- 10
##' V0 <- 100
##' V1 <- 1
##' V2 <- 10
##' V3 <- 1
##' n_sim <- 25000
##' alpha <- 1
##' beta <- 5
##' prob_detect_dilution_2_betabinom_pois(S, lambda, V0, V1, V2, V3, n_sim, alpha, beta)
##' @usage  prob_detect_dilution_2_betabinom_pois(S, lambda, V0, V1, V2, V3, n_sim, alpha, beta)
##' @export
prob_detect_dilution_2_betabinom_pois <- function(S, lambda, V0, V1, V2, V3, n_sim, alpha, beta) {
  N <- round(mean(rpois(50000, S * lambda)))
  n1 <- matrix(V0 * extraDistr::rbbinom(n_sim, N, alpha, beta), nrow = n_sim, ncol = 1)
  pd <- matrix(NA, nrow = n_sim, ncol = 1)
  for (j in 1:n_sim) {
    pd[j, ] <- 1 - (beta(alpha, beta + n1[j, ])) / (beta(alpha, beta))
  }
  result <- mean(pd)
  return(result)
}
