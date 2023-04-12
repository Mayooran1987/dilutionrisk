##' \code{\link{prob_detect_dilution_2_binom_pois}} provides a probability of detection for a dilution scheme in the second dilution stage based on the binomial distribution, while the count of microorganisms is modelled by Poisson distribution.
##' @title Binomial distribution-based probability of detection for a dilution scheme in the second dilution stage testing.
##' @param S amount of sample (in grams) used for diluted solution preparation
##' @param lambda the expected number of microorganisms count per gram
##' @param V0 dilution volume in the first dilution stage testing.
##' @param V1 the volume of the diluted solution used for testing (which is equal to the plated amount in this study)
##' @param V2 dilution volume in the first dilution stage testing.
##' @param V3 the volume of the diluted solution used for testing (which is equal to the plated amount in this study)
##' @param n_sim number of simulations (large simulations provide a more precise estimation).
##' @details \code{\link{prob_detect_dilution_2_binom_pois}} provides a probability of detection for a dilution scheme in the second dilution stage based on the binomial distribution, while the count of microorganisms is modelled by Poisson distribution.(This section will be updated later on.)
##' @return Probability of detection for a dilution scheme in the second dilution stage testing.
##' @examples
##' S <- 20
##' lambda  <- 10
##' V0 <- 100
##' V1 <- 1
##' V2 <- 10
##' V3 <- 1
##' n_sim <- 25000
##' prob_detect_dilution_2_binom_pois(S,lambda, V0, V1, V2, V3, n_sim)
##' @usage  prob_detect_dilution_2_binom_pois(S,lambda, V0, V1, V2, V3, n_sim)
##' @export
prob_detect_dilution_2_binom_pois <- function(S, lambda, V0, V1, V2, V3, n_sim) {
  # N <- (C_cfu / V1) * ((S + V0) / S)
  N <- round(mean(rpois(5000, S * lambda)))
  # N <- rpois(1,S*lambda)
  # number of successes must be at least one.
  rcond.binomial <- function(n_sim, N, p) {
    probs <- stats::dbinom(1:N, N, p) / (1 - p)
    sample(1:N, n_sim, replace = TRUE, prob = probs)
  }
  n1 <- matrix(V0 * rcond.binomial(n_sim, N, V1 / V0), nrow = n_sim, ncol = 1)
  pd <- matrix(NA, nrow = n_sim, ncol = 1)
  for (j in 1:n_sim) {
    pd[j, ] <- 1 - dbinom(0, n1[j, ], V3 / V2)
  }
  result <- mean(pd)
  return(result)
}
