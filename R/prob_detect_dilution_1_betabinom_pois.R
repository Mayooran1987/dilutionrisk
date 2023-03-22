##' \code{\link{prob_detect_dilution_1_betabinom_pois}} provides a probability of detection for a dilution scheme in the first dilution stage based on the beta binomial distribution, while the count of microorganisms is modelled by Poisson distribution.
##' @title Beta binomial distribution-based probability of detection for a dilution scheme in the first dilution stage testing.
##' @param S amount of sample (in grams) used for diluted solution preparation
##' @param lambda the expected number of microorganisms count per gram
##' @param V0 dilution volume in the first dilution stage testing.
##' @param V1 the volume of the diluted solution used for testing (which is equal to the plated amount in this study)
##' @param alpha non-negative parameters of the beta distribution.
##' @param beta non-negative parameters of the beta distribution.
##' @details \code{\link{prob_detect_dilution_1_betabinom_pois}} provides a probability of detection for a dilution scheme in the first dilution stage based on the beta binomial distribution, while the count of microorganisms is modelled by Poisson distribution.(This section will be updated later on.)
##' @return Probability of detection for a dilution scheme in the first dilution stage testing.
##' @examples
##' S <- 25
##' lambda  <- 10
##' V0 <- 100
##' V1 <- 1
##' alpha <- 1
##' beta <- 5
##' prob_detect_dilution_1_betabinom_pois(S, lambda, V0, V1, alpha, beta)
##' @usage  prob_detect_dilution_1_betabinom_pois(S, lambda, V0, V1, alpha, beta)
##' @export
prob_detect_dilution_1_betabinom_pois <- function(S, lambda, V0, V1, alpha, beta) {
  pd <- matrix(NA, nrow = 1, ncol = length(S))
  for (j in 1:length(S)) {
    # pd[, j] <- prob_detect_dilution_1_binom(S[j], lambda, V0, V1)
    N[j] <- round(mean(rpois(50000,S[j]*lambda) ))
    pd[, j] <- 1 - (beta(alpha, beta + N[j])) / (beta(alpha, beta))
  }
  return(pd)
}
