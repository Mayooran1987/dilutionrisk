##' \code{\link{prob_detect_dilution_1_binom_pois}} provides a probability of detection for a dilution scheme in the first dilution stage based on the binomial distribution, while the count of microorganisms is modelled by Poisson distribution.
##' @title Binomial distribution-based probability of detection for a dilution scheme in the first dilution stage testing.
##' @param S amount of sample (in grams) used for diluted solution preparation
##' @param lambda the expected number of microorganisms count per gram
##' @param V0 dilution volume in the first dilution stage testing.
##' @param V1 the volume of the diluted solution used for testing (which is equal to the plated amount in this study)
##' @details \code{\link{prob_detect_dilution_1_binom_pois}} provides a probability of detection for a dilution scheme in the first dilution stage based on the binomial distribution, while the count of microorganisms is modelled by Poisson distribution.(This section will be updated later on.)
##' @return Probability of detection for a dilution scheme in the first dilution stage testing.
##' @examples
##' S <- c(1,10,20)
##' lambda  <- 10
##' V0 <- 100
##' V1 <- 1
##' prob_detect_dilution_1_binom_pois(S,lambda, V0, V1)
##' @usage  prob_detect_dilution_1_binom_pois(S,lambda, V0, V1)
##' @export
prob_detect_dilution_1_binom_pois <- function(S,lambda, V0, V1) {
  pd <- matrix(NA, nrow = 1, ncol = length(S))
  for (j in 1:length(S)) {
    # pd[, j] <- prob_detect_dilution_1_binom(S[j], lambda, V0, V1)
    N[j] <-  round(mean(rpois(50000,S[j]*lambda) ))
    pd[, j] <- 1 - stats::dbinom(0, N[j], V1 / V0)
  }
  return(pd)
}
