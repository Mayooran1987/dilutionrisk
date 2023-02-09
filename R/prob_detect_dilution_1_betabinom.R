##' \code{\link{prob_detect_dilution_1_betabinom}} provides the probability of detection for a dilution scheme in the first dilution stage testing.
##' @title Probability of detection for a dilution scheme in the first dilution stage testing.
##' @param S amount of sample (in grams) used for diluted solution preparation
##' @param C_cfu number of colonies counted on the two dishes retained from two successive dilutions.
##' @param V0 dilution volume in the first dilution stage testing.
##' @param V1 the volume of the diluted solution used for testing (which is equal to the plated amount in this study)
##' @param alpha non-negative parameters of the beta distribution.
##' @param beta non-negative parameters of the beta distribution.
##' @details \code{\link{prob_detect_dilution_1_betabinom}} provides the probability of detection for a dilution scheme in the first dilution stage testing.
##' @return Probability of detection for a dilution scheme in the first dilution stage testing.
##' @examples
##' S <- 25
##' C_cfu  <- 25
##' V0 <- 100
##' V1 <- 1
##' alpha <- 1
##' beta <- 5
##' prob_detect_dilution_1_betabinom(S, V0, V1, C_cfu, alpha, beta)
##' @usage  prob_detect_dilution_1_betabinom(S, V0, V1, C_cfu, alpha, beta)
##' @export
prob_detect_dilution_1_betabinom <- function(S, V0, V1, C_cfu, alpha, beta) {
  N <- (C_cfu / V1) * ((S + V0) / S)
  Pd <- 1 - (beta(alpha, beta + N)) / (beta(alpha, beta))
  return(Pd)
}
