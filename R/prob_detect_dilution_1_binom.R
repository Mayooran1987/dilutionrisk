##' \code{\link{prob_detect_dilution_1_binom}} provides the probability of detection for a dilution scheme in the first dilution stage testing.
##' @title Probability of detection for a dilution scheme in the first dilution stage testing.
##' @param S amount of sample (in grams) used for diluted solution preparation
##' @param C_cfu number of colonies counted on the two dishes retained from two successive dilutions.
##' @param V0 dilution volume in the first dilution stage testing.
##' @param V1 the volume of the diluted solution used for testing (which is equal to the plated amount in this study)
##' @details \code{\link{prob_detect_dilution_1_binom}} provides the probability of detection for a dilution scheme in the first dilution stage testing.
##' @return Probability of detection for a dilution scheme in the first dilution stage testing.
##' @examples
##' S <- c(1,10,20)
##' C_cfu  <- 25
##' V0 <- 100
##' V1 <- 1
##' prob_detect_dilution_1_binom(S, V0, V1, C_cfu )
##' @usage  prob_detect_dilution_1_binom(S, V0, V1, C_cfu)
##' @export
prob_detect_dilution_1_binom <- function(S, V0, V1, C_cfu) {
  N <- (C_cfu / V1) * ((S + V0) / S)
  Pd <- 1 - stats::dbinom(0, N, V1 / V0)
  return(Pd)
}
