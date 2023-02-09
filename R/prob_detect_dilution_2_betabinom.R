##' \code{\link{prob_detect_dilution_2_betabinom}} provides the probability of detection for a dilution scheme in the second dilution stage testing.
##' @title Probability of detection for a dilution scheme in the second dilution stage testing.
##' @param S amount of sample (in grams) used for diluted solution preparation
##' @param V0 dilution volume in the first dilution stage testing.
##' @param V1 the volume of the diluted solution used for testing (which is equal to the plated amount in this study)
##' @param V2 dilution volume in the first dilution stage testing.
##' @param V3 the volume of the diluted solution used for testing (which is equal to the plated amount in this study)
##' @param C_cfu number of colonies counted on the two dishes retained from two successive dilutions.
##' @param n_sim number of simulations (large simulations provide a more precise estimation).
##' @param alpha non-negative parameters of the beta distribution.
##' @param beta non-negative parameters of the beta distribution.
##' @details \code{\link{prob_detect_dilution_2_betabinom}} provides the probability of detection for a dilution scheme in the second dilution stage testing.
##' @return Probability of detection for a dilution scheme in the second dilution stage testing.
##' @examples
##' S <- 20
##' C_cfu  <- 25
##' V0 <- 100
##' V1 <- 1
##' V2 <- 10
##' V3 <- 1
##' C_cfu  <- 25
##' n_sim <- 25000
##' alpha <- 1
##' beta <- 5
##' prob_detect_dilution_2_betabinom(S, V0,V1,V2,V3,C_cfu, n_sim, alpha, beta)
##' @usage  prob_detect_dilution_2_betabinom(S, V0,V1,V2,V3,C_cfu, n_sim, alpha, beta)
##' @export
prob_detect_dilution_2_betabinom <- function(S, V0, V1, V2, V3, C_cfu, n_sim, alpha, beta) {
  N <- (C_cfu / V1) * ((S + V0) / S)
  n1 <- matrix(extraDistr::rbbinom(n_sim, N, alpha, beta), nrow = n_sim, ncol = 1)
  pd <- matrix(NA, nrow = n_sim, ncol = 1)
  for (j in 1:n_sim) {
    pd[j, ] <- 1 - (beta(alpha, beta + n1[j, ])) / (beta(alpha, beta))
  }
  result <- mean(pd)
  return(result)
}
