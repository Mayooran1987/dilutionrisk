##' \code{\link{prob_acceptance_heterogeneous}} provides a probability of acceptance in the original sample when heterogeneous contamination.
##' @title Probability of acceptance estimation when diluted sample has heterogeneous contaminations.
##' @param c acceptance number
##' @param meanlog the mean concentration (on the log scale).
##' @param sdlog the standard deviation of the normal distribution (on the log scale).
##' @param a lower domain of the number of cell counts.
##' @param b upper domain of the number of cell counts.
##' @param FDF final dilution factor.
##' @param USL upper specification limit.
##' @param n number of samples which are used for inspection.
##' @param n_sim number of simulations (large simulations provide more precise estimation).
##' @details \code{\link{prob_detection_heterogeneous}} provides a probability of acceptance when the diluted sample has heterogeneous contaminants (this section will be updated later on).
##' @return Probability of acceptance when the diluted sample has heterogeneous contaminants.
##' @examples
##' c <- 0
##' meanlog <- -2
##' sdlog <- 0.8
##' a <- 0
##' b <- 300
##' FDF <- 0.001
##' USL <- 1000
##' n <- 5
##' n_sim <- 10000
##' prob_acceptance_heterogeneous(c, meanlog, sdlog, a, b, FDF, USL, n, n_sim)
##' @usage  prob_acceptance_heterogeneous(c, meanlog, sdlog, a, b, FDF, USL, n, n_sim)
##' @export
prob_acceptance_heterogeneous <- function(c, meanlog, sdlog = 0.8, a, b, FDF, USL, n, n_sim){
  pd <- prob_detection_heterogeneous(meanlog, sdlog, a, b, FDF, USL, n_sim)
  pa <- stats::pbinom(c, n, pd)
  return(pa)
}

