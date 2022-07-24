##' \code{\link{prob_acceptance_homogeneous}} provides a probability of acceptance in the original sample when homogeneous contamination.
##' @title Probability of acceptance estimation when diluted sample has homogeneous contaminations.
##' @param c acceptance number
##' @param lambda the expected cell count (\eqn{\lambda}).
##' @param a lower domain of the number of cell counts.
##' @param b upper domain of the number of cell counts.
##' @param pf plating factor (pf = 1/final dilution factor).
##' @param USL upper specification limit.
##' @param n number of samples which are used for inspection.
##' @param n_sim number of simulations (large simulations provide more precise estimation).
##' @details \code{\link{prob_detection_homogeneous}} provides a probability of acceptance when the diluted sample has homogeneous contaminants (this section will be updated later on).
##' @return Probability of acceptance when the diluted sample has homogeneous contaminants.
##' @examples
##' c <- 0
##' lambda <- 10
##' a <- 0
##' b <- 300
##' pf <- 1000
##' USL <- 1000
##' n <- 5
##' n_sim <- 10000
##' prob_acceptance_homogeneous(c, lambda, a, b, pf, USL, n, n_sim)
##' @usage  prob_acceptance_homogeneous(c, lambda, a, b, pf, USL, n, n_sim)
##' @export
prob_acceptance_homogeneous <- function(c, lambda, a, b, pf, USL, n, n_sim){
  pd <- prob_detection_homogeneous(lambda, a, b, pf, USL, n_sim)
  pa <- stats::pbinom(c, n, pd)
  return(pa)
}

