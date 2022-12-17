##' \code{\link{prob_acceptance_homogeneous}} provides a probability of acceptance in the original sample when samples collected from a homogeneous batch.
##' @title Probability of acceptance estimation when diluted sample collected from a homogeneous batch.
##' @param c acceptance number
##' @param lambda the expected cell count (\eqn{\lambda}).
##' @param a lower domain of the number of cell counts.
##' @param b upper domain of the number of cell counts.
##' @param f final dilution factor.
##' @param u amount put on the plate.
##' @param USL upper specification limit.
##' @param n number of samples which are used for inspection.
##' @param type what type of the results you would like to consider such as "theory" or "simulation" (default "theory").
##' @param n_sim number of simulations (large simulations provide a more precise estimation).
##' @details \code{\link{prob_detection_homogeneous}} provides a probability of acceptance when samples collected from a homogeneous batch (this section will be updated later on).
##' @return Probability of acceptance when the diluted sample collected from a homogeneous batch.
##' @examples
##' c <- 2
##' lambda <- 2000
##' a <- 0
##' b <- 300
##' f <- 0.01
##' u <- 0.1
##' USL <- 1000
##' n <- 5
##' prob_acceptance_homogeneous(c, lambda, a, b, f, u, USL, n)
##' @usage  prob_acceptance_homogeneous(c, lambda, a, b, f, u, USL, n, type, n_sim)
##' @export
prob_acceptance_homogeneous <- function(c, lambda, a, b, f, u, USL, n, type = "theory", n_sim = NA){
  pd <- prob_detection_homogeneous(lambda, a, b, f, u, USL, type, n_sim)
  pa <- stats::pbinom(c, n, pd)
  return(pa)
}




