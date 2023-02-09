##' \code{\link{prob_acceptance_heterogeneous}} provides a probability of acceptance in the original sample when samples collected from a heterogeneous batch.
##' @title Probability of acceptance estimation when diluted sample collected from a heterogeneous batch.
##' @param c acceptance number
##' @param mu the mean microbial count (on the log scale).
##' @param sd the standard deviation of the normal distribution (on the log scale).
##' @param a lower domain of the number of cell counts.
##' @param b upper domain of the number of cell counts.
##' @param f final dilution factor.
##' @param u amount put on the plate.
##' @param USL upper specification limit.
##' @param n number of samples which are used for inspection.
##' @param type what type of the results you would like to consider such as "theory" or "simulation" (default "theory").
##' @param n_sim number of simulations (large simulations provide a more precise estimation).
##' @details \code{\link{prob_detection_heterogeneous}} provides a probability of acceptance when diluted sample collected from a heterogeneous batch (this section will be updated later on).
##' @return Probability of acceptance when sample collected from a heterogeneous batch.
##' @examples
##' c <- 2
##' mu <- 7
##' sd <- 0.2
##' a <- 0
##' b <- 300
##' f <- 0.01
##' u <- 0.1
##' USL <- 1000
##' n <- 5
##' prob_acceptance_heterogeneous(c, mu, sd, a, b, f, u, USL, n)
##' @usage  prob_acceptance_heterogeneous(c, mu, sd, a, b, f, u, USL, n, type, n_sim)
##' @export
prob_acceptance_heterogeneous <- function(c, mu, sd, a, b, f, u, USL, n, type = "theory", n_sim = NA) {
  pd <- prob_detection_heterogeneous(mu, sd, a, b, f, u, USL, type, n_sim)
  pa <- stats::pbinom(c, n, pd)
  return(pa)
}
