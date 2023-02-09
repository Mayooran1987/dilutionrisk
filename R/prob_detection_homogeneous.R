##' \code{\link{prob_detection_homogeneous}} provides a probability of detection in the original sample when samples collected from a homogeneous batch.
##' @title Probability of detection estimation when diluted sample collected from a homogeneous batch.
##' @param lambda the expected microbial count (\eqn{\lambda}).
##' @param a lower domain of the number of microbial count.
##' @param b upper domain of the number of microbial count.
##' @param f final dilution factor.
##' @param u amount put on the plate.
##' @param USL upper specification limit.
##' @param type what type of the results you would like to consider such as "theory" or "simulation" (default "theory").
##' @param n_sim number of simulations (large simulations provide a more precise estimation).
##' @details \code{\link{prob_detection_homogeneous}} provides a probability of detection when the diluted sample has homogeneous contaminants. We define the random variable \eqn{X_{i}} is the number of colonies on the \eqn{i^{th}} plate.
##' In practice, the acceptance for countable numbers of colonies on a plate must be between 30 and 300. Therefore, we can utilise bounded distributions to model the number of colonies on a plate. In the homogeneous case, we employed truncated Poisson distribution to model (this section will be updated later on).
##' @return Probability of detection when the diluted sample collected from a homogeneous batch.
##' @examples
##' lambda <- 2000
##' a <- 0
##' b <- 300
##' f <- 0.01
##' u <- 0.1
##' USL <- 1000
##' n_sim <- 50000
##' prob_detection_homogeneous(lambda, a, b, f, u, USL)
##' prob_detection_homogeneous(lambda, a, b, f, u, USL, type = "simulation", n_sim)
##' @usage  prob_detection_homogeneous(lambda, a, b, f, u, USL, type, n_sim)
##' @export
prob_detection_homogeneous <- function(lambda, a, b, f, u, USL, type = "theory", n_sim = NA) {
  if (type == "theory") {
    pd_theory <- function(lambda, f, u, USL) {
      USL1 <- USL * f * u
      pd <- matrix(NA, nrow = USL1, ncol = 1)
      for (i in 1:USL1) {
        pd[i, 1] <- (u * f * lambda)^i / (factorial(i) * (exp(u * f * lambda) - 1))
      }
      Pd <- 1 - sum(pd)
      return(Pd)
    }
    result <- pd_theory(lambda, f, u, USL)
    return(result)
  } else if (type == "simulation") {
    if (is.na(n_sim) == TRUE) {
      stop("please set the number of simualtions, include n_sim =")
    } else {
      rtpois <- function(n, lambda, a = -Inf, b = Inf) {
        if (length(n) > 1) n <- length(n)
        cpp_rtpois(n, lambda, lower = a, upper = b)
        # .Call('_dilutionrisk_cpp_rtpois', PACKAGE = 'dilutionrisk', n, lambda, lower = a, upper = b)
      }
      sim1 <- matrix(NA, nrow = n_sim, ncol = 1)
      # lambda <- 10^(mu + (sd^2/2) * log(10, exp(1)))
      for (j in 1:n_sim) {
        sim1[j, ] <- (rtpois(1, lambda * f * u, a, b)) * (1 / (f * u))
      }
      pd <- length(which(sim1[, 1] > USL)) / n_sim
      return(pd)
    }
  } else {
    print("please include what type (theory/ simulation) you would like to consider")
  }
}
