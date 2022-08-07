##' \code{\link{prob_detection_heterogeneous}} provides a probability of detection in the original sample when samples collected from a heterogeneous batch.
##' @title Probability of detection estimation when diluted sample collected from a heterogeneous batch.
##' @param meanlog the mean concentration (on the log scale).
##' @param sdlog the standard deviation of the normal distribution (on the log scale).
##' @param a lower domain of the number of cell counts.
##' @param b upper domain of the number of cell counts.
##' @param f final dilution factor.
##' @param u amount put on the plate.
##' @param USL upper specification limit.
##' @param n_sim number of simulations (large simulations provide a more precise estimation).
##' @details \code{\link{prob_detection_heterogeneous}} provides a probability of detection when the diluted sample has heterogeneous contaminants. We define the random variable \eqn{X_{i}} is the number of colonies on the \eqn{i^{th}} plate.
##' In practice, the acceptance for countable numbers of colonies on a plate must be between 30 and 300. Therefore, we can utilise bounded distributions to model the number of colonies on a plate. In the heterogeneous case, we employed truncated Poisson lognormal distribution to model (this section will be updated later on).
##' @return Probability of detection when sample collected from a heterogeneous batch.
##' @examples
##' meanlog <- 2
##' sdlog <- 0.2
##' a <- 0
##' b <- 300
##' f <- 0.01
##' u <- 0.1
##' USL <- 1000
##' n_sim <- 50000
##' prob_detection_heterogeneous(meanlog, sdlog, a, b, f, u, USL, n_sim)
##' @usage  prob_detection_heterogeneous(meanlog, sdlog, a, b, f, u, USL, n_sim)
##' @export
prob_detection_heterogeneous <- function(meanlog, sdlog, a, b, f, u, USL, n_sim){
  if (meanlog > (log(b,exp(1)) - log((f*u),exp(1))) | meanlog < log(a,exp(1)) - log(f*u,exp(1)))
    stop("The truncated poisson lognormal (TPLN) random variable must be bounded by a and b,
         which means that the meanlog must be less than or equal to log(b)-log(fu) and meanlog must be greater than or equal to log(a)-log(fu).")
  sim1 <- matrix(NA, nrow =  n_sim, ncol = 1)
# lambda <- 10^(mu + (sd^2/2) * log(10, exp(1)))
  for (j in 1:n_sim) {
    sim1[j,] <-   (rtrunpoilog(1, (meanlog + log(f*u,exp(1))), sdlog, a, b))*(1/(f*u))
    # sim1[j,] <-   (rtrunpoilog(1, (mean_con * f*u), sdlog, a, b))*(1/(f*u))
}
  pd <- length(which(sim1[,1] > USL))/n_sim
return(pd)
}

