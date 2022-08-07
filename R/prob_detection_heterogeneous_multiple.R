##' \code{\link{prob_detection_heterogeneous_multiple}} provides a probability of detection for multiple dilution schemes in the original sample when samples collected from a heterogeneous batch.
##' @title Probability of detection estimation for multiple dilution schemes when diluted samples are collected from a heterogeneous batch.
##' @param meanlog the mean concentration (on the log scale).
##' @param sdlog the standard deviation of the normal distribution (on the log scale).
##' @param n_sim number of simulations (large simulations provide more precise estimations).
##' @param a lower domain of the number of cell counts.
##' @param b upper domain of the number of cell counts.
##' @param f vector of final dilution factor.
##' @param u vector of amount put on the plate.
##' @param USL upper specification limit.
##' @details \code{\link{prob_detection_heterogeneous_multiple}} provides a probability of detection when diluted samples are collected from a heterogeneous batch. We define the random variable \eqn{X_{i}} is the number of colonies on the \eqn{i^{th}} plate.
##' In practice, the acceptance for countable numbers of colonies on a plate must be between 30 and 300. Therefore, we can utilise bounded distributions to model the number of colonies on a plate. In the heterogeneous case, we employed truncated Poisson lognormal distribution to the model. (this section will be updated later on).
##' @return Probability of detection when samples collected from a heterogeneous batch.
##' @examples
##' meanlog <- 7
##' sdlog <- 0.2
##' a <- 0
##' b <- 300
##' f <- c(0.01,0.1,1)
##' u <- c(0.1,0.1,0.1)
##' USL <- 1000
##' n_sim <- 50000
##' prob_detection_heterogeneous_multiple(meanlog, sdlog, a, b, f, u, USL, n_sim)
##' @usage  prob_detection_heterogeneous_multiple(meanlog, sdlog, a, b, f, u, USL, n_sim)
##' @export
prob_detection_heterogeneous_multiple <- function(meanlog, sdlog, a, b, f, u, USL, n_sim){
  if (length(f) != length(u)) stop("please use equal length of f and u", call. = FALSE)
  sim1 <- matrix(NA, nrow =  1, ncol = length(f))
  for (i in 1:length(f)) {
    sim1[,i] <-   prob_detection_heterogeneous(meanlog, sdlog, a, b, f[i], u[i], USL, n_sim)
  }
  results <- as.matrix.data.frame(sim1)
  return(results)
}
