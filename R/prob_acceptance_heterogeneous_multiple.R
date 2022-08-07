##' \code{\link{prob_acceptance_heterogeneous_multiple}} provides a probability of acceptance in the original sample when samples collected from a heterogeneous batch.
##' @title Probability of acceptance estimation when diluted samples are collected from a heterogeneous batch.
##' @param c acceptance number
##' @param meanlog the mean concentration (on the log scale).
##' @param sdlog the standard deviation of the normal distribution (on the log scale).
##' @param a lower domain of the number of cell counts.
##' @param b upper domain of the number of cell counts.
##' @param f final dilution factor.
##' @param u amount put on the plate.
##' @param USL upper specification limit.
##' @param n number of samples which are used for inspection.
##' @param n_sim number of simulations (large simulations provide more precise estimations).
##' @details \code{\link{prob_acceptance_heterogeneous_multiple}} provides a probability of acceptance when diluted samples are collected from a heterogeneous batch (this section will be updated later on).
##' @return Probability of acceptance when samples collected from a heterogeneous batch.
##' @examples
##' c <- 2
##' meanlog <- 7
##' sdlog <- 0.2
##' a <- 0
##' b <- 300
##' f <- c(0.01,0.1,1)
##' u <- c(0.1,0.1,0.1)
##' USL <- 1000
##' n <- 5
##' n_sim <- 50000
##' prob_acceptance_heterogeneous_multiple (c, meanlog, sdlog, a, b, f, u, USL, n, n_sim)
##' @usage  prob_acceptance_heterogeneous_multiple (c, meanlog, sdlog, a, b, f, u, USL, n, n_sim)
##' @export
prob_acceptance_heterogeneous_multiple <- function(c, meanlog, sdlog, a, b, f, u, USL, n, n_sim){
  sim2 <- NULL
  if (length(f) != length(u)) stop("please use equal length of f and u", call. = FALSE)
  sim1 <- matrix(NA, nrow =  1, ncol = length(f))
  for (i in 1:length(f)) {
    sim2[i] <-  prob_detection_heterogeneous_multiple(meanlog, sdlog, a, b, f[i], u[i], USL, n_sim)
    sim1[,i] <-  stats::pbinom(c, n, sim2[i])
  }
  results <- as.matrix.data.frame(sim1)
  return(results)
}


