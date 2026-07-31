##' \code{\link{prob_detection_homogeneous_multiple}} provides a probability of detection for multiple dilution schemes in the original sample when samples collected from a homogeneous batch.
##' @title Probability of detection estimation for multiple dilution schemes when diluted samples are collected from a homogeneous batch.
##' @param lambda the expected microbial count (\eqn{\lambda}).
##' @param a lower domain of the number of microbial count.
##' @param b upper domain of the number of microbial count.
##' @param f final dilution factor.
##' @param u amount put on the plate.
##' @param USL upper specification limit.
##' @param type what type of the results you would like to consider such as "theory" or "simulation" (default "theory").
##' @param n_sim number of simulations (large simulations provide more precise estimations).
##' @details \code{\link{prob_detection_homogeneous_multiple}} provides a probability of detection when the diluted sample has homogeneous contaminants. We define the random variable \eqn{X_{i}} is the number of colonies on the \eqn{i^{th}} plate.
##' In practice, the acceptance for countable numbers of colonies on a plate must be between 30 and 300. Therefore, we can utilise bounded distributions to model the number of colonies on a plate. In the homogeneous case, we employed truncated Poisson distribution to model (this section will be updated later on).
##' @return Probability of detection when diluted samples are collected from a homogeneous batch.
##' @examples
##' lambda <- 1000
##' a <- 0
##' b <- 300
##' f <- c(0.01,0.1,1)
##' u <- c(0.1,0.1,0.1)
##' USL <- 1000
##' n_sim <- 50000
##' prob_detection_homogeneous_multiple(lambda, a, b, f, u, USL)
##' @usage  prob_detection_homogeneous_multiple(lambda, a, b, f, u, USL, type, n_sim)
##' @export
prob_detection_homogeneous_multiple <- function(lambda, a, b, f, u, USL, type = "theory", n_sim = NA){
  if (length(f) != length(u)) stop("please use equal length of f and u", call. = FALSE)
  pd <- matrix(NA, nrow =  1, ncol = length(f))
  for (i in 1:length(f)) {
    pd[,i] <-   prob_detection_homogeneous(lambda, a, b, f[i], u[i], USL, type, n_sim)
  }
  results <- as.matrix.data.frame(pd)
  return(results)
}
