##' \code{\link{prob_acceptance_homogeneous_multiple}} provides a probability of acceptance for multiple dilution schemes in the original sample when samples collected from a homogeneous batch
##' @title Probability of acceptance estimation for multiple dilution schemes when diluted samples are collected from a homogeneous batch.
##' @param c acceptance number
##' @param lambda the expected microbial count (\eqn{\lambda}).
##' @param a lower domain of the number of microbial count.
##' @param b upper domain of the number of microbial count.
##' @param f final dilution factor.
##' @param u amount put on the plate.
##' @param USL upper specification limit.
##' @param n number of samples which are used for inspection.
##' @param n_sim number of simulations (large simulations provide more precise estimations).
##' @details \code{\link{prob_detection_homogeneous_multiple}} provides a probability of acceptance for multiple dilution schemes in the original sample when samples collected from a homogeneous batch (this section will be updated later on).
##' @return Probability of acceptance  when diluted samples are collected from a homogeneous batch.
##' @examples
##' c <- 2
##' lambda <- 1000
##' a <- 0
##' b <- 300
##' f <- c(0.01,0.1,1)
##' u <- c(0.1,0.1,0.1)
##' USL <- 1000
##' n <- 5
##' n_sim <- 50000
##' prob_acceptance_homogeneous_multiple(c, lambda, a, b, f, u, USL, n, n_sim)
##' @usage  prob_acceptance_homogeneous_multiple(c, lambda, a, b, f, u, USL, n, n_sim)
##' @export
prob_acceptance_homogeneous_multiple <- function(c, lambda, a, b, f, u, USL, n, n_sim){
pd <- NULL
if (length(f) != length(u)) stop("please use equal length of f and u", call. = FALSE)
pa <- matrix(NA, nrow =  1, ncol = length(f))
for (i in 1:length(f)) {
  pd[i] <-  prob_detection_homogeneous(lambda, a, b, f[i], u[i], USL, n_sim)
  pa[,i] <-  stats::pbinom(c, n, pd[i])
}
results <- as.matrix.data.frame(pa)
return(results)
}

