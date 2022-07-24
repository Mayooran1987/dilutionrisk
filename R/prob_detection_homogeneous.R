##' \code{\link{prob_detection_homogeneous}} provides a probability of detection in the original sample when homogeneous contamination.
##' @title Probability of detection estimation when diluted sample has homogeneous contaminations.
##' @param lambda the expected cell count (\eqn{\lambda}).
##' @param n_sim number of simulations (large simulations provide more precise estimation).
##' @param a lower domain of the number of cell counts.
##' @param b upper domain of the number of cell counts.
##' @param pf plating factor (pf = 1/final dilution factor).
##' @param USL upper specification limit.
##' @details \code{\link{prob_detection_homogeneous}} provides a probability of detection when the diluted sample has homogeneous contaminants. We define the random variable \eqn{X_{i}} is the number of colonies on the \eqn{i^{th}} plate.
##' In practice, the acceptance for countable numbers of colonies on a plate must be between 30 and 300. Therefore, we can utilise bounded distributions to model the number of colonies on a plate. Homogeneous case, we employed truncated Poisson distribution to model (this section will be updated later on).
##' @return Probability of detection when the diluted sample has homogeneous contaminants.
##' @examples
##' lambda <- 20
##' a <- 0
##' b <- 300
##' pf <- 1000
##' USL <- 1000
##' n_sim <- 10000
##' prob_detection_homogeneous(lambda, a, b, pf, USL, n_sim)
##' @usage  prob_detection_homogeneous(lambda, a, b, pf, USL, n_sim)
##' @export
prob_detection_homogeneous <- function(lambda, a, b, pf, USL, n_sim){
  sim1 <- matrix(NA, nrow =  n_sim, ncol = pf)
  # lambda <- 10^(mu + (sd^2/2) * log(10, exp(1)))
  for (j in 1:n_sim) {
    # options(width = 80)
    # extra <- nchar('||100%')
    # width <- options()$width
    # step <- round((j / n_sim) * (width - extra))
    # text <- sprintf('|%s%s|% 3s%%', strrep('-', step),
    #                 strrep(' ', width - step - extra), round((j / n_sim) * 100))
    # cat(text)
    # sim1[j,] <-   (truncnorm::rtruncnorm(f, a, b, mu, sd)) * f
    sim1[j,] <-   (extraDistr::rtpois(pf, lambda, a, b))*pf
    # if (j == n_sim) cat(': simulations done!')
    # else cat('\014')
  }
  # pd <- length(which(sim1[,1] > USL))/100
  sim2 <- matrix(NA, nrow =  n_sim, ncol = 1)
  for (j in 1:n_sim) {
    sim2[j,] <-   length(which(sim1[j,] > USL))/pf
  }
  pd <- apply(sim2, 2, mean)
  return(pd)
}
