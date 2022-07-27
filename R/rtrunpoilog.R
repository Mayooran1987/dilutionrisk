##' \code{\link{rtrunpoilog}} provides generated random numbers from truncated poisson lognormal distribution with given parameters.
##' @title Generates random deviates from truncated poisson lognormal distribution.
##' @param n number of observations. If length(n) > 1 then the length is taken to be the number required.
##' @param meanlog the mean concentration (on the log scale).
##' @param sdlog the standard deviation of the normal distribution (on the log scale).
##' @param a lower truncation points.
##' @param b upper truncation points (a < x <= b).
##' @details \code{\link{rtrunpoilog}} provides generated random numbers from truncated poisson lognormal distribution with given parameters. (this section will be updated later on).
##' @return  rtrunpoilog generates random numbers from truncated poisson lognormal distribution.
##' @examples
##' n <- 100
##' meanlog <-  0
##' sdlog <-  1
##' a <- 0
##' b <- 300
##' rtrunpoilog(n, meanlog, sdlog, a, b)
##' @usage  rtrunpoilog(n, meanlog, sdlog, a, b)
##' @export
rtrunpoilog <- function(n, meanlog = 0, sdlog = 1, a, b){
  lambda <- stats::rlnorm(n = n, meanlog = meanlog, sdlog = sdlog)
  result <- extraDistr::rtpois(n = n, lambda = lambda,a,b)
  return(result)
}

