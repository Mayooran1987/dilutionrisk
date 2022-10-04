##' \code{\link{rtrunpoilog}} provides generated random numbers from truncated Poisson lognormal distribution with given parameters.
##' @title Generates random deviates from truncated Poisson lognormal distribution.
##' @param n number of observations. If length(n) > 1 then the length is taken to be the number required.
##' @param mu the mean microbial count (on the log scale).
##' @param sd the standard deviation of the normal distribution (on the log scale).
##' @param a lower truncation points (lower domain of the number of microbial count).
##' @param b upper truncation points (upper domain of the number of microbial count).
##' @details \code{\link{rtrunpoilog}} provides generated random numbers from truncated Poisson lognormal distribution with given parameters. (this section will be updated later on).
##' @return  \code{\link{rtrunpoilog}} generates random numbers from truncated Poisson lognormal distribution.
##' @examples
##' n <- 100
##' mu <-  0
##' sd <-  1
##' a <- 0
##' b <- 300
##' rtrunpoilog(n, mu, sd, a, b)
##' @usage  rtrunpoilog(n, mu, sd, a, b)
##' @export
rtrunpoilog <- function(n, mu, sd, a, b){
  if (mu > log(b,exp(1)) | mu < log(a,exp(1)))
    stop("The truncated Poisson lognormal (TPLN) random variable must be bounded by a and b, which means that the mu must be less than or equal to log(b) and mu must be greater than or equal to log(a)")
  lambda <- stats::rlnorm(n = n, meanlog = mu, sdlog = sd)
  rtpois <- function(n, lambda, a = -Inf, b = Inf) {
    if (length(n) > 1) n <- length(n)
    cpp_rtpois(n, lambda, lower = a, upper = b)
    # .Call('_dilutionrisk_cpp_rtpois', PACKAGE = 'dilutionrisk', n, lambda, lower = a, upper = b)
  }
  result <- rtpois(n = n, lambda, a, b)
  # sim1 <- matrix(NA, nrow =  n, ncol = 1)
  # for (j in 1:n) {
  #   sim1[j,] <-   extraDistr::rtpois(1, lambda = lambda1[j], a, b)
  # }
  # result <- as.numeric(sim1)
  return(result)
}

