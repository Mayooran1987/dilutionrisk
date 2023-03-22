##' \code{\link{prob_detect_dilution_2_binom_pln}}  provides a probability of detection for a dilution scheme in the second dilution stage based on the binomial distribution, while the count of microorganisms is modelled by Poisson lognormal distribution.
##' @title Binomial distribution-based probability of detection for a dilution scheme in the second dilution stage testing.
##' @param S amount of sample (in grams) used for diluted solution preparation
##' @param sd the standard deviation of the normal distribution (on the log scale).
##' @param mu the mean microbial count (on the log scale).
##' @param V0 dilution volume in the first dilution stage testing.
##' @param V1 the volume of the diluted solution used for testing (which is equal to the plated amount in this study)
##' @param V2 dilution volume in the first dilution stage testing.
##' @param V3 the volume of the diluted solution used for testing (which is equal to the plated amount in this study)
##' @param n_sim number of simulations (large simulations provide a more precise estimation).
##' @details \code{\link{prob_detect_dilution_2_binom_pln}} provides a probability of detection for a dilution scheme in the second dilution stage based on the binomial distribution, while the count of microorganisms is modelled by Poisson lognormal distribution.(This section will be updated later on.)
##' @return Probability of detection for a dilution scheme in the second dilution stage testing.
##' @examples
##' S <- 20
##' mu  <- 25
##' sd <- 0.8
##' V0 <- 100
##' V1 <- 1
##' V2 <- 10
##' V3 <- 1
##' n_sim <- 25000
##' prob_detect_dilution_2_binom_pln (S,mu,sd, V0, V1, V2, V3, n_sim)
##' @usage  prob_detect_dilution_2_binom_pln (S,mu,sd, V0, V1, V2, V3, n_sim)
##' @export
prob_detect_dilution_2_binom_pln <- function(S,mu,sd, V0, V1, V2, V3, n_sim) {
  rpoislog <- function(T, mu, sig, nu = 1, condS = FALSE, keep0 = FALSE){
    sim <- function(nr) {
      lamx <- rnorm(nr)
      x <- rpois(nr, exp(sig * lamx + mu + log(nu)))
      if (!keep0)
        x <- x[x > 0]
      return(x)
    }
    if (T < 1)
      stop("S is not positive")
    if (!is.finite(T))
      stop("T is not finite")
    if ((T/trunc(T)) != 1)
      stop("T is not an integer")
    if (sig < 0)
      stop("sig is not positive")
    if (nu < 0)
      stop("nu is not positive")
    if (condS) {
      simVec <- vector("numeric", 0)
      fac <- 2
      nr <- T
      while (length(simVec) < T) {
        simvals <- sim(nr * fac)
        simVec <- c(simVec, simvals)
        fac <- (1/(length(simvals)/(nr * fac))) * 2
        fac <- ifelse(is.finite(fac), fac, 1000)
        nr <- T - length(simvals)
      }
      simVec <- simVec[1:T]
    }
    else simVec <- sim(T)
    return(simVec)
  }
  # N <- (C_cfu / V1) * ((S + V0) / S)
  lambda <- 10^(mu + (sd^2/2) * log(10, exp(1)))
  N <- round(mean(rpoislog(50000, S*lambda, sd, keep0 = TRUE)))
  # N <- rpois(1,S*lambda)
  # number of successes must be at least one.
  rcond.binomial <- function(n_sim, N, p) {
    probs <- stats::dbinom(1:N, N, p) / (1 - p)
    sample(1:N, n_sim, replace = TRUE, prob = probs)
  }
  n1 <- matrix(100*rcond.binomial(n_sim, N, V1 / V0), nrow = n_sim, ncol = 1)
  pd <- matrix(NA, nrow = n_sim, ncol = 1)
  for (j in 1:n_sim) {
    pd[j, ] <- 1 - dbinom(0, n1[j, ], V3 / V2)
  }
  result <- mean(pd)
  return(result)
}
