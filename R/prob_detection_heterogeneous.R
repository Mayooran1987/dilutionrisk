##' \code{\link{prob_detection_heterogeneous}} provides a probability of detection in the original sample when samples collected from a heterogeneous batch.
##' @title Probability of detection estimation when diluted sample collected from a heterogeneous batch.
##' @param mu the mean microbial count (on the log scale).
##' @param sd the standard deviation of the normal distribution (on the log scale).
##' @param a lower domain of the number of cell counts.
##' @param b upper domain of the number of cell counts.
##' @param f final dilution factor.
##' @param u amount put on the plate.
##' @param USL upper specification limit.
##' @param type what type of the results you would like to consider such as "theory" or "simulation" (default "theory").
##' @param n_sim number of simulations (large simulations provide a more precise estimation).
##' @details \code{\link{prob_detection_heterogeneous}} provides a probability of detection when the diluted sample has heterogeneous contaminants. We define the random variable \eqn{X_{i}} is the number of colonies on the \eqn{i^{th}} plate.
##' In practice, the acceptance for countable numbers of colonies on a plate must be between 30 and 300. Therefore, we can utilise bounded distributions to model the number of colonies on a plate. In the heterogeneous case, we employed truncated Poisson lognormal distribution to model (this section will be updated later on).
##' @return Probability of detection when sample collected from a heterogeneous batch.
##' @examples
##' mu <- 2
##' sd <- 0.2
##' a <- 0
##' b <- 300
##' f <- 0.01
##' u <- 0.1
##' USL <- 1000
##' n_sim <- 50000
##' prob_detection_heterogeneous(mu, sd, a, b, f, u, USL)
##' @usage  prob_detection_heterogeneous(mu, sd, a, b, f, u, USL, type, n_sim)
##' @export
prob_detection_heterogeneous <- function(mu, sd, a, b, f, u, USL, type = "theory", n_sim = NA){
  if (mu > (log(b,exp(1)) - log((f*u),exp(1))) | mu < log(a,exp(1)) - log(f*u,exp(1)))
    stop("The truncated poisson lognormal (TPLN) random variable must be bounded by a and b,
         which means that the mu must be less than or equal to log(b)-log(fu) and mu must be greater than or equal to log(a)-log(fu).")
  if (type == "theory") {
    pd_theory_heterogeneous <- function(mu, sd, f, u, USL){
      USL1 <- USL*f*u
      pd <- matrix(NA, nrow = USL1, ncol = 1)
      for (i in 1:USL1) {
        mu_d <- mu + log(f*u, exp(1))
        # pd[i,1] <- integrate(function(t) exp(t*i - 0.5 * ((t - mu_d)/sd)^2)/(exp(exp(t))-1), lower = 0, upper = 500)$value/(factorial(i)*sqrt(2 * pi) * sd)
        pd[i,1] <- (integrate(function(t) t^(i - 1)*exp( -0.5 * ((log(t,exp(1)) - mu_d)/sd)^2)/(exp(t) - 1), lower = 0, upper = Inf )$value)/(exp(lgamma(i + 1))*sqrt(2 * pi) * sd)

      }
      P_d <- 1 - sum(pd)
      return(P_d)
    }
    result <- pd_theory_heterogeneous(mu, sd, f, u, USL)
    return(result)
  } else if (type == "simulation") {
    if (is.na(n_sim) == TRUE) {
      stop("please set the number of simualtions")
    } else {
      sim1 <- matrix(NA, nrow =  n_sim, ncol = 1)
      # lambda <- 10^(mu + (sd^2/2) * log(10, exp(1)))
      for (j in 1:n_sim) {
        sim1[j,] <-   (rtrunpoilog(1, (mu + log(f*u,exp(1))), sd, a, b))*(1/(f*u))
        # sim1[j,] <-   (rtrunpoilog(1, (mean_con * f*u), sd, a, b))*(1/(f*u))
      }
      result <- length(which(sim1[,1] > USL))/n_sim
      return(result)
    }
  }else {
    print("please include what type (theory/ simulation) you would like to consider")
  }
}
