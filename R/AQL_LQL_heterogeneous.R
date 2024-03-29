##' \code{\link{AQL_LQL_heterogeneous}} provides estimated AQL and LQL values for given dilution schemes when samples are collected from a heterogeneous batch.
##' @title AQL and LQL estimations for given dilution schemes when diluted samples are collected from a heterogeneous batch.
##' @param c acceptance number
##' @param mu_low the lower value of the mean microbial count(\eqn{\mu}) for use in the graphical display's x-axis.
##' @param mu_high the upper value of the mean microbial count(\eqn{\mu}) for use in the graphical display's x-axis.
##' @param sd the standard deviation of the normal distribution (on the log scale).
##' @param a lower domain of the number of microbial count.
##' @param b upper domain of the number of microbial count.
##' @param f final dilution factor.
##' @param u amount put on the plate.
##' @param USL upper specification limit.
##' @param n number of samples which are used for inspection.
##' @param type what type of the results you would like to consider such as "theory" or "simulation" (default "theory").
##' @param n_sim number of simulations (large simulations provide more precise estimations).
##' @param alpha producer’s risk
##' @param beta consumer’s risk
##' @param OC if we need AQL and LQL displayed with the OC curve, set \code{OC = "TRUE"}; otherwise, the output only provides the estimated values.
##' @details \code{\link{AQL_LQL_heterogeneous}} provides estimated AQL and LQL values for given dilution schemes when samples are collected from a heterogeneous batch.
##' Acceptable Quality Level (AQL) is the acceptable or good quality level at which the probability of acceptance is kept at a high level,
##' which is associated with the producer's risk. Conversely, the limiting Quality Level (LQL) refers to the rejectable or poor quality level
##' at which the probability of acceptance is kept at a low level, which is associated with the consumer's risk.
##' @return AQL and LQL  when diluted samples are collected from a heterogeneous batch.
##' @examples
##' c <- 2
##' mu_low <- 4
##' mu_high <- 9
##' sd <- 0.2
##' a <- 0
##' b <- 300
##' f <- 0.01
##' u <- 0.1
##' USL <- 1000
##' alpha <- 0.05
##' beta <- 0.10
##' n <- 5
##' AQL_LQL_heterogeneous(c, mu_low, mu_high, sd, a, b, f, u, USL, n, type = "theory",
##'                       alpha, beta, OC= "TRUE")
##' AQL_LQL_heterogeneous(c, mu_low, mu_high, sd, a, b, f, u, USL, n, type = "theory",
##'                       alpha, beta, OC= "FALSE")
##' @usage  AQL_LQL_heterogeneous(c,mu_low,mu_high,sd,a,b,f,u,USL,n,type,alpha,beta,OC,n_sim)
##' @export
AQL_LQL_heterogeneous <- function(c, mu_low, mu_high, sd, a, b, f, u, USL, n, type = "theory", alpha, beta, OC, n_sim = NA) {
  mu <- seq(mu_low, mu_high, 0.001)
  Pa <- matrix(NA, nrow = length(mu), ncol = length(f))
  for (i in 1:length(mu)) {
    Pa[i, ] <- round(cbind(prob_acceptance_heterogeneous_multiple(c, mu[i], sd, a, b, f, u, USL, n, type, n_sim)), 4)
  }
  Prob_df <- data.frame(mu, Pa)
  # Prob_df <- data.frame(p, Pa)
  if (OC == TRUE) {
    plot_AQL_LQL <- ggplot2::ggplot(Prob_df) +
      ggplot2::geom_line(ggplot2::aes(x = mu, y = Pa)) +
      ggplot2::xlab(expression("log mean microbial count  (" ~ mu * ~")")) +
      ggplot2::ylab(expression("Probability of acceptance" ~ (P[a]))) +
      ggplot2::theme_classic() +
      ggplot2::geom_hline(yintercept = c(Pa[which(mu == mu[which(abs(Pa - (1 - alpha)) == min(abs(Pa - (1 - alpha))))])]), linetype = "dashed") +
      ggplot2::geom_vline(xintercept = c(mu[which(abs(Pa - (1 - alpha)) == min(abs(Pa - (1 - alpha))))]), linetype = "dashed") +
      ggplot2::geom_hline(yintercept = c(Pa[which(mu == mu[which(abs(Pa - beta) == min(abs(Pa - beta)))])]), linetype = "dashed") +
      ggplot2::geom_vline(xintercept = c(mu[which(abs(Pa - beta) == min(abs(Pa - beta)))]), linetype = "dashed") +
      # ggplot2::geom_segment(ggplot2::aes(x=0,xend=400,y=1-alpha,yend=1-alpha),linetype = "dashed") +
      # ggplot2::geom_segment(ggplot2::aes(x=400,xend=400,y=0,yend=1-alpha),linetype = "dashed") +
      # ggplot2::geom_segment(ggplot2::aes(x=0,xend=2356,y=beta,yend=beta),linetype = "dashed") +
      # ggplot2::geom_segment(ggplot2::aes(x=2356,xend=2356,y=0,yend=beta),linetype = "dashed") +
      ggplot2::annotate("text",
        x = c(mu[which(abs(Pa - (1 - alpha)) == min(abs(Pa - (1 - alpha))))], mu[which(abs(Pa - beta) == min(abs(Pa - beta)))]),
        y = c(0, 0),
        label = sprintf(
          c("\n AQL = %0.4f", "\n LQL = %0.4f"),
          c(mu[which(abs(Pa - (1 - alpha)) == min(abs(Pa - (1 - alpha))))], mu[which(abs(Pa - beta) == min(abs(Pa - beta)))])
        ), size = 3
      ) +
      ggplot2::annotate("text",
        x = c(mu_low + 0.5, mu_low + 0.5),
        y = c(Pa[which(mu == mu[which(abs(Pa - (1 - alpha)) == min(abs(Pa - (1 - alpha))))])], Pa[which(mu == mu[which(abs(Pa - beta) == min(abs(Pa - beta)))])]),
        label = sprintf(c("\n 1-\u03B1 = %0.2f", "\n \u03B2 = %0.2f"), c(Pa[which(mu == mu[which(abs(Pa - (1 - alpha)) == min(abs(Pa - (1 - alpha))))])], Pa[which(mu == mu[which(abs(Pa - beta) == min(abs(Pa - beta)))])])), size = 4
      )

    results <- plot_AQL_LQL
  } else {
    # AQL <- mu[which(Pa == 1 - alpha)]
    # LQL <- mu[which(Pa == beta)]
    AQL <- mu[which(abs(Pa - (1 - alpha)) == min(abs(Pa - (1 - alpha))))]
    LQL <- mu[which(abs(Pa - beta) == min(abs(Pa - beta)))]

    results <- paste("AQL = ", AQL, ", LQL = ", LQL)
  }
  return(results)
}
