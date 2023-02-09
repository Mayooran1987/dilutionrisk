##' These functions provides true concentration level in the original sample when diluted samples collected from a heterogeneous batch.
##' @title True concentration level estimation when diluted sample collected from a heterogeneous batch.
##' @param mu the mean microbial count (on the log scale).
##' @param mu_low the lower value of the mean microbial count (\eqn{\mu}) for use in the graphical display's x-axis (on the log scale).
##' @param mu_high the upper value of the mean microbial count (\eqn{\mu}) for use in the graphical display's x-axis (on the log scale).
##' @param sd the standard deviation of the normal distribution (on the log scale).
##' @param a lower domain of the number of cell counts.
##' @param b upper domain of the number of cell counts.
##' @param f final dilution factor.
##' @param u amount put on the plate.
##' @param USL upper specification limit.
##' @param n_sim number of simulations (large simulations provide a more precise estimation).
##' @details Let Y be the count of microorganisms and C be the true concentration level (in counts per ml).
##' When diluted sample collected from heterogeneous (non-homogeneous) batch, Y can be modelled by Poisson lognormal distribution
##' with parameter \eqn{\mu, \sigma}. Let X be the count of microorganisms on a plate,
##' and it can be modelled by truncated Poisson lognormal distribution with parameters \eqn{\mu_d,\sigma, a, b}.
##' Also, \eqn{\lambda_d} can be written in terms of \eqn{\mu},f and u. It is given by
##' \deqn{\mu_d = \mu + \log(f)+ \log(u)}
##' And the true concentration level  is given  by
##' \deqn{C = \frac{X}{f*u}}
##' where \eqn{f} is final dilution factor and \eqn{u} is amount of diluted sample on plate. Based on the literatures, we used \eqn{\sigma = 0.2} in these dilution process; see Gonzales-Barron et al. (2013, p. 370) and Schothorst et al. (2009).
##' @references
##' \itemize{
##' \item Gonzales-Barron, U.A., Pilão Cadavez, V.A., Butler, F., 2013. Statistical approaches for the design of sampling plans for microbiological monitoring of foods, in:
##' Mathematical and Statistical Methods in Food Science and Technology. Wiley, Chichester, UK, pp.363–384.
##' \item Schothorst, M. van, Zwietering, M.H., Ross, T., Buchanan, R.L., Cole, M.B., 2009. Relating microbiological criteria to food safety objectives and performance objectives. Food Control 20, \href{https://doi.org/10.1016/j.foodcont.2008.11.005}{967–979}.
##' }
##' @return true concentration level  when sample collected from a heterogeneous batch.
##' @examples
##' mu_low <- 0
##' mu_high <- 10
##' sd <- 0.2
##' a <- 0
##' b <- 300
##' f <- c(0.01,0.1)
##' u <- c(0.1,0.1)
##' USL <- 1000
##' n_sim <- 5000
##' true_concentration_curves_heterogeneous(mu_low, mu_high, sd, a, b, f, u, USL, n_sim)
##' @name true_concentration_heterogeneous
##' @aliases true_concentration_heterogeneous
##' @rdname true_concentration_heterogeneous
##' @export
true_concentration_heterogeneous <- function(mu, sd, a, b, f, u, USL, n_sim) {
  sim1 <- matrix(NA, nrow = n_sim, ncol = 1)
  for (j in 1:n_sim) {
    sim1[j, ] <- (rtrunpoilog(1, (mu + log(f * u, exp(1))), sd, a, b)) * (1 / (f * u))
  }
  mean_concentr <- apply(sim1, 2, mean)
  return(mean_concentr)
}

##' @rdname true_concentration_heterogeneous
##' @export
true_concentration_heterogeneous_multiple <- function(mu, sd, a, b, f, u, USL, n_sim) {
  if (length(f) != length(u)) stop("please use equal length of f and u", call. = FALSE)
  sim1 <- matrix(NA, nrow = 1, ncol = length(f))
  for (i in 1:length(f)) {
    sim1[, i] <- true_concentration_heterogeneous(mu, sd, a, b, f[i], u[i], USL, n_sim)
  }
  results <- as.matrix.data.frame(sim1)
  return(results)
}

##' @rdname true_concentration_heterogeneous
##' @export
true_concentration_curves_heterogeneous <- function(mu_low, mu_high, sd, a, b, f, u, USL, n_sim) {
  C <- NULL
  Dilution_scheme <- NULL
  f_spr <- function(f, u) {
    sprintf("Scheme (f=%.3f, u=%.1f)", f, u)
  }
  mu <- seq(mu_low, mu_high, 0.1)
  C <- matrix(NA, nrow = length(mu), ncol = length(f))
  for (i in 1:length(mu)) {
    C[i, ] <- cbind(true_concentration_heterogeneous_multiple(mu[i], sd, a, b, f, u, USL, n_sim))
  }
  Prob <- data.frame(mu, C)
  colnames(Prob) <- c("mu", f_spr(f, u))
  melten.Prob <- reshape2::melt(Prob, id = "mu", variable.name = "Dilution_scheme", value.name = "C")
  plot_sam <- ggplot2::ggplot(melten.Prob) +
    ggplot2::geom_line(ggplot2::aes(x = mu, y = C, group = Dilution_scheme, colour = Dilution_scheme)) +
    ggplot2::theme_classic() +
    ggplot2::xlab(expression(" log mean microbial count  (" ~ mu * ~")")) +
    ggplot2::ylab(expression("Estimated true concentration(" ~ widehat(C) ~ ")")) +
    ggthemes::scale_colour_colorblind() +
    ggplot2::theme(
      plot.title = ggplot2::element_text(hjust = 0.5), legend.position = c(0.25, 0.85), axis.line.x.top = ggplot2::element_line(color = "red"),
      axis.ticks.x.top = ggplot2::element_line(color = "red"), axis.text.x.top = ggplot2::element_text(color = "red"), axis.title.x.top = ggplot2::element_text(color = "red")
    ) +
    ggplot2::scale_x_continuous(sec.axis = ggplot2::sec_axis(~.,
      name = "mean microbial count", breaks = seq(min(mu), max(mu), 1),
      labels = c(sprintf("%0.2f", exp(seq(min(mu), max(mu), 1))))
    ))
  # plot_sam
  return(plot_sam)
}
