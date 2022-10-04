##' These functions provides coefficient of variation in the original sample when diluted samples collected from a heterogeneous batch.
##' @title coefficient of variation estimation when diluted sample collected from a heterogeneous batch.
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
##' @details These functions provides coefficient of variation in the original sample when diluted samples collected from a heterogeneous batch.
##' @return coefficient of variation  when sample collected from a heterogeneous batch.
##' @examples
##' mu_low <- -5
##' mu_high <- 10
##' sd <- 0.2
##' a <- 0
##' b <- 300
##' f <- c(0.01,0.1)
##' u <- c(0.1,0.1)
##' USL <- 1000
##' n_sim <- 50000
##' cv_curves_heterogeneous(mu_low, mu_high, sd, a, b, f, u, USL, n_sim)
##' @name cv_heterogeneous
##' @aliases cv_heterogeneous
##' @rdname cv_heterogeneous
##' @export
cv_heterogeneous <- function(mu, sd, a, b, f, u, USL, n_sim){
  if (mu > (log(b,exp(1)) - log((f*u),exp(1))) | mu < log(a,exp(1)) - log(f*u,exp(1)))
    stop("The truncated poisson lognormal (TPLN) random variable must be bounded by a and b,
         which means that the mu must be less than or equal to log(b)-log(fu) and mu must be greater than or equal to log(a)-log(fu).")
  sim1 <- matrix(NA, nrow =  n_sim, ncol = 1/(f*u))
  # lambda <- 10^(mu + (sd^2/2) * log(10, exp(1)))
  for (j in 1:n_sim) {
    sim1[j,] <-   (rtrunpoilog(1/(f*u), (mu + log(f*u,exp(1))), sd, a, b))*(1/(f*u))
    # sim1[j,] <-   (rtrunpoilog(1, (mean_con * f*u), sd, a, b))*(1/(f*u))
  }
  sim2 <- apply(sim1,2,mean)
  cv <- sqrt(var(sim2)) / mean(sim2)
  # cv <- goeveg::cv(sim2, na.rm = FALSE)
  return(cv)
}
# cv_heterogeneous(mu, sd, a, b, f, u, USL, n_sim)

##' @rdname cv_heterogeneous
##' @export
cv_heterogeneous_multiple <- function(mu, sd, a, b, f, u, USL, n_sim){
  if (length(f) != length(u)) stop("please use equal length of f and u", call. = FALSE)
  sim1 <- matrix(NA, nrow =  1, ncol = length(f))
  for (i in 1:length(f)) {
    sim1[,i] <-   cv_heterogeneous(mu, sd, a, b, f[i], u[i], USL, n_sim)
  }
  results <- as.matrix.data.frame(sim1)
  return(results)
}
##' @rdname cv_heterogeneous
##' @export
cv_curves_heterogeneous <- function(mu_low, mu_high, sd, a, b, f, u, USL, n_sim){
  cv <- NULL
  Dilution_scheme <- NULL
  f_spr <- function(f, u) {
    sprintf("Scheme (f=%.3f, u=%.1f)", f, u)
  }
  mu <- seq(mu_low, mu_high, 0.01)
  # lambda <- 10^(mu + (sd^2/2) * log(10, exp(1)))
  cv <- matrix(NA, nrow = length(mu), ncol = length(f))
  for (i in 1:length(mu)) {
    cv[i,] <-  cbind(cv_heterogeneous_multiple(mu[i], sd, a, b, f, u, USL, n_sim))
  }
  Prob <- data.frame(mu, cv)
  colnames(Prob ) <- c("mu", f_spr(f,u))
  melten.Prob <- reshape2::melt(Prob, id = "mu", variable.name = "Dilution_scheme", value.name = "cv")
  plot_sam <- ggplot2::ggplot(melten.Prob) + ggplot2::geom_line(ggplot2::aes(x = mu, y = cv, group = Dilution_scheme, colour = Dilution_scheme)) +
    # ggplot2::ggtitle("OC curve based on Lognormal distribution") +
    ggplot2::theme_classic() + ggplot2::xlab(expression(" log mean microbial count  (" ~ mu*~")")) + ggplot2::ylab(expression("coefficient of variation"~(CV))) + ggthemes::scale_colour_colorblind() +
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5), legend.position = c(0.25, 0.85), axis.line.x.top = ggplot2::element_line(color = "red"),
                   axis.ticks.x.top = ggplot2::element_line(color = "red"), axis.text.x.top = ggplot2::element_text(color = "red"), axis.title.x.top = ggplot2::element_text(color = "red")) +
    ggplot2::scale_x_continuous(sec.axis = ggplot2::sec_axis(~., name = "mean microbial count", breaks = seq(min(mu),max(mu),1),
                                                             labels = c(sprintf("%0.2f", exp(seq(min(mu),max(mu),1))))))
  # plot_sam
  return(plot_sam)
}



