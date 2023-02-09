##' These functions provides coefficient of variation in the original sample  when diluted samples collected from a homogeneous batch.
##' @title coefficient of variation estimation when diluted sample collected from a homogeneous batch.
##' @param lambda the expected microbial count (\eqn{\lambda}).
##' @param lambda_low the lower value of the expected microbial count (\eqn{\lambda}) for use in the graphical display's x-axis.
##' @param lambda_high the upper value of the expected microbial count (\eqn{\lambda}) for use in the graphical display's x-axis.
##' @param a lower domain of the number of microbial count.
##' @param b upper domain of the number of microbial count.
##' @param f final dilution factor.
##' @param u amount put on the plate.
##' @param USL upper specification limit.
##' @param n_sim number of simulations (large simulations provide a more precise estimation).
##' @details These functions provides coefficient of variation in the original sample  when diluted samples collected from a homogeneous batch.
##' @return coefficient of variation when the diluted sample collected from a homogeneous batch.
##' @examples
##' lambda_low <- 1000
##' lambda_high <- 8000
##' a <- 0
##' b <- 300
##' f <- c(0.01,0.1)
##' u <- c(0.1,0.1)
##' USL <- 1000
##' n_sim <- 500
##' cv_curves_homogeneous(lambda_low, lambda_high, a, b, f, u, USL, n_sim)
##' @name cv_homogeneous
##' @aliases cv_homogeneous
##' @rdname cv_homogeneous
##' @export
cv_homogeneous <- function(lambda, a, b, f, u, USL, n_sim) {
  rtpois <- function(n, lambda, a = -Inf, b = Inf) {
    if (length(n) > 1) n <- length(n)
    cpp_rtpois(n, lambda, lower = a, upper = b)
  }
  sim1 <- matrix(NA, nrow = n_sim, ncol = 1 / (f * u))
  # lambda <- 10^(mu + (sd^2/2) * log(10, exp(1)))
  for (j in 1:n_sim) {
    sim1[j, ] <- (rtpois(1 / (f * u), lambda * f * u, a, b)) * (1 / (f * u))
    # sim1[j,] <-   (rtrunpoilog(1, (mean_con * f*u), sd, a, b))*(1/(f*u))
  }
  sim2 <- apply(sim1, 2, mean)
  cv <- sqrt(var(sim2)) / mean(sim2)
  # cv <- goeveg::cv(sim2, na.rm = FALSE)
  return(cv)
}

# cv_homogeneous(lambda, a, b, f, u, USL, n_sim)
##' @rdname cv_homogeneous
##' @export
cv_homogeneous_multiple <- function(lambda, a, b, f, u, USL, n_sim) {
  if (length(f) != length(u)) stop("please use equal length of f and u", call. = FALSE)
  sim1 <- matrix(NA, nrow = 1, ncol = length(f))
  for (i in 1:length(f)) {
    sim1[, i] <- cv_homogeneous(lambda, a, b, f[i], u[i], USL, n_sim)
  }
  results <- as.matrix.data.frame(sim1)
  return(results)
}
##' @rdname cv_homogeneous
##' @export
cv_curves_homogeneous <- function(lambda_low, lambda_high, a, b, f, u, USL, n_sim) {
  cv <- NULL
  Dilution_scheme <- NULL
  f_spr <- function(f, u) {
    sprintf("Scheme (f=%.3f, u=%.1f)", f, u)
  }
  lambda <- seq(lambda_low, lambda_high, 0.1)
  # lambda <- 10^(mu + (sd^2/2) * log(10, exp(1)))
  Pd <- matrix(NA, nrow = length(lambda), ncol = length(f))
  for (i in 1:length(lambda)) {
    Pd[i, ] <- cbind(cv_homogeneous_multiple(lambda[i], a, b, f, u, USL, n_sim))
  }
  Prob <- data.frame(lambda, Pd)
  colnames(Prob) <- c("lambda", f_spr(f, u))
  melten.Prob <- reshape2::melt(Prob, id = "lambda", variable.name = "Dilution_scheme", value.name = "cv")
  plot_sam <- ggplot2::ggplot(melten.Prob) +
    ggplot2::geom_line(ggplot2::aes(x = lambda, y = cv, group = Dilution_scheme, colour = Dilution_scheme)) +
    ggplot2::theme_classic() +
    ggplot2::xlab(expression("expected microbial count  (" ~ lambda * ~")")) +
    ggplot2::ylab(expression("coefficient of variation" ~ (CV))) +
    ggthemes::scale_colour_colorblind() +
    ggplot2::theme(
      plot.title = ggplot2::element_text(hjust = 0.5), legend.position = c(0.85, 0.25), axis.line.x.top = ggplot2::element_line(color = "red"),
      axis.ticks.x.top = ggplot2::element_line(color = "red"), axis.text.x.top = ggplot2::element_text(color = "red"), axis.title.x.top = ggplot2::element_text(color = "red")
    )
  return(plot_sam)
}
