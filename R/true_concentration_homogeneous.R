##' These functions provides true concentration level in the original sample  when diluted samples collected from a homogeneous batch.
##' @title True concentration level estimation when diluted sample collected from a homogeneous batch.
##' @param lambda the expected microbial count (\eqn{\lambda}).
##' @param lambda_low the lower value of the expected microbial count (\eqn{\lambda}) for use in the graphical display's x-axis.
##' @param lambda_high the upper value of the expected microbial count (\eqn{\lambda}) for use in the graphical display's x-axis.
##' @param a lower domain of the number of microbial count.
##' @param b upper domain of the number of microbial count.
##' @param f final dilution factor.
##' @param u amount put on the plate.
##' @param USL upper specification limit.
##' @param n_sim number of simulations (large simulations provide a more precise estimation).
##' @details Let Y be the count of microorganisms and C be the true concentration level (in counts per ml).
##' When diluted sample collected from homogeneous batch, Y can be modelled by Poisson distribution
##' with parameter \eqn{\lambda}. Let X be the count of microorganisms on a plate,
##' and it can be modelled by truncated Poisson distribution with parameters \eqn{\lambda_d, a, b}.
##' Also, \eqn{\lambda_d} can be written in terms of \eqn{\lambda},f and u. It is given by
##'\deqn{\lambda_d = \lambda * f *u}
##' And the true concentration level  is given  by
##' \deqn{C = \frac{X}{f*u}}
##' @return true concentration level when the diluted sample collected from a homogeneous batch.
##' @examples
##' lambda_low <- 0
##' lambda_high <- 5000
##' a <- 0
##' b <- 300
##' f <- c(0.01,0.1)
##' u <- c(0.1,0.1)
##' USL <- 1000
##' n_sim <- 50000
##' true_concentration_curves_homogeneous(lambda_low, lambda_high, a, b, f, u, USL, n_sim)
##' @name true_concentration_homogeneous
##' @aliases true_concentration_homogeneous
##' @rdname true_concentration_homogeneous
##' @export
true_concentration_homogeneous <- function(lambda, a, b, f, u, USL, n_sim){
  rtpois <- function(n, lambda, a = -Inf, b = Inf) {
    if (length(n) > 1) n <- length(n)
    cpp_rtpois(n, lambda, lower = a, upper = b)
  }
  sim1 <- matrix(NA, nrow =  n_sim, ncol = 1)
  for (j in 1:n_sim) {
    sim1[j,] <-   (rtpois(1, lambda*f*u, a, b))*(1/(f*u))
  }
  mean_concentr <- apply(sim1,2,mean)
  return(mean_concentr)
}

##' @rdname true_concentration_homogeneous
##' @export
true_concentration_homogeneous_multiple <- function(lambda, a, b, f, u, USL, n_sim){
  if (length(f) != length(u)) stop("please use equal length of f and u", call. = FALSE)
  pd <- matrix(NA, nrow =  1, ncol = length(f))
  for (i in 1:length(f)) {
    pd[,i] <-   true_concentration_homogeneous(lambda, a, b, f[i], u[i], USL, n_sim)
  }
  results <- as.matrix.data.frame(pd)
  return(results)
}

##' @rdname true_concentration_homogeneous
##' @export
true_concentration_curves_homogeneous <- function(lambda_low, lambda_high, a, b, f, u, USL, n_sim){
  C <- NULL
  Dilution_scheme <- NULL
  f_spr <- function(f, u) {
    sprintf("Scheme (f=%.3f, u=%.1f)", f, u)
  }
  lambda <- seq(lambda_low, lambda_high, 0.1)
  C <- matrix(NA, nrow = length(lambda), ncol = length(f))
  for (i in 1:length(lambda)) {
    C[i,] <-  cbind(true_concentration_homogeneous_multiple(lambda[i], a, b, f, u, USL, n_sim))
  }
  Prob <- data.frame(lambda, C)
  colnames(Prob ) <- c("lambda", f_spr(f,u))
  melten.Prob <- reshape2::melt(Prob, id = "lambda", variable.name = "Dilution_scheme", value.name = "C")
  plot_sam <- ggplot2::ggplot(melten.Prob) + ggplot2::geom_line(ggplot2::aes(x = lambda, y = C, group = Dilution_scheme, colour = Dilution_scheme)) +
    ggplot2::theme_classic() + ggplot2::xlab(expression("expected microbial count (" ~ lambda*~")")) + ggplot2::ylab(expression("Estimated true concentration("~widehat(C)~")")) + ggthemes::scale_colour_colorblind() +
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5), legend.position = c(0.85, 0.25), axis.line.x.top = ggplot2::element_line(color = "red"),
                   axis.ticks.x.top = ggplot2::element_line(color = "red"), axis.text.x.top = ggplot2::element_text(color = "red"), axis.title.x.top = ggplot2::element_text(color = "red"))
  return(plot_sam)
}



