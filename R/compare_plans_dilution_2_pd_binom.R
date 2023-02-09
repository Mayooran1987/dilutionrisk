##' \code{\link{compare_plans_dilution_2_pd_binom}} provides graphical displays of the probability of the detection curves for dilution schemes in the first dilution stage.
##' @title Comparison based on probability of detection curves for different dilution schemes in the first dilution stage.
##' @param S amount of sample (in grams) used for diluted solution preparation
##' @param V0 dilution volume in the first dilution stage testing.
##' @param V1 the volume of the diluted solution used for testing (which is equal to the plated amount in this study)
##' @param V2 dilution volume in the first dilution stage testing.
##' @param V3 the volume of the diluted solution used for testing (which is equal to the plated amount in this study)
##' @param n_sim number of simulations (large simulations provide a more precise estimation).
##' @param Clim_lower an lower value of the number of colonies counted on the two dishes retained from two successive dilutions for comparison purposes.
##' @param Clim_upper an upper value of the number of colonies counted on the two dishes retained from two successive dilutions for comparison purposes.
##' @details \code{\link{compare_plans_dilution_2_pd_binom}} provides graphical displays of the probability of the detection curves for dilution schemes in the first dilution stage.
##' @return Comparison based on probability of detection curves for different dilution schemes in the first dilution stage.
##' @examples
##' S <- c(1,10,20)
##' V0 <- 100
##' V1 <- 1
##' V2 <- 10
##' V3 <- 1
##' Clim_lower <- 1
##' Clim_upper <- 1000
##' n_sim <- 20000
##' compare_plans_dilution_2_pd_binom(S, V0, V1, V2, V3, n_sim, Clim_lower, Clim_upper)
##' @usage  compare_plans_dilution_2_pd_binom(S, V0, V1, V2, V3, n_sim, Clim_lower, Clim_upper)
##' @export
compare_plans_dilution_2_pd_binom <- function(S, V0, V1, V2, V3, n_sim, Clim_lower, Clim_upper) {
  Sampling_scheme <- NULL # Initalizing
  p_d <- NULL
  C <- NULL
  C_cfu <- seq(Clim_lower, Clim_upper, by = 1)
  f_spr <- function(S) {
    sprintf("sample weight (%.0f gram)", S)
  }
  prob_detect_dilution_2_multi_binom <- function(S, V0, V1, V2, V3, C_cfu, n_sim) {
    pd <- matrix(NA, nrow = 1, ncol = length(C_cfu))
    for (j in 1:length(C_cfu)) {
      pd[1, j] <- prob_detect_dilution_2_binom(S, V0, V1, V2, V3, C_cfu[j], n_sim)
    }
    result <- as.numeric(pd)
    return(result)
  }
  pd <- matrix(NA, nrow = length(C_cfu), ncol = length(S))
  for (j in 1:length(S)) {
    pd[, j] <- prob_detect_dilution_2_multi_binom(S[j], V0, V1, V2, V3, C_cfu, n_sim)
  }
  Prob <- data.frame(C_cfu, pd)
  colnames(Prob) <- c("C", f_spr(S))
  melten.Prob <- reshape2::melt(Prob, id = "C", variable.name = "Sampling_scheme", value.name = "p_d")
  plot_sam <- ggplot2::ggplot(melten.Prob) +
    ggplot2::geom_line(ggplot2::aes(x = C, y = p_d, group = Sampling_scheme, colour = Sampling_scheme)) +
    ggplot2::ylab(expression(P[D])) +
    ggplot2::xlab(expression(C[cfu])) +
    ggplot2::theme_classic() +
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5, size = 10), legend.position = c(0.75, 0.25)) +
    ggthemes::scale_colour_colorblind()
  return(plot_sam)
}
