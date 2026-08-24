#' dilutionrisk: Risk Assessment for Dilution Testing
#'
#' Provides probability estimations and graphical displays for modelling and
#' assessment of risk based on aerobic plate count (APC) on diluted testing.
#'
#' @description
#' The \code{dilutionrisk} package implements statistical methods for
#' microbiological risk assessment in food safety, focusing on:
#'
#' \itemize{
#'   \item \strong{Probability of Detection (PD)}: Probability that the
#'     estimated microbial count exceeds the upper specification limit.
#'   \item \strong{Probability of Acceptance (PA)}: Probability that a batch
#'     is accepted based on sampling plans.
#'   \item \strong{Operating Characteristic (OC) Curves}: Visual comparison
#'     of different dilution schemes.
#'   \item \strong{Coefficient of Variation (CV)}: Variability estimation
#'     for different dilution schemes.
#'   \item \strong{True Concentration Estimation}: Estimating the true
#'     microbial concentration in the original sample.
#' }
#'
#' @details
#' The package supports two scenarios:
#'
#' \strong{Homogeneous Batches:}
#' Assumes the microbial distribution is uniform throughout the batch.
#' The count follows a truncated Poisson distribution.
#'
#' \strong{Heterogeneous Batches:}
#' Assumes the microbial distribution is non-uniform.
#' The count follows a truncated Poisson-lognormal distribution.
#'
#' Both scenarios incorporate dilution factors and plate amounts to model
#' the dilution process in microbiological testing.
#'
#' @section Key Functions:
#' \describe{
#'   \item{Probability Detection}{\code{\link{prob_detection_homogeneous}},
#'     \code{\link{prob_detection_heterogeneous}}}
#'   \item{Probability Acceptance}{\code{\link{prob_acceptance_homogeneous}},
#'     \code{\link{prob_acceptance_heterogeneous}}}
#'   \item{OC Curves}{\code{\link{OC_curves_homogeneous}},
#'     \code{\link{OC_curves_heterogeneous}}}
#'   \item{PD Curves}{\code{\link{pd_curves_homogeneous}},
#'     \code{\link{pd_curves_heterogeneous}}}
#'   \item{CV Estimation}{\code{\link{cv_homogeneous}},
#'     \code{\link{cv_heterogeneous}}}
#'   \item{True Concentration}{\code{\link{true_concentration_homogeneous}},
#'     \code{\link{true_concentration_heterogeneous}}}
#' }
#'
#' @references
#' \itemize{
#'   \item Gonzales-Barron, U.A., Pilão Cadavez, V.A., Butler, F., 2013.
#'     Statistical approaches for the design of sampling plans for
#'     microbiological monitoring of foods. In: Mathematical and Statistical
#'     Methods in Food Science and Technology. Wiley, pp. 363–384.
#'   \item Schothorst, M. van, Zwietering, M.H., Ross, T., Buchanan, R.L.,
#'     Cole, M.B., 2009. Relating microbiological criteria to food safety
#'     objectives and performance objectives. Food Control 20, 967–979.
#' }
#'
#' @keywords package
#' @name dilutionrisk
#' @keywords internal
"_PACKAGE"
#' @useDynLib dilutionrisk, .registration = TRUE
#' @importFrom stats pbinom rlnorm
#' @importFrom methods setLoadAction
#' @importFrom ggplot2 ggplot aes geom_line geom_vline annotate theme_classic
#' @importFrom ggplot2 xlab ylab theme element_text element_rect element_line
#' @importFrom ggplot2 scale_x_continuous sec_axis
#' @importFrom ggthemes scale_colour_colorblind
#' @importFrom reshape2 melt
#' @importFrom utils setTxtProgressBar txtProgressBar
NULL

#' Register Rcpp Export Functions
#'
#' Internal function to register Rcpp exported functions.
#'
#' @keywords internal
#' @noRd
.onLoad <- function(libname, pkgname) {
  # Rcpp registration happens automatically via Rcpp attributes
  # This is a placeholder if needed
  invisible()
}
