% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/compare_plans_dilution_1_pd_betabinom_pln.R
\name{compare_plans_dilution_1_pd_betabinom_pln}
\alias{compare_plans_dilution_1_pd_betabinom_pln}
\title{Comparison based on beta binomial distribution-based  probability of detection curves for different dilution schemes in the first dilution stage.}
\usage{
compare_plans_dilution_1_pd_betabinom_pln(S, sd, V0, V1, mu_lower, mu_upper, alpha, beta)
}
\arguments{
\item{S}{amount of sample (in grams) used for diluted solution preparation}

\item{sd}{the standard deviation of the normal distribution (on the log scale).}

\item{V0}{dilution volume in the first dilution stage testing.}

\item{V1}{the volume of the diluted solution used for testing (which is equal to the plated amount in this study)}

\item{mu_lower}{the lower value of the mean microbial count(\eqn{\mu}) for use in the graphical display's x-axis.}

\item{mu_upper}{the upper value of the mean microbial count(\eqn{\mu}) for use in the graphical display's x-axis.}

\item{alpha}{non-negative parameters of the beta distribution.}

\item{beta}{non-negative parameters of the beta distribution.}
}
\value{
Comparison based on probability of detection curves for different dilution schemes in the first dilution stage.
}
\description{
\code{\link{compare_plans_dilution_1_pd_betabinom_pln}} provides graphical displays of the probability of the detection curves for dilution schemes in the first dilution stage based on the binomial distribution, while the count of microorganisms is modelled by Poisson lognormal distribution.
}
\details{
\code{\link{compare_plans_dilution_1_pd_betabinom_pln}} provides graphical displays of the probability of the detection curves for dilution schemes in the first dilution stage based on the binomial distribution, while the count of microorganisms is modelled by Poisson lognormal distribution. (This section will be updated later on.)
}
\examples{
S <- c(25,125,250)
sd <- 0.8
V0 <- 100
V1 <- 1
mu_lower <- -8
mu_upper <- 8
alpha <- 1
beta <- 5
compare_plans_dilution_1_pd_betabinom_pln(S, sd, V0, V1, mu_lower, mu_upper, alpha, beta)
}
