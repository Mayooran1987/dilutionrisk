% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/prob_detection_homogeneous_multiple.R
\name{prob_detection_homogeneous_multiple}
\alias{prob_detection_homogeneous_multiple}
\title{Probability of detection estimation for multiple dilution schemes when diluted samples are collected from a homogeneous batch.}
\usage{
prob_detection_homogeneous_multiple(lambda, a, b, f, u, USL, type, n_sim)
}
\arguments{
\item{lambda}{the expected microbial count (\eqn{\lambda}).}

\item{a}{lower domain of the number of microbial count.}

\item{b}{upper domain of the number of microbial count.}

\item{f}{final dilution factor.}

\item{u}{amount put on the plate.}

\item{USL}{upper specification limit.}

\item{type}{what type of the results you would like to consider such as "theory" or "simulation" (default "theory").}

\item{n_sim}{number of simulations (large simulations provide more precise estimations).}
}
\value{
Probability of detection when diluted samples are collected from a homogeneous batch.
}
\description{
\code{\link{prob_detection_homogeneous_multiple}} provides a probability of detection for multiple dilution schemes in the original sample when samples collected from a homogeneous batch.
}
\details{
\code{\link{prob_detection_homogeneous_multiple}} provides a probability of detection when the diluted sample has homogeneous contaminants. We define the random variable \eqn{X_{i}} is the number of colonies on the \eqn{i^{th}} plate.
In practice, the acceptance for countable numbers of colonies on a plate must be between 30 and 300. Therefore, we can utilise bounded distributions to model the number of colonies on a plate. In the homogeneous case, we employed truncated Poisson distribution to model (this section will be updated later on).
}
\examples{
lambda <- 1000
a <- 0
b <- 300
f <- c(0.01,0.1,1)
u <- c(0.1,0.1,0.1)
USL <- 1000
n_sim <- 50000
prob_detection_homogeneous_multiple(lambda, a, b, f, u, USL)
}
