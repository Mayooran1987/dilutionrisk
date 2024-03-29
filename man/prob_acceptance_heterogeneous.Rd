% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/prob_acceptance_heterogeneous.R
\name{prob_acceptance_heterogeneous}
\alias{prob_acceptance_heterogeneous}
\title{Probability of acceptance estimation when diluted sample collected from a heterogeneous batch.}
\usage{
prob_acceptance_heterogeneous(c, mu, sd, a, b, f, u, USL, n, type, n_sim)
}
\arguments{
\item{c}{acceptance number}

\item{mu}{the mean microbial count (on the log scale).}

\item{sd}{the standard deviation of the normal distribution (on the log scale).}

\item{a}{lower domain of the number of cell counts.}

\item{b}{upper domain of the number of cell counts.}

\item{f}{final dilution factor.}

\item{u}{amount put on the plate.}

\item{USL}{upper specification limit.}

\item{n}{number of samples which are used for inspection.}

\item{type}{what type of the results you would like to consider such as "theory" or "simulation" (default "theory").}

\item{n_sim}{number of simulations (large simulations provide a more precise estimation).}
}
\value{
Probability of acceptance when sample collected from a heterogeneous batch.
}
\description{
\code{\link{prob_acceptance_heterogeneous}} provides a probability of acceptance in the original sample when samples collected from a heterogeneous batch.
}
\details{
\code{\link{prob_detection_heterogeneous}} provides a probability of acceptance when diluted sample collected from a heterogeneous batch (this section will be updated later on).
}
\examples{
c <- 2
mu <- 7
sd <- 0.2
a <- 0
b <- 300
f <- 0.01
u <- 0.1
USL <- 1000
n <- 5
prob_acceptance_heterogeneous(c, mu, sd, a, b, f, u, USL, n)
}
