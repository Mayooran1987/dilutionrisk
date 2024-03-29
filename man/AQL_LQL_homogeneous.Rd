% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/AQL_LQL_homogeneous.R
\name{AQL_LQL_homogeneous}
\alias{AQL_LQL_homogeneous}
\title{AQL and LQL estimations for given dilution schemes when diluted samples are collected from a homogeneous batch.}
\usage{
AQL_LQL_homogeneous(c,lambda_low,lambda_high,a,b,f,u,USL,n,type,alpha,beta,OC,n_sim)
}
\arguments{
\item{c}{acceptance number}

\item{lambda_low}{the lower value of the expected microbial count(\eqn{\lambda}) for use in the graphical display's x-axis.}

\item{lambda_high}{the upper value of the expected microbial count(\eqn{\lambda}) for use in the graphical display's x-axis.}

\item{a}{lower domain of the number of microbial count.}

\item{b}{upper domain of the number of microbial count.}

\item{f}{final dilution factor.}

\item{u}{amount put on the plate.}

\item{USL}{upper specification limit.}

\item{n}{number of samples which are used for inspection.}

\item{type}{what type of the results you would like to consider such as "theory" or "simulation" (default "theory").}

\item{alpha}{producer’s risk}

\item{beta}{consumer’s risk}

\item{OC}{if we need AQL and LQL displayed with the OC curve, set \code{OC = "TRUE"}; otherwise, the output only provides the estimated values.}

\item{n_sim}{number of simulations (large simulations provide more precise estimations).}
}
\value{
AQL and LQL  when diluted samples are collected from a homogeneous batch.
}
\description{
\code{\link{AQL_LQL_homogeneous}} provides estimated AQL and LQL values for given dilution schemes when samples are collected from a homogeneous batch.
}
\details{
\code{\link{AQL_LQL_homogeneous}} provides estimated AQL and LQL values for given dilution schemes when samples are collected from a homogeneous batch.
Acceptable Quality Level (AQL) is the acceptable or good quality level at which the probability of acceptance is kept at a high level,
which is associated with the producer's risk. Conversely, the limiting Quality Level (LQL) refers to the rejectable or poor quality level
at which the probability of acceptance is kept at a low level, which is associated with the consumer's risk.
}
\examples{
c <- 2
n <- 5
lambda_low <- 2
lambda_high <- 5000
a <- 0
b <- 300
f <- 0.01
u <- 0.1
USL <- 1000
alpha <- 0.05
beta <- 0.10
AQL_LQL_homogeneous(c, lambda_low,lambda_high, a, b, f, u, USL, n, type = "theory",
                    alpha, beta, OC = "FALSE")
AQL_LQL_homogeneous(c, lambda_low,lambda_high, a, b, f, u, USL, n, type = "theory",
                    alpha, beta, OC = "TRUE")
}
