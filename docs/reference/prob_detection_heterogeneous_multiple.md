# Probability of detection estimation for multiple dilution schemes when diluted samples are collected from a heterogeneous batch.

`prob_detection_heterogeneous_multiple` provides a probability of
detection for multiple dilution schemes in the original sample when
samples collected from a heterogeneous batch.

## Usage

``` r
prob_detection_heterogeneous_multiple(mu, sd, a, b, f, u, USL, type , n_sim)
```

## Arguments

- mu:

  the mean microbial count (on the log scale).

- sd:

  the standard deviation of the normal distribution (on the log scale).

- a:

  lower domain of the number of cell counts.

- b:

  upper domain of the number of cell counts.

- f:

  vector of final dilution factor.

- u:

  vector of amount put on the plate.

- USL:

  upper specification limit.

- type:

  what type of the results you would like to consider such as "theory"
  or "simulation" (default "theory").

- n_sim:

  number of simulations (large simulations provide a more precise
  estimation).

## Value

Probability of detection when samples collected from a heterogeneous
batch.

## Details

`prob_detection_heterogeneous_multiple` provides a probability of
detection when diluted samples are collected from a heterogeneous batch.
We define the random variable \\X\_{i}\\ is the number of colonies on
the \\i^{th}\\ plate. In practice, the acceptance for countable numbers
of colonies on a plate must be between 30 and 300. Therefore, we can
utilise bounded distributions to model the number of colonies on a
plate. In the heterogeneous case, we employed truncated Poisson
lognormal distribution to the model. (this section will be updated later
on).

## Examples

``` r
mu <- 7
sd <- 0.2
a <- 0
b <- 300
f <- c(0.01,0.1,1)
u <- c(0.1,0.1,0.1)
USL <- 1000
n_sim <- 50000
prob_detection_heterogeneous_multiple(mu, sd, a, b, f, u, USL, n_sim)
#> Error in match.arg(type, c("theory", "simulation")): 'arg' must be NULL or a character vector
```
