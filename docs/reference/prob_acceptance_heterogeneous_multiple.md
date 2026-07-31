# Probability of acceptance estimation when diluted samples are collected from a heterogeneous batch.

`prob_acceptance_heterogeneous_multiple` provides a probability of
acceptance in the original sample when samples collected from a
heterogeneous batch.

## Usage

``` r
prob_acceptance_heterogeneous_multiple (c, mu, sd, a, b, f, u, USL, n, type, n_sim)
```

## Arguments

- c:

  acceptance number

- mu:

  the mean microbial count (on the log scale).

- sd:

  the standard deviation of the normal distribution (on the log scale).

- a:

  lower domain of the number of cell counts.

- b:

  upper domain of the number of cell counts.

- f:

  final dilution factor.

- u:

  amount put on the plate.

- USL:

  upper specification limit.

- n:

  number of samples which are used for inspection.

- type:

  what type of the results you would like to consider such as "theory"
  or "simulation" (default "theory").

- n_sim:

  number of simulations (large simulations provide more precise
  estimations).

## Value

Probability of acceptance when samples collected from a heterogeneous
batch.

## Details

`prob_acceptance_heterogeneous_multiple` provides a probability of
acceptance when diluted samples are collected from a heterogeneous batch
(this section will be updated later on).

## Examples

``` r
c <- 2
mu <- 7
sd <- 0.2
a <- 0
b <- 300
f <- c(0.01,0.1,1)
u <- c(0.1,0.1,0.1)
USL <- 1000
n <- 5
n_sim <- 50000
prob_acceptance_heterogeneous_multiple (c, mu, sd, a, b, f, u, USL, n)
#> Error: 'f' must be a numeric scalar between 0 and 1
```
