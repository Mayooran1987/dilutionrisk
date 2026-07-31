# Probability of acceptance estimation for multiple dilution schemes when diluted samples are collected from a homogeneous batch.

`prob_acceptance_homogeneous_multiple` provides a probability of
acceptance for multiple dilution schemes in the original sample when
samples collected from a homogeneous batch

## Usage

``` r
prob_acceptance_homogeneous_multiple(c, lambda, a, b, f, u, USL, n, type, n_sim)
```

## Arguments

- c:

  acceptance number

- lambda:

  the expected microbial count (\\\lambda\\).

- a:

  lower domain of the number of microbial count.

- b:

  upper domain of the number of microbial count.

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

Probability of acceptance when diluted samples are collected from a
homogeneous batch.

## Details

`prob_detection_homogeneous_multiple` provides a probability of
acceptance for multiple dilution schemes in the original sample when
samples collected from a homogeneous batch (this section will be updated
later on).

## Examples

``` r
c <- 2
lambda <- 1000
a <- 0
b <- 300
f <- c(0.01,0.1,1)
u <- c(0.1,0.1,0.1)
USL <- 1000
n <- 5
n_sim <- 50000
prob_acceptance_homogeneous_multiple(c, lambda, a, b, f, u, USL, n)
#> Error: 'f' must be a numeric scalar between 0 and 1
```
