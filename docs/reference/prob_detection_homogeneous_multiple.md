# Probability of detection estimation for multiple dilution schemes when diluted samples are collected from a homogeneous batch.

`prob_detection_homogeneous_multiple` provides a probability of
detection for multiple dilution schemes in the original sample when
samples collected from a homogeneous batch.

## Usage

``` r
prob_detection_homogeneous_multiple(lambda, a, b, f, u, USL, type, n_sim)
```

## Arguments

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

- type:

  what type of the results you would like to consider such as "theory"
  or "simulation" (default "theory").

- n_sim:

  number of simulations (large simulations provide more precise
  estimations).

## Value

Probability of detection when diluted samples are collected from a
homogeneous batch.

## Details

`prob_detection_homogeneous_multiple` provides a probability of
detection when the diluted sample has homogeneous contaminants. We
define the random variable \\X\_{i}\\ is the number of colonies on the
\\i^{th}\\ plate. In practice, the acceptance for countable numbers of
colonies on a plate must be between 30 and 300. Therefore, we can
utilise bounded distributions to model the number of colonies on a
plate. In the homogeneous case, we employed truncated Poisson
distribution to model (this section will be updated later on).

## Examples

``` r
lambda <- 1000
a <- 0
b <- 300
f <- c(0.01,0.1,1)
u <- c(0.1,0.1,0.1)
USL <- 1000
n_sim <- 50000
prob_detection_homogeneous_multiple(lambda, a, b, f, u, USL)
#> Error: 'f' must be a numeric scalar between 0 and 1
```
