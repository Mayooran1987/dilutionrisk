# Comparison based on probability of detection curves for different dilution schemes when diluted samples collected from a homogeneous batch.

`pd_validation_homogeneous` provides the probability of detection curves
for validate the results when samples collected from a homogeneous
batch.

## Usage

``` r
pd_validation_homogeneous(lambda_low, lambda_high, a, b, f, u, USL, n_sim)
```

## Arguments

- lambda_low:

  the lower value of the expected microbial count (\\\lambda\\) for use
  in the graphical display's x-axis.

- lambda_high:

  the upper value of the expected microbial count (\\\lambda\\) for use
  in the graphical display's x-axis.

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

- n_sim:

  number of simulations (large simulations provide more precise
  estimations).

## Value

Probability of detection curves when diluted samples collected from a
homogeneous batch.

## Details

[`pd_curves_homogeneous`](https://Mayooran1987.github.io/dilutionrisk/reference/pd_curves_homogeneous.md)
provides probability of detection curves for different dilution schemes
when samples collected from a homogeneous batch (this section will be
updated later on).

## Examples

``` r
lambda_low <- 0
lambda_high <- 5000
a <- 0
b <- 300
f <- 0.01
u <- 0.1
USL <- 1000
n_sim <- 50000
pd_validation_homogeneous(lambda_low, lambda_high, a, b, f, u, USL, n_sim = 500)
#> Error: 'lambda' must be a positive numeric scalar
pd_validation_homogeneous(lambda_low, lambda_high, a, b, f, u, USL, n_sim = 50000)
#> Error: 'lambda' must be a positive numeric scalar
```
