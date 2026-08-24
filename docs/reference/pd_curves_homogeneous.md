# Probability of Detection Curves for Homogeneous Batches

Generates probability of detection curves comparing different dilution
schemes when samples are collected from a homogeneous batch.

## Usage

``` r
pd_curves_homogeneous(
  lambda_low,
  lambda_high,
  a,
  b,
  f,
  u,
  USL,
  type = "theory",
  n_sim = NA
)
```

## Arguments

- lambda_low:

  Lower bound of expected microbial count for x-axis.

- lambda_high:

  Upper bound of expected microbial count for x-axis.

- a:

  Lower bound of cell count domain. Must be non-negative.

- b:

  Upper bound of cell count domain. Must be greater than `a`.

- f:

  Vector of final dilution factors.

- u:

  Vector of amounts placed on the plate.

- USL:

  Upper specification limit for microbial count.

- type:

  Type of calculation: "theory" (default) or "simulation".

- n_sim:

  Number of simulations. Required when `type = "simulation"`.

## Value

A ggplot object showing PD curves for different dilution schemes.

## Examples

``` r
pd_curves_homogeneous(lambda_low = 0, lambda_high = 5000,
                      a = 0, b = 300, f = c(0.01, 0.1), u = c(0.1, 0.1),
                      USL = 1000)
#>   |                                                                              |                                                                      |   0%
#> Error: 'lambda' must be a positive numeric scalar
```
