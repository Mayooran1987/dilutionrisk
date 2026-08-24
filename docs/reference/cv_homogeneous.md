# Coefficient of Variation for Homogeneous Batches

Estimates the coefficient of variation (CV) in the original sample when
diluted samples are collected from a homogeneous batch.

## Usage

``` r
cv_homogeneous(lambda, a, b, f, u, USL, n_sim)

cv_homogeneous_multiple(lambda, a, b, f, u, USL, n_sim)

cv_curves_homogeneous(lambda_low, lambda_high, a, b, f, u, USL, n_sim)
```

## Arguments

- lambda:

  Expected microbial count (\\\lambda\\). Must be non-negative.

- a:

  Lower bound of cell count domain. Must be non-negative.

- b:

  Upper bound of cell count domain. Must be greater than `a`.

- f:

  Final dilution factor. Must be between 0 and 1.

- u:

  Amount placed on the plate. Must be positive.

- USL:

  Upper specification limit for microbial count.

- n_sim:

  Number of simulations. Larger values provide more precise estimates.
  Recommended minimum: 10,000.

- lambda_low:

  Lower bound of expected microbial count for x-axis in graphical
  displays.

- lambda_high:

  Upper bound of expected microbial count for x-axis in graphical
  displays.

## Value

Numeric value representing the coefficient of variation for the
homogeneous batch.

## Examples

``` r
# Basic usage
cv_homogeneous(lambda = 2000, a = 0, b = 300,
               f = 0.01, u = 0.1, USL = 1000, n_sim = 5000)
#> [1] 0.01022323
```
