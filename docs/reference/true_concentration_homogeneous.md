# True Concentration Estimation for Homogeneous Batches

Estimates the true microbial concentration in the original sample when
diluted samples are collected from a homogeneous batch.

## Usage

``` r
true_concentration_homogeneous(lambda, a, b, f, u, USL, n_sim)

true_concentration_homogeneous_multiple(lambda, a, b, f, u, USL, n_sim)

true_concentration_curves_homogeneous(
  lambda_low,
  lambda_high,
  a,
  b,
  f,
  u,
  USL,
  n_sim
)
```

## Arguments

- lambda:

  Expected microbial count (\\\lambda\\). Must be positive.

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

  Lower bound of expected microbial count for x-axis.

- lambda_high:

  Upper bound of expected microbial count for x-axis.

## Value

Numeric value representing the estimated true concentration.

## Details

The true concentration is estimated as: \$\$C = \frac{X}{f \times u}\$\$
where X is the count of microorganisms on a plate following a truncated
Poisson distribution, f is the final dilution factor, and u is the
amount placed on the plate.

## Examples

``` r
if (FALSE) { # \dontrun{
# Basic usage
true_concentration_homogeneous(lambda = 2000, a = 0, b = 300,
                               f = 0.01, u = 0.1, USL = 1000, n_sim = 50000)

# Multiple dilution schemes
true_concentration_homogeneous_multiple(lambda = 2000, a = 0, b = 300,
                                        f = c(0.01, 0.1), u = c(0.1, 0.1),
                                        USL = 1000, n_sim = 50000)

# Plot concentration curves
true_concentration_curves_homogeneous(lambda_low = 0, lambda_high = 5000,
                                      a = 0, b = 300, f = c(0.01, 0.1),
                                      u = c(0.1, 0.1), USL = 1000,
                                      n_sim = 50000)
} # }
```
