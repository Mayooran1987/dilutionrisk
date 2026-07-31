# Coefficient of Variation for Heterogeneous Batches

Estimates the coefficient of variation (CV) in the original sample when
diluted samples are collected from a heterogeneous (non-homogeneous)
batch.

## Usage

``` r
cv_heterogeneous(mu, sd, a, b, f, u, USL, n_sim)

cv_heterogeneous_multiple(mu, sd, a, b, f, u, USL, n_sim)

cv_curves_heterogeneous(mu_low, mu_high, sd, a, b, f, u, USL, n_sim)
```

## Arguments

- mu:

  Mean microbial count on the log scale. Must satisfy
  `log(a) - log(f*u) <= mu <= log(b) - log(f*u)`.

- sd:

  Standard deviation of the normal distribution on the log scale. Must
  be positive.

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

- mu_low:

  Lower bound of mean microbial count (\\\mu\\) for x-axis in graphical
  displays (log scale).

- mu_high:

  Upper bound of mean microbial count (\\\mu\\) for x-axis in graphical
  displays (log scale).

## Value

Numeric value representing the coefficient of variation for the
heterogeneous batch.

## References

- Gonzales-Barron, U.A., Pilão Cadavez, V.A., Butler, F., 2013.
  Statistical approaches for the design of sampling plans for
  microbiological monitoring of foods. In: Mathematical and Statistical
  Methods in Food Science and Technology. Wiley, pp. 363–384.

## Examples

``` r
if (FALSE) { # \dontrun{
# Basic usage
cv_heterogeneous(mu = 2, sd = 0.2, a = 0, b = 300,
                 f = 0.01, u = 0.1, USL = 1000, n_sim = 50000)

# Multiple dilution schemes
cv_heterogeneous_multiple(mu = 2, sd = 0.2, a = 0, b = 300,
                          f = c(0.01, 0.1), u = c(0.1, 0.1),
                          USL = 1000, n_sim = 50000)

# Plot CV curves
cv_curves_heterogeneous(mu_low = -5, mu_high = 10, sd = 0.2,
                        a = 0, b = 300, f = c(0.01, 0.1),
                        u = c(0.1, 0.1), USL = 1000, n_sim = 50000)
} # }
```
