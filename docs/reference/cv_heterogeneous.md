# Coefficient of Variation for Heterogeneous Batches

Estimates the coefficient of variation (CV) in the original sample when
diluted samples are collected from a heterogeneous (non-homogeneous)
batch. The CV measures the relative variability of the estimated
microbial concentration.

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

## Value

Numeric value representing the coefficient of variation for the
heterogeneous batch.

## Details

The coefficient of variation is defined as: \$\$CV =
\frac{\sigma}{\mu}\$\$ where \\\sigma\\ is the standard deviation and
\\\mu\\ is the mean of the estimated microbial concentration.

For heterogeneous batches, the microbial count follows a truncated
Poisson-lognormal distribution, accounting for variability in
contamination levels across the batch.

## References

- Gonzales-Barron, U.A., Pilão Cadavez, V.A., Butler, F., 2013.
  Statistical approaches for the design of sampling plans for
  microbiological monitoring of foods. In: Mathematical and Statistical
  Methods in Food Science and Technology. Wiley, pp. 363–384.

- Schothorst, M. van, Zwietering, M.H., Ross, T., Buchanan, R.L., Cole,
  M.B., 2009. Relating microbiological criteria to food safety
  objectives and performance objectives. Food Control 20, 967–979.

## Examples

``` r
if (FALSE) { # \dontrun{
# Basic usage - estimate CV for a single dilution scheme
cv_heterogeneous(mu = 2, sd = 0.2, a = 0, b = 300,
                 f = 0.01, u = 0.1, USL = 1000, n_sim = 10000)

# Compare multiple dilution schemes
cv_heterogeneous_multiple(mu = 2, sd = 0.2, a = 0, b = 300,
                          f = c(0.01, 0.1), u = c(0.1, 0.1),
                          USL = 1000, n_sim = 10000)

# Generate CV curves across a range of mean microbial counts
cv_curves_heterogeneous(mu_low = -5, mu_high = 10, sd = 0.2,
                        a = 0, b = 300, f = c(0.01, 0.1),
                        u = c(0.1, 0.1), USL = 1000, n_sim = 10000)
} # }

# Quick example with small simulation size for testing
cv_heterogeneous(mu = 2, sd = 0.2, a = 0, b = 300,
                 f = 0.01, u = 0.1, USL = 1000, n_sim = 500)
#> [1] 0.5100018
```
