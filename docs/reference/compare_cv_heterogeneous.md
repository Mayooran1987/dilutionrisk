# Helper function to understand CV behavior

Compares CV for different dilution schemes and parameters. Useful for
sensitivity analysis.

## Usage

``` r
compare_cv_heterogeneous(
  mu_seq,
  sd_values,
  a,
  b,
  f,
  u,
  USL,
  n_sim = 10000,
  ...
)
```

## Arguments

- mu_seq:

  Vector of mu values to evaluate.

- sd_values:

  Vector of sd values to evaluate.

- a:

  Lower bound.

- b:

  Upper bound.

- f:

  Dilution factor.

- u:

  Plate amount.

- USL:

  Upper specification limit.

- n_sim:

  Number of simulations.

- ...:

  Additional arguments passed to cv_heterogeneous.

## Value

A data frame with CV values.

## Examples

``` r
if (FALSE) { # \dontrun{
# Compare CV across different mu and sd values
cv_sensitivity <- compare_cv_heterogeneous(
  mu_seq = seq(0, 5, by = 0.5),
  sd_values = c(0.1, 0.3, 0.5),
  a = 0, b = 300,
  f = 0.01, u = 0.1,
  USL = 1000, n_sim = 10000
)
print(cv_sensitivity)
} # }
```
