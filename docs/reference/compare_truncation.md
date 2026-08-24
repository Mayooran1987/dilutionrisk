# Compare Truncated and Regular Distributions

Helper function to visualize the effect of truncation on Poisson and
Poisson-lognormal distributions.

## Usage

``` r
compare_truncation(
  n,
  lambda = NULL,
  mu = NULL,
  sd = NULL,
  a = 0,
  b = 300,
  distribution = c("poisson", "poilog"),
  ...
)
```

## Arguments

- n:

  Number of observations for simulation.

- lambda:

  Poisson parameter for comparison.

- mu:

  Mean parameter for Poisson-lognormal.

- sd:

  Standard deviation parameter for Poisson-lognormal.

- a:

  Lower bound for truncation.

- b:

  Upper bound for truncation.

- distribution:

  Type of distribution: "poisson" or "poilog".

- ...:

  Additional arguments passed to plot.

## Value

Invisibly returns a list of generated samples.

## Examples

``` r
# Compare truncated vs regular Poisson
if (FALSE) { # \dontrun{
compare_truncation(n = 1000, lambda = 10,
                   a = 0, b = 20,
                   distribution = "poisson")

# Compare truncated vs regular Poisson-lognormal
compare_truncation(n = 1000, mu = 2, sd = 0.5,
                   a = 0, b = 50,
                   distribution = "poilog")
} # }
```
