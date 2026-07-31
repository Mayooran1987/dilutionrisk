# Operating Characteristic (OC) Curves for Heterogeneous Batches

Generates operating characteristic curves comparing different dilution
schemes when samples are collected from a heterogeneous
(non-homogeneous) batch.

## Usage

``` r
OC_curves_heterogeneous(
  c,
  mu_low,
  mu_high,
  sd,
  a,
  b,
  f,
  u,
  USL,
  n,
  type = "theory",
  n_sim = NA
)
```

## Arguments

- c:

  Acceptance number (maximum number of defective units allowed).

- mu_low:

  Lower bound of mean microbial count for x-axis (log scale).

- mu_high:

  Upper bound of mean microbial count for x-axis (log scale).

- sd:

  Standard deviation of the normal distribution on the log scale.

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

- n:

  Number of samples inspected. Must be positive integer.

- type:

  Type of calculation: "theory" (default) or "simulation".

- n_sim:

  Number of simulations. Required when `type = "simulation"`.

## Value

A ggplot object showing OC curves for different dilution schemes.

## Examples

``` r
if (FALSE) { # \dontrun{
OC_curves_heterogeneous(c = 2, mu_low = 4, mu_high = 9, sd = 0.2,
                        a = 0, b = 300, f = c(0.01, 0.1), u = c(0.1, 0.1),
                        USL = 1000, n = 5)
} # }
```
