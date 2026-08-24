# Operating Characteristic (OC) Curves for Homogeneous Batches

Generates operating characteristic curves comparing different dilution
schemes when samples are collected from a homogeneous batch.

## Usage

``` r
OC_curves_homogeneous(
  c,
  lambda_low,
  lambda_high,
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

- n:

  Number of samples inspected. Must be positive integer.

- type:

  Type of calculation: "theory" (default) or "simulation".

- n_sim:

  Number of simulations. Required when `type = "simulation"`.

## Value

A ggplot object showing OC curves for different dilution schemes.

## Details

Operating characteristic curves show the probability of acceptance as a
function of the true microbial concentration. They are useful for
comparing the performance of different dilution schemes and sampling
plans.

## Examples

``` r
OC_curves_homogeneous(c = 2, lambda_low = 0, lambda_high = 5000,
                      a = 0, b = 300, f = c(0.01, 0.1), u = c(0.1, 0.1),
                      USL = 1000, n = 5)
#>   |                                                                              |                                                                      |   0%
#> Error: 'lambda' must be a positive numeric scalar

```
