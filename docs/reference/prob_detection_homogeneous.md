# Probability of Detection for Homogeneous Batches

Calculates the probability of detection (PD) in the original sample when
diluted samples are collected from a homogeneous batch.

## Usage

``` r
prob_detection_homogeneous(
  lambda,
  a,
  b,
  f,
  u,
  USL,
  type = "theory",
  n_sim = NA
)

prob_detection_homogeneous_multiple(
  lambda,
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

- type:

  Type of calculation: "theory" (default) or "simulation".

- n_sim:

  Number of simulations. Required when `type = "simulation"`.

## Value

Numeric value representing the probability of detection.

## Details

The probability of detection is defined as the probability that the
estimated microbial count exceeds the upper specification limit (USL).
For homogeneous batches, the count follows a truncated Poisson
distribution.

The theoretical probability is calculated as: \$\$P_d = 1 -
\sum\_{i=1}^{USL \times f \times u} \frac{(\lambda f u)^i}{i!(e^{\lambda
f u} - 1)}\$\$

## References

- Schothorst, M. van, Zwietering, M.H., Ross, T., Buchanan, R.L., Cole,
  M.B., 2009. Relating microbiological criteria to food safety
  objectives and performance objectives. Food Control 20, 967–979.

## Examples

``` r
if (FALSE) { # \dontrun{
# Theoretical calculation
prob_detection_homogeneous(lambda = 2000, a = 0, b = 300,
                           f = 0.01, u = 0.1, USL = 1000)

# Simulation-based calculation
prob_detection_homogeneous(lambda = 2000, a = 0, b = 300,
                           f = 0.01, u = 0.1, USL = 1000,
                           type = "simulation", n_sim = 50000)
} # }
```
