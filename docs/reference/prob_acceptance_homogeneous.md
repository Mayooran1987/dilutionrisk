# Probability of Acceptance for Homogeneous Batches

Calculates the probability of acceptance (PA) in the original sample
when diluted samples are collected from a homogeneous batch.

## Usage

``` r
prob_acceptance_homogeneous(
  c,
  lambda,
  a,
  b,
  f,
  u,
  USL,
  n,
  type = "theory",
  n_sim = NA
)

prob_acceptance_homogeneous_multiple(
  c,
  lambda,
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

- n:

  Number of samples inspected. Must be positive integer.

- type:

  Type of calculation: "theory" (default) or "simulation".

- n_sim:

  Number of simulations. Required when `type = "simulation"`.

## Value

Numeric value representing the probability of acceptance.

## Details

The probability of acceptance is calculated using the binomial
distribution: \$\$P_a = P(X \le c) = \sum\_{k=0}^{c} \binom{n}{k} p_d^k
(1-p_d)^{n-k}\$\$ where \\p_d\\ is the probability of detection and
\\n\\ is the number of samples inspected.

## Examples

``` r
if (FALSE) { # \dontrun{
# Basic usage
prob_acceptance_homogeneous(c = 2, lambda = 2000, a = 0, b = 300,
                            f = 0.01, u = 0.1, USL = 1000, n = 5)

# Multiple dilution schemes
prob_acceptance_homogeneous_multiple(c = 2, lambda = 2000, a = 0, b = 300,
                                     f = c(0.01, 0.1), u = c(0.1, 0.1),
                                     USL = 1000, n = 5)
} # }
```
