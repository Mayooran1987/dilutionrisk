# Probability of Detection for Heterogeneous Batches

Calculates the probability of detection (PD) in the original sample when
diluted samples are collected from a heterogeneous (non-homogeneous)
batch.

## Usage

``` r
prob_detection_heterogeneous(
  mu,
  sd,
  a,
  b,
  f,
  u,
  USL,
  type = "theory",
  n_sim = NA
)

prob_detection_heterogeneous_multiple(
  mu,
  sd,
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

- type:

  Type of calculation: "theory" (default) or "simulation".

- n_sim:

  Number of simulations. Required when `type = "simulation"`.

## Value

Numeric value representing the probability of detection.

## Details

The probability of detection is defined as the probability that the
estimated microbial count exceeds the upper specification limit (USL).
For heterogeneous batches, the count follows a truncated
Poisson-lognormal distribution.

The theoretical probability is calculated by integrating over the
lognormal mixing distribution: \$\$P_d = 1 - \sum\_{i=1}^{USL \times f
\times u} \frac{1}{i!\sqrt{2\pi}\sigma} \int_0^\infty \frac{t^{i-1}
\exp(-0.5((\log(t)-\mu_d)/\sigma)^2)}{e^t - 1} dt\$\$

## References

- Gonzales-Barron, U.A., Pilão Cadavez, V.A., Butler, F., 2013.
  Statistical approaches for the design of sampling plans for
  microbiological monitoring of foods. In: Mathematical and Statistical
  Methods in Food Science and Technology. Wiley, pp. 363–384.

## Examples

``` r
# Theoretical calculation
prob_detection_heterogeneous(mu = 2, sd = 0.2, a = 0, b = 300,
                             f = 0.01, u = 0.1, USL = 1000)
#> [1] 0.003764234

# Simulation-based calculation
prob_detection_heterogeneous(mu = 2, sd = 0.2, a = 0, b = 300,
                             f = 0.01, u = 0.1, USL = 1000,
                             type = "simulation", n_sim = 5000)
#> [1] 0
```
