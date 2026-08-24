# Truncated Poisson-Lognormal Random Number Generation

Generates random deviates from a truncated Poisson-lognormal
distribution. This is a compound distribution where the Poisson mean
parameter follows a lognormal distribution, conditioned on the Poisson
outcome being within the interval a, b.

Generates random deviates from a truncated Poisson-lognormal
distribution.

## Usage

``` r
rtrunpoilog(n, mu, sd, a, b)

rtrunpoilog(n, mu, sd, a, b)
```

## Arguments

- n:

  Number of observations. If `length(n) > 1`, the length is taken to be
  the number required.

- mu:

  Mean of the lognormal distribution on the log scale.

- sd:

  Standard deviation of the lognormal distribution on the log scale.

- a:

  Lower truncation point (lower bound of cell count).

- b:

  Upper truncation point (upper bound of cell count).

## Value

A vector of random numbers from the truncated Poisson-lognormal
distribution. Returns `NA` if rejection sampling fails to converge.

A vector of random numbers from the truncated Poisson-lognormal
distribution.

## Details

The Poisson-lognormal distribution is a mixture distribution where:

- \\\Lambda \sim \text{Lognormal}(\mu, \sigma^2)\\

- \\X \| \Lambda = \lambda \sim \text{Poisson}(\lambda)\\

The truncation limits the support to \\a \le X \le b\\.

The probability mass function is given by: \$\$P(X = x) = \frac{1}{x!
\sqrt{2\pi}\sigma} \int_0^\infty \lambda^{x-1} e^{-\lambda}
\exp\left(-\frac{(\log(\lambda)-\mu)^2}{2\sigma^2}\right) d\lambda\$\$

This distribution is particularly useful for modeling microbial counts
in heterogeneous (non-homogeneous) batches where the contamination level
varies across the batch.

The Poisson-lognormal distribution is a mixture distribution where:

- \\\Lambda \sim \text{Lognormal}(\mu, \sigma^2)\\

- \\X \| \Lambda = \lambda \sim \text{Poisson}(\lambda)\\

The truncation limits the support to \\a \le X \le b\\.

## References

- Gonzales-Barron, U.A., Pilão Cadavez, V.A., Butler, F., 2013.
  Statistical approaches for the design of sampling plans for
  microbiological monitoring of foods. In: Mathematical and Statistical
  Methods in Food Science and Technology. Wiley, pp. 363–384.

- Bulmer, M. G. (1974). On fitting the Poisson lognormal distribution to
  species-abundance data. Biometrics, 30(1), 101-110.

## Examples

``` r
# Generate 10 observations from truncated Poisson-lognormal
# with mu = 2, sd = 0.5, bounded between 0 and 300
set.seed(123)
rtrunpoilog(n = 10, mu = 2, sd = 0.5, a = 0, b = 300)
#>  [1]  9  8 17  9  9 17  7  2 10  9

# Generate with different parameters
set.seed(456)
rtrunpoilog(n = 5, mu = 1, sd = 0.3, a = 0, b = 100)
#> [1] 1 2 5 3 2

# Compare with regular Poisson-lognormal (when bounds are wide)
set.seed(789)
n_sim <- 1000
trunc_sample <- rtrunpoilog(n = n_sim, mu = 2, sd = 0.5, a = 0, b = 300)
regular_sample <- rpois(n = n_sim, lambda = rlnorm(n_sim, 2, 0.5))

cat("Truncated mean:", mean(trunc_sample), "\n")
#> Truncated mean: 8.363 
cat("Regular mean:", mean(regular_sample), "\n")
#> Regular mean: 8.599 
cat("Truncated variance:", var(trunc_sample), "\n")
#> Truncated variance: 27.75298 
cat("Regular variance:", var(regular_sample), "\n")
#> Regular variance: 27.72592 

# When bounds are restrictive (0-10 only)
set.seed(321)
trunc_restricted <- rtrunpoilog(n = n_sim, mu = 3, sd = 0.8, a = 0, b = 10)
#> Error: The truncated Poisson lognormal (TPLN) random variable must be bounded by a and b,
#> which means that mu must be between log(a) and log(b).
#> Current mu = 3.0000, bounds: [-Inf, 2.3026]
cat("Restricted truncated mean:", mean(trunc_restricted), "\n")
#> Error: object 'trunc_restricted' not found

# Generate samples with different sd values
mu_val <- 2
sd_values <- c(0.2, 0.5, 1.0)
set.seed(654)
for (sd_val in sd_values) {
  sample <- rtrunpoilog(n = 10, mu = mu_val, sd = sd_val, a = 0, b = 300)
  cat("mu =", mu_val, "sd =", sd_val, ":",
      summary(sample), "\n")
}
#> mu = 2 sd = 0.2 : 4 6 6 7.1 7.75 13 
#> mu = 2 sd = 0.5 : 3 6.25 7.5 7.4 8.75 12 
#> mu = 2 sd = 1 : 1 2.25 9 12.8 16.25 36 

# Note: Rejection sampling may be slow for very restrictive bounds
# or when sd is very large

# Generate 100 random values
set.seed(123)
rtrunpoilog(n = 100, mu = 0, sd = 1, a = 0, b = 300)
#>   [1]  0  3  5  1  1  8  1  0  0  0  3  1  0  1  0  7  1  0  3  2  0  3  1  1  0
#>  [26]  0  2  1  0  6  2  0  2  0  1  1  1  2  0  1  0  1  0 10  2  1  0  2  5  1
#>  [51]  0  0  1  3  1  6  0  1  1  2  3  2  1  2  0  1  1  0  0  8  1  0  1  0  1
#>  [76]  4  3  0  0  1  2  0  0  1  0  1  1  1  0  4  2  0  0  1  5  1 16  2  0  1

# Generate a single value
rtrunpoilog(n = 1, mu = 0.5, sd = 0.8, a = 10, b = 200)
#> Error: The truncated Poisson lognormal (TPLN) random variable must be bounded by a and b,
#> which means that mu must be between log(a) and log(b).
#> Current mu = 0.5000, bounds: [2.3026, 5.2983]
```
