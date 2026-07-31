# Truncated Poisson-Lognormal Random Number Generation

Generates random deviates from a truncated Poisson-lognormal
distribution.

## Usage

``` r
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
distribution.

## Details

The Poisson-lognormal distribution is a mixture distribution where:

- \\\Lambda \sim \text{Lognormal}(\mu, \sigma^2)\\

- \\X \| \Lambda = \lambda \sim \text{Poisson}(\lambda)\\

The truncation limits the support to \\a \le X \le b\\.

## Examples

``` r
if (FALSE) { # \dontrun{
# Generate 100 random values
set.seed(123)
rtrunpoilog(n = 100, mu = 0, sd = 1, a = 0, b = 300)

# Generate a single value
rtrunpoilog(n = 1, mu = 0.5, sd = 0.8, a = 10, b = 200)
} # }
```
