# Truncated Poisson Random Number Generation

Generates random deviates from a truncated Poisson distribution. The
function uses rejection sampling to generate observations from a Poisson
distribution conditioned to lie within the interval a, b.

## Usage

``` r
rtpois(n, lambda, a, b)
```

## Arguments

- n:

  Number of observations. Must be a positive integer.

- lambda:

  Expected count parameter. Must be non-negative.

- a:

  Lower truncation point. Must be non-negative.

- b:

  Upper truncation point. Must be greater than `a`.

## Value

A vector of random numbers from the truncated Poisson distribution.
Returns `NA` if rejection sampling fails to converge.

## Details

The truncated Poisson distribution has probability mass function:
\$\$P(X = x) = \frac{\lambda^x e^{-\lambda}}{x! \sum\_{k=a}^{b}
\frac{\lambda^k e^{-\lambda}}{k!}}\$\$ for \\a \le x \le b\\.

The function uses rejection sampling with a maximum of 1,000,000
iterations to ensure convergence. For lambda = 0, all observations are 0
if within bounds.

## References

- Johnson, N. L., Kemp, A. W., & Kotz, S. (2005). Univariate Discrete
  Distributions (3rd ed.). Wiley.

## Examples

``` r
# Generate 10 observations from truncated Poisson with lambda = 5
# bounded between 0 and 300
set.seed(123)
rtpois(n = 10, lambda = 5, a = 0, b = 300)
#>  [1] 4 7 4 8 9 2 5 8 5 5

# Generate 5 observations with lower bound 10 and upper bound 50
rtpois(n = 5, lambda = 20, a = 10, b = 50)
#> [1] 27 22 14 12 25

# Lambda = 0 case (all zeros)
rtpois(n = 5, lambda = 0, a = 0, b = 300)
#> [1] 0 0 0 0 0

# Generate for different lambda values
lambdas <- c(1, 5, 10, 20)
set.seed(456)
for (lam in lambdas) {
  cat("Lambda =", lam, ":",
      rtpois(n = 3, lambda = lam, a = 0, b = 100), "\n")
}
#> Lambda = 1 : 0 0 1 
#> Lambda = 5 : 7 7 4 
#> Lambda = 10 : 5 9 7 
#> Lambda = 20 : 21 24 22 

# Compare with regular Poisson (when bounds are wide)
set.seed(789)
trunc_sample <- rtpois(n = 1000, lambda = 10, a = 0, b = 300)
regular_sample <- rpois(n = 1000, lambda = 10)

# Both should be similar when bounds are not restrictive
cat("Truncated mean:", mean(trunc_sample), "\n")
#> Truncated mean: 9.963 
cat("Regular mean:", mean(regular_sample), "\n")
#> Regular mean: 9.927 

# When bounds are restrictive (0-10 only)
set.seed(321)
trunc_restricted <- rtpois(n = 1000, lambda = 20, a = 0, b = 10)
cat("Restricted truncated mean:", mean(trunc_restricted), "\n")
#> Restricted truncated mean: 9.253 

# Note: Rejection sampling may be slow for very restrictive bounds
# or very large lambda values
```
