# True Concentration Estimation for Heterogeneous Batches

Estimates the true microbial concentration in the original sample when
diluted samples are collected from a heterogeneous (non-homogeneous)
batch.

## Usage

``` r
true_concentration_heterogeneous(mu, sd, a, b, f, u, USL, n_sim)

true_concentration_heterogeneous_multiple(mu, sd, a, b, f, u, USL, n_sim)

true_concentration_curves_heterogeneous(
  mu_low,
  mu_high,
  sd,
  a,
  b,
  f,
  u,
  USL,
  n_sim
)
```

## Arguments

- mu:

  Mean microbial count on the log scale. Must satisfy
  `log(a) - log(f*u) <= mu <= log(b) - log(f*u)`.

- sd:

  Standard deviation of the normal distribution on the log scale.

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

- n_sim:

  Number of simulations. Larger values provide more precise estimates.
  Recommended minimum: 10,000.

- mu_low:

  Lower bound of mean microbial count for x-axis (log scale).

- mu_high:

  Upper bound of mean microbial count for x-axis (log scale).

## Value

Numeric value representing the estimated true concentration.

## Details

The true concentration is estimated as: \$\$C = \frac{X}{f \times u}\$\$
where X is the count of microorganisms on a plate following a truncated
Poisson-lognormal distribution, f is the final dilution factor, and u is
the amount placed on the plate.

## References

- Gonzales-Barron, U.A., Pilão Cadavez, V.A., Butler, F., 2013.
  Statistical approaches for the design of sampling plans for
  microbiological monitoring of foods. In: Mathematical and Statistical
  Methods in Food Science and Technology. Wiley, pp. 363–384.

## Examples

``` r
# Basic usage
true_concentration_heterogeneous(mu = 2, sd = 0.2, a = 0, b = 300,
                                 f = 0.01, u = 0.1, USL = 1000, n_sim = 5000)
#> [1] 5.8

# Multiple dilution schemes
true_concentration_heterogeneous_multiple(mu = 2, sd = 0.2, a = 0, b = 300,
                                          f = c(0.01, 0.1), u = c(0.1, 0.1),
                                          USL = 1000, n_sim = 5000)
#>      [,1] [,2]
#> [1,]  5.4  7.1

# Plot concentration curves
true_concentration_curves_heterogeneous(mu_low = 0, mu_high = 10, sd = 0.2,
                                        a = 0, b = 300, f = c(0.01, 0.1),
                                        u = c(0.1, 0.1), USL = 1000,
                                        n_sim = 5000)
#>   |                                                                              |                                                                      |   0%  |                                                                              |=                                                                     |   1%  |                                                                              |=                                                                     |   2%  |                                                                              |==                                                                    |   3%  |                                                                              |===                                                                   |   4%  |                                                                              |===                                                                   |   5%  |                                                                              |====                                                                  |   6%  |                                                                              |=====                                                                 |   7%  |                                                                              |======                                                                |   8%  |                                                                              |======                                                                |   9%  |                                                                              |=======                                                               |  10%  |                                                                              |========                                                              |  11%  |                                                                              |========                                                              |  12%  |                                                                              |=========                                                             |  13%  |                                                                              |==========                                                            |  14%  |                                                                              |==========                                                            |  15%  |                                                                              |===========                                                           |  16%  |                                                                              |============                                                          |  17%  |                                                                              |============                                                          |  18%  |                                                                              |=============                                                         |  19%  |                                                                              |==============                                                        |  20%  |                                                                              |===============                                                       |  21%  |                                                                              |===============                                                       |  22%  |                                                                              |================                                                      |  23%  |                                                                              |=================                                                     |  24%  |                                                                              |=================                                                     |  25%  |                                                                              |==================                                                    |  26%  |                                                                              |===================                                                   |  27%  |                                                                              |===================                                                   |  28%  |                                                                              |====================                                                  |  29%  |                                                                              |=====================                                                 |  30%  |                                                                              |=====================                                                 |  31%  |                                                                              |======================                                                |  32%  |                                                                              |=======================                                               |  33%  |                                                                              |========================                                              |  34%  |                                                                              |========================                                              |  35%  |                                                                              |=========================                                             |  36%  |                                                                              |==========================                                            |  37%  |                                                                              |==========================                                            |  38%  |                                                                              |===========================                                           |  39%  |                                                                              |============================                                          |  40%  |                                                                              |============================                                          |  41%  |                                                                              |=============================                                         |  42%  |                                                                              |==============================                                        |  43%  |                                                                              |==============================                                        |  44%  |                                                                              |===============================                                       |  45%  |                                                                              |================================                                      |  46%  |                                                                              |=================================                                     |  47%  |                                                                              |=================================                                     |  48%  |                                                                              |==================================                                    |  49%  |                                                                              |===================================                                   |  50%  |                                                                              |====================================                                  |  51%  |                                                                              |=====================================                                 |  52%  |                                                                              |=====================================                                 |  53%  |                                                                              |======================================                                |  54%  |                                                                              |=======================================                               |  55%  |                                                                              |========================================                              |  56%  |                                                                              |========================================                              |  57%  |                                                                              |=========================================                             |  58%  |                                                                              |==========================================                            |  59%  |                                                                              |==========================================                            |  60%  |                                                                              |===========================================                           |  61%  |                                                                              |============================================                          |  62%  |                                                                              |============================================                          |  63%  |                                                                              |=============================================                         |  64%  |                                                                              |==============================================                        |  65%  |                                                                              |==============================================                        |  66%  |                                                                              |===============================================                       |  67%  |                                                                              |================================================                      |  68%  |                                                                              |=================================================                     |  69%  |                                                                              |=================================================                     |  70%  |                                                                              |==================================================                    |  71%  |                                                                              |===================================================                   |  72%  |                                                                              |===================================================                   |  73%  |                                                                              |====================================================                  |  74%  |                                                                              |=====================================================                 |  75%  |                                                                              |=====================================================                 |  76%  |                                                                              |======================================================                |  77%  |                                                                              |=======================================================               |  78%  |                                                                              |=======================================================               |  79%  |                                                                              |========================================================              |  80%  |                                                                              |=========================================================             |  81%  |                                                                              |==========================================================            |  82%  |                                                                              |==========================================================            |  83%  |                                                                              |===========================================================           |  84%  |                                                                              |============================================================          |  85%  |                                                                              |============================================================          |  86%  |                                                                              |=============================================================         |  87%  |                                                                              |==============================================================        |  88%  |                                                                              |==============================================================        |  89%  |                                                                              |===============================================================       |  90%  |                                                                              |================================================================      |  91%  |                                                                              |================================================================      |  92%  |                                                                              |=================================================================     |  93%  |                                                                              |==================================================================    |  94%  |                                                                              |===================================================================   |  95%  |                                                                              |===================================================================   |  96%  |                                                                              |====================================================================  |  97%  |                                                                              |===================================================================== |  98%
#> Warning: Rejection sampling reached maximum iterations for lambda = 386.687468797014, bounds = [0, 300]
#> Warning: Rejection sampling reached maximum iterations for lambda = 411.552798858616, bounds = [0, 300]
#>   |                                                                              |===================================================================== |  99%
#> Warning: Rejection sampling reached maximum iterations for lambda = 466.873504005675, bounds = [0, 300]
#> Warning: Rejection sampling reached maximum iterations for lambda = 411.321107905397, bounds = [0, 300]
#> Warning: Rejection sampling reached maximum iterations for lambda = 412.335562568732, bounds = [0, 300]
#> Warning: Rejection sampling reached maximum iterations for lambda = 429.049679137279, bounds = [0, 300]
#> Warning: Rejection sampling reached maximum iterations for lambda = 416.249069598012, bounds = [0, 300]
#> Warning: Rejection sampling reached maximum iterations for lambda = 409.44442449395, bounds = [0, 300]
#> Warning: Rejection sampling reached maximum iterations for lambda = 392.055203120516, bounds = [0, 300]
#> Warning: Rejection sampling reached maximum iterations for lambda = 398.420621193089, bounds = [0, 300]
#> Warning: Rejection sampling reached maximum iterations for lambda = 424.489181036379, bounds = [0, 300]
#> Warning: Rejection sampling reached maximum iterations for lambda = 433.343547063478, bounds = [0, 300]
#>   |                                                                              |======================================================================| 100%
#> Warning: Removed 2 rows containing missing values or values outside the scale range
#> (`geom_line()`).


```
