# Visualize CV Comparison

Creates a heatmap or contour plot of CV values across mu and sd
parameters.

## Usage

``` r
plot_cv_comparison(cv_data, type = "heatmap")
```

## Arguments

- cv_data:

  Data frame from compare_cv_heterogeneous.

- type:

  Type of plot: "heatmap" or "contour".

## Value

A ggplot object.

## Examples

``` r
if (FALSE) { # \dontrun{
cv_data <- compare_cv_heterogeneous(
  mu_seq = seq(0, 5, by = 0.25),
  sd_values = seq(0.1, 0.8, by = 0.05),
  a = 0, b = 300,
  f = 0.01, u = 0.1,
  USL = 1000, n_sim = 5000
)
plot_cv_comparison(cv_data, type = "heatmap")
} # }
```
