# Package index

## Random Number Generation

Functions for generating random numbers from truncated distributions

- [`rtpois()`](https://mayooran1987.github.io/dilutionrisk/reference/rtpois.md)
  : Truncated Poisson Random Number Generation
- [`rtrunpoilog()`](https://mayooran1987.github.io/dilutionrisk/reference/rtrunpoilog.md)
  : Truncated Poisson-Lognormal Random Number Generation

## Coefficient of Variation

Estimate coefficient of variation for homogeneous and heterogeneous
batches

- [`cv_homogeneous()`](https://mayooran1987.github.io/dilutionrisk/reference/cv_homogeneous.md)
  [`cv_homogeneous_multiple()`](https://mayooran1987.github.io/dilutionrisk/reference/cv_homogeneous.md)
  [`cv_curves_homogeneous()`](https://mayooran1987.github.io/dilutionrisk/reference/cv_homogeneous.md)
  : Coefficient of Variation for Homogeneous Batches
- [`cv_heterogeneous()`](https://mayooran1987.github.io/dilutionrisk/reference/cv_heterogeneous.md)
  [`cv_heterogeneous_multiple()`](https://mayooran1987.github.io/dilutionrisk/reference/cv_heterogeneous.md)
  [`cv_curves_heterogeneous()`](https://mayooran1987.github.io/dilutionrisk/reference/cv_heterogeneous.md)
  : Coefficient of Variation for Heterogeneous Batches

## Probability of Detection

Calculate probability of detection for microbial contamination

- [`prob_detection_homogeneous()`](https://mayooran1987.github.io/dilutionrisk/reference/prob_detection_homogeneous.md)
  [`prob_detection_homogeneous_multiple()`](https://mayooran1987.github.io/dilutionrisk/reference/prob_detection_homogeneous.md)
  : Probability of Detection for Homogeneous Batches
- [`prob_detection_homogeneous_multiple()`](https://mayooran1987.github.io/dilutionrisk/reference/prob_detection_homogeneous_multiple.md)
  : Probability of detection estimation for multiple dilution schemes
  when diluted samples are collected from a homogeneous batch.
- [`prob_detection_heterogeneous()`](https://mayooran1987.github.io/dilutionrisk/reference/prob_detection_heterogeneous.md)
  [`prob_detection_heterogeneous_multiple()`](https://mayooran1987.github.io/dilutionrisk/reference/prob_detection_heterogeneous.md)
  : Probability of Detection for Heterogeneous Batches
- [`prob_detection_heterogeneous_multiple()`](https://mayooran1987.github.io/dilutionrisk/reference/prob_detection_heterogeneous_multiple.md)
  : Probability of detection estimation for multiple dilution schemes
  when diluted samples are collected from a heterogeneous batch.

## Probability of Acceptance

Calculate probability of acceptance for sampling plans

- [`prob_acceptance_homogeneous()`](https://mayooran1987.github.io/dilutionrisk/reference/prob_acceptance_homogeneous.md)
  [`prob_acceptance_homogeneous_multiple()`](https://mayooran1987.github.io/dilutionrisk/reference/prob_acceptance_homogeneous.md)
  : Probability of Acceptance for Homogeneous Batches
- [`prob_acceptance_homogeneous_multiple()`](https://mayooran1987.github.io/dilutionrisk/reference/prob_acceptance_homogeneous_multiple.md)
  : Probability of acceptance estimation for multiple dilution schemes
  when diluted samples are collected from a homogeneous batch.
- [`prob_acceptance_heterogeneous()`](https://mayooran1987.github.io/dilutionrisk/reference/prob_acceptance_heterogeneous.md)
  [`prob_acceptance_heterogeneous_multiple()`](https://mayooran1987.github.io/dilutionrisk/reference/prob_acceptance_heterogeneous.md)
  : Probability of Acceptance for Heterogeneous Batches
- [`prob_acceptance_heterogeneous_multiple()`](https://mayooran1987.github.io/dilutionrisk/reference/prob_acceptance_heterogeneous_multiple.md)
  : Probability of acceptance estimation when diluted samples are
  collected from a heterogeneous batch.

## Operating Characteristic Curves

Generate operating characteristic curves for dilution schemes

- [`OC_curves_homogeneous()`](https://mayooran1987.github.io/dilutionrisk/reference/OC_curves_homogeneous.md)
  : Operating Characteristic (OC) Curves for Homogeneous Batches
- [`OC_curves_heterogeneous()`](https://mayooran1987.github.io/dilutionrisk/reference/OC_curves_heterogeneous.md)
  : Operating Characteristic (OC) Curves for Heterogeneous Batches

## Probability of Detection Curves

Generate probability of detection curves for dilution schemes

- [`pd_curves_homogeneous()`](https://mayooran1987.github.io/dilutionrisk/reference/pd_curves_homogeneous.md)
  : Probability of Detection Curves for Homogeneous Batches
- [`pd_curves_heterogeneous()`](https://mayooran1987.github.io/dilutionrisk/reference/pd_curves_heterogeneous.md)
  : Probability of Detection Curves for Heterogeneous Batches

## True Concentration Estimation

Estimate true microbial concentration from diluted samples

- [`true_concentration_homogeneous()`](https://mayooran1987.github.io/dilutionrisk/reference/true_concentration_homogeneous.md)
  [`true_concentration_homogeneous_multiple()`](https://mayooran1987.github.io/dilutionrisk/reference/true_concentration_homogeneous.md)
  [`true_concentration_curves_homogeneous()`](https://mayooran1987.github.io/dilutionrisk/reference/true_concentration_homogeneous.md)
  : True Concentration Estimation for Homogeneous Batches
- [`true_concentration_heterogeneous()`](https://mayooran1987.github.io/dilutionrisk/reference/true_concentration_heterogeneous.md)
  [`true_concentration_heterogeneous_multiple()`](https://mayooran1987.github.io/dilutionrisk/reference/true_concentration_heterogeneous.md)
  [`true_concentration_curves_heterogeneous()`](https://mayooran1987.github.io/dilutionrisk/reference/true_concentration_heterogeneous.md)
  : True Concentration Estimation for Heterogeneous Batches

## Validation Tools

Compare theoretical and simulation-based results

- [`pd_validation_homogeneous()`](https://mayooran1987.github.io/dilutionrisk/reference/pd_validation_homogeneous.md)
  : Comparison based on probability of detection curves for different
  dilution schemes when diluted samples collected from a homogeneous
  batch.
- [`pd_validation_heterogeneous()`](https://mayooran1987.github.io/dilutionrisk/reference/pd_validation_heterogeneous.md)
  : Comparison based on probability of detection curves for different
  dilution schemes when diluted samples collected from a heterogeneous
  batch.

## Helper Functions

Utility functions for sensitivity analysis and comparison

- [`compare_cv_heterogeneous()`](https://mayooran1987.github.io/dilutionrisk/reference/compare_cv_heterogeneous.md)
  : Helper function to understand CV behavior
- [`plot_cv_comparison()`](https://mayooran1987.github.io/dilutionrisk/reference/plot_cv_comparison.md)
  : Visualize CV Comparison
