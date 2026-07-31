# Changelog

## dilutionrisk 0.1.0 (Development)

### Major Changes

#### New Features

##### Homogeneous Batch Models

- Added
  [`prob_detection_homogeneous()`](https://Mayooran1987.github.io/dilutionrisk/reference/prob_detection_homogeneous.md) -
  Calculate probability of detection for homogeneous batches
- Added `prob_detection_homogeneous_multiple()` - Calculate PD for
  multiple dilution schemes
- Added
  [`prob_acceptance_homogeneous()`](https://Mayooran1987.github.io/dilutionrisk/reference/prob_acceptance_homogeneous.md) -
  Calculate probability of acceptance for sampling plans
- Added `prob_acceptance_homogeneous_multiple()` - Calculate PA for
  multiple dilution schemes
- Added
  [`cv_homogeneous()`](https://Mayooran1987.github.io/dilutionrisk/reference/cv_homogeneous.md) -
  Estimate coefficient of variation for homogeneous batches
- Added
  [`cv_homogeneous_multiple()`](https://Mayooran1987.github.io/dilutionrisk/reference/cv_homogeneous.md) -
  Estimate CV for multiple dilution schemes
- Added
  [`cv_curves_homogeneous()`](https://Mayooran1987.github.io/dilutionrisk/reference/cv_homogeneous.md) -
  Generate CV curves for homogeneous batches
- Added
  [`OC_curves_homogeneous()`](https://Mayooran1987.github.io/dilutionrisk/reference/OC_curves_homogeneous.md) -
  Generate operating characteristic curves
- Added
  [`pd_curves_homogeneous()`](https://Mayooran1987.github.io/dilutionrisk/reference/pd_curves_homogeneous.md) -
  Generate probability of detection curves
- Added
  [`true_concentration_homogeneous()`](https://Mayooran1987.github.io/dilutionrisk/reference/true_concentration_homogeneous.md) -
  Estimate true concentration
- Added
  [`true_concentration_homogeneous_multiple()`](https://Mayooran1987.github.io/dilutionrisk/reference/true_concentration_homogeneous.md) -
  Estimate concentration for multiple schemes
- Added
  [`true_concentration_curves_homogeneous()`](https://Mayooran1987.github.io/dilutionrisk/reference/true_concentration_homogeneous.md) -
  Generate concentration curves

##### Heterogeneous Batch Models

- Added
  [`prob_detection_heterogeneous()`](https://Mayooran1987.github.io/dilutionrisk/reference/prob_detection_heterogeneous.md) -
  Calculate probability of detection for heterogeneous batches
- Added `prob_detection_heterogeneous_multiple()` - Calculate PD for
  multiple dilution schemes
- Added
  [`prob_acceptance_heterogeneous()`](https://Mayooran1987.github.io/dilutionrisk/reference/prob_acceptance_heterogeneous.md) -
  Calculate probability of acceptance for sampling plans
- Added `prob_acceptance_heterogeneous_multiple()` - Calculate PA for
  multiple dilution schemes
- Added
  [`cv_heterogeneous()`](https://Mayooran1987.github.io/dilutionrisk/reference/cv_heterogeneous.md) -
  Estimate coefficient of variation for heterogeneous batches
- Added
  [`cv_heterogeneous_multiple()`](https://Mayooran1987.github.io/dilutionrisk/reference/cv_heterogeneous.md) -
  Estimate CV for multiple dilution schemes
- Added
  [`cv_curves_heterogeneous()`](https://Mayooran1987.github.io/dilutionrisk/reference/cv_heterogeneous.md) -
  Generate CV curves for heterogeneous batches
- Added
  [`OC_curves_heterogeneous()`](https://Mayooran1987.github.io/dilutionrisk/reference/OC_curves_heterogeneous.md) -
  Generate operating characteristic curves
- Added
  [`pd_curves_heterogeneous()`](https://Mayooran1987.github.io/dilutionrisk/reference/pd_curves_heterogeneous.md) -
  Generate probability of detection curves
- Added
  [`true_concentration_heterogeneous()`](https://Mayooran1987.github.io/dilutionrisk/reference/true_concentration_heterogeneous.md) -
  Estimate true concentration
- Added
  [`true_concentration_heterogeneous_multiple()`](https://Mayooran1987.github.io/dilutionrisk/reference/true_concentration_heterogeneous.md) -
  Estimate concentration for multiple schemes
- Added
  [`true_concentration_curves_heterogeneous()`](https://Mayooran1987.github.io/dilutionrisk/reference/true_concentration_heterogeneous.md) -
  Generate concentration curves

##### Distribution Functions

- Added
  [`rtrunpoilog()`](https://Mayooran1987.github.io/dilutionrisk/reference/rtrunpoilog.md) -
  Generate random numbers from truncated Poisson-lognormal distribution

##### Validation Functions

- Added
  [`pd_validation_homogeneous()`](https://Mayooran1987.github.io/dilutionrisk/reference/pd_validation_homogeneous.md) -
  Validate homogeneous model with simulations
- Added
  [`pd_validation_heterogeneous()`](https://Mayooran1987.github.io/dilutionrisk/reference/pd_validation_heterogeneous.md) -
  Validate heterogeneous model with simulations

#### Performance Improvements

- Implemented C++ code using Rcpp for efficient random number generation
- Optimized probability calculations using log-sum-exp for numerical
  stability
- Added progress bars for long simulations
- Improved vectorization for multiple dilution scheme calculations

#### Documentation

- Added comprehensive package documentation
- Added vignettes:
  - Introduction to dilutionrisk
  - Homogeneous Batch Models
  - Heterogeneous Batch Models
  - Validation and Simulation
  - Case Studies
- Created pkgdown website for online documentation
- Added PDF manual
- Added README with examples and installation instructions

#### Testing

- Added unit tests using testthat
- Added test coverage tracking with covr
- Added continuous integration with GitHub Actions
- Added CI/CD workflows:
  - R-CMD-check
  - Test coverage
  - pkgdown build and deploy

#### Bug Fixes

- Fixed numerical stability issues in probability calculations
- Fixed edge cases when USL \* f \* u \< 1
- Fixed input validation for all functions
- Fixed memory issues in large simulations

#### Deprecated Functions

- None (initial release)

#### Breaking Changes

- None (initial release)

------------------------------------------------------------------------

## dilutionrisk 0.0.1 (Pre-release)

### Initial Development

#### Features

- Basic probability of detection calculations
- Basic probability of acceptance calculations
- Simulation-based validation
- Truncated Poisson-lognormal distribution

#### Known Issues

- Numerical stability issues in extreme cases
- No formal testing framework
- Limited documentation
- No vignettes

------------------------------------------------------------------------

### Version History

#### Version 0.1.0 - Planned Release

Complete all core functionality

Add comprehensive documentation

Implement C++ optimization

Add test suite

Create vignettes

Build pkgdown site

Add CI/CD workflows

Fix all known bugs

Achieve \>90% code coverage

#### Version 1.0.0 - Future Release

Submit to CRAN

Add Bayesian estimation methods

Add more dilution schemes

Add more validation methods

Add interactive Shiny application

Add more case studies

Add support for other microbial count distributions

------------------------------------------------------------------------

### Contributors

- **Mayooran Thevaraja** - Author, Maintainer
- **Mark Bebbington** - Author

------------------------------------------------------------------------

### How to Contribute

We welcome contributions to `dilutionrisk`! Please see our contributing
guidelines.

#### Reporting Issues

If you find a bug or have a feature request, please open an issue on
GitHub: <https://github.com/Mayooran1987/dilutionrisk/issues>

#### Submitting Changes

1.  Fork the repository
2.  Create a new branch for your changes
3.  Make your changes
4.  Submit a pull request

------------------------------------------------------------------------

### References

For more information about the package and its development, see:

- Package website: <https://Mayooran1987.github.io/dilutionrisk/>
- GitHub repository: <https://github.com/Mayooran1987/dilutionrisk>
- Package manual:
  [dilutionrisk_0.0.1.pdf](https://github.com/Mayooran1987/dilutionrisk/blob/main/dilutionrisk_0.0.1.pdf)

------------------------------------------------------------------------

### License

This package is licensed under GPL-2.

------------------------------------------------------------------------

### Acknowledgments

We thank the R community and our respective institutions for their
support.

------------------------------------------------------------------------

**Last updated**: 2026-08-01
