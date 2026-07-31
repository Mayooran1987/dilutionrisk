# dilutionrisk: Risk Assessment for Dilution Testing

The `dilutionrisk` package implements statistical methods for
microbiological risk assessment in food safety, focusing on:

- **Probability of Detection (PD)**: Probability that the estimated
  microbial count exceeds the upper specification limit.

- **Probability of Acceptance (PA)**: Probability that a batch is
  accepted based on sampling plans.

- **Operating Characteristic (OC) Curves**: Visual comparison of
  different dilution schemes.

- **Coefficient of Variation (CV)**: Variability estimation for
  different dilution schemes.

- **True Concentration Estimation**: Estimating the true microbial
  concentration in the original sample.

## Details

Provides probability estimations and graphical displays for modelling
and assessment of risk based on aerobic plate count (APC) on diluted
testing.

The package supports two scenarios:

**Homogeneous Batches:** Assumes the microbial distribution is uniform
throughout the batch. The count follows a truncated Poisson
distribution.

**Heterogeneous Batches:** Assumes the microbial distribution is
non-uniform. The count follows a truncated Poisson-lognormal
distribution.

Both scenarios incorporate dilution factors and plate amounts to model
the dilution process in microbiological testing.

## Key Functions

- Probability Detection:

  [`prob_detection_homogeneous`](https://Mayooran1987.github.io/dilutionrisk/reference/prob_detection_homogeneous.md),
  [`prob_detection_heterogeneous`](https://Mayooran1987.github.io/dilutionrisk/reference/prob_detection_heterogeneous.md)

- Probability Acceptance:

  [`prob_acceptance_homogeneous`](https://Mayooran1987.github.io/dilutionrisk/reference/prob_acceptance_homogeneous.md),
  [`prob_acceptance_heterogeneous`](https://Mayooran1987.github.io/dilutionrisk/reference/prob_acceptance_heterogeneous.md)

- OC Curves:

  [`OC_curves_homogeneous`](https://Mayooran1987.github.io/dilutionrisk/reference/OC_curves_homogeneous.md),
  [`OC_curves_heterogeneous`](https://Mayooran1987.github.io/dilutionrisk/reference/OC_curves_heterogeneous.md)

- PD Curves:

  [`pd_curves_homogeneous`](https://Mayooran1987.github.io/dilutionrisk/reference/pd_curves_homogeneous.md),
  [`pd_curves_heterogeneous`](https://Mayooran1987.github.io/dilutionrisk/reference/pd_curves_heterogeneous.md)

- CV Estimation:

  [`cv_homogeneous`](https://Mayooran1987.github.io/dilutionrisk/reference/cv_homogeneous.md),
  [`cv_heterogeneous`](https://Mayooran1987.github.io/dilutionrisk/reference/cv_heterogeneous.md)

- True Concentration:

  [`true_concentration_homogeneous`](https://Mayooran1987.github.io/dilutionrisk/reference/true_concentration_homogeneous.md),
  [`true_concentration_heterogeneous`](https://Mayooran1987.github.io/dilutionrisk/reference/true_concentration_heterogeneous.md)

## References

- Gonzales-Barron, U.A., Pilão Cadavez, V.A., Butler, F., 2013.
  Statistical approaches for the design of sampling plans for
  microbiological monitoring of foods. In: Mathematical and Statistical
  Methods in Food Science and Technology. Wiley, pp. 363–384.

- Schothorst, M. van, Zwietering, M.H., Ross, T., Buchanan, R.L., Cole,
  M.B., 2009. Relating microbiological criteria to food safety
  objectives and performance objectives. Food Control 20, 967–979.

## See also

Useful links:

- <https://github.com/Mayooran1987/dilutionrisk>

- Report bugs at <https://github.com/Mayooran1987/dilutionrisk/issues>

## Author

**Maintainer**: Mayooran Thevaraja <mayooran@eng.jfn.ac.lk>
([ORCID](https://orcid.org/0000-0002-4786-7248)) (University of Jaffna)

Authors:

- Mark Bebbington <m.bebbington@massey.ac.nz>
  ([ORCID](https://orcid.org/0000-0003-3504-7418)) (Massey University)
