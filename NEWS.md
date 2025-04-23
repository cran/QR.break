# QR.break (version 1.0.2)

* Minor improvements to output: additional estimation results are now properly saved and returned by `rq.break()`
* Small updates to documentation for clarity

# QR.break (version 1.0.1)

* Initial CRAN release
* Implemented methods for structural break detection in linear quantile regression based on:
  * Qu, Z. (2008). Testing for structural Change in Regression Quantiles. Journal of Econometrics, 146(1), 170-184.
  * Oka, T., & Qu, Z. (2011). Estimating Structural Changes in Regression Quantiles. Journal of Econometrics, 162(2), 248-267.

* Core functionality includes:
  * Detection of structural breaks in quantile regression models
  * Determination of optimal number of breaks
  * Estimation of break point locations
  * Support for both single and multiple quantile analysis
  * Compatibility with time series and repeated cross-sectional data
