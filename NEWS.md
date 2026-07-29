# QR.break (version 1.0.3)

* A vignette has been added, outlining the methods and working through both
  example datasets. Open it with `vignette("QR-break", package = "QR.break")`.

* `rq.break()`, `sq()`, `dq()`, `sq.test.0vs1()`, `dq.test.0vs1()`, `sq.test.lvsl_1()` and
  `dq.test.lvsl_1()` include an optional argument `norm.method`, which selects how the subgradient
  process underlying the SQ and DQ tests is normalized. The argument is added at the end of each
  argument list.

* `norm.method = "cholesky"` is the default and reproduces the results of earlier versions
  exactly. The normalization is the Cholesky factor of `x'x`. Because that factor is constructed
  sequentially, the test statistics depend on the order
  in which the columns of `x` are listed.

* `norm.method = "spectral"` uses the symmetric square root of the correlation matrix of the
  regressors, rescaled by their standard deviations. The normalization is then free of any
  dependence on the ordering of the columns of `x`, on the units in which they are measured, and
  on their signs. Both choices give the same limiting null distribution, therefore the same critical
  values apply and both tests are valid.

* Estimated break dates are unaffected by `norm.method`, since they are obtained by minimizing
  the check function, which does not depend on the normalization.

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
