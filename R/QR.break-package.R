#' Structural Breaks in Quantile Regression
#'
#' Methods for detecting structural breaks, determining
#' the number of breaks, and estimating break locations in linear quantile regression,
#' using a single or multiple quantiles, based on Qu (2008) and Oka and Qu (2011).
#' Applicable to both time series and repeated cross-sectional data.
#'
#' @section Main Functions:
#' \itemize{
#'   \item \code{\link{rq.break}}: Main function for detecting structural breaks in quantile regression
#'   \item \code{\link{sq}}: Performs structural break testing using single quantile approach
#'   \item \code{\link{dq}}: Tests for breaks using multiple quantiles
#'   \item \code{\link{brdate}}: Estimates potential break dates
#' }
#'
#' @references
#' Qu, Z. (2008). Testing for Structural Change in Regression Quantiles.
#' Journal of Econometrics, 146(1), 170-184
#' <doi:10.1016/j.jeconom.2008.08.006>
#'
#' Oka, T., and Qu, Z. (2011).
#' Estimating Structural Changes in Regression Quantiles.
#' Journal of Econometrics, 162(2), 248-267
#' <doi:10.1016/j.jeconom.2011.01.005>
#'
#' @docType package
#' @name QR.break
#' @keywords internal
"_PACKAGE"

NULL
