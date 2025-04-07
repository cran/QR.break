#' US Real GDP Growth Data
#'
#' This dataset contains quarterly real GDP growth rates in the US from 1947 Q4 to 2009 Q2. 
#' It is used in Oka and Qu (2011) to examine whether and where distributional breaks occur in US GDP.
#' 
#'
#' @docType data
#' @keywords datasets
#' @name gdp
#' @usage data(gdp)
#' @format A data frame with four variables:
#' \itemize{
#'   \item{\code{yq}} A character vector representing the year and quarter (e.g., "1947 Q2").
#'   \item{\code{gdp}} A numeric vector of real GDP growth rates (annualized by multiplying by 4).
#'   \item{\code{lag1}} A numeric vector representing the first-order lagged value of \code{gdp}.
#'   \item{\code{lag2}} A numeric vector representing the second-order lagged value of \code{gdp}.
#' }
#'
#' @references
#' Oka, T. and Z. Qu (2011). Estimating Structural Changes in Regression Quantiles.
#' \emph{Journal of Econometrics}, \emph{162}(2), 248–267.
#'
#' @examples
#' data(gdp)
#' names(gdp)
#' summary(gdp)
#'
"gdp"

