#' The Dataset for Young Drivers
#'
#' This is a repeated cross-sectional data set on young drivers 
#' (under 21 years old) involved in motor vehicle accidents in the state of California 
#' from 1983 to 2007 (quarterly data). The data are obtained from the National Highway 
#' Traffic Safety Administration (NHTSA), which include the blood alcohol concentration (BAC) 
#' of the driver, their age, gender, and whether the crash was fatal.
#' 
#' Motor vehicle crashes are the leading cause of death among youth aged 15–20, with a high 
#' proportion involving drunk driving. The BAC level is an important measure of alcohol impairment.
#' Oka and Qu (2011) used this data to examine whether and how young drivers' drinking behaviors 
#' have changed over time.
#'
#' @docType data
#' @keywords datasets
#' @name driver
#' @usage data(driver)
#' @format A data frame with five variables:
#' \itemize{
#'   \item \code{yq}: Year and quarter ("Year Quarter" format, e.g., "1983 Q2").
#'   \item \code{bac}: The blood alcohol concentration (BAC) level of the driver.
#'   \item \code{age}: The driver's age.
#'   \item \code{gender}: A gender dummy, with 1 for male and 0 for female.
#'   \item \code{winter}: A dummy variable for the fourth quarter, with 1 for Q4 and 0 otherwise.
#' }
#'
#' @references
#' Oka, T. and Z. Qu (2011). Estimating Structural Changes in Regression Quantiles. \emph{Journal of Econometrics}, 162(2), 248–267.
#'
#' @examples
#' data(driver)
#' names(driver)
#' summary(driver)
#' 
#' 
"driver"

