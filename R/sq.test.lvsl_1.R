#' Sequential Test for an Additional Break in a Conditional Quantile
#'
#' This function tests the null hypothesis of \eqn{L} breaks against the alternative hypothesis of \eqn{L+1} breaks
#' in a single conditional quantile.
#'
#'
#' @param y A numeric vector of dependent variables (\eqn{NT \times 1}).
#' @param x A numeric matrix of regressors (\eqn{NT \times p}).
#' @param v.tau A numeric value representing the quantile of interest.
#' @param vec.date A numeric vector of break dates estimated under the null hypothesis.
#' @param n.size An integer specifying the size of the cross-section (\eqn{N}).
#'
#' @return A numeric value representing the test statistic.
#'
#' @details The function sequentially tests for breaks by splitting the sample conditional on the
#' break dates under the null hypothesis. At each step, it applies \code{sq.test.0vs1()} to compare
#' the hypothesis of no additional break against one more break.
#'
#'
#' @references
#' Qu, Z. (2008). Testing for Structural Change in Regression Quantiles.
#' *Journal of Econometrics*, 146(1), 170-184.
#'
#' Oka, T. and Z. Qu (2011). Estimating Structural Changes in Regression Quantiles.
#' *Journal of Econometrics*, 162(2), 248-267.
#'
#'
#' @examples
#' ## data
#' data(gdp)
#' y = gdp$gdp
#' x = gdp[,c("lag1", "lag2")]
#'
#' ## quantile
#' v.tau = 0.8
#'
#' # cross-sectional size
#' n.size = 1
#'
#' ## break date
#' vec.date = 146
#'
#' ## sq-test: 1 vs 2
#' result = sq.test.lvsl_1(y, x, v.tau, n.size, vec.date)
#' print(result)
#'
#'
#' @export

sq.test.lvsl_1 = function(y, x, v.tau, n.size = 1, vec.date) #order n.size, vec.date adjusted
{
    ## the number of breaks
    n.break  = length(vec.date)

    ## sequential test
    vec.test = matrix(0, (n.break+1), 1)
    rem.y    = y
    rem.x    = x
    pre.date = c(0, vec.date)
    for(j in 1:n.break){

        v.date = vec.date[j] - pre.date[j]

        ## split sample given a break date
        temp   = sample.split(rem.y, rem.x, v.date, n.size)
        rem.y  = temp$y2
        rem.x  = temp$x2

        ## test statistics
       # vec.test[j] = SQtest(temp$y1, temp$x1, v.tau, n.size) switched 11/2023
        vec.test[j] = sq.test.0vs1(temp$y1, temp$x1, v.tau, n.size)
    }

    ## the last regime
    #vec.test[(n.break+1)] = SQtest(rem.y, rem.x, v.tau, n.size) switched 11/2023
    vec.test[(n.break+1)] = sq.test.0vs1(rem.y, rem.x, v.tau, n.size)

    ## return: maximum over regimes
    return( max(vec.test) )
}



