#' Sequential Test for Additional Breaks within a Range of Quantiles
#'
#' This procedure tests for the existence of \eqn{L} breaks against \eqn{L+1} breaks
#' based on multiple quantiles:
#' \eqn{H_0: L} breaks vs. \eqn{H_1: L+1} breaks.
#'
#' @description
#' This function performs a sequential test to determine whether the number of
#' breaks in a quantile regression model should be increased from \eqn{L} to \eqn{L+1} using
#' multiple quantiles.
#'
#' @param y A numeric vector of dependent variables (\eqn{NT \times 1}).
#' @param x A numeric matrix of regressors (\eqn{NT \times p}).
#' @param q.L A numeric value specifying the lower bound of the quantile range.
#' @param q.R A numeric value specifying the upper bound of the quantile range.
#' @param vec.date A numeric vector (\eqn{L \times 1}) of estimated break dates under the null hypothesis.
#' @param n.size An integer specifying the size of cross-sections (\eqn{N}).
#'
#' @return A numeric value representing the DQ test statistic.
#'
#' @references
#' Qu, Z. (2008). Testing for Structural Breaks in Regression Quantiles.
#' \emph{Journal of Econometrics}, 146(1), 170-184.
#'
#'
#' @examples
#' # Load data
#' data(gdp)
#' y = gdp$gdp
#' x = gdp[,c("lag1", "lag2")]
#'
#' # Set quantile range (left and right limits)
#' q.L = 0.2
#' q.R = 0.8
#'
#' # Set N parameter
#' n.size = 1
#'
#' # Specify break date under H_0
#' vec.date = 146
#'
#' # Run the test
#' result = dq.test.lvsl_1(y, x, q.L, q.R, n.size, vec.date)
#' print(result)
#'
#' @export

dq.test.lvsl_1 = function(y, x, q.L, q.R, n.size = 1, vec.date) #order of vec.date and n.size adjusted 11/2023
{
    ## the number of breaks
    n.break  = length(vec.date)

    ## test for n.break+1 breaks
    vec.test = matrix(0, (n.break+1), 1)
    rem.y    = y
    rem.x    = x
    pre.date = c(0, vec.date)
    for(j in 1:n.break){

        v.date = vec.date[j] - pre.date[j]

        ## split sample
        temp   = sample.split(rem.y, rem.x, v.date, n.size)
        rem.y  = temp$y2
        rem.x  = temp$x2

        ## test statistics
        # vec.test[j] = DQtest(temp$y1, temp$x1, q.L, q.R, n.size) switched 11/23
        vec.test[j] = dq.test.0vs1(temp$y1, temp$x1, q.L, q.R, n.size)

    }
    ## the last regime
   # vec.test[(n.break+1)] = DQtest(rem.y, rem.x, q.L, q.R, n.size) switched 11/2023
    vec.test[(n.break+1)] = dq.test.0vs1(rem.y, rem.x, q.L, q.R, n.size)

    ## return: maximum over tests
    return( max(vec.test) )
}

