#' Test for a Structural Break in a Conditional Quantile
#'
#' @description The function implements a break test to evaluate whether a single structural break exists at a given quantile.
#'
#' @param y A numeric vector of dependent variables (\eqn{NT \times 1}).
#' @param x A numeric matrix of regressors (\eqn{NT \times p}).
#' @param v.tau A numeric value representing the quantile level.
#' @param n.size An integer specifying the size of the cross-section (\eqn{N}).
#'
#' @return A numeric value representing the test statistic for the presence of a structural break.
#'
#'
#' @references
#' Koenker, R. and G. Bassett Jr, (1978).
#' Regression Quantiles. \emph{Econometrica}, \emph{46}(1), 33–50.
#'
#' Qu, Z. (2008).
#' Testing for Structural Change in Regression Quantiles.
#' \emph{Journal of Econometrics}, \emph{146}(1), 170–184.
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
#' # sq test: 0 vs 1
#' result = sq.test.0vs1(y, x, v.tau, n.size)
#' print(result)
#'
#' @import quantreg
#' @export


sq.test.0vs1 = function(y, x, v.tau, n.size = 1)
{
    ## the number of time periods
    t.size = length(y) / n.size
    x      = as.matrix(x)

    ## regressors
    bigX    = as.matrix( cbind(1, x) )
    sqXX    = chol( t(bigX) %*% bigX)
    invsqXX = solve( t(sqXX) )

    ## quantile regression using the full sample
    res = rq(y ~ x, tau = v.tau)$res

    ## Test
    temp = (res <= 0.0) - v.tau
    H1n  = invsqXX %*% t(bigX) %*% temp
    difH = matrix(0, 1, t.size)

    ## loop
    for (j in 2: t.size){

        end1    = n.size * j
        HH      = t(bigX)[,1:end1] %*% temp[1:end1]
        difH[j] = max( abs(invsqXX %*% HH - (j / t.size) * H1n) )

    }

    v.max = max(difH) / sqrt( v.tau * (1 - v.tau) )

    ## Result
    return( v.max )
}

