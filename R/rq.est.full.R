#' Estimating Break Sizes and Confidence Intervals Given Break Dates
#'
#' This procedure estimates a linear quantile regression
#' given a set of break dates. It is structured to compute break sizes between adjacent regimes and their confidence intervals.
#'
#' @param y A numeric vector of dependent variables (\eqn{NT \times 1}).
#' @param x A numeric matrix of regressors (\eqn{NT \times p}).
#' @param v.tau A numeric value representing the quantile of interest.
#' @param vec.date A numeric vector of break dates, specified by the user.
#' @param n.size An integer specifying the size of the cross-section (\eqn{N}).
#'
#' @return An object from the quantile regression estimates, \code{rq()}, with structural breaks.
#'
#' @references
#' Koenker, R. and G. Bassett Jr. (1978).
#' Regression Quantiles. \emph{Econometrica}, 46(1), 33–50.
#'
#' Oka, T. and Z. Qu (2011).
#' Estimating Structural Changes in Regression Quantiles.
#' \emph{Journal of Econometrics}, 162(2), 248–267.
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
#' ## break date
#' vec.date = 146
#'
#' # cross-sectional size
#' n.size = 1
#'
#' ## estimation
#' rq.est.full(y, x, v.tau, vec.date, n.size)
#'
#'
#'
#' @import quantreg
#' @export

rq.est.full = function(y, x, v.tau, vec.date, n.size=1)
{
    ## add an intercept
    x   = as.matrix(x)  ## matrix format
    x.1 = cbind(1, x)

    ## sample sizes and number of breaks
    nt.size = length(y)
    t.size  = nt.size / n.size
    p.size  = ncol(x.1)
    n.break = length(vec.date)

    ## note: matrix looks like stairs to obtain break sizes
    b.size = matrix(0, nt.size, (n.break*p.size))
    for(i in 1:n.break){
        r.beg = n.size * vec.date[i] + 1
        c.beg = p.size * (i -1) + 1
        c.end = p.size *  i
        b.size[r.beg:nt.size,c.beg:c.end] = x.1[r.beg:nt.size,]
    }
    big.x = cbind(x, b.size)

    ## estimation
    fit = rq(y ~ big.x, tau = v.tau)   ## w/o intercepts b/c rq() adds it

    ## return
    return(fit)
}

