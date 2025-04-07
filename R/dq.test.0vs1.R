#' Test for the Presence of a Break within a Range of Quantiles
#'
#' This procedure computes test statistics for detecting a single structural break
#' within a range of quantiles. The null hypothesis is that there is no break
#' in any quantile in the specified range; the alternative is that at least one
#' quantile in the range is affected by a break.
#'
#' @param y A numeric vector of dependent variables (\eqn{NT \times 1}).
#' @param x A numeric matrix of regressors (\eqn{NT \times p}), excluding the constant term.
#' @param q.L A numeric value specifying the lower bound of the quantile range.
#' @param q.R A numeric value specifying the upper bound of the quantile range.
#' @param n.size An integer specifying the size of the cross-section (\eqn{N}).
#'
#' @return A numeric scalar representing the DQ test statistic.
#'
#' @references
#' Koenker, R. and Bassett Jr, G. (1978).
#' Regression quantiles. \emph{Econometrica}, \emph{46}(1), 33–50.
#'
#' Qu, Z. (2008).
#' Testing for structural change in regression quantiles.
#' \emph{Journal of Econometrics}, \emph{146}(1), 170–184.
#'
#' @examples
#' ## data
#' data(gdp)
#' y = gdp$gdp
#' x = gdp[,c("lag1", "lag2")]
#'
#' ## qunatile range (left and right limits)
#' q.L = 0.2
#' q.R = 0.8
#'
#' ## N
#' n.size = 1
#'
#' # dq-test
#' result = dq.test.0vs1(y, x, q.L, q.R, n.size)
#' print(result)
#'
#' @import quantreg
#' @export

dq.test.0vs1 = function(y, x, q.L, q.R, n.size = 1)
{

    t.size = length(y) / n.size
    x      = as.matrix(x)  # matrix

    ## regressors
    bigX    = cbind(1, x)
    p.size  = ncol(bigX)
    sqXX    = chol( t(bigX) %*% bigX)
    invsqXX = solve( t(sqXX) )

    ## trimming
    seq.tau = seq(q.L, q.R, by = 1 / t.size)
    n.tau   = length(seq.tau)
    Qstat   = matrix(0, n.tau, 1)

    for(k in 1:n.tau){

        ## tau
        v.tau = seq.tau[k]

        ## fit
        res  = rq(y ~ 1 + x, tau = v.tau)$residuals
        temp = (res <= 0.0) - v.tau

        ## test
        H1n     = invsqXX %*% t(bigX) %*% temp
        H1n     = t(H1n)
        Hlambda = matrix(0, t.size, p.size)
        difH    = Hlambda
        for (t in 2:t.size) {

            HH          = t(bigX)[,1:(t*n.size)] %*% temp[1:(t*n.size)]
            Hlambda[t,] = invsqXX %*% HH
            difH[t,]    = Hlambda[t,] - (t / t.size) * H1n

        }

        ## maximum over v.tau
        Qstat[k] = max(abs(difH))
    }

    ## return
    return(max(Qstat))
}

