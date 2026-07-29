#' Test for the Presence of a Break over a Range of Quantiles
#'
#' This procedure computes test statistics for detecting a single structural break
#' over a range of quantiles. The null hypothesis is that there is no break
#' in any quantile in the specified range; the alternative is that at least one
#' quantile in the range is affected by a break.
#'
#' @usage dq.test.0vs1(y, x, q.L, q.R, n.size, norm.method)
#'
#' @param y A numeric vector of dependent variables (\eqn{NT \times 1}).
#' @param x A numeric matrix of regressors (\eqn{NT \times p}), excluding the constant term.
#' @param q.L A numeric value specifying the lower bound of the quantile range.
#' @param q.R A numeric value specifying the upper bound of the quantile range.
#' @param n.size An integer specifying the size of the cross-section (\eqn{N}).
#' @param norm.method A character string specifying how the subgradient process is normalized;
#' either \code{"cholesky"} (the default) or \code{"spectral"}. See Details.
#'
#' @return A numeric scalar representing the DQ test statistic.
#'
#' @details
#' Let \eqn{W = \sum_{t} \sum_{i} x_{it} x_{it}'}{W = sum_t sum_i x_it x_it'}, where for a
#' single time series the sums over \eqn{i} are absent. The test normalizes the subgradient
#' process by an inverse square root of \eqn{W}. See Qu (2008). Such a matrix is not unique,
#' and \code{norm.method} selects which one is used.
#'
#' With \code{"cholesky"}, the normalization is \eqn{(R')^{-1}}, where \eqn{R} is the Cholesky
#' factor of \eqn{W} returned by \code{chol()}, with \eqn{R'R = W}. This is the normalization
#' used in versions 1.0.2 and earlier, and is the default setting. This option implicitly
#' gives more weight to regressors appearing earlier in the regression. It is invariant to the
#' units in which the regressors are measured, and to their signs.
#'
#' With \code{"spectral"}, the normalization is \eqn{C^{-1/2} D^{-1}}, where
#' \eqn{D = \mathrm{diag}(\sqrt{\mathrm{diag}(W)})}{D = diag(sqrt(diag(W)))} and
#' \eqn{C = D^{-1} W D^{-1}} is the correlation matrix of the regressors, and \eqn{C^{-1/2}}
#' is its symmetric square root obtained from the spectral decomposition. It yields invariance
#' to the order of the regressors, to their units of measurement, and to their signs.
#'
#' These invariance properties refer to the normalization. Because the quantile regression is
#' re-estimated, the computed statistics can differ slightly when the regressors are reordered
#' or rescaled, by an amount that decreases with the sample size.
#'
#' Both choices give the same limiting null distribution, so the same critical values apply.
#' The second option might yield lower power. See \code{\link{rq.break}} for further discussion.
#'
#' @references
#' Kessy, A., Lewin, A. and Strimmer, K. (2018).
#' Optimal whitening and decorrelation.
#' \emph{The American Statistician}, \emph{72}(4), 309–314.
#'
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

dq.test.0vs1 = function(y, x, q.L, q.R, n.size = 1, norm.method = c("cholesky", "spectral"))
{
    norm.method = match.arg(norm.method)

    t.size = length(y) / n.size
    x      = as.matrix(x)  # matrix

    ## regressors
    bigX   = cbind(1, x)
    p.size = ncol(bigX)
    MM     = t(bigX) %*% bigX

    ## normalization: invsqXX satisfies invsqXX %*% MM %*% t(invsqXX) = I
    if (norm.method == "cholesky"){

        ## Cholesky factor, as in versions 1.0.2 and earlier
        sqXX    = chol( MM )
        invsqXX = solve( t(sqXX) )

    } else {

        ## symmetric square root of the correlation matrix, rescaled:
        ## invsqXX = Mtilde^(-1/2) D^(-1), with D = diag( sqrt(diag(MM)) )
        d.scale = sqrt( diag(MM) )
        if ( any(d.scale <= 0) ){
            stop("A column of the regressor matrix is identically zero.")
        }
        eig = eigen( MM / outer(d.scale, d.scale), symmetric = TRUE )
        if ( min(eig$values) <= 0 | min(eig$values) / max(eig$values) < 1e-12 ){
            stop("The regressor matrix is (nearly) collinear; the test statistic is not ",
                 "well defined. This can occur when a regressor is (nearly) constant. ")
        }
        invsqXX = sweep( eig$vectors %*% ( t(eig$vectors) / sqrt(eig$values) ),
                         2, d.scale, '/' )

    }

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

