#' Test for a Structural Break in a Conditional Quantile
#'
#' @description The function implements a break test to evaluate whether a single structural break exists at a given quantile.
#'
#' @usage sq.test.0vs1(y, x, v.tau, n.size, norm.method)
#'
#' @param y A numeric vector of dependent variables (\eqn{NT \times 1}).
#' @param x A numeric matrix of regressors (\eqn{NT \times p}).
#' @param v.tau A numeric value representing the quantile level.
#' @param n.size An integer specifying the size of the cross-section (\eqn{N}).
#' @param norm.method A character string specifying how the subgradient process is normalized;
#' either \code{"cholesky"} (the default) or \code{"spectral"}. See Details.
#'
#' @return A numeric value representing the test statistic for the presence of a structural break.
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
#'
#' @references
#' Kessy, A., A. Lewin and K. Strimmer (2018).
#' Optimal Whitening and Decorrelation.
#' \emph{The American Statistician}, \emph{72}(4), 309–314.
#'
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


sq.test.0vs1 = function(y, x, v.tau, n.size = 1, norm.method = c("cholesky", "spectral"))
{
    norm.method = match.arg(norm.method)

    ## the number of time periods
    t.size = length(y) / n.size
    x      = as.matrix(x)

    ## regressors
    bigX = as.matrix( cbind(1, x) )
    MM   = t(bigX) %*% bigX

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

