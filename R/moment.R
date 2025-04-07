#' Moments for Statistical Inference
#'
#' This function computes sample moments used for statistical inference
#' in quantile regression.
#'
#' @param y A numeric vector of dependent variables (NT x 1).
#' @param x A numeric matrix of regressors (NT x p), a column of ones will be added automatically for the intercept.
#' @param v.tau A numeric value specifying the quantile of interest.
#'
#' @return A list with components:
#'   \describe{
#'     \item{H}{Hessian matrix.}
#'     \item{J}{Jacobian matrix.}
#'     \item{mean.f}{Estimated density.}
#'   }
#'
#' @references
#' Koenker, R. and G. Bassett Jr. (1978). 
#' Regression quantiles. \emph{Econometrica}, 46(1), 33-50.
#'
#' Oka, T. and Z. Qu (2011).  
#' Estimating structural changes in regression quantiles.  
#' \emph{Journal of Econometrics}, 162(2), 248-267.
#'
#' @importFrom quantreg rq.fit.fnb bandwidth.rq
#' 
#' @noRd

moment = function(y, x, v.tau)
{

  eps = .Machine$double.eps ^ (2/3) #variable to control for estimation precision

  ## size
  nt.size = length(y)

  ## regressor
  x.1 = cbind(1, x)

  ## bandwidth
  h = bandwidth.rq(v.tau, nt.size, hs = FALSE)
  if (v.tau + h > 1)
    stop("v.tau + h > 1:  error in summary.rq")
  if (v.tau - h < 0)
    stop("v.tau - h < 0:  error in summary.rq")

  ## quotient
  bhi = rq.fit.fnb(x.1, y, tau = v.tau + h)$coef
  blo = rq.fit.fnb(x.1, y, tau = v.tau - h)$coef
  dyhat = x.1 %*% (bhi - blo)
  if (any(dyhat <= 0)) {
    warning(paste(sum(dyhat <= 0), "non-positive fis"))
  }
  f = pmax(0, (2 * h)/(dyhat - eps))
  H = (1 / nt.size) * crossprod( (f * x.1), x.1)
  J = (1 / nt.size) * crossprod(x.1)
  mean.f = mean(f)

  ## return
  list(H = H, J = J, mean.f = mean.f)
}

