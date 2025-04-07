#' Confidence Intervals for Break Dates
#'
#' This function constructs confidence intervals for break dates
#' based on a single quantile or multiple quantiles (specified by the user).
#'
#' @param y A numeric vector of dependent variables (\eqn{NT \times 1}).
#' @param x A numeric matrix of regressors (\eqn{NT \times p}).
#' @param vec.tau A numeric vector of quantiles of interest.
#' @param vec.date A numeric vector of estimated break dates.
#' @param n.size An integer specifying the size of the cross section (\eqn{N}).
#' @param v.b A numeric value specifying the confidence level:
#'  \itemize{
#'    \item \code{v.b = 1} for the 90% confidence interval.
#'    \item \code{v.b = 2} for the 95% confidence interval.
#' }
#'
#' @return A numeric matrix where:
#' \itemize{
#'    \item The 1st column contains the break dates.
#'    \item The 2nd and 3rd columns contain the lower and upper bounds of the confidence intervals, respectively.
#' }
#'
#' @references
#' Oka, T. and Z. Qu (2011).
#' Estimating Structural Changes in Regression Quantiles.
#' \emph{Journal of Econometrics}, 162(2), 248–267.
#' 
#' @examples
#' 
#' # data 
#' data(gdp)
#' y = gdp$gdp
#' x = gdp[,c("lag1", "lag2")] 
#' 
#' # quantiles 
#' vec.tau  = 0.8
#' 
#' # break dates (point estimates)
#' vec.date = c(146, 200)
#' 
#' # Calculate confidence intervals for break dates
#' res = ci.date.m(y, x, vec.tau, vec.date, n.size = 1, v.b = 2)
#' print(res)
#'
#' @export

ci.date.m = function(y, x, vec.tau, vec.date, n.size=1, v.b=2)
{
    # matrix 
    x = as.matrix(x)  # <=== Formatting 
  
    ## size
    n.break = length(vec.date)
    n.tau   = length(vec.tau)
    p.size  = ncol(x) + 1    ## add intercept

    ## storage
    mat.ci     = matrix(0, n.break, 3)
    mat.ci[,1] = matrix(vec.date, n.break, 1)

    mat.size = matrix(0, (n.break*p.size), n.tau)
    for(k in 1:n.tau){
        v.tau = vec.tau[k]
        vec.b = rq.est.full(y, x, v.tau, vec.date, n.size)$coef[]
        mat.size[,k] = vec.b[(p.size+1):((n.break+1)*p.size)]
    }

    ## confidence interval
    ## 1st regime
    sub01  = sample.split(y, x, vec.date[1], n.size)
    y.L    = sub01$y1
    x.L    = sub01$x1
    y.rem  = sub01$y2
    x.rem  = sub01$x2

    ## If the number of breaks is more than two.
    if(n.break >= 2){
        for(i in 1:(n.break-1)){

            ## split sample
            date02 = vec.date[(i+1)] - vec.date[i]
            sub02  = sample.split(y.rem, x.rem, date02, n.size)
            y.R    = sub02$y1
            x.R    = sub02$x1
            y.rem  = sub02$y2
            x.rem  = sub02$x2

            ## confidence interval for 1st - (m-1)th break
            v.date    = vec.date[i] # <== Put a real break date
            beg.row   = p.size * (i -1) + 1
            end.row   = p.size * i
            temp.size = as.matrix(mat.size[beg.row:end.row,], ncol = n.tau)
            mat.ci[i,] = ci.date.m.sub(y.L, x.L, y.R, x.R, v.date, temp.size, vec.tau, n.size, v.b)

            ## replace
            y.L = y.R
            x.L = x.R
        }
    }

    ## For the last break
    v.date    = vec.date[n.break] # <== Put a real break date
    beg.row   = p.size * (n.break -1) + 1
    end.row   = p.size *  n.break
    temp.size = as.matrix(mat.size[beg.row:end.row,], ncol = n.tau)
    mat.ci[n.break,] = ci.date.m.sub(y.L, x.L, y.rem, x.rem, v.date, temp.size, vec.tau, n.size, v.b)

    ## return
    return(mat.ci)
}
