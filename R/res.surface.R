#' Compute Critical Values for the DQ test using a Response Surface
#'
#' This function returns critical values obtained from a response surface analysis.
#' Note that this procedure only applies when the trimming is symmetric, i.e., \eqn{q.R=1-q.L}.
#'
#' @param p The number of parameters in the model.
#' @param l The number of breaks under the null \eqn{H_0} (i.e., \eqn{l+1} under \eqn{H_1}).
#' @param q.L The lower bound of the quantile range.
#' @param q.R The upper bound of the quantile range (not used in the function because \eqn{q.R=1-q.L}).
#' @param d.Sym A logical value. Must be TRUE, as this method applies only to symmetric trimming (\eqn{q.R=1-q.L}).
#'
#' @return A numeric vector of length 3 containing critical values at the 10%, 5%, and 1% significance levels.
#'
#' @export
#' 
#' @examples 
#' # The number of regerssors 
#' p = 5 
#' ## The number of breaks under the null 
#' l = 2 
#' 
#' # qunatile range (left and right limits)
#' q.L = 0.2
#' q.R = 0.8 
#' 
#' # symmetric quantile trimming is true 
#' d.Sym = TRUE 
#' 
#' ## critical values from response surface 
#' cvs = res.surface(p, l, q.L, q.R, d.Sym)
#' 
#' print(cvs)
#' 

res.surface = function(p, l, q.L, q.R, d.Sym)
{
    ## We need symmetry
    if (d.Sym != TRUE){
        stop('This program is only for DQ test with a symmetric quantile range, i.e., [tau, 1 - tau]')
    }

    ## storage
    vec.cv = matrix(0, 3, 1)

    ## linear part 1
    vec.cv[1] = 0.9481 + 0.0062 * p + 0.0166 * (l + 1) -0.1386 * (1 / p)
    vec.cv[2] = 0.9944 + 0.0058 * p + 0.0157 * (l + 1) -0.1284 * (1 / p)
    vec.cv[3] = 1.0929 + 0.0050 * p + 0.0134 * (l + 1) -0.1134 * (1 / p)


    ## linear part 2
    vec.cv[1] = vec.cv[1] -0.0004 * (l+1) * p + 0.0018 * (l+1) * q.L
    vec.cv[2] = vec.cv[2] -0.0004 * (l+1) * p + 0.0017 * (l+1) * q.L
    vec.cv[3] = vec.cv[3] -0.0002 * (l+1) * p + 0.0010 * (l+1) * q.L

    ## exponential part
    vec.cv[1] = vec.cv[1]*exp(-0.0801 * (1/(l+1)) -0.0004 * (1 /(q.L*(l+1))) -0.0254 * q.L)
    vec.cv[2] = vec.cv[2]*exp(-0.0716 * (1/(l+1)) -0.0005 * (1 /(q.L*(l+1))) -0.0203 * q.L)
    vec.cv[3] = vec.cv[3]*exp(-0.0565 * (1/(l+1)) -0.0000 * (1 /(q.L*(l+1))) -0.0062 * q.L)

    ## return
    return(vec.cv)
}


