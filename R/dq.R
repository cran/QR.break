#' Sequential Determination of the Number of Breaks Using the DQ Test
#'
#' This function determines the number of breaks in a model by
#' sequentially applying the DQ(\eqn{l | l+1}) test. It tests for additional breaks
#' by comparing the test statistic to critical values at various significance levels.
#'
#' @param y A numeric vector of dependent variables (\eqn{NT \times 1}).
#' @param x A numeric matrix of regressors (\eqn{NT \times p}).
#' @param vec.tau A numeric vector of quantiles of interest.
#' @param q.L A numeric value specifying the left-end quantile range for the DQ test.
#' @param q.R A numeric value specifying the right-end quantile range for the DQ test.
#' @param n.size An integer specifying the size of the cross-section (\eqn{N}).
#' @param m.max An integer indicating the maximum number of breaks allowed.
#' @param trim.size A numeric trimming value (the minimum length of a regime).
#' @param mat.date A numeric matrix of break dates.
#' @param d.Sym A logical value indicating whether the quantile range is symmetric satisfying \eqn{q.R=1-q.L}.
#' @param table.cv A matrix of simulated critical values for cases not covered by the response surface.
#'
#' @return A list containing:
#' \describe{
#'   \item{\code{test}}{A numeric vector of DQ test statistics.}
#'   \item{\code{cv}}{A numeric matrix of critical values for the DQ test at 10, 5, and 1 percent significance levels.}
#'   \item{\code{date}}{A numeric matrix of estimated break dates.}
#'   \item{\code{nbreak}}{A numeric vector indicating the number of detected breaks at different significance levels.}
#' }
#'
#' @examples
#'
#' \donttest{
#' # This example may take substantial time for automated package 
#' # checks since it involves dynamic programming
#' data(gdp)
#' y = gdp$gdp
#' x = gdp[,c("lag1", "lag2")]
#' n.size = 1
#' T.size = length(y) # number of time periods
#'
#' # setting
#' vec.tau   = seq(0.20, 0.80, by = 0.150)
#' trim.e    = 0.2
#' trim.size = round(T.size * trim.e)  #minimum length of a regime
#' m.max     = 3
#'
#' # dynamic program algorithm to compute the objective function
#' out.long   = gen.long(y, x, vec.tau, n.size, trim.size)
#' mat.long.s = out.long$mat.long  ## for break estimation using individual quantile
#' vec.long.m = out.long$vec.long  ## for break estimation using multiple quantiles jointly
#'
#' # break date
#' mat.date = brdate(y, x, n.size, m.max, trim.size, vec.long.m)
#'
#' ## qunatile ranges (left and right)
#' q.L   = 0.2
#' q.R   = 0.8
#' d.Sym = TRUE ## symmetric trimming of quantiles
#' table.cv = NULL ##covered by the response surface because d.Sym = TRUE
#'
#' # determine the number of breaks
#' out.m = dq(y, x, vec.tau, q.L, q.R, n.size, m.max, trim.size, mat.date, d.Sym, table.cv)
#'
#' # result
#' print(out.m)
#' }
#'
#'
#' @export

dq = function(y, x, vec.tau, q.L, q.R, n.size=1, m.max, trim.size, mat.date,d.Sym, table.cv)
{
    ## the number of regressors
    p.size = ncol(x) + 1  # add intercept

    vec.nb   = matrix(0, 3, 1)
    vec.test = matrix(0, 1, m.max)
    mat.cv   = matrix(0, 3, m.max)

    ## 0 vs 1 break
    vec.test[1] = dq.test.0vs1(y, x, q.L, q.R, n.size)
    mat.cv[,1]  = get.cv.dq(0, p.size, q.L, q.R,d.Sym,table.cv)
    for (a in 1:3){
        if (mat.cv[a,1] < vec.test[1]){
            vec.nb[a] = 1
        }
    }

    ## k vs k+1 break
    if (m.max >= 2){
        for (k in 1:(m.max-1) ) { # k is the number of breaks under the null hypothesis
            if (max(vec.nb) == k){

                ## the 1st to kth break dates are stored in the kth column
                vec.loc = mat.date[1:k,k]

                ## test: k vs k+1
                vec.test[(k+1)] = dq.test.lvsl_1(y, x, q.L, q.R, n.size, vec.loc)
                mat.cv[,(k+1)]  = get.cv.dq(k, p.size, q.L, q.R,d.Sym, table.cv)

                ## for each significance level
                for (a in 1:3){
                    if (vec.nb[a] == k  &  mat.cv[a,(k+1)] < vec.test[(k+1)]){
                        vec.nb[a] = vec.nb[a] + 1
                    }
                }
            }
        }
    }

    ## store optimal break dates
    mat.date.opt = matrix(0, 3, m.max)
    for(a in 1:3){
        nb = vec.nb[a]
        if (nb >=1){
            mat.date.opt[a,1:nb] = mat.date[1:nb,nb]
        }
    }

    ## return
    list(test = vec.test, cv = mat.cv, date = mat.date.opt, nbreak = vec.nb)
}
