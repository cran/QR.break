#' Determine the Number of Breaks Using the SQ(l|l+1) Test
#'
#' This procedure sequentially applies the SQ test to determine the number of breaks, based on a single quantile.
#'
#' @param y A numeric vector of dependent variables (\eqn{NT \times 1}).
#' @param x A numeric matrix of regressors (\eqn{NT \times p}).
#' @param v.tau A numeric value representing the quantile of interest.
#' @param n.size An integer specifying the size of the cross-section (\eqn{N}).
#' @param m.max An integer specifying the maximum number of breaks allowed.
#' @param trim.size A numeric value specifying the trimming size (the minimum length of a segment).
#' @param mat.date A numeric matrix of break dates.
#'
#' @return A list with the following components:
#' \describe{
#'   \item{\code{test}}{A numeric vector of SQ test statistics.}
#'   \item{\code{cv}}{A numeric matrix of critical values for the SQ test, with the 1st, 2nd, and 3rd rows corresponding to the 10%, 5%, and 1% significance levels.}
#'   \item{\code{date}}{A numeric matrix of break dates, with the 1st, 2nd, and 3rd rows corresponding to the 10%, 5%, and 1% significance levels.}
#'   \item{\code{nbreak}}{A numeric vector indicating the number of breaks at the 10%, 5%, and 1% significance levels.}
#' }
#'
#' @examples
#' \donttest{
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
#' # the maximum number of breaks 
#' m.max = 3
#' 
#' ## trim 
#' T.size    = length(y)
#' trim.e    = 0.2
#' trim.size = round(T.size * trim.e)  #minimum length of a regime
#' 
#' # get.long 
#' out.long   = gen.long(y, x, v.tau, n.size, trim.size)
#' mat.long.s = out.long$mat.long  ## for individual quantile
#' 
#' # mat.date 
#' mat.date = brdate(y, x, n.size, m.max, trim.size, mat.long.s)
#' 
#' # sq 
#' result = sq(y, x, v.tau, n.size, m.max, trim.size, mat.date)
#' print(result)
#' }
#'
#' @export

sq = function(y, x, v.tau, n.size=1, m.max, trim.size, mat.date)
{
    ## matrix 
    x = as.matrix(x)
  
  
    ## the number of regressors
    p.size = ncol(x) + 1  # plus intercept

    ## stroage
    vec.nb       = matrix(0, 3, 1)
    vec.test     = matrix(0, 1, m.max)
    mat.cv       = matrix(0, 3, m.max)

    ## 0 vs 1
    vec.test[1] = sq.test.0vs1(y, x, v.tau, n.size)
    mat.cv[,1]  = get.cv.sq(0, p.size) # critical value
    for (a in 1:3){
        if (mat.cv[a,1] < vec.test[1]){
            vec.nb[a] = 1
        }
    }

    ## k vs k+1
    if (m.max >= 2){
        for(k in 1:(m.max-1)){
            if (max(vec.nb) == k){

                ## the 1st to kth break dates are stored in the kth column
                vec.loc = mat.date[1:k,k]

                ## test: k vs k+1
                vec.test[(k+1)] = sq.test.lvsl_1(y, x, v.tau, n.size, vec.loc)
                mat.cv[,(k+1)]  = get.cv.sq(k, p.size) # critical value

                ## for each significance level
                for (a in 1:3){
                    if (vec.nb[a] == k & mat.cv[a,(k+1)] < vec.test[(k+1)]){
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
        if(nb >= 1){
            mat.date.opt[a,1:nb] = mat.date[1:nb,nb]
        }
    }

    ## return
    list(test = vec.test, cv = mat.cv, date = mat.date.opt, nbreak = vec.nb)
}

