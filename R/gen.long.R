#' Dynamic Programming Algorithm
#'
#' This function computes the objective function values for all possible segments of the sample.
#'
#'
#' @param y A numeric vector of dependent variables (\code{NT × 1}).
#' @param x A numeric matrix of regressors (\code{NT × p}).
#' @param vec.tau A vector of quantiles of interest.
#' @param n.size The size of the cross-section; default is set to 1.
#' @param trim.size The minimum length of a regime (integer).
#'
#' @return A list containing:
#'   \describe{
#'     \item{mat.long}{A matrix of objective function values for separate quantiles.}
#'     \item{vec.long}{A matrix of objective function values for combined quantiles.}
#'   }
#'
#' @references
#' Bai, J and P. Perron (2003).
#' Computation and Analysis of Multiple Structural Change Models.
#' *Journal of Applied Econometrics*, 18(1), 1-22.
#'
#' @examples
#'
#' \donttest{
#' # This example may take substantial time for automated package checks
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
#'
#' out.long   = gen.long(y, x, vec.tau, n.size, trim.size)
#'
#' }
#'
#' @export

gen.long = function(y, x, vec.tau, n.size=1, trim.size)
{

    t.size = length(y) / n.size
    n.tau  = length(vec.tau)
    x = as.matrix(x)

    ## The following matrix, mat.long, and vector, vec.long, contain the objective function
    ## corresponding to all possible segments. The objective function value for the
    ## segment starting at j and lasting for k periods is stored as the
    ## T(j-1) - (j-1)(j-2)/2 + k-th element of the vector/matrix. The matrix and vector
    ## correspond to separate quantiles and the combined quantiles, respectively.

    mat.long = matrix(0, (t.size*(t.size+1)/2), n.tau)  ## separate quantiles
    vec.long = matrix(0, (t.size*(t.size+1)/2), 1)      ## combined multiple quantiles

    for(i in 1:(t.size-trim.size+1)){
        mat.out = gen.mat.rho(y, x, vec.tau, n.size, trim.size, i, t.size)

        ## rows
        row01 = (i-1) * t.size + i - (i - 1) * i / 2
        row02 =  i    * t.size     - (i - 1) * i / 2

        ## store
        mat.long[row01:row02,] = mat.out[i:t.size,]
        #vec.long[row01:row02]  = rowSums(mat.out[i:t.size,]) #for segments starting at i, ending at i, i+1,...t.size; the first  trim.size - 1 values are zeros --Qu
        if (n.tau > 1){
          vec.long[row01:row02]  = rowSums(mat.out[i:t.size,]) #for segments starting at i, ending at i, i+1,...t.size; the first  trim.size - 1 values are zeros --Qu

        }
    }

    ## return
    list(mat.long = mat.long, vec.long = vec.long)
}
