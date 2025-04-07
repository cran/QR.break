#' Sample Splitting Subroutine
#'
#' This procedure splits the sample into two subsamples based on a specified break date.
#'
#'
#' @param y A numeric vector of dependent variables (\eqn{NT \times 1}).
#' @param x A numeric matrix of regressors (\eqn{NT \times p}).
#' @param v.date A numeric value representing the break date for splitting the sample.
#' @param n.size An integer specifying the size of the cross-section (\eqn{N}).
#'
#' @return A list with the following components:
#' \describe{
#'   \item{\code{y1}}{A numeric vector of dependent variables in the first subsample.}
#'   \item{\code{x1}}{A numeric matrix of regressors in the first subsample.}
#'   \item{\code{y2}}{A numeric vector of dependent variables in the second subsample.}
#'   \item{\code{x2}}{A numeric matrix of regressors in the second subsample.}
#' }
#'
#'
#' @noRd

sample.split = function(y, x, v.date, n.size=1)
{
    
    nt.01   = v.date * n.size ## the end of a regime
    nt.size = length(y)

    ## split 
    y1 = as.matrix( y[1:nt.01]  )
    x1 = as.matrix( x[1:nt.01,] )  
    y2 = as.matrix( y[(nt.01+1):nt.size]  )
    x2 = as.matrix( x[(nt.01+1):nt.size,] ) 

    ## return
    list(y1 = y1, x1 = x1, y2 = y2, x2 = x2)
}


