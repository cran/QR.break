#' Values of the Objective Function for All Segments within the Date Range
#'
#' This procedure calculates the values of the objective function for all segments
#' within a specified range of dates. 
#'
#' @param y A numeric vector of dependent variables (\eqn{NT \times 1}).
#' @param x A numeric matrix of regressors (\eqn{NT \times p}).
#' @param vec.tau A numeric vector of quantiles of interest.
#' @param n.size An integer specifying the size of the cross-section. Default is `1`.
#' @param trim.size A numeric value specifying the trimming size (the minimum length of a segment).
#' @param start An integer representing the starting date of the range.
#' @param last An integer representing the ending date of the range.
#'
#' @return A matrix of objective function values for all segments within the date range.
#'
#' @noRd

gen.mat.rho = function(y, x, vec.tau, n.size=1, trim.size, start, last)
{
    ## number of quantiles
    n.tau = length(vec.tau)

    ## storage
    mat.rho = matrix(0, last, n.tau)

    ## quantile regression
    col01 = n.size * (start - 1) + 1
    loc01 = start + trim.size - 1      ## the start
    loc02 = last                       ## the end
    for(j in loc01:loc02){
        col02 = n.size * j

        ## calcuate rho for each quntile
        temp.rho = matrix(0, n.tau, 1)
        for(k in 1:n.tau){
            y.temp = y[col01:col02]
            x.temp = x[col01:col02,]
            mat.rho[j,k] = rq(y.temp ~ x.temp, tau = vec.tau[k])$rho
        }
    }

    ## return
    return(mat.rho)
}

