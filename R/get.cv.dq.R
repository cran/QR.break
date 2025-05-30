#' Critical Values for the DQ Test
#'
#' This function calculates critical values for the sequential DQ test,
#' which tests the null hypothesis of L breaks against the alternative
#' hypothesis of L+1 breaks based on multiple quantiles.
#'
#'
#' @param n.break The number of breaks under the null hypothesis.
#' @param p.size The number of regressors (including the intercept).
#' @param q.L The lower end of the quantile range.
#' @param q.R The upper end of the quantile range.
#' @param d.Sym A logical value indicating whether the response surface method is used; if d.Sym=FALSE, the critical values are generated by simulations.
#' @param table.cv A matrix containing simulated critical values for cases not covered by the response surface.
#'
#' @return A numeric vector of length 3 containing critical values at the 10%, 5%, and 1% significance levels.
#'
#'
#' @noRd
#'
get.cv.dq = function(n.break, p.size, q.L, q.R, d.Sym, table.cv)
{
    ## row.beg = 4 * (p.size-1) + 4 # the first three lines include q.L, q.R &  #test
    ## row.end = 4 * (p.size-1) + 6
    ## vec.cv  = mat.cv[row.beg:row.end, n.break]

    if (d.Sym == TRUE){

        ## response surface
        vec.cv = res.surface(p.size, n.break, q.L, q.R, d.Sym)

    } else {

        #file.name = paste('table.cv.DQ', q.L, '_', q.R , '_', p.size, '.txt', sep='') #can save to file
        #cat('Note: DQ test uses simulated critical values, saved in', file.name, ' \n')
        #mat.cv = as.matrix(read.table(file.name, header=T), ncol = 3)
        vec.cv = table.cv[4:6,(n.break+1)]

    }

    ## return
    return(vec.cv)
}

