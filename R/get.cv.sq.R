#' Critical Values for the SQ Test
#'
#' @description 
#' This function extracts critical values for the SQ test based on a single quantile, 
#' given the number of breaks under the null hypothesis and the number of regressors.
#'
#' @param n.break An integer specifying the number of breaks under the null hypothesis.
#' @param p.size An integer representing the number of regressors.
#'
#' @return A numeric vector of length 3 containing critical values at the 
#' 10\%, 5\%, and 1\% significance levels.
#'
#' @details 
#' The function retrieves critical values from an internally stored 
#' matrix of critical values, computed using \code{input.all.cv()}. 
#' The rows selected correspond to the given number of regressors \code{p.size}, 
#' and the column corresponds to \code{n.break + 1}.
#'
#' @noRd
#' 
get.cv.sq = function(n.break, p.size){
    
    mat.cv = as.matrix(input.all.cv())

    # mat.cv = as.matrix(read.table("table.cv.SQtest010910.txt", header=T), ncol = 5)
    
    row.beg = 4 * (p.size-1) + 1 + 1 # the first line includes the number of test 
    row.end = 4 * (p.size-1) + 4

    vec.cv  = mat.cv[row.beg:row.end, (n.break+1)]

    ## return
    return(vec.cv)
}

