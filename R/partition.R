#' Optimal Single-Break Partition
#'
#' @description This function finds the optimal single-break partition within a segment 
#' by evaluating the objective function at each possible break point.
#'
#' @param vec.long A numeric vector containing values of the objective function.
#' @param start An integer specifying the start of the segment under consideration.
#' @param b1 An integer indicating the first possible break date.
#' @param b2 An integer indicating the last possible break date.
#' @param last An integer specifying the end of the segment under consideration.
#' @param t.size An integer representing the size of the time series.
#'
#' @return A list with the following components:
#' \describe{
#'   \item{\code{loc.min}}{The optimal break point.}
#'   \item{\code{rho.min}}{The value of the objective function at the optimal break point.}
#' }
#'
#' @details The function evaluates potential break points within the given range \code{b1:b2}. 
#' It calculates the sum of two objective function values corresponding to the segments 
#' before and after the break. The break point that minimizes this sum is returned as loc.min.
#'
#'
#' @noRd

partition = function(vec.long, start, b1, b2, last, t.size)
{
    ## initialization 
    ini = (start - 1) * t.size - (start - 2) * (start - 1) / 2 + 1

    ## storage
    temp.rho = matrix(0, t.size, 1)

    for(j in b1:b2){
        
        loc01 = j - start + ini
        loc02 = j * t.size - (j - 1) * j / 2 + last - j

        ## 1st term: rho of [start,j] --- the "break" is the last day of the previous regime (not the first day of the new regime)
        ## 2nd term: rho of [(j+1),last-j] -- typo, should be [j+1,last] - Qu
        temp.rho[j] = vec.long[loc01] + vec.long[loc02]
    }

    ## minimization
    loc.min = (b1 - 1) + which.min(temp.rho[b1:b2])
    rho.min = temp.rho[loc.min]

    ## return
    list(loc.min = loc.min, rho.min = rho.min)
}
