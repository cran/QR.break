#' Subroutine for Confidence Interval Construction for Break Dates
#'
#' This subroutine constructs confidence intervals for break dates, based on
#' specified quantiles of interest. It is used within the `ci.date.m` function.
#' Related routines include \code{sample.split()} and \code{moment()}.
#'
#'
#' @param y.L A numeric vector of dependent variable observations on the left side of the break date.
#' @param x.L A numeric matrix of regressor observations on the left side of the break date.
#' @param y.R A numeric vector of dependent variable observations on the right side of the break date.
#' @param x.R A numeric matrix of regressor observations on the right side of the break date.
#' @param v.date A numeric value representing the estimated break date.
#' @param mat.size Estimated regression coefficients from the \code{rq.est.full()} function.
#' @param vec.tau A numeric vector of quantiles of interest for the confidence interval calculation.
#' @param n.size An integer specifying the size of the cross-sections (default is 1 for time series).
#' @param v.b A numeric value specifying the confidence interval level:
#'  \itemize{
#'    \item \code{v.b = 1} for the 90% confidence interval.
#'    \item \code{v.b = 2} for the 95% confidence interval.
#' }
#'
#' @return A numeric vector containing three values:
#' \describe{
#'    \item{1st column}{The estimated break date.}
#'    \item{2nd column}{The lower bound of the confidence interval.}
#'    \item{3rd column}{The upper bound of the confidence interval.}
#' }
#'
#'
#' @noRd


ci.date.m.sub = function(y.L, x.L, y.R, x.R, v.date, mat.size, vec.tau, n.size=1, v.b)
{
    ## size
    n.tau = length(vec.tau)

    ## Adjustment to use the shrinking-break asymptotic framework
    mat.size = mat.size * sqrt(n.size)

    ## vec.size
    mat.L = matrix(0, n.tau, 2)
    mat.R = matrix(0, n.tau, 2)

    ## taus
    for(i in 1:n.tau){
        v.tau    = vec.tau[i]
        vec.size = mat.size[,i]
        ## moment
        mL  = moment(y.L, x.L, v.tau)
        mR  = moment(y.R, x.R, v.tau)
        ## pi
        mat.L[i,1] = t(vec.size) %*% mL$H %*% vec.size
        mat.R[i,1] = t(vec.size) %*% mR$H %*% vec.size
        ## sigma
        for(s in 1:n.tau){
            ## size for s
            vec.size.S = mat.size[,s]
            ## mix: Note, J does not depdent on tau
            temp = min(v.tau, vec.tau[s]) - v.tau * vec.tau[s]
            ## - - difference m1$J and m2$J
            mat.L[i,2] = mat.L[i,2] + temp * t(vec.size) %*% mL$J %*% vec.size.S
            mat.R[i,2] = mat.R[i,2] + temp * t(vec.size) %*% mR$J %*% vec.size.S
        }
    }

    ## summation
    pi.L = sum(mat.L[,1])
    pi.R = sum(mat.R[,1])
    s2.L = sum(mat.L[,2])
    s2.R = sum(mat.R[,2])

    ## confidence interval
    vec.Q = c(7.7, 11.0)  ## 90% and 95% confidence levels
    vec.ci = matrix(0, 1, 3)
    vec.ci[1] = v.date
    vec.ci[2] = v.date - round( (vec.Q[v.b] * s2.L / (pi.L ^ 2) )) - 1
    vec.ci[3] = v.date + round( (vec.Q[v.b] * s2.R / (pi.R ^ 2) )) + 1

    return(vec.ci)
}


