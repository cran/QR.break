#' Critical Value Simulation for the DQ Test
#'
#' @param x A numeric matrix of regressors (excluding the intercept) (\eqn{NT \times p}).
#' @param m.max An integer specifying the maximum number of breaks allowed.
#' @param vec.tau A numeric vector of quantiles of interest.
#'
#' @return
#' A numeric matrix containing the simulated critical values for the DQ test. The matrix has \code{4 + 3} rows and \code{m.max} columns, structured as follows:
#' \itemize{
#'   \item Row 1: Minimum quantile value (\code{min(vec.tau)}), repeated across columns.
#'   \item Row 2: Maximum quantile value (\code{max(vec.tau)}), repeated across columns.
#'   \item Row 3: Number of regressors (including intercept), repeated across columns.
#'   \item Rows 4-6: Simulated critical values corresponding to significance levels 0.90, 0.95, and 0.99, respectively.
#' }
#'
#' @references
#' Qu, Z. (2008).
#' Testing for Structural Change in Regression Quantiles.
#' \emph{Journal of Econometrics}, \emph{146}(1), 170–184.
#
#'
#' @noRd

critical.DQtest_specific = function(x, m.max, vec.tau)
{
    ## keep track of computation time
    time = proc.time()

    ## the number of regression coefficients allowed to change
    p.size = ncol(x) + 1 ## include intercept

    ## creates a grid over quantiles
    n.grid   = 500
    n.sim    = 50000 ## number of simulation replications to compute the critical value

    ## significance level
    vec.a = c(0.90, 0.95, 0.99)
    n.a   = length(vec.a)

    vec.nn = seq(1, m.max, by = 1)
    mat.prob = matrix(0, 3, m.max)
    mat.prob[1,] = 0.90 ^ (1 / vec.nn)
    mat.prob[2,] = 0.95 ^ (1 / vec.nn)
    mat.prob[3,] = 0.99 ^ (1 / vec.nn)
    mat.index = round(mat.prob * n.sim)

    ## discretization of the quantile range
    beg.tau   = min(vec.tau)
    end.tau   = max(vec.tau)
    cont.tau = seq(beg.tau, end.tau, by = 1 / n.grid)
    n.tau    = length(cont.tau)

    ## table to save critical values
    table.cv     = matrix(0, 2, m.max)
    table.cv[1,] = beg.tau
    table.cv[2,] = end.tau

    ## Set random seed for reproducibility
    ## set.seed(7, kind = NULL) #removed

    ## simulation
    vec.p  = vec.prob.DQ(n.grid, p.size, n.sim, cont.tau)
    vec.p  = sort(vec.p)
    mat.cv   = matrix(0, 3, m.max)
    vec.name = matrix(p.size, 1, m.max)

    for(k in 1:3){
        mat.cv[k,] = vec.p[mat.index[k,]]
    }

    table.cv = rbind(table.cv, vec.name, mat.cv)

    ## save the values in a file for easy access
    #file.name = paste('table.cv.DQ', beg.tau, '_', end.tau , '_', p.size, '.txt', sep='')
    #write.table(table.cv,  file = file.name)
    ## end of computation
    ## cat('Time used for simulating the critical values',   proc.time() - time, '\n')
    return(table.cv)
}

