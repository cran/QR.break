#' Critical Value Simulation for the SQ Test
#'
#' This function simulates critical values for the SQ test based on the specified number
#' of breaks and the number of coefficients allowed to change.
#'
#'
#' @param p.size An integer specifying the total number of coefficients allowed to change.
#' @param m.max An integer specifying the maximum number of breaks allowed under the alternative hypothesis. The tests evaluates l breaks against l+1 breaks for l between 0 and m.max-1.
#'
#' @details
#' This function calculates the critical values at three significance levels:
#' 10%, 5%, and 1%.
#'
#' @return
#' A numeric matrix (`table.cv`) containing the critical values for the SQ test. The matrix consists of:
#' \itemize{
#'   \item The first row contains the number of coefficients (\eqn{p.size}).
#'   \item The second row contains the critical values for the 10% significance level.
#'   \item The third row contains the critical values for the 5% significance level.
#'   \item The fourth row contains the critical values for the 1% significance level.
#'   \item The column index corresponds to the number of breaks under the alternative hypothesis.
#' }
#'
#'
#' @examples
#' \donttest{
#'   # Simulating critical values (time-consuming)
#'   p.size = 5  # Example number of regressors
#'   m.max  = 3  # Example number of breaks
#'   result = cv.sq(p.size, m.max)
#'   print(result)
#' }
#'
#' @references
#' Qu, Z. (2008). Testing for Structural Change in Regression Quantiles.
#' \emph{Journal of Econometrics}, 146(1), 170-184.
#'
#' Oka, T. and Z. Qu (2011).
#' Estimating Structural Changes in Regression Quantiles.
#' \emph{Journal of Econometrics}, 162(2), 248–267.
#'
#'
#' @noRd


cv.sq = function(p.size, m.max)
{

    tau = 0.5 #set to the median. The choice of quantile does not affect the critical value

    ## start time
    time = proc.time()

    ## setting
    n.grid   = max(1000,50*p.size)
    n.sim    = 500000

    ## significance level
    vec.a = c(0.90, 0.95, 0.99)
    n.a   = length(vec.a)

    vec.nn = seq(1, m.max, by = 1)
    mat.prob = matrix(0, 3, m.max)
    mat.prob[1,] = 0.90 ^ (1 / vec.nn)
    mat.prob[2,] = 0.95 ^ (1 / vec.nn)
    mat.prob[3,] = 0.99 ^ (1 / vec.nn)
    mat.index = round(mat.prob * n.sim)

    n.tau    = length(tau)

    ## seed
    ## set.seed(7, kind = NULL) #removed

    ## simulation
    vec.p  = vec.prob.DQ(n.grid, p.size, n.sim, tau) # this subroutine works for both sq and dq.
    vec.p  = sort(vec.p)
    mat.cv   = matrix(0, 3, m.max)
    vec.name = matrix(p.size, 1, m.max)

    for(k in 1:3){
        mat.cv[k,] = vec.p[mat.index[k,]]
    }

    if (length(tau)==1){
    ## this is for the SQ test
    table.cv = rbind( vec.name, mat.cv/sqrt(tau*(1-tau)) ) #this differs from DQ
    #file.name = paste('table.cv.SQ.new', p.size, '.txt', sep='')
    #write.table(table.cv,  file = file.name)
    }
    else {
    stop("The number of quantiles must be equal to 1")

    }

    ## end time
    ##cat('Critical values for SQ tests for up to', m.max, 'breaks, allowing', p.size, "coefficients to change. \n")
    ## cat("rows: the number of coefficients, followed by critical values at significant levels", vec.a, "\n");
    ## cat("columns: the number of breaks under the alternative hypothesis\n")
    ## cat('Time used for simulating the critical values',   proc.time() - time, '\n')
   return(table.cv)
}

