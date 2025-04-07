#' Simulate Critical Values for the DQ Test
#'
#' This subroutine, used in \code{critical.DQtest_specific()}, simulates critical values
#' for the DQ test.
#'
#' @param n.grid An integer specifying the number of grid points over the quantile range.
#' @param n.dim An integer indicating the number of parameters allowed to change.
#' @param n.sim An integer specifying the number of simulation replications.
#' @param vec.tau An vector of quantiles of interest.
#'
#' @return A numeric vector (`vec.p`) containing the simulated values needed to compute the critical values for the DQ test.
#'
#'
#' @noRd

vec.prob.DQ = function(n.grid, n.dim, n.sim, vec.tau)
{
    ## number of quantiles
    n.tau = length(vec.tau)

    ## simulation
    vec.lamda = seq(0, 1, length = n.grid)
    vec.p     = matrix(0, n.sim, 1)
    for(s in 1:n.sim){

        ## Over dimension
        vec.max = matrix(0, n.dim, 1)
        for(d in 1:n.dim){

            ## e ~ U[0,1]
            vec.e = runif(n.grid, 0, 1)

            vec.limit = matrix(0, n.tau, 1)
            for(i in 1:n.tau){

                ## value tau
                v.tau = vec.tau[i]

                ## brownian motion
                vec.bm = cumsum( (vec.e <= v.tau) )

                ## brownian bridge
                vec.BB = (vec.bm - vec.lamda * vec.bm[n.grid]) / sqrt(n.grid)

                ## max over lamda
                vec.limit[i] = max(abs(vec.BB))
            }

            ## max over tau
            vec.max[d] = max(vec.limit)
        }

        ## maximization over dimension
        vec.p[s] = max(vec.max)
    }

    ## return
    return(vec.p)
}

