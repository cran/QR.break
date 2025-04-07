#' Estimating Break Dates in a Quantile Regression
#'
#' This function estimates break dates in a quantile regression model,
#' allowing for up to `m` breaks.
#' When `m = 1`, this function finds the optimal single-break partition within
#' a segment by evaluating the objective function at each possible break point
#' to determine the break date.
#' When `m > 1`, a dynamic programming algorithm is used to efficiently determine
#' the break dates.
#'
#' @param y A numeric vector of dependent variables (\eqn{NT \times 1}).
#' @param x A numeric matrix of regressors (\eqn{NT \times p}); a column of ones will be automatically added to \eqn{x}.
#' @param n.size An integer specifying the number of cross-sectional units (\eqn{N}), equal to 1 for time series data.
#' @param m An integer specifying the maximum number of breaks allowed.
#' @param trim.size An integer representing the minimum length of a regime, which ensures break dates are not too short or too close to the sample's boundaries.
#' @param vec.long A numeric vector/matrix used in dynamic programming, storing pre-computed objective function values for all possible segments of the sample for optimization.
#'                 Produced by the function `gen.long()`.
#'
#' @return An upper triangular matrix (`m × m`) containing estimated break dates. The row index represents break dates, and the column index corresponds to the permitted number of breaks before the ending date.
#'
#' @details
#' The function first determines the optimal one-break partition.
#' If multiple breaks are allowed (`m > 1`), it applies a dynamic programming
#' algorithm as in Bai and Perron (2003) to minimize the objective function.
#' The algorithm finds break dates by iterating over all possible partitions,
#' returning optimal break locations and associated objective function values.
#'
#'
#' @references
#' Bai, J. and P. Perron (2003).
#' Computation and analysis of multiple structural change models.
#' *Journal of Applied Econometrics*, 18(1), 1-22.
#'
#' Oka, T. and Z. Qu (2011). Estimating structural changes in regression quantiles,
#' *Journal of Econometrics*, 162(2), 248-267.
#'
#' @export
#'
#'
#' @examples
#' \donttest{
#'   # data
#'   data(gdp)
#'   y = gdp$gdp
#'   x = gdp[,c("lag1", "lag2")]
#'   n.size = 1
#'   T.size = length(y) # number of time periods
#'
#'   # setting
#'   vec.tau   = seq(0.20, 0.80, by = 0.150)
#'   trim.e    = 0.2
#'   trim.size = round(T.size * trim.e)  #minimum length of a regime
#'   m.max     = 3
#'
#'   # compute the objective function under all possible partitions
#'   out.long   = gen.long(y, x, vec.tau, n.size, trim.size)
#'   vec.long.m = out.long$vec.long  ## for break estimation using multiple quantiles
#'
#'   # break date
#'   mat.date = brdate(y, x, n.size, m.max, trim.size, vec.long.m)
#'   print(mat.date)
#' }

brdate = function(y, x, n.size=1, m, trim.size, vec.long)
{
    ## number of time periods
    t.size = length(y) / n.size

    mat.date = matrix(0, m, m)
    ## an upper-triangular matrix containing the estimated break dates
    ## from one to m

    mat.opt.loc = matrix(0, t.size, m)
    ## a row index corresponding to the ending dates, column index corresponds
    ## to the number of breaks permitted before the ending date the cell
    ## contains the break date

    mat.opt.rho = matrix(0, t.size ,m)
    ## same as above, the cell contains the minimal value of the objective
    ## function corresponding to that break dates

    dvec = matrix(0, t.size,1)
    ## the value is the date after which we inserting the break point. The
    ## cell contains the corresponding value of the objective function

    vec.global = matrix(0, m, 1)
    ## Global minimum when i breaks are allowed

    ## m = 1
    if (m == 1){
        res01 = partition(vec.long, 1, trim.size, t.size-trim.size, t.size, t.size)
        mat.date[1,1] = res01$loc.min
        vec.global[1] = res01$rho.min
    } else {
        ## when m > 1, a dynamic programming algorithm is used to find the global minimum.
        ## The first step is to obtain the optimal one-break partitions for all
        ## possible ending dates from 2h to T-mh+1.
        ## The optimal dates are stored in a vector optdat.
        ## The associated objective functions are stored in a vector optssr.

        ## First loop. Looking for the optimal one-break partitions for break
        ## dates between h and T-h; j1 is the last date of the segment.
        for (j1 in (2*trim.size):t.size) {    #starting date,  last possible break, ending date.
            res02 = partition(vec.long, 1, trim.size, (j1-trim.size), j1, t.size)
            mat.opt.rho[j1,1] = res02$rho.min         # not a typo
            mat.opt.loc[j1,1] = res02$loc.min

            vec.global[1] = mat.opt.rho[t.size,1]     # not a typo
            mat.date[1,1] = mat.opt.loc[t.size,1]


            ## Next, the algorithm looks for optimal 2,3,... breaks partitions
            ## The index used is ib.
            for (ib in 2:m){
                if (ib == m){
                    ## if we have reached the maximum number of breaks allowed,
                    ## then do the following
                    jlast = t.size

                    beg01 = ib * trim.size
                    end01 = t.size - trim.size
                    for (jb in beg01:end01){
                        dvec[jb] = mat.opt.rho[jb,ib-1] +
                            vec.long[(jb+1)*t.size-jb*(jb+1)/2]
                    }
                    min.loc = (beg01 - 1) + which.min(dvec[beg01:end01])
                    mat.opt.rho[jlast,ib] = dvec[min.loc]
                    mat.opt.loc[jlast,ib] = min.loc
                } else {
                    ## if we have not reached the maximum number of breaks allowed,
                    ## we loop over the possible last dates of the segments,
                    ## between (ib+1)*h and T.
                    for(jlast in ((ib+1)*trim.size):t.size){
                        beg01 = ib * trim.size
                        end01 = jlast - trim.size
                        for(jb in beg01:end01){
                            dvec[jb] = mat.opt.rho[jb,ib-1] +
                                vec.long[jb*t.size-jb*(jb-1)/2+jlast-jb]
                        }
                        min.loc = (beg01 - 1) + which.min(dvec[beg01:end01])
                        mat.opt.rho[jlast,ib] = dvec[min.loc]
                        mat.opt.loc[jlast,ib] = min.loc
                    }
                }

                mat.date[ib,ib] = mat.opt.loc[t.size,ib]
                for(i in 1:(ib-1)){
                    xx = ib - i
                    mat.date[xx,ib] = mat.opt.loc[mat.date[xx+1,ib],xx]
                }
                vec.global[ib] = mat.opt.rho[t.size,ib]
            }
        }         # closing the 'if' for the case m > 1

    }

    ## return
    return(mat.date)
}


