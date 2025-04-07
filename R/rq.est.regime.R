#' Regime-Specific Coefficients and Confidence Intervals Given Break Dates
#'
#' This function estimates the coefficients for each regime, given the break dates.
#'
#'
#' @param y A vector of dependent variables (\eqn{NT \times 1}).
#' @param x A matrix of regressors (\eqn{NT \times p}).
#' @param v.tau The quantile of interest.
#' @param vec.date A vector of estimated break dates.
#' @param n.size The cross-sectional sample size (\eqn{N}).
#'
#' @return A list containing the estimated coefficients for each regime.
#'
#' @examples
#' ## data
#' data(gdp)
#' y        = gdp$gdp
#' x        = gdp[,c("lag1", "lag2")]
#'
#' ## quantile
#' v.tau = 0.8
#'
#' ## break date
#' vec.date = 146
#'
#' # cross-sectional size
#' n.size = 1
#'
#' ## estimation
#' result = rq.est.regime(y, x, v.tau, vec.date, n.size)
#' print(result)
#'
#' @export

rq.est.regime = function(y, x, v.tau, vec.date, n.size=1)
{

    ## the number of breaks
    n.break  = length(vec.date)

    ## Initialize an empty list to store results
    out <- vector("list", n.break + 1)  # Ensure enough space for all regimes
    names(out) <- paste0("Regime_", seq_len(n.break + 1))  # Name list elements

    ## sequentially split sample
    rem.y  = y
    rem.x  = x
    pre.date = c(0, vec.date)
    for(j in 1:n.break){
        #cat('  < Regime:', j,'>\n')
        v.date = vec.date[j] - pre.date[j]
        ## split sample given a break date
        temp   = sample.split(rem.y, rem.x, v.date, n.size)
        rem.y  = temp$y2
        rem.x  = temp$x2
        ## estimation
        fit.regime  = rq(temp$y1 ~ temp$x1, v.tau)
        sum.regime  = summary.rq(fit.regime, se = "nid", covariance = TRUE)

        ## print results; suppressed
        rownames(sum.regime$coef) <- c("Intercept", paste0("x", seq_len(nrow(sum.regime$coef) - 1)))
        #print(format(sum.regime$coef, digits = 4), quote = FALSE)
        out[[paste0("Regime_", j)]] <- sum.regime$coef
        #cat('\n')
    }

    ## the final regime
    #cat('  < Regime:', (n.break + 1), '>\n')
    fit.regime  = rq(rem.y~rem.x, v.tau)
    sum.regime  = summary.rq(fit.regime, se = "nid", covariance = TRUE)
    rownames(sum.regime$coef) <- c("Intercept", paste0("x", seq_len(nrow(sum.regime$coef) - 1)))

    #print(format(sum.regime$coef, digits = 4), quote = FALSE)

    out[[paste0("Regime_", n.break + 1)]] <- sum.regime$coef
    #print(format(sum.regime$coef, digits = 4), quote = FALSE)
    #cat('\n')
    return(out)
}

