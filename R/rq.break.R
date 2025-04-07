#' Testing for Breaks and Estimating Break Dates and Sizes with Confidence Intervals
#'
#' @description
#' This is the main function of this package for testing breaks in quantile regression models
#' and estimating break dates and break sizes with corresponding confidence intervals.
#'
#' @usage rq.break(y, x, vec.tau, N, trim.e, vec.time, m.max, v.a, v.b, verbose)
#'
#' @param y a numeric vector, the outcome variable (NT x 1), the first N units are from the first period,
#' the next N from the second period, and so forth.
#' @param x A matrix of regressors (NT x p), structured in the same way as y, a column of ones will be automatically added to x.
#' @param vec.tau a numeric vector, quantiles used for break estimation, for example \code{vec.tau = seq(0.20, 0.80, by = 0.10)}
#' @param N a numeric vector, the number of cross-sectional units. Set to 1 for a time series quantile regression.
#' @param trim.e a scalar between 0 and 1, the trimming proportion.
#' For example, if \code{trim.e=0.1}, the minimum regime length is 0.1 times the data span.
#' @param vec.time a vector of dates, needed for reporting the estimated break dates, in the format of (starting date...ending date);
#' If set to NULL, the break dates will be reported as indices (e.g., 55 for the 55th observation in the sample).
#' @param m.max the maximum number of breaks allowed.
#' @param v.a the significance level used for determining the number of breaks; 1, 2 or 3 for 10%, 5% or 1%, respectively
#' @param v.b the coverage level for constructing the confidence intervals of break dates; 1 or 2 for 90% and 95%, respectively.
#' @param verbose Logical; set to TRUE to print estimates to the console. Default is FALSE.
#'
#' @return A list containing:
#' - `$s.out`: A list with break testing results, estimated break dates, confidence intervals, and coefficient estimates
#'   based on individual quantiles.
#' - `$m.out`: A list with break testing results, estimated break dates, confidence intervals, and coefficient estimates
#'   obtained by testing and estimating breaks using multiple quantiles simultaneously.
#'
#' Each list (`s.out` or `m.out`) contains:
#' - `test_tau`: A matrix of test statistics and critical values for break detection at quantile `tau`.
#' - `nbreak_tau`: The number of detected breaks at quantile `tau`.
#' - `br_est_tau`: A matrix of estimated break dates and their confidence intervals at quantile `tau`.
#' - `br_est_time_tau`: The same as `br_est_tau`, but with break dates reported in calendar format (if `vec.time` is provided and is not NULL).
#' - `coef_tau`: Estimated regression coefficients for each regime at quantile `tau`.
#' - `bsize_tau`: Break size estimates for each transition between regimes at quantile `tau`.
#'
#'
#' @import quantreg
#' @importFrom stats runif
#' @importFrom utils read.table write.table
#'
#' @references
#' Koenker, R. and G. Bassett Jr. (1978).
#' Regression quantiles. \emph{Econometrica}, 46(1), 33-50.
#'
#' Oka, T. and Z. Qu (2011). Estimating Structural Changes in Regression Quantiles.
#' \emph{Journal of Econometrics}, 162(2), 248-267.
#'
#' Qu, Z. (2008). Testing for Structural Change in Regression Quantiles.
#' \emph{Journal of Econometrics}, 146(1), 170-184.
#'
#' @examples
#' \donttest{
#' ## Example 1
#' ## Time series example, using GDP data
#' ## data
#' data(gdp)
#' y        = gdp$gdp
#' x        = gdp[,c("lag1", "lag2")]
#' vec.time = gdp$yq
#'
#' ## the maximum number of breaks allowed
#' m.max = 3
#'
#' ## the signifance level for sequenatial testing
#' ## 1, 2 or 3 for 10%, 5% or 1%, respectively
#' v.a = 2
#'
#' ## the significance level for the confidence intervals of estimated break dates.
#' ## 1 or 2 for 90% and 95%, respectively.
#' v.b = 2
#'
#' ## the size of the cross-section
#' N = 1
#'
#' ## the trimming proportion for estimating the break dates
#' ## (represents the minimum length of a regime; used to exclude
#' ## the boundaries of the sample)
#' trim.e = 0.15
#'
#' ## quantiles
#' vec.tau = seq(0.20, 0.80, by = 0.150)
#'
#' verbose = FALSE #do not print
#'
#' ## main estimation
#' res = rq.break(y, x, vec.tau, N, trim.e, vec.time, m.max, v.a, v.b, verbose)
#'
#' }
#'
#' @export

rq.break = function(y, x, vec.tau, N=1, trim.e, vec.time, m.max, v.a, v.b, verbose=FALSE)
{

    x <- as.matrix(x)

    if (nrow(x) != length(y)) {
      stop("Error: The number of rows in x must match the length of y.")
    }

    if (length(vec.tau) == 0) {
      stop("Error: vec.tau must contain at least one quantile for estimation.")
    }

    if (N < 1) {
      stop("Error: N must be at least 1.")
    }

    if (trim.e < 0 || trim.e > 0.5) {
      stop("Error: trim.e must be between 0 and 0.5 (inclusive).")
    }

    if (m.max*trim.e > 1) {
      stop("Error: m.max*trim.e exceeds 1. This occurs because too many regimes are allowed or the minimum length of a regime is too large. Consider decreasing m.max, trim.e, or both.")
    }

    if (!(v.a %in% c(1, 2, 3))) {
      stop("Error: v.a must be 1, 2, or 3 to select a significance level of 10%, 5%, or 1%.")
    }

    if (!(v.b %in% c(1, 2))) {
      stop("Error: v.b must be 1 or 2 to select a significance level of 10% or 5%.")
    }

    n.size    = N                       # cross section size, 1 for time series data
    N.tau     = length(vec.tau)         #number of quantiles used in the estimation
    T.size    = length(y) / n.size      # number of time periods
    trim.size = round(T.size * trim.e)  #minimum length of a regime
    p.size = ncol(x) + 1                # user specified regressors plus an intercept

    if (!is.null(vec.time) && length(vec.time) != T.size) {
      stop("Error: The length of vec.time must match the length of y.")
    }

    if (trim.e * nrow(x) < p.size) {
      stop("Error: trim.e * nrow(x) must be at least the number of regressors; otherwise, the estimation results are not unique. Consider increasing trim.e.")
    }

    if (p.size > 100) {
      stop("The number of regressors cannot exceed 100.")
    }

    if (m.max > 10) {
      stop("The maximum number of breaks can not exceed 10.")
    }

    q.L = min(vec.tau) #lower bound of quantile range
    q.R = max(vec.tau) #upper bound of quantile range

    ##===========================================================
    ## Analysis begins
    ##===========================================================

    time = proc.time()  ## starting time

    s.out <- list()     #output
    m.out <- list()     #output

    ## obtain value of objective function
    out.long   = gen.long(y, x, vec.tau, n.size, trim.size)
    mat.long.s = out.long$mat.long  ## for individual quantile
    vec.long.m = out.long$vec.long  ## for multiple quantiles

    ## Single Quantile
    if (verbose){
    cat('================================================================== \n')
    cat('===== Analysis based on a single conditional quantile ===== \n')
    }

    vec.nb      = matrix(0, (N.tau*3), 1)
    mat.loc.opt = matrix(0, (N.tau*3), m.max)

    for(i in 1:N.tau){

        v.tau  = vec.tau[i]

        if (verbose){
        cat('================================================================== \n')
        cat('For quantile level: ', v.tau, '\n')
        }

        ## break dates given the maximum number of breaks
        vec.long.s = mat.long.s[,i]
        mat.date   = brdate(y, x, n.size, m.max, trim.size, vec.long.s)

        ## determine the number of breaks
        out.s = sq(y, x, v.tau, n.size, m.max, trim.size, mat.date)

        ## estimation and inference
        vec.level = c(10, 5, 1)

        if (verbose){
        cat('--- Break testing results at the', vec.level[v.a],'% significance level---\n' )
        cat('(Note: The k-th column and beyond may display zero if the (k-1)-th break is already insignificant)\n')
        }
        # Determine the number of columns dynamically
        num_cols <- length(out.s$test)  # Assuming test and cv have the same number of columns
        # Generate column names dynamically
        col_names <- paste0(seq_len(num_cols), " Breaks")

	      table_matrix <- matrix(
          	c(out.s$test, out.s$cv[v.a, ]),  # Combine SQ test and Critical values
          	nrow = 2, byrow = TRUE,
          	dimnames = list(c("SQ test", "Critical values"), col_names)  # Set row and column names
        )

	      name_test <- paste0("test_", v.tau)
	      s.out[[name_test]]<-table_matrix

	      if (verbose){
        print(format(table_matrix, digits = 4), quote = FALSE)
        cat('The number of breaks detected (SQ above the critical value): ', out.s$nbreak[v.a], '\n')
	      }

        name_break <- paste0("nbreak_", v.tau)
        s.out[[name_break]]<-out.s$nbreak[v.a]

        ## estimation
        n.break      = out.s$nbreak[v.a]
        ind01        = i + N.tau * (v.a - 1)
        vec.nb[ind01]= n.break

        if (n.break >= 1){
            vec.date = out.s$date[v.a,1:n.break]
            fit      = rq.est.full(y, x, v.tau, vec.date, n.size)
            result   = summary.rq(fit, se = "nid", covariance = TRUE)
            mat.ci   = ci.date.m(y, x, v.tau, vec.date, n.size, v.b)
            if (verbose){
	          cat('\n')
            cat('--- Estimated Break Dates and their', 100-vec.level[v.b], '% Confidence Intervals ---\n')
            }

            # number of rows = number of breaks
            p_num_breaks <- nrow(mat.ci)  # Get the number of rows dynamically
            # Generate row names dynamically
            rownames(mat.ci) <- paste0("Break ", seq_len(p_num_breaks))
            # Assign fixed column names
            colnames(mat.ci) <- c("Estimate", "CI_Lower_Bound", "CI_Upper_Bound")
            if (verbose){
            cat('\n')
            cat('(a) Breaks as indices (e.g., 50 for the 50th observation in the sample):')
            cat('\n')
            }
            name_ci <- paste0("br_est_", v.tau)
            s.out[[name_ci]]<-mat.ci
            if (verbose){
            print( mat.ci)
            }

            if (min(mat.ci[,2]) < 1 | T.size < max(mat.ci[,3])){
                warning('  confidence interval is out of the range\n')
            } else if (is.null(vec.time)) {
              # Do nothing (explicitly left empty)
            }
            else
            {
                 mat.ci.time <- matrix( vec.time[mat.ci], ncol = 3)
                # Generate row names dynamically
                rownames(mat.ci.time) <- paste0("Break ", seq_len(p_num_breaks))
                # Assign fixed column names
                colnames(mat.ci.time) <- c("Estimate", "CI_Lower_Bound", "CI_Upper_Bound")
                name_ci <- paste0("br_est_time_", v.tau)
                s.out[[name_ci]]<-mat.ci.time
                if (verbose){
                cat('\n')
                cat ('(b) Breaks in date format:')
                cat ('\n')
                print( mat.ci.time)
                }
            }

            if (verbose){
            cat ('\n')
            cat('- - - - - - Coeffcient Estimation Results - - - - - - \n')
            ## Estimates for each regime
            cat ('\n')
            cat(' (a) Coefficients estimates for each regime \n')
            cat ('[rows: intercept, first regressor, second regressor ...] \n')
            cat ('\n')
            }

            coef.est<-rq.est.regime(y, x, v.tau, vec.date, n.size)
            if (verbose){
              # Properly formatted print
              for (name in names(coef.est)) {
                cat("\n", name, "\n")
                print(format(coef.est[[name]], digits = 4), quote = FALSE)
              }
            }
            name_coef <- paste0("coef_", v.tau)
            s.out[[name_coef]]<-coef.est

            ## Estimates for the break sizes
            if (verbose){
            cat(' (b) Break sizes \n')
            cat ('[rows: break in intercept, break in first coef, break in second coef ...] \n ')
            cat ('\n')
            }

            for (j in 1:n.break){
                beg01 = 1 + p.size * j
                end01 =     p.size * (j+1)
                coef_subset <- result$coef[beg01:end01,]
                # Rename the rows
                rownames(coef_subset) <- c("Intercept", paste0("x", seq_len(nrow(coef_subset) - 1)))
                # Print the updated coefficient matrix
                if (verbose){
                cat('  < Regime', (j + 1), '  minus Regime', j, '>\n')
                print(format(coef_subset, digits = 4), quote = FALSE)
                cat('  \n')
                }
                name_coef <- paste0("bsize_", v.tau,'_Regime_', (j + 1), '_minus Regime_', j)
                s.out[[name_coef]]<-result$coef[beg01:end01,]
                rownames(s.out[[name_coef]]) <- c("Intercept", paste0("x", seq_len(nrow(result$coef[beg01:end01,]) - 1)))

            }

            mat.loc.opt[ind01,1:n.break] = vec.date
        }
    }

    ## figure for single quantile
    ## fig.s(y, x, vec.tau, vec.nb, mat.date, vec.time, n.size)

    if (length(vec.tau) > 1) {

      ## For the DQ test; d.Sym=TRUE if vec.tau specifies a symmetric quantile range
      ## i.e, q.L = (1 - q.R); Otherwise, d.Sym=FALSE; for the latter case, the critical values of the DQ test is generated via simulations
      d.Sym <- abs(q.L - (1 - q.R)) <= 1e-5
      d.Sim <- !d.Sym
      ## simulation to generate critical values, not needed if the analysis is for a single quantile only
        if ((d.Sym == FALSE & d.Sim == TRUE) | p.size>20 | m.max>5){
          if (verbose){
            cat("The critical values of the DQ test is obtained via simulations\n")
          }
          table.cv<-critical.DQtest_specific(x, m.max, vec.tau) #simulate critical value for the DQ test;
          #if not true, compute critical values from a response surface
        } else {
        table.cv<-NULL
        }
    ## Multiple quantiles
    if (verbose){
    cat('================================================================== \n')
    cat('===== Analysis based on multiple conditional quantiles  ===== \n')
    cat('================================================================== \n')
    }
    mat.date = brdate(y, x, n.size, m.max, trim.size, vec.long.m)
    ## (i) determine the number of breaks
    out.m = dq(y, x, vec.tau, q.L, q.R, n.size, m.max, trim.size, mat.date,d.Sym,table.cv)
    vec.level = c(10, 5, 1) # significance level

    ## the number of breaks
    n.break = out.m$nbreak[v.a]
    if (verbose){
        cat('----- Break testing results at the', vec.level[v.a], '% significance level: ----\n' )
    }
        # Determine the number of columns dynamically
        num_cols <- length(out.m$test)  # Assuming test and cv have the same number of columns
        # Generate column names dynamically
        col_names <- paste0(seq_len(num_cols), " Breaks")
        table_matrix <- matrix(
          	c(out.m$test, out.m$cv[v.a, ]),  # Combine test and Critical values
          	nrow = 2, byrow = TRUE,
          	dimnames = list(c("DQ test", "Critical values"), col_names)  # Set row and column names
        		       )
        if (verbose){
        print(format(table_matrix, digits = 4), quote = FALSE)
        }
        name_test <- paste0("test_joint")
        m.out[[name_test]] <- table_matrix
        if (verbose){
        cat('The number of breaks detected (DQ above the critical value): ', out.m$nbreak[v.a], '\n')
        }
        name_break <- paste0("nbreak_joint")
        m.out[[name_break]] <- out.m$nbreak[v.a]

    if (n.break >= 1){
        vec.date = out.m$date[v.a,1:n.break]
        ## confidence interval
        mat.ci = ci.date.m(y, x, vec.tau, vec.date, n.size, v.b)
            if (verbose){
              cat('\n')
              cat('--- Estimated Break Dates and their', 100-vec.level[v.b], '% Confidence Intervals:---\n')
            }
            p_num_breaks <- nrow(mat.ci)  # Get the number of rows dynamically
            # Generate row names dynamically
            rownames(mat.ci) <- paste0("Break ", seq_len(p_num_breaks))
            # Assign fixed column names
            colnames(mat.ci) <- c("Estimate", "CI_Lower_Bound", "CI_Upper_Bound")
            name_ci <- paste0("br_est_joint")
            m.out[[name_ci]] <- mat.ci
            if (verbose){
            cat ('\n')
            cat('(a) Breaks as indices (e.g., 50 for the 50th observation in the sample):')
            cat ('\n')
            print( mat.ci)
            }

        if (min(mat.ci[,2]) < 1 | T.size < max(mat.ci[,3])){
            warning('  confidence interval is out of the range\n')
        } else if (is.null(vec.time)) {
            # Do nothing (explicitly left empty)
          }
         else {

              mat.ci.time <- matrix( vec.time[mat.ci], ncol = 3)

                # Generate row names dynamically
                rownames(mat.ci.time) <- paste0("Break ", seq_len(p_num_breaks))
                # Assign fixed column names
                colnames(mat.ci.time) <- c("Estimate", "CI_Lower_Bound", "CI_Upper_Bound")

                name_ci <- paste0("br_est_joint_time")
                m.out[[name_ci]] <- mat.ci.time
                if (verbose){
                cat ('\n')
                cat('(b) Breaks in date format: \n')
                print( mat.ci.time)
                }

          }

        if (verbose){
        cat ('\n')
        cat('------------ Coefficient Estimation Results ------------\n')
        }
        for(k in 1:N.tau){
            v.tau  = vec.tau[k]
            if (verbose){
            cat('========================================================== \n')
            cat('- For quantile level: ', v.tau ,'-\n')
            cat ('\n')
            ## Estimates for each regime by splitting sample
            cat(' (a) Coefficients estimates for each regime \n')
            cat ('[rows: intercept, first regressor, second regressor...] \n')
            }
            coef.est<-rq.est.regime(y, x, v.tau, vec.date, n.size)
            name_coef <- paste0("coef")
            m.out[[name_coef]] <- coef.est

            if (verbose){
              # Properly formatted print
              for (name in names(coef.est)) {
                cat("\n", name, "\n")
                print(format(coef.est[[name]], digits = 4), quote = FALSE)
              }

            }

            if (verbose){
            ## Estimates of the break sizes
            cat(' (b) Break sizes \n')
            cat ('[break in intercept, break in first coef, break in second coef...] \n')
            }
            fit    = rq.est.full(y, x, v.tau, vec.date, n.size)
            result = summary.rq(fit, se = 'nid', covariance = TRUE)

            for (j in 1:n.break){
                beg01 = 1 + p.size * j
                end01 =     p.size * (j+1)

                coef_subset <- result$coef[beg01:end01,]

                # Rename the rows
                rownames(coef_subset) <- c("Intercept", paste0("x", seq_len(nrow(coef_subset) - 1)))

                # Print the updated coefficient matrix

                if (verbose){
                print(format(coef_subset, digits = 4), quote = FALSE)
                cat('  \n')
                }
                name_coef <- paste0("bsize_", v.tau,'_Regime_', (j + 1), '_minus Regime_', j)
                m.out[[name_coef]]<-result$coef[beg01:end01,]
                rownames(m.out[[name_coef]]) <- c("Intercept", paste0("x", seq_len(nrow(result$coef[beg01:end01,]) - 1)))
            }

        }

    }
}

    return(list(s.out = s.out, m.out= m.out))
    if (verbose){
    cat('================================================================== \n')
    cat('Estimation Time (min) \n')
    print( (proc.time() - time) / 60)
    cat('================================================================== \n')
    }
}



