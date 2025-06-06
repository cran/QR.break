% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/rq.est.regime.R
\name{rq.est.regime}
\alias{rq.est.regime}
\title{Regime-Specific Coefficients and Confidence Intervals Given Break Dates}
\usage{
rq.est.regime(y, x, v.tau, vec.date, n.size = 1)
}
\arguments{
\item{y}{A vector of dependent variables (\eqn{NT \times 1}).}

\item{x}{A matrix of regressors (\eqn{NT \times p}).}

\item{v.tau}{The quantile of interest.}

\item{vec.date}{A vector of estimated break dates.}

\item{n.size}{The cross-sectional sample size (\eqn{N}).}
}
\value{
A list containing the estimated coefficients for each regime.
}
\description{
This function estimates the coefficients for each regime, given the break dates.
}
\examples{
## data
data(gdp)
y        = gdp$gdp
x        = gdp[,c("lag1", "lag2")]

## quantile
v.tau = 0.8

## break date
vec.date = 146

# cross-sectional size
n.size = 1

## estimation
result = rq.est.regime(y, x, v.tau, vec.date, n.size)
print(result)

}
