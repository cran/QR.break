# QR.break
## Structural Breaks in Quantile Regression

## Overview

Methods for detecting structural breaks, determining the number of breaks, and estimating break locations in linear quantile regression. The package implements approaches based on Qu (2008) and Oka and Qu (2011), supporting both single and multiple quantiles analysis for time series and repeated cross-sectional data.

## Installation

You can install the development version from GitHub with:

```r
install.packages("QR.break")
```

## Example

This basic example shows how to detect structural breaks in a quantile regression model:

```r
library(QR.break)

# Example 1
# Time series example, using US GDP data
# Data
data(gdp)

y = gdp[,"gdp"]
x = gdp[, c("lag1", "lag2")]

# Quantiles
vec.tau = seq(0.20, 0.80, by = 0.150)
N = 1
trim.e = 0.15
vec.time = gdp[,"yq"]
m.max = 3
v.a = 2
v.b = 2
verbose = TRUE #print
options(warn=-1) #sometimes fit is non-unique

## Structural breaks in quantile regression
result = rq.break(y, x, vec.tau, N, trim.e, vec.time, m.max, v.a, v.b, verbose)

print(result)



## Example 2
## Repeated cross-section example, using youth driving data
data(driver)
Driving_data<-driver
y <- Driving_data[,"bac"]
x <- Driving_data[, c("age", "gender", "winter")]
vec.tau = seq(0.70, 0.85, 0.05)
N <- 108
trim.e <- 0.05
vec.time <- unique(Driving_data[,"yq"])
m.max <- 3
v.a <-2
v.b <-2
verbose = TRUE #print
options(warn=-1) #sometimes fit is non-unique
result <- rq.break(y, x, vec.tau, N, trim.e, vec.time, m.max, v.a, v.b, verbose)


```

## Features

- Structural break detection for quantile regression models
- Support for single or multiple quantiles
- Methods to determine the optimal number of breaks
- Estimation of break locations
- Compatible with time series and cross-sectional data


## References

- Qu, Z. (2008). "Testing for Structural Change in Regression Quantiles". Journal of Econometrics, 146(1), 170-184.
- Oka, T., & Qu, Z. (2011). "Estimating Structural Changes in Regression Quantiles". Journal of Econometrics, 162(2), 248-267.

