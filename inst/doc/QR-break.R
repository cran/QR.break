## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment  = "#>",
  fig.width  = 7,
  fig.height = 3.6,
  fig.align  = "center"
)
options(warn = -1)   # the quantile regression fit is occasionally non-unique

## ----load---------------------------------------------------------------------
library(QR.break)

## ----gdp-data-----------------------------------------------------------------
data(gdp)
str(gdp)
head(gdp, 3)

## ----gdp-setup----------------------------------------------------------------
y <- gdp[, "gdp"]
x <- gdp[, c("lag1", "lag2")]

## ----gdp-tau------------------------------------------------------------------
vec.tau <- seq(0.20, 0.80, by = 0.15)
vec.tau

## ----gdp-run, eval = FALSE----------------------------------------------------
# res <- rq.break(y, x,
#                 vec.tau     = vec.tau,
#                 N           = 1,
#                 trim.e      = 0.15,
#                 vec.time    = gdp[, "yq"],
#                 m.max       = 3,
#                 v.a         = 2,
#                 v.b         = 2,
#                 verbose     = TRUE,
#                 norm.method = "cholesky")

## ----gdp-plot-----------------------------------------------------------------
tt  <- seq_len(nrow(gdp))
brk <- 146; lo <- 120; hi <- 147     # DQ estimate and its 95% interval

op <- par(mar = c(3.5, 4, 2.5, 1))
plot(tt, gdp$gdp, type = "n", xaxt = "n", bty = "n", ylim = c(-12, 18),
     xlab = "", ylab = "Real GDP growth (%, annualized)")
rect(lo, -12, hi, 18, col = "#EDE7F6", border = NA)
abline(h = 0, col = "#CFCFD4")
lines(tt, gdp$gdp, col = "#5A5A66", lwd = 1.4)
segments(brk, -12, brk, 15, col = "#6C4FB8", lwd = 2)
text(brk, 16.5, "  break: 1984 Q1", adj = c(0, 0.5), col = "#6C4FB8", cex = 0.85)
text(lo,  16.5, "95% CI  ",         adj = c(1, 0.5), col = "#8E7BC6", cex = 0.8)
at <- seq(2, nrow(gdp), by = 40)
axis(1, at = at, labels = sub(" Q[1-4]$", "", gdp$yq[at]), col = "#CFCFD4")
par(op)

## ----err1, error = TRUE-------------------------------------------------------
try({
rq.break(y, x, vec.tau, N = 1, trim.e = 0.2, vec.time = gdp[, "yq"],
         m.max = 6, v.a = 2, v.b = 2)
})

## ----err2, error = TRUE-------------------------------------------------------
try({
rq.break(y, x, vec.tau, N = 1, trim.e = 0.01, vec.time = gdp[, "yq"],
         m.max = 3, v.a = 2, v.b = 2)
})

## ----driver-data--------------------------------------------------------------
data(driver)
str(driver)

## ----driver-setup-------------------------------------------------------------
y        <- driver[, "bac"]
x        <- driver[, c("age", "gender", "winter")]
vec.time <- unique(driver[, "yq"])    # length T = 100, one label per quarter
length(vec.time)

## ----driver-zeros-------------------------------------------------------------
mean(driver$bac == 0)
quantile(driver$bac, c(0.50, 0.60, 0.65, 0.70, 0.80, 0.85))

## ----driver-run, eval = FALSE-------------------------------------------------
# res.d <- rq.break(y, x,
#                   vec.tau     = seq(0.70, 0.85, by = 0.05),
#                   N           = 108,
#                   trim.e      = 0.05,
#                   vec.time    = vec.time,
#                   m.max       = 3,
#                   v.a         = 2,
#                   v.b         = 2,
#                   verbose     = TRUE,
#                   norm.method = "cholesky")

## ----sq-test------------------------------------------------------------------
y <- gdp[, "gdp"]
x <- gdp[, c("lag1", "lag2")]

sq.test.0vs1(y, x, v.tau = 0.8, n.size = 1)

## ----dq-test------------------------------------------------------------------
dq.test.0vs1(y, x, q.L = 0.2, q.R = 0.8, n.size = 1)

## ----res-surface--------------------------------------------------------------
res.surface(p = 3, l = 0, q.L = 0.2, q.R = 0.8, d.Sym = TRUE)   # 10%, 5%, 1%

## ----seq-tests----------------------------------------------------------------
sq.test.lvsl_1(y, x, v.tau = 0.8, n.size = 1, vec.date = 146)
dq.test.lvsl_1(y, x, q.L = 0.2, q.R = 0.8, n.size = 1, vec.date = 146)

## ----est-regime---------------------------------------------------------------
rq.est.regime(y, x, v.tau = 0.8, vec.date = 146, n.size = 1)

## ----ci-date------------------------------------------------------------------
ci.date.m(y, x, vec.tau = 0.8, vec.date = 146, n.size = 1, v.b = 2)

