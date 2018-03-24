# evalstats.R
# copied from ../SWOT/stanOutput.R
# stats on prediction, actual series --------------------------------------

RRMSE <- function(pred, meas) 
  sqrt(mean((pred - meas)^2 / meas^2))

MRR <- function(pred, meas)
  mean((meas - pred) / meas)

SDRR <- function(pred, meas)
  sd((meas - pred) / meas)

NSE <- function(pred, meas)
  1 - var(meas - pred) / var(meas)

NRMSE <- function(pred, meas)
  sqrt(mean((meas - pred)^2)) / mean(meas)

rBIAS <- function(pred, meas)
  mean(pred - meas) / mean(meas)

CoV <- function(pred, meas)
  sd(pred - meas) / mean(meas)

Ej <- function(pred, meas, j = 1, bench = mean(meas))
  1 - mean(abs(meas - pred)) / mean(abs(pred - bench))

logNSE <- function(pred, meas)
  1 - var(log(meas / pred)) / var(log(meas / mean(meas)))
