# bamr.R
# Mark Hagemann
# 04/03/2018
# Functions to put into the bamr package, when I'm done here

estimate_logA0 <- function(Wobs) {
  lwbar <- apply(log(Wobs), 1, mean)
  lwsd <- apply(log(Wobs), 1, sd)
  logA0hat <- -1.782 + 1.438 * lwbar - 2.268 * lwsd
  logA0hat
}