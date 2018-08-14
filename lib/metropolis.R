# metropolis.R
# Mark Hagemann
# 8/8/2018
# Based on work done in reports/A0_likelihood.Rmd.


#' Log likelihood for mass-conserved Manning equation
llfun_mm <- function(swotlist, zero = "same") {
  logW <- log(swotlist$W)
  logS <- log(swotlist$S)
  x <- 1 / 2 * logS - 2 / 3 * logW
  
  if (zero != "same") {
    swotlist$dA <- rezero_dA(swotlist$dA, zero = zero)
  }
  
  out <- function(A0, logqn, sigma) {
    logA <- log(swot_vec2mat(A0, x) + swotlist$dA)
    qnmat <- swot_vec2mat(logqn, x)
    vec <- as.vector(x + 5 / 3 * logA - qnmat)
    ll <- sum(dnorm(vec, mean = 0, sd = sigma, log = TRUE))
    ll
  }
  out
}

#' Conditional log likelihood for logA0
llfun_A0 <- function(llfun, logqn, sigma) {
  out <- function(A0) {
    llfun(A0 = A0, 
          logqn = logqn,
          sigma = sigma)
  }
}

#' Conditional log likelihood for A0 (not log-transformed)
llfun_logA0 <- function(llfun, logqn, sigma) {
  out <- function(logA0) {
    llfun(A0 = exp(logA0), 
          logqn = logqn,
          sigma = sigma)
  }
}


# Sample ------------------------------------------------------------------
ll_mm <- llfun_mm(reachdata$Ganges)

sample_logA0 <- function(inputs, state) {
  prop <- rnorm(inputs$nx, state$logA0, state$stepsize)
  ratnum <- ll_mm(exp(prop), logqn = state$qn, sigma = 1 / sqrt(state$prec)) +
    sum(dnorm(prop, inputs$logA0_hat, inputs$logA0_sd, log = TRUE))
  ratdenom <- ll_mm(state$A0, logqn = state$qn, sigma = 1 / sqrt(state$prec)) +
    sum(dnorm(log(state$A0), inputs$logA0_hat, inputs$logA0_sd, log = TRUE))
  rat <- ratnum - ratdenom
  
  accept <- log(runif(1)) < rat
  
  if (accept) {
    # cat("|")
    state <- prop
  } else {
    # cat(".")
  }
  
}





