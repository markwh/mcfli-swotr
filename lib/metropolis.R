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

#' Function that decays to 1 for optimizing step size
#' 
#' @param iter chain iteration number. Later iterations will yield values closer to 1.
#' @param scale scale factor: larger values will yield early values farther from 1. 
decfun <- function(iter, scale = 5) {
  exp(scale / sqrt(iter))
}


# Sample ------------------------------------------------------------------
# ll_mm <- llfun_mm(reachdata$Ganges)

met_sample_A0 <- function(inputs, state) {
  logA0 <- log(state$A0)
  stepsize <- state$stepsize_logA0
  llfun <- inputs$llfun
  bds <- log(inputs$A0_min)
  prop <- truncnorm::rtruncnorm(inputs$nx, a = bds, b = Inf, 
                                mean = logA0, sd = state$stepsize)
  
  ratnum <- llfun(exp(prop), logqn = state$qn, sigma = sqrt(state$sigsq_epsilon)) +
    sum(dnorm(prop, inputs$logA0_hat, inputs$logA0_sd, log = TRUE))
  ratdenom <- llfun(state$A0, logqn = state$qn, sigma = sqrt(state$sigsq_epsilon)) +
    sum(dnorm(log(state$A0), inputs$logA0_hat, inputs$logA0_sd, log = TRUE))
  rat <- ratnum - ratdenom
  
  accept <- log(runif(1)) < rat
  
  if (accept) {
    # cat("|")
    logA0 <- prop
    if (runif(1) > inputs$optar)
      stepsize <- stepsize * decfun(i, scale = 5)
  } else if (runif(1) < inputs$optar) {
    stepsize <- stepsize / decfun(i, scale = 5)
    # cat(".")
  }
  
  out <- list(A0 = exp(logA0), stepsize = stepsize, accept = accept)
  attr(out$A0, "llik") <- ifelse(accept, ratnum, ratdenom)
  out
}





