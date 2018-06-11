# BATMAN models
# Extend estA0 to other transformations, using various optimization and Bayesian tools

# This version was me overthinking MLE on a lognormal. I don't need to estimate sigma. 
defunct_batman_log <- function(swotlist, ...) {
  W <- swotlist$W
  dA <- swotlist$dA
  S <- swotlist$S
  
  X <- -2/3 * log(W) + 1/2 * log(S)
  
  times <- 1:ncol(X)
  reaches <- 1:nrow(X)
  
  objfun <- function(pars) {
    
    A0 <- pars[reaches]
    logQn <- pars[max(reaches) + times]
    sigsq <- pars[max(reaches) + max(times) + 1]
    
    A <- dA + swot_vec2mat(A0, dA)
    logQnmat <- swot_vec2mat(logQn, dA)
    
    objmat1 <- 1 / (2 * sigsq) * (X + 5/3 * log(A) - logQnmat)^2
    objmat <- objmat1 #+ sigsq / 8 #+ log(A)
    obj <- sum(objmat) + max(times) * max(reaches) / 2 * log(2 * pi * sigsq)
    obj
  }
  
  initA0 <- apply(dA, 1, function(x) max(x) - min(x))
  initlogQn <- apply(X + 5/3 *  swot_vec2mat(initA0, X), 2, mean)
  initsigsq <- 5^2
  
  inits <- c(initA0, initlogQn, initsigsq)
  
  # browser()
  
  # optres <- nlm(objfun, p = inits, ...)
  optres <- nlminb(start = inits, objfun, 
                   lower = c(rep(0, length(initA0)),
                             rep(-5, length(initlogQn)), 0), 
                   ...)
  
  ests <- optres$par
  code <- optres$convergence
  minval <- optres$objective
  
  # ests <- optres$estimate
  # code <- optres$code 
  # minval <- optres$minimum
  # 
  A0ests <- ests[reaches]
  logQnests <- ests[max(reaches) + times]
  sigsqest <- ests[max(reaches) + max(times) + 1]
  
  out <- list(A0 = A0ests, logQn = logQnests, sigsq = sigsqest, code = code, 
              obj = minval)
  out
}

batman_A0priorfun <- function(Wmat) {
  lwbar <- apply(log(Wmat), 1, mean)
  lwsd <- apply(log(Wmat), 1, sd)
  logA0hat <- -1.4058 + 1.4931 * lwbar - 0.2293 * lwsd
  logA0hat
  
  outfun <- function(A0) {
    log()
  }
}

batman_log <- function(swotlist, steptol = 1e-12, gradtol = 1e-12, ...) {
  W <- swotlist$W
  dA <- swotlist$dA
  S <- swotlist$S
  
  X <- -2/3 * log(W) + 1/2 * log(S)
  
  times <- 1:ncol(X)
  reaches <- 1:nrow(X)
  
  objfun <- function(pars) {
    
    A0 <- pars[reaches]
    logQn <- pars[max(reaches) + times]
    # sigsq <- pars[max(reaches) + max(times) + 1]
    
    A <- dA + swot_vec2mat(A0, dA)
    logQnmat <- swot_vec2mat(logQn, dA)
    
    objmat <- (X + 5/3 * log(A) - logQnmat)^2
    obj <- sum(objmat)
    
    obj
  }
  
  gradfun <- function(pars) {
    A0 <- pars[reaches]
    logQn <- pars[max(reaches) + times]
    A <- dA + swot_vec2mat(A0, dA)
    logQnmat <- swot_vec2mat(logQn, dA)
    
    dlogAmat <- 2 * (X + 5/3 * log(A) - logQnmat) * (5 / (3 * A))
    dlogAvec <- apply(dlogAmat, 1, sum)
    
    dlogQnmat <- -2 * (X + 5/3 * log(A) - logQnmat)
    dlogQnvec <- apply(dlogQnmat, 2, sum)
    
    out <- c(dlogAvec, dlogQnvec)
  }
  attr(objfun, "gradient") <- gradfun  
  
  initA0 <- apply(dA, 1, function(x) max(x) - min(x))
  initlogQn <- apply(X + 5/3 *  swot_vec2mat(initA0, X), 2, mean)

  inits <- c(initA0, initlogQn)
  
  optres <- nlm(objfun, p = inits, steptol = steptol, gradtol = gradtol, ...)
  ests <- optres$estimate
  code <- optres$code
  minval <- optres$minimum

  A0ests <- ests[reaches]
  logQnests <- ests[max(reaches) + times]

  out <- list(A0 = A0ests, logQn = logQnests, code = code, 
              obj = minval)
  out
}

batman_linA <- function(swotlist, ...) {
  W <- swotlist$W
  S <- swotlist$S
  X <- W^(-2/5) * S^(3/10)
  dA <- swotlist$dA
  minA0 <- -apply(dA, 1, min) + 1
  
  times <- 1:ncol(X)
  reaches <- 1:nrow(X)
  
  objfun <- function(pars) {
    
    # if (sum(A0 < minA0) > 0) return(1e15)
    # if (sum(Qn <= 0) > 0) return(1e15)
    
    Qn35 <- pars[times]
    A0 <- pars[max(times) + reaches]
    
    A <- dA + swot_vec2mat(A0, dA)
    Qn35mat <- swot_vec2mat(Qn35, dA)
    objmat <- (A * X - Qn35mat)^2
    obj <- sum(objmat)
    obj
    
  }
  
  initA0 <- apply(dA, 1, function(x) max(x) - min(x))
  
  initQn35 <- apply(X * (dA + swot_vec2mat(initA0, dA)), 2, mean)
  inits <- c(initQn35, initA0)
  
  optres <- nlm(objfun, p = inits, ...)
  ests <- optres$estimate
  code <- optres$code 
  
  A0ests <- ests[max(times) + reaches]
  Qnests <- ests[times] ^ (5/3)

  out <- list(A0 = A0ests, Qn = Qnests, code = code, 
              obj = optres$minimum)
  out
}


# batman_linA(sscase)
# batman_linA(reachdata$Po)


batman_linQ <- function(swotlist, steptol = 1e-12, gradtol = 1e-12, ...) {
  W <- swotlist$W
  S <- swotlist$S
  X <- W^(-2/3) * S^(1/2)
  dA <- swotlist$dA
  minA0 <- -apply(dA, 1, min) + 1
  
  times <- 1:ncol(X)
  reaches <- 1:nrow(X)
  
  objfun <- function(pars) {
    Qn <- pars[times]
    A0 <- pars[max(times) + reaches]
    
    A <- dA + swot_vec2mat(A0, dA)
    Qnmat <- swot_vec2mat(Qn, dA)
    objmat <- (A^(5/3) * X - Qnmat)^2
    obj <- sum(objmat)
    obj
    
  }
  
  initA0 <- apply(dA, 1, function(x) max(x) - min(x))
  
  initQn35 <- apply(X * (dA + swot_vec2mat(initA0, dA)), 2, mean)
  inits <- c(initQn35, initA0)
  
  optres <- nlm(objfun, p = inits, steptol = steptol, gradtol = gradtol, ...)
  ests <- optres$estimate
  code <- optres$code 
  
  A0ests <- ests[max(times) + reaches]
  Qnests <- ests[times]
  
  out <- list(A0 = A0ests, Qn = Qnests, code = code, 
              obj = optres$minimum)
  out
}


