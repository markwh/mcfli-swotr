
#' Create a log-likelihood function for a bamdata object. 
#' 
#' @param bamdata A bamdata object, as returned by \code{bam_data}
#' @param negative If TRUE, return negative log likelihood. 
bam_llfun <- function(bamdata, negative = FALSE) {
  
  errmatfun <- bam_errmatfun(bamdata) 
  
  llfun <- function(A0, logn, Q, sigma = 0.25) {
    errs <- errmatfun(A0, logn, Q)
    outmat <- dnorm(errs, mean = 0, sd = sigma, log = TRUE)
    out <- sum(outmat)
    if (negative) {
      out <- -out
    }
    out
  }
}

#' Same as bam_llfun, but takes a single vector argument.
#' 
#' vectorizes in the order of bam_llfun arguments: A0, logn, Q, sigma.
#' @param bamdata A bamdata object, as returned by \code{bam_data}
#' @param negative If TRUE, return negative log likelihood. 
bam_llfun_vector <- function(bamdata, negative = FALSE) {
  ns <- nrow(bamdata$Wobs)
  nt <- ncol(bamdata$Wobs)
  
  fun1 <- bam_llfun(bamdata = bamdata, negative = negative)
  
  out <- function(params) {
    stopifnot(length(params) == (ns + nt + 2))
    A0vec <- params[1:ns]
    logn <- params[ns + 1]
    Qvec <- params[(ns + 1) + 1:nt]
    sigma <- params[ns + nt + 2]
    
    return(fun1(A0 = A0vec, logn = logn, Q = Qvec, sigma = sigma))
  }
  out
}



#' Another likelihood function for bamdata.
#'
#' Same as above, but with the following changes:
#' 
#' - log-transform Q and A0
#' - separate logQ into logqbar and logqdot
#' 
#' vectorizes in the order of bam_llfun arguments: 
#'   logA0, logn, logqbar, sigma, logqdot.
#' @param bamdata A bamdata object, as returned by \code{bam_data}
#' @param negative If TRUE, return negative log likelihood. 
bam2_llfun_vector <- function(bamdata, negative = FALSE) {
  ns <- nrow(bamdata$Wobs)
  nt <- ncol(bamdata$Wobs)
  
  errmatfun <- bam_errmatfun(bamdata) 
  
  llfun <- function(params) {
    stopifnot(length(params) == (ns + nt + 3))
    A0vec <- exp(params[1:ns])
    logn <- params[ns + 1]
    logqbar <- params[ns + 2]
    sigma <- params[ns + 3]
    logqdotvec <- params[(ns + 3) + 1:nt]
    
    Q <- exp(logqbar + logqdotvec)

    errs <- errmatfun(A0vec, logn, Q)
    outmat <- dnorm(errs, mean = 0, sd = sigma, log = TRUE)
    out <- sum(outmat)
    if (negative) {
      out <- -out
    }
    out
  }
  
  llfun
}


#' Yet another likelihood function for bamdata.
#'
#' Same as above, but with the following changes:
#' 
#' - Adds a hierarchical parameter for logA0_bar
#' 
#' vectorizes in the order of bam_llfun arguments: 
#'   logA0bar, logn, logqbar, sigma, logA0dot, logqdot.
#' @param bamdata A bamdata object, as returned by \code{bam_data}
#' @param negative If TRUE, return negative log likelihood. 
bam3_llfun_vector <- function(bamdata, negative = FALSE) {
  ns <- nrow(bamdata$Wobs)
  nt <- ncol(bamdata$Wobs)
  
  errmatfun <- bam_errmatfun(bamdata) 
  
  llfun <- function(params) {
    stopifnot(length(params) == (ns + nt + 4))
    logA0bar <- params[1]
    logn <- params[2]
    logqbar <- params[3]
    sigma <- params[4]
    logA0dotvec <- params[4 + 1:ns]
    logqdotvec <- params[(ns + 4) + 1:nt]
    
    Q <- exp(logqbar + logqdotvec)
    A0vec <- exp(logA0bar + logA0dotvec)
    
    errs <- errmatfun(A0vec, logn, Q)
    outmat <- dnorm(errs, mean = 0, sd = sigma, log = TRUE)
    out <- sum(outmat)
    if (negative) {
      out <- -out
    }
    out
  }
  
  llfun
}

#' A0 and sigma only log likelihood function factory
#' 
#' Based on BAM likelihood, with q and n eliminated
#' @param bamdata A bamdata object, as returned by \code{bam_data}
#' @param negative If TRUE, return negative log likelihood. 
#' @param gradient If TRUE, returned function returns gradient as an attribute.
A0_llfun_vector <- function(bamdata, negative = FALSE, 
                            gradient = FALSE) {
  ns <- nrow(bamdata$Wobs)
  nt <- ncol(bamdata$Wobs)
  
  xmat <- 1/2 * log(bamdata$Sobs) - 2/3 * log(bamdata$Wobs)
  xbarvec <- apply(xmat, 2, mean)
  ymat <- xmat - swot_vec2mat(xbarvec, xmat)
  
  llfun <- function(params) {
    stopifnot(length(params) == (ns + 1))
    
    sigma <- params[ns + 1]
    A0vec <- params[1:ns]
    Amat <- bamdata$dAobs + swot_vec2mat(A0vec, ymat)
    logAmat <- log(Amat)
    logAbarvec <- apply(logAmat, 2, mean)
    logAmat_ctr <- logAmat - swot_vec2mat(logAbarvec, ymat)
    mumat <- ymat + 5/3 * logAmat_ctr
    
    ldens1 <- dnorm(as.vector(mumat), mean = 0, sd = sigma, log = TRUE)
    ldensvec <- ldens1 - as.vector(logAmat) # Jacobian adjustment
    llik <- sum(ldensvec)
    
    # Begin graident computation
    Mmat <- 10 / 3 * (ymat + 5/3 * logAmat_ctr)
    Vmat <- -1 / (ns * Amat)
    Pmat <- -ns * Mmat * Vmat
    
    bit1 <- as.vector(rep(1, ns) %*% Mmat %*% t(Vmat))
    bit2 <- apply(Pmat, 1, sum)
    
    A0grad1 <- -1 / (2 * sigma^2) * (bit1 + bit2)
    A0jacadj <- apply(- 1 / Amat, 1, sum) # derivative of log Jacobian
    A0grad <- A0grad1 + A0jacadj
    
    sigmagrad1 <- -ns * nt / (sigma)
    # sigmagrad1 <- 0 # because I transformed already and use z score. 
    sigmagrad2 <- sum(mumat^2) / sigma^3
    grad <- c(A0grad, sigmagrad1 + sigmagrad2)
    
    if (negative) {
      llik <- -llik
      grad <- -grad
    }
    
    if (gradient) {
      attr(llik, "gradient") <- grad
    }
    
    llik
  }
  
  llfun
}


bam_errmatfun <- function(bamdata) {
  logW <- log(bamdata$Wobs)
  logS <- log(bamdata$Sobs)
  
  lhs <- 4 * logW - 3 * logS
  
  nx <- nrow(logW)
  nt <- ncol(logW)
  
  dA <- bamdata$dAobs
  dA_adj <- dA - swot_vec2mat(apply(dA, 1, min), dA)
  

  
  out1 <- function(A0, logn, Q) {
    
    Qvec <- rep(Q, length.out = nx * nt)
    logQ <- log(matrix(Qvec, nrow = nx, byrow = TRUE))

    logA <- log(dA_adj + swot_vec2mat(A0, dA_adj))    
    rhs <- 10 * logA - 6 * logn - 6 * logQ
    
    out2 <- lhs - rhs
    out2
  }
  out1
}



