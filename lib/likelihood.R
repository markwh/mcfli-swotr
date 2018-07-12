
#' Create a log-likelihood function for a bamdata object. 
#' 
#' @param bamdata A bamdata object, as returned by \code{bam_data}
bam_llfun <- function(bamdata) {
  
  errmatfun <- bam_errmatfun(bamdata) 
  
  llfun <- function(A0, logn, Q, sigma = 0.25) {
    errs <- errmatfun(A0, logn, Q)
    outmat <- dnorm(errs, mean = 0, sd = sigma, log = TRUE)
    out <- sum(outmat)
    out
  }
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



