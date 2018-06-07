# Functions for working with "simplest" McFli a la Mike's DAWG blog posts.
# 5/1/2018
# Mark Hagemann
# Some modified from 4/30 and just previous notebooks. 
# Prefix for these is "sm" as in "simplest McFli" or "simplest Mannings"

#' Generate a function for reach response (a function of discharge)
#' 
#' @param n Manning's n
#' @param bw bottom width (m)
#' @param ss side slope (unitless)
#' @param bs bottom slope (unitless)
#' 
sm_reachfun <- function(n, bw, ss, bs) {
  scalfun <- function(Q) {
    
    lhs <- Q * n * bs ^ (-1/2)
    obfun <- function(W) {
      rhs <- (ss / 4 * (W^2 - bw^2)) ^ (5/3) * W^(-2/3)
      out <- (lhs - rhs)^2
    }
    
    optw <- optimize(f = obfun, interval = c(bw, bw + 10 * bw / ss + 10000 / ss))
    W <- optw$minimum
    A <- (Q * n * bs^(-1/2) * W^(2/3))^(3/5)
    H <- (W - bw) * ss / 2
    
    out <- data.frame(A = A, W = W, H = H)
    out
  }
  vecfun <- Vectorize(scalfun, SIMPLIFY = FALSE)
  out <- function(Q) {
    dplyr::bind_rows(vecfun(Q))
  } 
}


#' Regression model matrix for simplest McFli
#' 
#' @param reachdata1 a data.frame as returned by a sm_reachfun function. 
#' @param reachdata2 a data.frame as returned by a sm_reachfun function.
#' @param A0_ref Where should A0 correspond to? Default is minimum area in each reach. 
sm_modmat <- function(reachdata1, reachdata2) {
  
  W1 <- reachdata1$W
  W2 <- reachdata2$W
  
  col1 <- W1^(-2/5)
  col2 <- -W2^(-2/5)
  
  out <- cbind(col1, col2)
  out
}

sm_respvec <- function(reachdata1, reachdata2, A0_ref = c("min", "first", "median")) {
  
  A0_ref <- match.arg(A0_ref)
  
  A0vec <- sm_A0(reachdata1, reachdata2, A0_ref = A0_ref)
  A01 <- A0vec[1]
  A02 <- A0vec[2]
  
  A1 <- reachdata1$A
  A2 <- reachdata2$A
  
  dA1 <- A1 - A01
  dA2 <- A2 - A02
  
  W1 <- reachdata1$W
  W2 <- reachdata2$W
  
  out <- dA2 * W2^(-2/5) - dA1 * W1^(-2/5)
  out
}

sm_A0 <- function(reachdata1, reachdata2, A0_ref = c("min", "first", "median")) {
  A0_ref <- match.arg(A0_ref)
  A1 <- reachdata1$A
  A2 <- reachdata2$A
  
  A01 <- ifelse(A0_ref == "min", min(A1), 
                ifelse(A0_ref == "first", A1[1], median(A1)))
  A02 <- ifelse(A0_ref == "min", min(A2), 
                ifelse(A0_ref == "first", A2[1], median(A2)))
  out <- c(A01, A02)
  out
}

#' Make a swotdata case out of two reaches' response functions.
#' 
sm_swotlist <- function(reachfun1, reachfun2, Q, 
                        A0_ref = c("min", "first", "median")) {
  A0_ref <- match.arg(A0_ref)
  
  df1 <- reachfun1(Q)
  df2 <- reachfun2(Q)
  
  Wmat <- rbind(df1$W, df2$W)
  Amat <- rbind(df1$A, df2$A)
  Hmat <- rbind(df1$H, df2$H)
  A0vec <- sm_A0(df1, df2, A0_ref = A0_ref)
  
  dAmat <- rbind(df1$A - A0vec[1], df2$A - A0vec[2])
  Qmat <- rbind(Q, Q)
  
  Svec <- c(environment(reachfun1)$bs, environment(reachfun2)$bs)
  Smat <- swot_vec2mat(Svec, Wmat)
  xmat <- swot_vec2mat(1:2, Wmat)
  out <- list(
    W = Wmat,
    S = Smat,
    H = Hmat,
    dA = dAmat,
    A = Amat,
    x = xmat,
    Q = Qmat
  )
  out
}
