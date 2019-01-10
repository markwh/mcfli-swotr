
#' Caclulate true Manning parameters and derivatives
#' 
#' Modified from airSWOT project. 
#' 
#' @param swotlist A list of SWOT observables
#' 
calcManParams <- function(swotlist, A0ref = c("minimum", "median")) {
  
  A0ref <- match.arg(A0ref)
  
  S <- swotlist$S
  W <- swotlist$W
  A <- swotlist$A
  Q <- swotlist$Q
  
  N <- 1 / Q * W^(-2/3) * A^(5/3) * S^(1/2)
  n <- geomMean(N, na.rm = TRUE)
  
  sigma <- sd(log(N) - log(n), na.rm = TRUE)
  
  if (A0ref == "minimum") {
    A0 <- apply(A, 1, min, na.rm = TRUE)
  } else if (A0ref == "median") {
    A0 <- apply(A, 1, median, na.rm = TRUE)
  }

  qvec <- apply(Q, 2, median)
  logQbar <- mean(log(Q))
  logQdot <- apply(log(Q), 2, mean) - logQbar
  
  out <- list(n = n, A0 = A0, logQbar = logQbar, logQdot = logQdot, Q = qvec, 
              sigma = sigma)
  out
}

manningN <- function(A, W, S, Q, Aexp = 5/3, Wexp = -2/3, Sexp = 1/2) {
  out <- A^Aexp * W^Wexp * S^Sexp / Q
  out
}

manningN_list <- function(datalist, Aexp = 5/3, Wexp = -2/3, Sexp = 1/2) {
  A <- datalist$A
  W <- datalist$W
  S <- datalist$S
  Q <- datalist$Q
  
  out <- A^Aexp * W^Wexp * S^Sexp / Q
  out
}

manningQ <- function(A, W, S, n, Aexp = 5/3, Wexp = -2/3, Sexp = 1/2) {
  out <- A^Aexp * W^Wexp * S^Sexp / n
  out
}

manningQ_list <- function(datalist, n, Aexp = 5/3, Wexp = -2/3, Sexp = 1/2) {
  A <- datalist$A
  W <- datalist$W
  S <- datalist$S
  
  out <- A^Aexp * W^Wexp * S^Sexp / n
  out
}
