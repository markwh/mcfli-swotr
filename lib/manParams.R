
#' Caclulate true Manning parameters and derivatives
#' 
#' Modified from airSWOT project. 
#' 
#' @param datalist A netcdf-derived list obtained from swotData::nc_list()
#' 
calcManParams <- function(datalist) {

  S <- datalist$S
  W <- datalist$W
  A <- datalist$A
  Q <- datalist$Q
  
  N <- 1 / Q * W^(-2/3) * A^(5/3) * S^(1/2)
  n <- median(N, na.rm = TRUE)
  
  A0 <- apply(A, 1, min, na.rm = TRUE)
  
  qvec <- apply(Q, 2, median)
  logQbar <- mean(log(Q))
  logQdot <- apply(log(Q), 2, mean) - logQbar
  
  out <- list(n = n, A0 = A0, logQbar = logQbar, logQdot = logQdot, Q = qvec)
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
